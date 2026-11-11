"""
    SMDFRHRM

Algebraic replica of the 12-node high-resolution model (HRM) of the SMDFR, the "plant" of
the digital-twin loop (`models/SMDFR_hrm_MCDET_2022b_control_CMEX.slx`). Stage F1b of
`docs/monolithic_formulation_plan.md` (§8.6). Every term is traced to its block of the
`.slx`; the structure was read without MATLAB with `tools/slx_dump.py`.

The HRM is the ROM of `SMDFRROM.jl` with the thermal-hydraulics split into 12 axially
chained nodes:

  * point kinetics with circulating precursors (6 groups)  `MATLAB Function`, `…1`, `…2`
  * reactivity from node-12 outlet / inlet temperatures      `MATLAB Function3`
  * power P = P₀N/N₀ distributed with Π₀                     `MATLAB Function13`, `…14`
  * fuel, wall and coolant energy balances per node          `MATLAB Function5`, `…11`, `…12`
  * outlet temperatures chained node to node (Tₑ = 2T − T_in) `MATLAB Function15`, `…16`
  * coolant inventory shared by mass, flow and h_wc ∝ ṁ^0.8   `scalar to vector3`, `Sum`,
                                                             `MATLAB Function18`, `…19`
  * temperature-dependent lead c_p (Gurvich, 1991)           `MATLAB Function7`
  * leak of Park et al. (2021), clipped at 0                  `Leak flow rate, Park 2021`, `MATLAB Function4`
  * PCT = maxᵢ T_w,i (no p₄ regression, unlike the ROM)       `tp_function`

State (50): x = [N; C₁…C₆; Cᵉ₁…Cᵉ₆; T_f,1…12; T_w,1…12; T_c,1…12; ΔM]
Simulink integrates with `ode1be` (backward Euler, Δt = 0.5 s, ONE Newton iteration).

Actuator chain (`SMDFR Units/Subsystem`): the ACS flow is CurrentToMassFlow(I) of the
current selected by three switches — the MPC current `u_ACS` when the root inport
`active_mpc` > 0 and the DET scheduler has activated the protection, otherwise the
output of a discrete PI controller on 1144 − PCT. The activation logic (C S-functions
`msfcn_limintm_hrm_c`, `msfcn_schedule_hrm_c`) is replicated in `DETLogic`; component
failures (`msfcn_indirectps_hrm_c`, rates ~1e-8 s⁻¹) are not modelled: availability = 1.
"""
module SMDFRHRM

using LinearAlgebra, ForwardDiff, MAT

export HRMParams, hrm_params, initial_state, rhs, algebraic, pct, leak_flow,
       current_to_massflow, step_ode1be, ExtInput, value, DETLogic, simulate,
       NX, NN, IDX

const NN = 12
const NX = 13 + 3NN + 1
const IDX = (N = 1, C = 2:7, Ce = 8:13, Tf = 14:25, Tw = 26:37, Tc = 38:49, dM = 50)

"""Model parameters: `matlab/parameters_dualfluid_hrm.m` plus the model workspace script
`matlab/simulation_core_dualfluid_FI3_sim_hrm_extract_MCDET.m` (overrides from the loop
script `DET_scheduler_v1_perc_online_LJ.m` passed as keywords to `hrm_params`)."""
Base.@kwdef struct HRMParams
    # neutronics
    Lambda::Float64 = 1.05e-6
    beta::Float64 = 3.02e-3
    beta_i::Vector{Float64} = [7.91e-5, 7.03e-4, 5.04e-4, 1.17e-3, 4.57e-4, 1.1e-4]
    lambda_i::Vector{Float64} = [1.27e-2, 3e-2, 1.1e-1, 3.19e-1, 1.18, 7.02]
    lambda_c::Float64 = 0.184
    lambda_e::Float64 = 0.1
    N_0::Float64 = 1.5e13
    P_0::Float64 = 1e8
    Pi_0::Vector{Float64}
    rho_0::Float64 = 1.11e-3
    alfa_f::Float64 = -2.08e-4
    alfa_c::Float64 = -8.28e-6
    rhoex::Float64 = 0.0                          # rhoex_value (model workspace)
    C_i0::Vector{Float64} = [4.29e14, 4.15e15, 3.83e15, 1.16e16, 5.7e15, 1.53e15]
    Ce_i0::Vector{Float64} = [7.01e14, 5.88e15, 3.36e15, 5.09e15, 8.20e14, 3.96e13]
    # thermal-hydraulics, per node
    A_fw::Vector{Float64} = [13.6; fill(9.0, 10); 13.6]
    A_wc::Vector{Float64} = [10.9; fill(10.3, 10); 10.9]
    mdot_f::Vector{Float64} = fill(327.7, NN)
    mdot_c::Vector{Float64} = fill(5550.5, NN)
    M_f::Vector{Float64} = [307.7; fill(116.4, 10); 307.7]
    M_w::Vector{Float64} = [39.3; fill(31.1, 10); 39.3]
    M_c::Vector{Float64} = [219.9; fill(770.5, 10); 219.9]
    c_pf::Vector{Float64} = fill(400.0, NN)
    c_pw::Vector{Float64} = fill(690.0, NN)
    h_fw::Vector{Float64} = [6546.0; fill(4305.1, 10); 7629.5]
    h_wc::Vector{Float64} = [27405.4; fill(13833.4, 10); 32177.3]
    # steady state of SA_DualFluid (initial condition and reference temperatures)
    Tf_0::Vector{Float64}; Tw_0::Vector{Float64}; Tc_0::Vector{Float64}
    Te_f0::Vector{Float64}; Te_c0::Vector{Float64}
    Tbar_f0::Float64; Tbar_c0::Float64
    Q0::Float64 = 8.2564e7 * 0.578                # `MATLAB Function15`
    # LOCA (Park et al., 2021) and loop settings
    fault_time::Float64 = 500.0
    mdot_c_drift::Float64 = 10.0
    A::Float64 = 0.3
    B::Float64 = 10.0
    noise::Float64 = 0.0                          # `Constant4` into `MATLAB Function4`
    # actuator chain
    pid_P::Float64 = -0.05; pid_I::Float64 = -0.01
    pid_umin::Float64 = 0.0; pid_umax::Float64 = 60.0
    pct_setpoint::Float64 = 1144.0                # `Constant9` of `Subsystem` (PI error)
    monitor_band::Tuple{Float64,Float64} = (1180.0, 1220.0)  # `Constant4` (limintm `a`)
    monitor_reset::Float64 = 1144.0               # `Constant3` (limintm lb_ps)
    monitor_start::Float64 = 200.0                # `Constant26`: monitor input is 1100 K before
end

"""
    hrm_params(; sa_file, kwargs...)

Build the parameters. The steady state comes from `models/SA_DualFluid.mat` (as
`parameters_dualfluid_hrm.m` does). `Tbar_f0 = mean(Te_f0(12), Tin_f0)` in MATLAB passes
`Tin_f0` as the DIMENSION argument, so it equals `Te_f0(12)` (same for the coolant);
this is reproduced, as for the ROM (plan §3.1, point 2).
"""
function hrm_params(; sa_file = joinpath(@__DIR__, "..", "models", "SA_DualFluid.mat"), kw...)
    sa = matread(sa_file)
    v(k) = vec(Float64.(sa[k]))
    Pi_0 = [1.58e-1, 6.42e-2, 6.98e-2, 7.1e-2, 7.14e-2, 7.12e-2, 7.07e-2, 7e-2,
            6.91e-2, 6.76e-2, 6.49e-2, 1.52e-1]
    Pi_0 = Pi_0 .+ (1 - sum(Pi_0)) / 12               # "Update 21-04-2024" normalisation
    Te_f0 = v("Tf_e"); Te_c0 = v("Tc_e")
    HRMParams(; Pi_0, Tf_0 = v("Tf"), Tw_0 = v("Tw"), Tc_0 = v("Tc"), Te_f0, Te_c0,
              Tbar_f0 = Te_f0[12], Tbar_c0 = Te_c0[12], kw...)
end

"""Initial state of the integrators (`Int-*` initial conditions; `Int-mdot_c` starts at 0)."""
function initial_state(p::HRMParams)
    x = zeros(NX)
    x[IDX.N] = p.N_0; x[IDX.C] = p.C_i0; x[IDX.Ce] = p.Ce_i0
    x[IDX.Tf] = p.Tf_0; x[IDX.Tw] = p.Tw_0; x[IDX.Tc] = p.Tc_0
    x
end

"""Lead heat capacity [J/(kg K)], `MATLAB Function7` (Gurvich, 1991)."""
cp_lead(T) = 175.1 - 4.961e-2T + 1.985e-5T^2 - 2.099e-9T^3 - 1.524e6 / T^2

"""Leak flow [kg/s] at time t: `Leak flow rate, Park 2021` then `MATLAB Function4`.
`Gain` multiplies A by `mdot_c_drift`, so ṁ_leak = −d + A·d·(t − t_f)^(1/B), clipped ≤ 0."""
function leak_flow(t, p::HRMParams)
    t >= p.fault_time || return 0.0
    leak = -p.mdot_c_drift + (t - p.fault_time)^(1 / p.B) * (p.A * p.mdot_c_drift)
    min(leak + p.noise, 0.0)
end

"""Algebraic signals of `SMDFR 1` as functions of the state."""
function algebraic(x::AbstractVector{T}, p::HRMParams) where {T}
    Tf, Tw, Tc = x[IDX.Tf], x[IDX.Tw], x[IDX.Tc]
    dM = x[IDX.dM]
    Mtot = sum(p.M_c)
    M_act = p.M_c .+ dM .* p.M_c ./ Mtot             # `scalar to vector3` + `Sum`
    mdot_c = p.mdot_c .* M_act ./ p.M_c               # `MATLAB Function18`
    h_wc = p.h_wc .* (mdot_c ./ 5550) .^ 0.8          # `MATLAB Function19`
    cp = cp_lead.(Tc)                                 # `MATLAB Function7`
    # `MATLAB Function15`: Tin_c uses Te_c0(12) (the vector is overwritten afterwards)
    Q = p.Q0 * mdot_c[1] / (5550 / 2) - (5500 - mdot_c[1])
    Tin_c = p.Te_c0[12] - Q / ((sum(mdot_c) / NN) * cp[1])
    Te_c = similar(Tc, promote_type(T, typeof(Tin_c)))
    Te_c[1] = 2Tc[1] - Tin_c
    for i in 2:NN
        Te_c[i] = 2Tc[i] - Te_c[i-1]
    end
    # `MATLAB Function16`: fuel inlet fixed at Te_f0(12)
    Tin_f = p.Te_f0[12]
    Te_f = similar(Tf)
    Te_f[1] = 2Tf[1] - Tin_f
    for i in 2:NN
        Te_f[i] = 2Tf[i] - Te_f[i-1]
    end
    # `MATLAB Function3`
    rho = p.rho_0 + p.alfa_f * ((Te_f[12] + Tin_f) / 2 - p.Tbar_f0) +
          p.alfa_c * ((Te_c[12] + Tin_c) / 2 - p.Tbar_c0) + p.rhoex
    P = p.P_0 / p.N_0 * x[IDX.N]                      # `MATLAB Function13`
    Pi = P .* p.Pi_0                                  # `MATLAB Function14`
    (; M_act, mdot_c, h_wc, cp, Tin_c, Te_c, Tin_f, Te_f, rho, P, Pi)
end

"""Peak cladding (wall) temperature, `tp_function`: maxᵢ T_w,i."""
pct(x) = maximum(@view x[IDX.Tw])

"""Time derivative ẋ = f(x, ṁ_leak, ṁ_ACS)."""
function rhs(x::AbstractVector{T}, leak, acs, p::HRMParams) where {T}
    a = algebraic(x, p)
    N, C, Ce = x[IDX.N], x[IDX.C], x[IDX.Ce]
    Tf, Tw, Tc = x[IDX.Tf], x[IDX.Tw], x[IDX.Tc]
    dx = similar(x, promote_type(T, typeof(leak), typeof(acs), eltype(a.Te_c)), NX)
    dx[IDX.N] = (a.rho - p.beta) * N / p.Lambda + sum(p.lambda_i .* C)          # Function2
    dx[IDX.C] = p.beta_i ./ p.Lambda .* N .- (p.lambda_i .+ p.lambda_c) .* C .+
                p.lambda_e .* Ce                                                # Function
    dx[IDX.Ce] = p.lambda_c .* C .- (p.lambda_i .+ p.lambda_e) .* Ce            # Function1
    for i in 1:NN
        Tf_up = i == 1 ? a.Tin_f : a.Te_f[i-1]
        Tc_up = i == 1 ? a.Tin_c : a.Te_c[i-1]
        # `MATLAB Function5` (node 1 uses the node-1 parameters, as written in the block)
        dx[IDX.Tf[i]] = (a.Pi[i] + p.mdot_f[i] * p.c_pf[i] * (Tf_up - a.Te_f[i]) -
                         p.h_fw[i] * p.A_fw[i] * (Tf[i] - Tw[i])) / (p.M_f[i] * p.c_pf[i])
        # `MATLAB Function11`
        dx[IDX.Tw[i]] = (p.h_fw[i] * p.A_fw[i] * (Tf[i] - Tw[i]) -
                         a.h_wc[i] * p.A_wc[i] * (Tw[i] - Tc[i])) / (p.M_w[i] * p.c_pw[i])
        # `MATLAB Function12` (actual inventory and c_p(T_c))
        dx[IDX.Tc[i]] = (a.mdot_c[i] * a.cp[i] * (Tc_up - a.Te_c[i]) +
                         a.h_wc[i] * p.A_wc[i] * (Tw[i] - Tc[i])) / (a.M_act[i] * a.cp[i])
    end
    dx[IDX.dM] = leak + acs                                                     # Sum1, Int-mdot_c
    dx
end

"""ACS pump map I [A] → ṁ [kg/s] (`CurrentToMassFlow.m`, ECM of Nowak et al., 2022)."""
function current_to_massflow(I)
    rho = 10449.0; Hd = 0.515; Wd = 0.515; S = Wd * Hd; L = 1.0; d = 0.001
    Bm = 1.0; DeltaP = 15000.0; f = 0.02; q0_LM = 0.118; q0_w = 1.2; K = 1.0
    R_LM = q0_LM * Wd / (Hd * L); R_w = q0_w * Wd / (2d * L); R_b = q0_LM * K^2 * Hd / d
    R_ver = R_w * R_b / (R_w + R_b); R_tot = R_LM + R_ver
    C = Bm * R_ver * Hd^2 / (R_tot * d); B = Bm^2 * R_tot * Hd^2 / d
    A = f * rho * L * (Wd + Hd) / (4S^3)
    disc = B^2 - 4A * (DeltaP - C * I)
    disc < 0 && return 0.0
    Q = max((-B + sqrt(disc)) / (2A), 0.0)
    max(rho * Q, 0.0)
end

"""One `ode1be` step: backward Euler with a single Newton iteration from x_n."""
function step_ode1be(x, leak, acs, h, p::HRMParams)
    f(z) = rhs(z, leak, acs, p)
    J = ForwardDiff.jacobian(f, x)
    x .+ (I - h * J) \ (h .* f(x))
end

"""
    ExtInput(t, cols)

Root-inport data (`setExternalInput`) with Simulink's default "Interpolate data": linear
interpolation inside the data and LINEAR EXTRAPOLATION beyond the last sample. The
extrapolation matters: each loop query re-simulates the HRM up to t_k with data only up
to t_{k-1} (`DET_scheduler_v1_perc_online_LJ.m`).
"""
struct ExtInput
    t::Vector{Float64}
    cols::Matrix{Float64}
end
function value(e::ExtInput, t, j)
    ts = e.t; n = length(ts)
    n == 1 && return e.cols[1, j]
    if t <= ts[1]
        return e.cols[1, j]
    elseif t >= ts[n]
        i = n - 1
    else
        i = searchsortedlast(ts, t)
        i == n && (i = n - 1)
    end
    w = (t - ts[i]) / (ts[i+1] - ts[i])
    (1 - w) * e.cols[i, j] + w * e.cols[i+1, j]
end

"""Discrete states of the activation logic and the PI controller (`SMDFR Units/Subsystem`)."""
Base.@kwdef mutable struct DETLogic
    xD::Float64 = 0.0        # limintm: current band state
    smax::Float64 = 0.0      # limintm: max band state since last reset (output y2)
    act::Float64 = 0.0       # schedule: act_PS1 (protection activated)
    xI::Float64 = 0.0        # PI integrator (Forward Euler)
    sw1_prev::Float64 = 0.0  # previous Switch1 output (PI external reset, falling edge)
end

"""
Major-step update at time t with the state x (Simulink order: outputs, then updates).
Returns the PI output held during the next step and the Switch1 value seen in minor steps.
"""
function major_step!(L::DETLogic, t, x, p::HRMParams, h)
    tw = pct(x)
    sw1 = L.act > 0 ? 1.0 : 0.0                        # Switch1: Step (=1) if act_PS1 > 0
    if L.sw1_prev > 0 && sw1 <= 0                      # PI ExternalReset = falling
        L.xI = 0.0
    end
    L.sw1_prev = sw1
    e = sw1 > 0 ? p.pct_setpoint - tw : 0.0            # Switch7 / Add1
    u_pid = clamp(p.pid_P * e + L.xI, p.pid_umin, p.pid_umax)
    # updates
    smax_out = L.smax                                  # limintm output y2 (before update)
    u_mon = (t - p.monitor_start) > 0 ? tw : 1100.0    # Switch2
    a1, a2 = p.monitor_band
    new = u_mon <= a1 ? 0 : (u_mon <= a2 ? 1 : 2)      # u > a2 would stop the simulation
    smax = max(new, Int(L.smax))
    if u_mon < p.monitor_reset && smax > 0
        smax = 0; new = 0
    end
    L.xD = new; L.smax = smax
    demand = smax_out == 1.0                           # schedule: tres_demand_ps = 1 (`Constant12`)
    L.act = demand ? 1.0 : 0.0                         # availability = 1 (no failures)
    L.xI += p.pid_I * h * e
    (; u_pid, sw1_minor = L.act > 0 ? 1.0 : 0.0, tw)
end

"""
    simulate(p, ext, t_end; h = 0.5, x0 = initial_state(p))

Simulate the HRM from t = 0 to `t_end` with the root inports `ext` (columns: u_ACS [A],
active_mpc). Timing of the inputs inside step n → n+1, identified against the archived
Simulink run (`julia/f1b_validate_hrm.jl`):
  * root inports (u_ACS, active_mpc) at t_{n+1}, linearly inter-/extrapolated — as found
    for the ROM in F1, and confirmed here by the extrapolation spike at 1160.5 s;
  * the leak, produced by Clock-driven MATLAB Function blocks, at t_n (with t_{n+1} the
    replica drifts by ½·ṁ_leak·Δt of inventory, ≈ 0.03 K of PCT);
  * discrete outputs (PI, activation logic) held from the major step at t_n.
Returns a NamedTuple of logs on the grid t0:h:t_end (`pct`, `leak`, `Te_c12`, `acs`,
`current`, `act`), the final state `x` and the final activation/PI state `logic`.
To continue a run, pass its `x` and `logic` with `t0` (the logic is copied, not mutated).
"""
function simulate(p::HRMParams, ext::ExtInput, t_end; h = 0.5, t0 = 0.0,
                  x0 = initial_state(p), logic::DETLogic = DETLogic())
    K = round(Int, (t_end - t0) / h)
    x = copy(x0); L = deepcopy(logic)
    t = t0 .+ h .* (0:K)
    pctv = zeros(K + 1); leak = zeros(K + 1); Tec12 = zeros(K + 1)
    acs = zeros(K + 1); cur = zeros(K + 1); actv = zeros(K + 1)
    for n in 0:K
        tn = t[n+1]
        pctv[n+1] = pct(x); leak[n+1] = leak_flow(tn, p)
        Tec12[n+1] = algebraic(x, p).Te_c[12]
        n == K && break
        ms = major_step!(L, tn, x, p, h)
        t1 = t[n+2]
        a = value(ext, t1, 2); u = value(ext, t1, 1)   # Switch4/5/6 in minor steps
        I_acs = a > 0 ? (ms.sw1_minor > 0 ? u : 0.0) : ms.u_pid
        m_acs = current_to_massflow(I_acs)
        acs[n+2] = m_acs; cur[n+2] = I_acs; actv[n+2] = ms.sw1_minor
        # leak (Clock-driven MATLAB Functions) at t_n; root inports at t_{n+1} (see docstring)
        x = step_ode1be(x, leak_flow(tn, p), m_acs, h, p)
    end
    (; t = collect(t), pct = pctv, leak, Te_c12 = Tec12, acs, current = cur, act = actv, x,
       logic = L)
end

end # module
