"""
    ClosedLoop

Closed-loop test bench of F5a (`docs/monolithic_formulation_plan.md` §8.6): the MPC
variants act on the HRM replica (`SMDFRHRM.jl`) as the plant, without Simulink.

Timing, as in the archived loop (`DET_scheduler_v1_perc_online_LJ.m`): at query t_q the
ROM is initialised from the HRM state at t_q − 10 s and warmed up to τ = t_q − Δt; the
MPC decides the flows ṁ₁…ṁ_K of the steps ending at t_q, …, t_q + (K−1)Δt; ṁ₁ is applied
to the plant over (τ, τ + Δt] (root inport value at t_q). Unlike the archived loop, the
plant is simulated once, forward, so the inport extrapolation artefact of §8.7 does not
appear.

Controllers (all on the ROM, horizon K = 26):
  :N3   monolithic NLP (F4Lib.solve_det, Ipopt), warm-started with the shifted solution;
  :N1   QP linearised once at the query point (condensed, Clarabel);
  :N2   real-time iteration: one QP linearised along the ROM trajectory of the shifted
        previous solution (ode1be);
  :LJ   replay of the archived Luus–Jaakola currents (only 1160–1179.5 s);
  :zero no ACS; :max constant ṁ_max.
`bias = true` adds the offset-free correction d = PCT_HRM(τ) − ŷ_ROM(τ) to the prediction.
"""
module ClosedLoop

include(joinpath(@__DIR__, "HRMQueries.jl"))
include(joinpath(@__DIR__, "F4Lib.jl"))
using .HRMQueries, .F4Lib, MAT, JuMP, Clarabel, ForwardDiff, LinearAlgebra, Statistics, Random

export Bench, setup, run_loop, metrics, thesis_metrics, time_to_limit, ACT_ARCH, ACT_PHYS, act_cap

const ROOT = joinpath(@__DIR__, "..")
const h = 0.5
const K = 26
const DELAY = 10.0
const R4 = F4Lib.SMDFRROM

const ACT_ARCH = (name = :arch, m_lo = NaN, m_hi = NaN, ramp = :soft)   # filled in setup
const ACT_PHYS = (name = :phys, m_lo = 0.0, m_hi = 5543.0, ramp = :hard)
"""Finite-capacity actuator of F5b.1a (plan §8.12): ṁ ∈ [0, m_max], hard ramp `dm` per step.
Use it with a plant built with `hrm_params(…; acs_max = m_max, acs_ramp = dm)`."""
act_cap(m_max, dm) = (name = :cap, m_lo = 0.0, m_hi = Float64(m_max), ramp = :hard, dm = Float64(dm))

"""Copy of the HRM parameters with some fields replaced."""
with_params(pH; kw...) = SMDFRHRM.HRMParams((haskey(kw, f) ? kw[f] : getfield(pH, f)
                                             for f in fieldnames(SMDFRHRM.HRMParams))...)

"""Everything a closed-loop run needs: plant parameters, the plant history before the
first query, the ROM parameters and the archived references."""
struct Bench
    pH; pR; pR4; y_ref::Float64; act_arch; cs::Matrix{Float64}
    t_start::Float64
    hist::Dict{Float64,Any}     # HRM (x, logic, pct) on the grid before t_start
    t_act::Float64              # ACS activation by the DET logic (NaN if unknown)
end

"""
    setup(; t_start = 1160.0, pH = nothing, lead = 17.0, t_max = 2500.0)

Archived realisation (default) or a new one with HRM parameters `pH` (leak A, B, drift).
For a new realisation the MPC starts `lead` s after the DET logic activates the ACS (the
archived loop: activation 1143 s, first query 1160 s); returns `nothing` if the ACS is never
demanded before `t_max`.
"""
function setup(; t_start = 1160.0, pH = nothing, lead = 17.0, t_max = 2500.0, cap = nothing)
    ref = matread(joinpath(ROOT, "results", "f1_rom_replica", "rom_reference_1180.mat"))
    ub = ref["u_bounds"]
    act_arch = (name = :arch, m_lo = current_to_massflow(minimum(ub[:, 1])),
                m_hi = current_to_massflow(maximum(ub[:, 2])), ramp = :soft)
    t_act = NaN
    if pH === nothing
        arc = matread(joinpath(ROOT, "results", "f1b_hrm_replica", "archive_1180.mat"))
        cs = arc["control_series"]
        pH = hrm_params(A = arc["sampled_A"], B = arc["sampled_B"], mdot_c_drift = arc["sampled_mdot_c_drift"])
        # F5b.1b0: archived realisation on a finite-capacity plant (`cap = (acs_max, acs_ramp)`)
        cap === nothing || (pH = with_params(pH; acs_max = cap.acs_max, acs_ramp = cap.acs_ramp))
    else
        # activation time: the ACS capacity acts only after the activation, so detect it on the
        # unlimited plant (a capped PI may exhaust the inventory before t_max)
        # (and without jumps of the break: they are set after the first query, hence after t_act)
        pA = with_params(pH; acs_max = Inf, acs_ramp = Inf, leak_jumps = Tuple{Float64,Float64}[])
        s = SMDFRHRM.simulate(pA, ExtInput([0.0, 1.0], zeros(2, 2)), t_max)
        ia = findfirst(>(0), s.act)
        ia === nothing && return nothing
        t_act = s.t[ia]
        t_start = ceil((t_act + lead) / h) * h
        tt = collect(0.0:h:(t_start - h))
        cs = vcat([0.0 0.0 0.0], hcat(tt, zeros(length(tt), 2)))
    end
    # plant history: no MPC before t_start (the archived series is zero there)
    last = findfirst(==(t_start - h), cs[:, 1])
    ext0 = ExtInput(cs[2:last, 1], cs[2:last, 2:3])
    t0 = t_start - DELAY - h
    s = SMDFRHRM.simulate(pH, ext0, t0)
    hist = Dict{Float64,Any}(t0 => (s.x, s.logic, s.pct[end]))
    for t in t0:h:(t_start - 2h)
        x, L, _ = hist[t]
        s2 = SMDFRHRM.simulate(pH, ext0, t + h; t0 = t, x0 = x, logic = L)
        hist[t+h] = (s2.x, s2.logic, s2.pct[end])
    end
    Bench(pH, rom_params(), R4.load_params(ref["par"]), Float64(ref["target"]), act_arch,
          cs, t_start, hist, t_act)
end

cost_step(y, y_ref) = (y - y_ref)^2

# ---------------------------------------------------------------- leak forecast (F5b.0)
"""
Local-linear-trend Kalman filter of the leak (plan §8.12, step 2), without Park's law:
level ℓ and slope s per step, ℓ⁺ = ℓ + s + w₁, s⁺ = s + w₂, y = ℓ + v, with
Q = diag(q_l, q_s), R = r (calibrated offline: `python/f5b0_leak_forecast.py`).
Intervention: if |ν|/√S > `thr`, a jump of the break is absorbed in the level (ℓ = y,
P_ℓℓ += ν²), the slope is kept. `thr = Inf` disables it.
"""
mutable struct LeakKF
    x::Vector{Float64}; P::Matrix{Float64}; ql::Float64; qs::Float64; r::Float64; thr::Float64
end
LeakKF(y0; q_l, q_s, r, thr = 5.0) = LeakKF([y0, 0.0], [r+1.0 0.0; 0.0 1.0], q_l, q_s, r, thr)

function kf_update!(f::LeakKF, y)
    F = [1.0 1.0; 0.0 1.0]
    f.x = F * f.x; f.P = F * f.P * F' + Diagonal([f.ql, f.qs])
    S = f.P[1, 1] + f.r; ν = y - f.x[1]
    Kg = f.P[:, 1] ./ S
    f.x = f.x .+ Kg .* ν; f.P = f.P - Kg * f.P[1:1, :]
    if abs(ν) / sqrt(S) > f.thr
        f.x[1] = y; f.P[1, 1] += ν^2
    end
    f
end

"""Forecast L̂(τ + jΔt), j = 1…K, and its covariance Σ_L (K×K)."""
function kf_forecast(f::LeakKF, K)
    m = [f.x[1] + f.x[2] * j for j in 1:K]
    Σ = [dot([1.0, i], f.P * [1.0, j]) + sum((f.ql + (i - l) * (j - l) * f.qs for l in 1:min(i, j)); init = 0.0)
         for i in 1:K, j in 1:K]
    (m = m, Σ = Σ)
end

"""Kalman parameters calibrated offline (results/f5b0_leak_forecast/kalman_params.json, extreme family)."""
const KF_PARAMS = Dict(0.1 => (q_l = 3.1623e-13, q_s = 3.1623e-10, r = 0.01),
                       0.3 => (q_l = 3.1623e-13, q_s = 1.0e-9, r = 0.09))
const KF_BURN = 60.0

# ---------------------------------------------------------------- N1 / N2 (condensed QP)
"""Linearisation of the ROM step at (x, d, u)."""
function lin(p, x, d, u)
    (f = R4.rhs(x, d, u, p), A = ForwardDiff.jacobian(z -> R4.rhs(z, d, u, p), x),
     B = ForwardDiff.derivative(v -> R4.rhs(x, d, v, p), u),
     E = ForwardDiff.derivative(v -> R4.rhs(x, v, u, p), d),
     D = ForwardDiff.derivative(n -> R4.polyval(p.poly, n), d + u + 5550), y = R4.output(x, d, u, p))
end

"""Condensed affine prediction ŷ = y0 + G ṁ of the backward-Euler ROM linearised at the
points `pts` = [(x̄_k, d̄_k, ū_k)], k = 1…K, with x̄_0 = `xbar0` (as f3_certificates_condensed.jl)."""
function condense(c, L, xbar0, pts)
    p = c.p; S = c.S; Sd = Diagonal(S); NX = R4.NX; Tw = R4.IDX.Tw
    xbar = vcat([xbar0], [pt[1] for pt in pts])
    Z0 = Sd \ (c.x_init - xbar0); Zu = zeros(NX, K); y0 = zeros(K); G = zeros(K, K)
    for k in 1:K
        x̄, d̄, ū = pts[k]; l = lin(p, x̄, d̄, ū)
        Mf = lu(Sd \ (I - h * l.A) * Sd)
        r = Sd \ (xbar[k] - xbar[k+1] + h * l.f); b = Sd \ (h * l.B); e = Sd \ (h * l.E) * (L[k] - d̄)
        Z0 = Mf \ (Z0 .+ r .- b .* ū .+ e)
        Zu = Mf \ Zu; Zu[:, k] .+= Mf \ b
        y0[k] = l.y + S[Tw] * Z0[Tw] + l.D * ((L[k] - d̄) - ū)
        G[k, :] = S[Tw] .* Zu[Tw, :]; G[k, k] += l.D
    end
    (y0 = y0, G = G)
end

function solve_qp(c, C, act)
    md = Model(Clarabel.Optimizer); set_silent(md)
    @variable(md, act.m_lo <= u[1:K] <= act.m_hi)
    @expression(md, y[k = 1:K], C.y0[k] + sum(C.G[k, j] * u[j] for j in 1:k))
    r = @variable(md, [1:K], lower_bound = 0)
    @constraint(md, [k = 1:K], y[k] - F4Lib.U_y <= r[k])
    obj = @expression(md, sum((y[k] - c.y_ref)^2 for k in 1:K) + F4Lib.lam_pct * sum(r .^ 2))
    Rm = get(c, :Rmove, 0.0)                      # F8: move suppression R Σ (u_k − u_{k−1})², u_0 = applied
    if Rm > 0
        uv = [c.m_prev; u]
        obj = @expression(md, obj + Rm * sum((uv[k+1] - uv[k])^2 for k in 1:K))
    end
    if get(c, :dsm, nothing) !== nothing         # F5b.1b: λ_s Σ [DSM_k − τ]₊², DSM on the bias-corrected ŷ
        z = c.dsm; den = z.U - (c.y_ref + c.d)    # U − y_ref (c.y_ref is shifted by −d)
        q = @variable(md, [1:K], lower_bound = 0)
        @constraint(md, [k = 1:K], (y[k] - c.y_ref) / den - z.thr <= q[k])
        obj = @expression(md, obj + z.lam * sum(q .^ 2))
    end
    if get(c, :ubound, nothing) !== nothing       # F5b.1c: ŷ_k + d ≤ U_α (soft, exact ℓ1 penalty:
        v = @variable(md, [1:K], lower_bound = 0)  # = hard bound when reachable, least total
        @constraint(md, [k = 1:K], y[k] + c.d - c.ubound <= v[k])   # violation when not; a
        obj = @expression(md, obj + c.lam_r * sum(v))                # quadratic slack with an
    end                                                              # unreachable U_α broke Clarabel
    if get(c, :exp_dT, nothing) !== nothing       # F5b.1c: ρ E_D[(ŷ_k + d + ΔT(D) − U)₊] (SAA on
        nq = length(c.exp_dT)                     # stratified quantiles of D; convex, not a bound)
        w = @variable(md, [1:K, 1:nq], lower_bound = 0)
        @constraint(md, [k = 1:K, i = 1:nq], y[k] + c.d + c.exp_dT[i] - c.U_lim <= w[k, i])
        obj = @expression(md, obj + c.rho / nq * sum(w))
    end
    if get(c, :scen, nothing) !== nothing          # F5b.1d: scenarios (j, η, D) inside the horizon
        sc = c.scen                                 # excess_s = a_s·u + b_s (affine in u, η multiplies u)
        A = zeros(0, K); bsc = Float64[]; wsc = Float64[]
        for j in sc.J, (ie, η) in enumerate(sc.etas), D in sc.Dq
            ks = min(K, j + ceil(Int, D / h))       # end of the failure, or the horizon end
            Drem = max(D - (ks - j + 1) * h, 0.0)  # failure time left after step ks
            e = [i < j ? 1.0 : η for i in 1:K]; e[ks+1:end] .= 0
            A = vcat(A, (C.G[ks, :] .* e)')
            push!(bsc, C.y0[ks] + c.d + (Drem > 0 ? dT_at(sc.curves[ie], Drem) : 0.0) - 1220.0)
            push!(wsc, sc.peta[ie] / (length(sc.J) * length(sc.Dq)))
        end
        nsc = length(bsc)
        z = @variable(md, [1:nsc], lower_bound = 0)
        if sc.kind == :cvar                         # CVaR_β(excess) ≤ 0, soft (ℓ1)
            tv = @variable(md); vv = @variable(md, lower_bound = 0)
            @constraint(md, [i = 1:nsc], z[i] >= sum(A[i, j] * u[j] for j in 1:K) + bsc[i] - tv)
            @constraint(md, tv + sum(wsc[i] * z[i] for i in 1:nsc) / (1 - sc.beta) <= vv)
            obj = @expression(md, obj + sc.lam * vv)
        else                                        # ρ E[(excess)₊]
            @constraint(md, [i = 1:nsc], z[i] >= sum(A[i, j] * u[j] for j in 1:K) + bsc[i])
            obj = @expression(md, obj + sc.rho * sum(wsc[i] * z[i] for i in 1:nsc))
        end
    end
    if act.ramp == :soft
        s = @variable(md, [1:K-1], lower_bound = 0)
        @constraint(md, [k = 1:K-1], u[k+1] - u[k] <= F4Lib.dm_lim + s[k])
        @constraint(md, [k = 1:K-1], u[k] - u[k+1] <= F4Lib.dm_lim + s[k])
        obj = @expression(md, obj + F4Lib.lam_m * sum(s .^ 2))
    else
        dm = get(act, :dm, F4Lib.dm_lim)            # hard ramp [kg/s per step]
        mm = [c.m_prev; u]
        @constraint(md, [k = 1:K], mm[k+1] - mm[k] <= dm)
        @constraint(md, [k = 1:K], mm[k] - mm[k+1] <= dm)
    end
    @objective(md, Min, obj)
    t = @elapsed optimize!(md)
    (m = clamp.(value.(u), act.m_lo, act.m_hi), status = termination_status(md), time = t)
end

"""ROM ode1be trajectory under flows m (rows 1…K+1: x at τ, t_q, …)."""
rom_traj(c, m, L) = R4.simulate(c.x_init, vcat(c.leak0, L), vcat(0.0, m), h, c.p; scheme = :ode1be, input_at = :next)

function decide(ctrl, c, L, act, prev)
    mid = (act.m_lo + min(act.m_hi, 50.0)) / 2
    guess = prev === nothing ? nothing : clamp.(vcat(prev[2:end], prev[end]), act.m_lo, act.m_hi)
    if ctrl == :N3
        t = @elapsed s = solve_det(c, L; m_lo = act.m_lo, m_hi = act.m_hi, ramp = act.ramp, guess = guess)
        return (m = s.m, status = string(s.status), time = t)
    elseif ctrl == :N1
        t = @elapsed begin
            C = condense(c, L, c.x_init, fill((c.x_init, c.leak0, mid), K))
            s = solve_qp(c, C, act)
        end
        return (m = s.m, status = string(s.status), time = t)
    elseif ctrl == :N2
        t = @elapsed begin
            m_inc = guess === nothing ? fill(mid, K) : guess
            X = rom_traj(c, m_inc, L)
            C = condense(c, L, X[1, :], [(X[k+1, :], L[k], m_inc[k]) for k in 1:K])
            s = solve_qp(c, C, act)
        end
        return (m = s.m, status = string(s.status), time = t)
    end
    error("unknown controller $ctrl")
end

# ---------------------------------------------------------------- risk bound (F5b.1c)
"""
    heatup_curve(c, yhat_now, tmax)

Heat-up of the PCT if the ACS failed now and stayed off (or, with `acs` > 0, delivered that
constant flow: partial failure, F5b.1d): the ROM (ACS = 0, leak held at its
current value) from the query state, ΔT(t) = max_{s ≤ t} ŷ_noACS(s) − ŷ_now on t = 0:Δt:tmax
(non-decreasing; it includes the instantaneous p₄ jump when the flow drops). If the ROM breaks
down (inventory exhausted) the last rate is extrapolated. Shared by every risk formulation of
F5b.1c (plan §8.14, §8.17): under the additive model PCT_peak = PCT(t_f) + ΔT(D), the duration
D enters only through ΔT.
"""
function heatup_curve(c, yhat_now, tmax; acs = 0.0)
    n = ceil(Int, tmax / h)
    X = nothing; nok = n
    while X === nothing && nok >= 4
        X = try
            R4.simulate(c.x_init, fill(c.leak0, nok + 1), fill(acs, nok + 1), h, c.p; scheme = :ode1be, input_at = :next)
        catch err
            err isa DomainError || rethrow(); nok ÷= 2; nothing
        end
    end
    y = [R4.output(X[i, :], c.leak0, acs, c.p) for i in axes(X, 1)]
    dT = accumulate(max, y) .- yhat_now
    if nok < n
        r = max((y[end] - y[end-1]) / h, 0.0)
        dT = vcat(dT, [dT[end] + r * h * k for k in 1:(n - nok)])
    end
    dT
end

"""ΔT at duration D from a heat-up curve (linear interpolation; the last slope beyond its end)."""
function dT_at(curve, D)
    x = D / h + 1; n = length(curve)
    x <= 1 && return curve[1]
    if x >= n
        return curve[n] + max(curve[n] - curve[n-1], 0.0) * (x - n)
    end
    k = floor(Int, x); w = x - k
    (1 - w) * curve[k] + w * curve[k+1]
end

"""Failure-length distribution of the DET-accelerated family: a mixture of exponentials
(weights `w`, repair rates `mu`). Survival, quantile (bisection) and sampling."""
failure_surv(t, dist) = sum(dist.w .* exp.(-dist.mu .* t))
function failure_quantile(p, dist)
    lo, hi = 0.0, 1e6
    for _ in 1:200; m = (lo + hi) / 2; 1 - failure_surv(m, dist) < p ? (lo = m) : (hi = m); end
    (lo + hi) / 2
end
function failure_sample(rng, dist)
    u = rand(rng); i = findfirst(>=(u), cumsum(dist.w)); i === nothing && (i = length(dist.w))
    randexp(rng) / dist.mu[i]
end

"""
    risk_bound(c, yhat_now, tau; U_lim = 1220.0)

Risk-adjusted upper bound of the PCT (plan §8.14): U_α = U_lim − ΔT(tau) with ΔT from
`heatup_curve`. With tau = q₁₋α(D), P(PCT_peak > U_lim | failure starting now) ≤ α holds if the
bias-corrected PCT stays below U_α. Returns (U = U_α, dT = ΔT).
"""
function risk_bound(c, yhat_now, tau; U_lim = 1220.0)
    dT = dT_at(heatup_curve(c, yhat_now, tau), tau)
    (U = U_lim - dT, dT = dT)
end

# ---------------------------------------------------------------- the loop
"""
    run_loop(B, ctrl; n = 61, act = B.act_arch, bias = false, eta = 1.0, fail = nothing,
             fbias = 0.0)

Closed loop from t_start over `n` queries. Returns a NamedTuple of per-query logs.
Plant perturbations the MPC does not know about (F5a.2):
  `eta`   ACS efficiency: the plant receives η·ṁ for a commanded ṁ (degraded ACS);
  `fail`  (t_a, t_b), or a vector of such windows (F5b.1c, from the DET component Markov
          chains): the ACS delivers nothing for t_a < t ≤ t_b (forced failure); a window
          (t_a, t_b, η) is a PARTIAL failure (F5b.1d): the ACS delivers η·ṁ (η = 0: total);
  `fbias` constant error of the leak forecast [kg/s] (ŷ uses ṁ_leak + fbias).
  `audit(c, L, act, prev, dec)` optional callback after each MPC decision (same query data).
`ctrl = :PI` leaves the plant to the DET's own PI (active_mpc = 0); `fail` applies to it
(the pump is forced off), `eta` only with `eta_pi = true` (F5b.1a; false reproduces F5a);
`pi_setpoint` replaces the PI setpoint (DET: 1144 K) after the first query. If the plant model breaks down (coolant inventory
exhausted, the HRM algebra fails) the loop stops and `log.breakdown[1]` is set.
When the PCT exceeds the upper edge of the DET monitor band (1220 K) the loop stops and
`log.limit[1]` is set: in Simulink the limintm block stops the simulation there (end state of
the event tree), so the replica's continuation (protection off) has no meaning.
By default the ROM warm-up uses the COMMANDED currents (what the twin knows), the plant the
real ones, and the MPC gets the true future leak (+ `fbias`).
What the twin knows (F5b.0, plan §8.12 step 2):
  `measured_acs = true`  the warm-up uses the flow the plant actually received;
  `leak_noise` σ [kg/s]  the twin sees y = ṁ_leak + σε (ε iid N(0,1), `noise_seed`), also in
                         the warm-up; the plant is unchanged;
  `forecast`             :oracle (true future leak + `fbias`), :persist (last measurement),
                         :kalman (local-linear-trend filter with jump intervention, `kf`
                         parameters, default `KF_PARAMS[leak_noise]`), :kalman_noint (no
                         intervention). The filter starts `KF_BURN` s before the first query.
`log.ferr` is the maximum forecast error over the horizon at each query.
`dsm = (lam, thr, U)` adds to the N1/N2 QP the margin penalty λ_s Σ [DSM_k − τ]₊² with
DSM_k = (ŷ_k − y_ref)/(U − y_ref) on the bias-corrected prediction (F5b.1b, plan §8.12 step 3.i;
the thesis form of `matlab/fitness_cost_dsm.m`).
`risk = (tau, lam)` adds the risk-adjusted soft bound ŷ_k + d ≤ U_α (`risk_bound`, exact ℓ1
penalty `lam` [1/K] on the slack). Other formulations (plan §8.17), all on `heatup_curve`:
  `(mode = :cvar, beta, Dq, lam)`     U = 1220 − CVaR_β[ΔT(D)] over the quantile grid `Dq` of D
                                      (the exact reduction of the Rockafellar–Uryasev SAA);
  `(mode = :wilks, dist, n = 59, lam)` U = 1220 − ΔT(max of n fresh draws of D) per query
                                      (online 95/95 PSM);
  `(mode = :expected, rho, Dq)`       cost ρ/N Σ_k Σ_i (ŷ_k + d + ΔT(D_i) − 1220)₊ (expected
                                      excess given a failure now; convex, no bound).
  `(mode = :scenario, fn, rho)`       F9: the expected excess with ONE scenario, margin
                                      δ̂ = fn(curve, c, yhat_now) (Homem-de-Mello et al., 2026).
  `Dcap` (400 s) caps D: a failure longer than that crosses 1220 K whatever the controller.
Partial failures (F5b.1d, plan §8.18; severities `etas` with probabilities `peta`, the ACS then
delivers η·`mmax` at most):
  `(mode = :var_eta, alpha, etas, peta, dist, mmax, lam)`  U = 1220 − m with
                                      Σ_η p_η P(ΔT_η(D) > m) = α (η-aware margin; still a bound);
  `(mode = :scen_cvar | :scen_exp, J, etas, peta, Dq, mmax, beta | rho, lam)`  scenarios
                                      (failure start j ∈ J inside the horizon, η, D): the
                                      delivered flow is η·u from j on, so the excess is affine in
                                      u with a scenario-dependent slope (not a bound); CVaR_β of
                                      the excess (Rockafellar–Uryasev) or ρ·E[excess₊]. for a failure of length `tau` = q₁₋α(D) (F5b.1c, plan §8.14).
`log.ubound` keeps U_α per query (NaN without `risk`).
`ctrl = :PID` is the thesis PID (`matlab/pid_controller.m`, used in
`DET_scheduler_v1_perc_online.m` and `run_phase2c_full_loop.m`): on e = y_ref − PCT(τ),
u_k = u_{k−1} + K_p e + I_k, I_k = I_{k−1} + K_i e Δt (K_p = −0.5, K_i = −0.4 A/K), clamped
to [56.68, 100] A and rate-limited to 60 A/s from u_{k−1} (u_0 = 0, as in the thesis); the
current goes to the pump map and the plant (capacity, η and failures apply as to the MPC).
`ctrl = :PI` is the DET's own PI (setpoint 1144 K), a different controller.
`Rmove = R` (F8) adds the move suppression R Σ_k (ṁ_k − ṁ_{k−1})² to N1/N2/N3, with ṁ_0 the
flow applied at the previous query (R = 0: the original formulation).
`keep_states = true` stores the plant state and DET logic after every query in `log.states`
(for `time_to_limit`).
"""
function run_loop(B::Bench, ctrl::Symbol; n = 61, act = B.act_arch, bias = false, eta = 1.0,
                  fail = nothing, fbias = 0.0, audit = nothing, measured_acs = false,
                  leak_noise = 0.0, noise_seed = 1, forecast = :oracle, kf = nothing, kf_thr = 5.0, eta_pi = false,
                  dsm = nothing, pi_setpoint = nothing, risk = nothing, Rmove = 0.0,
                  pid = (Kp = -0.5, Ki = -0.4, I_min = 56.68, I_max = 100.0, dI_max = 60.0),
                  keep_states = false)
    pH = B.pH; hist = copy(B.hist)
    pPI = with_params(pH; acs_eta = eta_pi ? eta : 1.0,
                      pct_setpoint = pi_setpoint === nothing ? pH.pct_setpoint : pi_setpoint)
    last = findfirst(==(B.t_start - h), B.cs[:, 1])
    T = B.cs[2:last, 1]; U = B.cs[2:last, 2:3]          # commanded inport series (grows)
    Up = copy(U)                                        # what the plant actually receives
    m_applied = 0.0; prev = nothing
    pid_u = 0.0; pid_I = 0.0                            # thesis PID state
    wrng = Random.Xoshiro(noise_seed + 7)               # Wilks draws (risk mode :wilks)
    log = (t_q = Float64[], pct = Float64[], m = Float64[], m_plant = Float64[], u = Float64[],
           yhat_now = Float64[], d = Float64[], time = Float64[], status = String[],
           m_plan = Vector{Float64}[], ferr = Float64[], ubound = Float64[], yhat_max = Float64[], states = Any[], breakdown = [false], limit = [false])
    # leak measurements and forecaster
    tb = B.t_start - KF_BURN - DELAY - h
    nz = randn(Random.Xoshiro(noise_seed), round(Int, (B.t_start + n * h - tb) / h) + 2)
    ymeas(t) = SMDFRHRM.leak_flow(t, pH) + leak_noise * nz[round(Int, (t - tb) / h) + 1]
    forecast in (:oracle, :persist, :kalman, :kalman_noint) || error("forecast $forecast")
    filt = nothing; t_kf = B.t_start - KF_BURN
    if forecast in (:kalman, :kalman_noint)
        par = kf === nothing ? KF_PARAMS[leak_noise] : kf
        filt = LeakKF(ymeas(t_kf); par..., thr = forecast == :kalman ? kf_thr : Inf)
    end
    acs_meas = measured_acs ? (t -> hist[t][2].m_last) : nothing
    fails = fail === nothing ? Tuple{Float64,Float64,Float64}[] : fail isa Tuple ? [fail] : fail
    fails = [length(f) == 2 ? (f[1], f[2], 0.0) : (f[1], f[2], f[3]) for f in fails]   # (t_a, t_b, η)
    leak_meas = leak_noise > 0 ? ymeas : nothing
    for k in 1:n
        t_q = B.t_start + (k - 1) * h; τ = t_q - h
        x, Lg, pct_now = hist[τ]
        xH10, L10, _ = hist[t_q-DELAY]
        eta_f = minimum((e for (a, b, e) in fails if a < t_q <= b); init = 1.0)   # 1: no failure
        failed = eta_f <= 0.0
        ctrl == :PI && 0 < eta_f < 1 && error("partial failures are not implemented for the DET PI")
        if ctrl == :PI                       # the DET's own PI drives the pump (no twin needed)
            T = vcat(T, t_q); U = vcat(U, [0.0 0.0]); Up = vcat(Up, failed ? [56.0 1.0] : [0.0 0.0])
            s = try
                SMDFRHRM.simulate(pPI, ExtInput(T, Up), t_q; t0 = τ, x0 = x, logic = Lg)
            catch err
                err isa DomainError || rethrow(); log.breakdown[1] = true; break
            end
            hist[t_q] = (s.x, s.logic, s.pct[end]); keep_states && push!(log.states, (s.x, s.logic))
            push!(log.t_q, t_q); push!(log.pct, s.pct[end]); push!(log.m, s.acs[end])
            push!(log.m_plant, s.acs[end]); push!(log.u, s.current[end]); push!(log.yhat_now, NaN)
            push!(log.d, NaN); push!(log.time, 0.0); push!(log.status, "-"); push!(log.m_plan, [s.acs[end]])
            push!(log.ferr, NaN); push!(log.ubound, NaN); push!(log.yhat_max, NaN)
            s.pct[end] > pH.monitor_band[2] && (log.limit[1] = true; break)
            continue
        end
        ext = ExtInput(T, U)
        w = try
            warm_rom(pH, B.pR, xH10, L10.xI, ext, t_q; delay = DELAY, h, leak_meas, acs_meas)
        catch err
            err isa DomainError || rethrow(); log.breakdown[1] = true; break
        end
        Ltrue = [SMDFRHRM.leak_flow(t_q + (j - 1) * h, pH) for j in 1:K]   # leak at τ + jΔt
        if forecast == :oracle
            L = Ltrue .+ fbias
        elseif forecast == :persist
            L = fill(ymeas(τ), K)
        else
            while t_kf < τ - 1e-9
                t_kf += h; kf_update!(filt, ymeas(t_kf))
            end
            L = kf_forecast(filt, K).m
        end
        yhat_now = R4.output(w.x_init, w.leak_warm[end], w.acs_warm[end], B.pR4)
        d = pct_now - yhat_now
        c = (; p = B.pR4, h, K, S = max.(abs.(w.x_init), 1.0), x_init = w.x_init,
             leak0 = w.leak_warm[end], y_ref = bias ? B.y_ref - d : B.y_ref, m_prev = m_applied,
             d = bias ? d : 0.0, dsm, Rmove = Float64(Rmove))
        ub = NaN
        if risk !== nothing
            mode = get(risk, :mode, :var)
            if mode == :var                                   # chance constraint, quantile q₁₋α
                ub = risk_bound(c, yhat_now, risk.tau).U
            else
                Dcap = get(risk, :Dcap, 400.0)                # longer failures cross anyway
                curve = heatup_curve(c, yhat_now, Dcap)
                if mode == :cvar                              # CVaR_β of ΔT(D) (stratified SAA)
                    dTs = sort([dT_at(curve, min(D, Dcap)) for D in risk.Dq])
                    nt = max(1, ceil(Int, (1 - risk.beta) * length(dTs)))
                    ub = 1220.0 - sum(dTs[end-nt+1:end]) / nt
                elseif mode == :wilks                         # 95/95 Wilks: max of 59 fresh draws
                    Dm = maximum(failure_sample(wrng, risk.dist) for _ in 1:get(risk, :n, 59))
                    ub = 1220.0 - dT_at(curve, min(Dm, Dcap))
                elseif mode == :var_eta                       # quantile of the joint (η, D) margin
                    curves = [heatup_curve(c, yhat_now, Dcap; acs = e * risk.mmax) for e in risk.etas]
                    pexc(m) = sum(pe * (let i = findfirst(>=(m), cv); i === nothing ? 0.0 :
                                  failure_surv((i - 1) * h, risk.dist) end) for (pe, cv) in zip(risk.peta, curves))
                    lo, hi = 0.0, maximum(maximum.(curves))
                    for _ in 1:60; m = (lo + hi) / 2; pexc(m) > risk.alpha ? (lo = m) : (hi = m); end
                    ub = 1220.0 - hi
                elseif mode in (:scen_cvar, :scen_exp)        # in-horizon partial-failure scenarios
                    curves = [heatup_curve(c, yhat_now, Dcap; acs = e * risk.mmax) for e in risk.etas]
                    c = (; c..., scen = (J = risk.J, etas = risk.etas, peta = risk.peta,
                         Dq = [min(D, Dcap) for D in risk.Dq], curves,
                         kind = mode == :scen_cvar ? :cvar : :exp, beta = get(risk, :beta, 0.95),
                         rho = get(risk, :rho, 0.0), lam = get(risk, :lam, 1e3)))
                elseif mode == :expected                      # ρ E[(peak − 1220)₊], no bound
                    c = (; c..., exp_dT = [dT_at(curve, min(D, Dcap)) for D in risk.Dq], rho = risk.rho,
                         U_lim = 1220.0)
                elseif mode == :scenario                      # F9: one-scenario MPC (Homem-de-Mello
                    c = (; c..., exp_dT = [risk.fn(curve, c, yhat_now)], rho = risk.rho,   # et al.):
                         U_lim = 1220.0)                                  # margin δ̂ = fn(context)
                else
                    error("risk mode $mode")
                end
            end
            mode in (:expected, :scenario, :scen_cvar, :scen_exp) || (c = (; c..., ubound = ub, lam_r = risk.lam))
        end
        if ctrl in (:N1, :N2, :N3)
            dec = decide(ctrl, c, L, act, prev)
            audit === nothing || audit(c, L, act, prev, dec)
            if !(dec.status in OK_STATUS) && prev !== nothing   # fallback: shifted previous plan
                dec = (m = clamp.(vcat(prev[2:end], prev[end]), act.m_lo, act.m_hi), status = dec.status, time = dec.time)
            end
            m1 = dec.m[1]; prev = dec.m
            # F8 check: largest PCT the ROM predicts for the applied plan (no bias correction, as
            # in the λ_pct[ŷ − U_y]₊ term of the cost)
            Xp = rom_traj(c, dec.m, L)
            yhmax = maximum(R4.output(Xp[j+1, :], L[j], dec.m[j], c.p) for j in 1:K)
        else
            if ctrl == :PID
                e = B.y_ref - pct_now
                pid_I += pid.Ki * e * h
                uu = clamp(pid_u + pid.Kp * e + pid_I, pid.I_min, pid.I_max)
                pid_u = pid_u + sign(uu - pid_u) * min(abs(uu - pid_u), pid.dI_max * h)
                m1 = current_to_massflow(pid_u)
            else
                m1 = ctrl == :zero ? 0.0 : ctrl == :max ? act.m_hi :
                     ctrl == :LJ ? current_to_massflow(B.cs[findfirst(==(t_q), B.cs[:, 1]), 2]) : error("ctrl")
            end
            dec = (m = fill(m1, K), status = "-", time = 0.0)
            yhmax = NaN
        end
        u1 = m1 > 0 ? F4Lib.massflow_to_current(m1) : 56.0                  # 56 A: pump off
        mp = failed ? 0.0 : eta_f * eta * m1                                 # real flow
        up = mp > 0 ? F4Lib.massflow_to_current(mp) : 56.0
        T = vcat(T, t_q); U = vcat(U, [u1 1.0]); Up = vcat(Up, [up 1.0])
        s = try
            SMDFRHRM.simulate(pH, ExtInput(T, Up), t_q; t0 = τ, x0 = x, logic = Lg)
        catch err
            err isa DomainError || rethrow(); log.breakdown[1] = true; break
        end
        # the DET logic passes the MPC flow only while the protection is active (e.g. not yet
        # in the first steps after the activation, scheduler delay), and a finite-capacity
        # plant (pH.acs_max, pH.acs_ramp) limits it (e.g. the ramp-up after a failure):
        # log what the plant got
        if s.act[end] > 0 && isinf(pH.acs_max) && isinf(pH.acs_ramp)
            @assert abs(s.acs[end] - mp) <= 1e-6 * max(1.0, mp) "plant did not receive the intended flow"
        end
        mp = s.acs[end]
        hist[t_q] = (s.x, s.logic, s.pct[end]); keep_states && push!(log.states, (s.x, s.logic))
        m_applied = m1
        push!(log.t_q, t_q); push!(log.pct, s.pct[end]); push!(log.m, m1); push!(log.m_plant, mp)
        push!(log.u, u1)
        push!(log.yhat_now, yhat_now); push!(log.d, d); push!(log.time, dec.time)
        push!(log.status, dec.status); push!(log.m_plan, dec.m); push!(log.ferr, maximum(abs.(L .- Ltrue)))
        push!(log.ubound, ub); push!(log.yhat_max, yhmax)
        s.pct[end] > pH.monitor_band[2] && (log.limit[1] = true; break)
    end
    log
end

const OK_STATUS = ("-", "LOCALLY_SOLVED", "OPTIMAL", "ALMOST_OPTIMAL", "ALMOST_LOCALLY_SOLVED")

"""Closed-loop metrics of a `run_loop` log (PCT logged at each t_q, after the applied step)."""
function metrics(lg, y_ref)
    e = lg.pct .- y_ref
    (; n = length(e), breakdown = lg.breakdown[1], limit = lg.limit[1],
     t_limit = lg.limit[1] ? lg.t_q[end] - lg.t_q[1] + h : NaN, pct_max = maximum(lg.pct), pct_final = lg.pct[end],
     t_above = h * count(>(0), e), ise = sum(e .^ 2), iae = h * sum(abs.(e)),
     err_last10 = sum(abs.(e[max(1, end - 19):end])) / min(20, length(e)),
     mass_cmd = h * sum(lg.m), mass_plant = h * sum(lg.m_plant),
     ramp_max = maximum(abs.(diff(vcat(0.0, lg.m)))),
     t_p50 = sort(lg.time)[cld(length(lg.time), 2)], t_p99 = sort(lg.time)[ceil(Int, 0.99 * length(lg.time))],
     t_max = maximum(lg.time), fails = count(s -> !(s in OK_STATUS), lg.status))
end

"""
    time_to_limit(pH, x, logic, t; U = 1220.0, t_cap = 600.0, chunk = 50.0)

Time τ for the PCT of the HRM (the plant) to exceed `U` if the ACS is off from `t` on (pump
failed: zero flow, whatever the PI or the MPC ask), starting from the plant state `x` and DET
logic `logic` at `t`. Simulated in chunks of `chunk` s up to `t_cap` (returns `t_cap` if not
reached) or until the inventory is exhausted (returns the time reached). Used by the
Rao–Blackwell risk estimator of F5b.1c (plan §8.16). `flow` > 0: partial failure, the ACS
delivers that constant flow (F5b.1d; the plant capacity and ramp still apply).
"""
function time_to_limit(pH, x, logic, t; U = 1220.0, t_cap = 600.0, chunk = 50.0, flow = 0.0)
    cur = flow > 0 ? F4Lib.massflow_to_current(flow) : 56.0          # 56 A: zero flow
    ext = ExtInput([t, t + t_cap + chunk], [cur 1.0; cur 1.0])
    t0 = t; xx = x; L = logic
    while t0 < t + t_cap
        s = try
            SMDFRHRM.simulate(pH, ext, t0 + chunk; t0 = t0, x0 = xx, logic = L)
        catch err
            err isa DomainError || rethrow(); return t0 - t
        end
        i = findfirst(>(U), s.pct)
        i === nothing || return s.t[i] - t
        # partial failure: if the PCT is already falling with a constant flow it will not cross
        # (the leak of Park's law only shrinks with time, so the net cooling only grows)
        flow > 0 && s.pct[end] < s.pct[1] - 0.5 && return t_cap
        t0 += chunk; xx = s.x; L = s.logic
    end
    t_cap
end

"""
    thesis_metrics(lg, y_ref; tol = 0.30, U = 1244.0, t0 = lg.t_q[1])

Performance indices of the thesis (`matlab/calculate_performance_metrics.m`,
`python/compile_thesis_kpis.py`), on the per-query PCT and the commanded flow [kg/s]
(the thesis used the current [A]; flows are comparable across actuators). e = y_ref − PCT.
  min_safety_margin  min(U − PCT)                   reach_time    first t with |e| ≤ tol
  min_overshoot      min(PCT − y_ref)               settling_time t after which |e| stays ≤ tol
  max_abs_dev_pt     max|PCT − y_ref|, without the first two samples (activation transient)
  steady_state_error mean e over the last 20 %      IAE, cumulative_error  trapz over t
  peak_PCT, final_PCT, total_variation Σ|Δṁ| of the commanded flow (tv_plant: of the flow the
  plant received, which differs when a PID asks beyond capacity), max_control_rate, coolant mass,
  runtime_median. Times are measured from `t0` (the first query); Inf if never reached.
"""
function thesis_metrics(lg, y_ref; tol = 0.30, U = 1244.0, t0 = lg.t_q[1])
    t = lg.t_q .- t0; pct = lg.pct; e = y_ref .- pct; m = lg.m
    trapz(x, y) = sum((x[2:end] .- x[1:end-1]) .* (y[2:end] .+ y[1:end-1]) ./ 2)
    inb = abs.(e) .<= tol
    reach = findfirst(inb); out = findlast(.!inb)
    settle = reach === nothing ? Inf : out === nothing ? t[reach] : out < length(t) ? t[out+1] : Inf
    du = diff(vcat(0.0, m))
    (; min_safety_margin = minimum(U .- pct), margin_1220 = 1220 - maximum(pct),
     min_overshoot = minimum(pct .- y_ref), max_abs_dev_pt = maximum(abs.(pct[min(3, end):end] .- y_ref)),
     reach_time = reach === nothing ? Inf : t[reach], settling_time = settle,
     steady_state_error = sum(e[round(Int, 0.8 * length(e)):end]) / length(e[round(Int, 0.8 * length(e)):end]),
     IAE = trapz(t, abs.(e)), cumulative_error = trapz(t, e), peak_PCT = maximum(pct), final_PCT = pct[end],
     total_variation = sum(abs.(du)), max_control_rate = maximum(abs.(du)) / h,
     tv_plant = sum(abs.(diff(vcat(0.0, lg.m_plant)))),
     coolant_mass = h * sum(lg.m_plant), runtime_median = sort(lg.time)[cld(length(lg.time), 2)])
end

end # module
