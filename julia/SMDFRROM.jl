"""
    SMDFRROM

Algebraic replica of the 0-D reduced-order model (ROM) of the SMDFR used as the MPC
prediction model (`models/rom_MCDET_v2_2022b_control_CMEX.slx`, subsystem
`SMDFR Units/SMDFR 1`). Every term is traced to a MATLAB Function block of the
`.slx`; see `docs/monolithic_formulation_plan.md` §3.

State (17): x = [N; C₁…C₆; Cᵉ₁…Cᵉ₆; T_f; T_w; T_c; ΔM]
  N      neutron population                       (Int_N)
  Cᵢ     in-core precursors                       (Int-C_i)
  Cᵉᵢ    out-of-core precursors                   (Int-Ce_i)
  T_f    fuel temperature            [K]          (Int-T_f)
  T_w    wall temperature            [K]          (Int-T_w)
  T_c    coolant temperature         [K]          (Int-T_c)
  ΔM     coolant inventory change    [kg]         (Int-mdot_c, IC 0 at t = 0)

Inputs per step: leak flow ṁ_leak [kg/s] (≤ 0 during the LOCA) and ACS flow ṁ_ACS ≥ 0.

Output (PCT estimate, `pred_tw_max_signal`):
  ŷ = T_w + p₄(ṁ_leak + ṁ_ACS + 5550) + 0.4
where p₄ is the degree-4 regression of ΔPCT = PCT_HRM − T_w,ROM on the net coolant
flow (`fit_poly`, constant term shifted by −2.5 K). The BiLSTM ΔPCT forecast that
reaches the ROM (`pred_delta_tw_max`) is NOT connected inside the ROM.

Simulink integrates with `ode1be` (backward Euler, Δt = 0.5 s, ONE Newton
iteration); `step_ode1be` reproduces that scheme.
"""
module SMDFRROM

using LinearAlgebra, ForwardDiff, MAT

export Params, load_params, rhs, output, current_to_massflow, massflow_to_current,
       step_ode1be, step_be,
       simulate, NX, IDX

const NX = 17
const IDX = (N = 1, C = 2:7, Ce = 8:13, Tf = 14, Tw = 15, Tc = 16, dM = 17)

"""Model parameters, as read from the ROM model workspace."""
Base.@kwdef struct Params
    Lambda::Float64; beta::Float64; beta_i::Vector{Float64}; lambda_i::Vector{Float64}
    lambda_c::Float64; lambda_e::Float64; N_0::Float64; P_0::Float64; Pi_0::Float64
    rho_0::Float64; rhoex::Float64; alfa_f::Float64; alfa_c::Float64
    Tbar_f0::Float64; Tbar_c0::Float64; Te_f0::Float64; Te_c0::Float64
    A_fw::Float64; A_wc::Float64; h_fw::Float64; h_wc::Float64
    M_f::Float64; M_w::Float64; M_c::Float64; c_pf::Float64; c_pw::Float64; c_pc::Float64
    mdot_f::Float64; mdot_c::Float64
    poly::Vector{Float64}          # descending powers, as MATLAB polyval
    Q0::Float64 = 8.2564e7 * 0.578 # `MATLAB Function15` (empirical inlet-temperature law)
end

"""Read parameters from the `par` struct written by `matlab/f1_export_rom_reference.m`."""
function load_params(par::AbstractDict)
    g(k) = par[k]
    s(k) = Float64(only(g(k)))
    v(k) = vec(Float64.(g(k)))
    Params(Lambda = s("Lambda"), beta = s("beta"), beta_i = v("beta_i"), lambda_i = v("lambda_i"),
        lambda_c = s("lambda_c"), lambda_e = s("lambda_e"), N_0 = s("N_0"), P_0 = s("P_0"),
        Pi_0 = s("Pi_0"), rho_0 = s("rho_0"), rhoex = s("rhoex_value"), alfa_f = s("alfa_f"),
        alfa_c = s("alfa_c"), Tbar_f0 = s("Tbar_f0"), Tbar_c0 = s("Tbar_c0"),
        Te_f0 = s("Te_f0"), Te_c0 = s("Te_c0"), A_fw = s("A_fw"), A_wc = s("A_wc"),
        h_fw = s("h_fw"), h_wc = s("h_wc"), M_f = s("M_f"), M_w = s("M_w"), M_c = s("M_c"),
        c_pf = s("c_pf"), c_pw = s("c_pw"), c_pc = s("c_pc"), mdot_f = s("mdot_f"),
        mdot_c = s("mdot_c"), poly = v("poly_vector"))
end

polyval(c, x) = foldl((acc, ci) -> acc * x + ci, c; init = zero(x))

"""Intermediate algebraic quantities (useful for checks against logged signals)."""
function algebraic(x, p::Params)
    Tf, Tc = x[IDX.Tf], x[IDX.Tc]
    M_act  = p.M_c + x[IDX.dM]                       # Sum: M_c + ∫ṁ_net
    mdot_c = p.mdot_c * M_act / p.M_c                # MATLAB Function18
    h_wc   = p.h_wc * (mdot_c / 5550)^0.8            # MATLAB Function19
    Q      = p.Q0 * mdot_c / (5550 / 2) - (5500 - mdot_c)
    Tin_c  = p.Te_c0 - Q / (mdot_c * p.c_pc)         # MATLAB Function15
    Te_c   = 2Tc - Tin_c
    Tin_f  = p.Te_f0                                 # MATLAB Function16
    Te_f   = 2Tf - Tin_f
    # MATLAB Function3: mean([Te,Tin]) is exactly the node temperature
    rho    = p.rho_0 + p.alfa_f * ((Te_f + Tin_f) / 2 - p.Tbar_f0) +
             p.alfa_c * ((Te_c + Tin_c) / 2 - p.Tbar_c0) + p.rhoex
    P      = p.P_0 / p.N_0 * x[IDX.N]                # MATLAB Function13
    (; M_act, mdot_c, h_wc, Tin_c, Te_c, Tin_f, Te_f, rho, P, Pi = P * p.Pi_0)
end

"""Time derivative ẋ = f(x, ṁ_leak, ṁ_ACS)."""
function rhs(x::AbstractVector{T}, leak, acs, p::Params) where {T}
    a = algebraic(x, p)
    N, C, Ce = x[IDX.N], x[IDX.C], x[IDX.Ce]
    Tf, Tw, Tc = x[IDX.Tf], x[IDX.Tw], x[IDX.Tc]
    dx = similar(x, promote_type(T, typeof(leak), typeof(acs)), NX)
    dx[IDX.N]  = (a.rho - p.beta) * N / p.Lambda + sum(p.lambda_i .* C)
    dx[IDX.C]  = p.beta_i ./ p.Lambda .* N .- (p.lambda_i .+ p.lambda_c) .* C .+ p.lambda_e .* Ce
    dx[IDX.Ce] = p.lambda_c .* C .- (p.lambda_i .+ p.lambda_e) .* Ce
    dx[IDX.Tf] = (a.Pi + p.mdot_f * p.c_pf * (a.Tin_f - a.Te_f) - p.h_fw * p.A_fw * (Tf - Tw)) /
                 (p.M_f * p.c_pf)
    dx[IDX.Tw] = (p.h_fw * p.A_fw * (Tf - Tw) - a.h_wc * p.A_wc * (Tw - Tc)) / (p.M_w * p.c_pw)
    dx[IDX.Tc] = (a.mdot_c * p.c_pc * (a.Tin_c - a.Te_c) + a.h_wc * p.A_wc * (Tw - Tc)) /
                 (a.M_act * p.c_pc)
    dx[IDX.dM] = leak + acs
    dx
end

"""PCT estimate ŷ (`pred_tw_max_signal`)."""
output(x, leak, acs, p::Params) = x[IDX.Tw] + polyval(p.poly, leak + acs + 5550) + 0.4

"""ECM coefficients of the ACS pump: ΔP = C·I − B·Q − A·Q² (Nowak et al., 2022)."""
function pump_coefficients()
    rho = 10449.0; Hd = 0.515; Wd = 0.515; S = Wd * Hd; L = 1.0; d = 0.001
    Bm = 1.0; DeltaP = 15000.0; f = 0.02; q0_LM = 0.118; q0_w = 1.2; K = 1.0
    R_LM = q0_LM * Wd / (Hd * L); R_w = q0_w * Wd / (2d * L); R_b = q0_LM * K^2 * Hd / d
    R_ver = R_w * R_b / (R_w + R_b); R_tot = R_LM + R_ver
    C = Bm * R_ver * Hd^2 / (R_tot * d); B = Bm^2 * R_tot * Hd^2 / d
    A = f * rho * L * (Wd + Hd) / (4S^3)
    (; rho, DeltaP, A, B, C)
end

"""ACS pump map I [A] → ṁ [kg/s], as `matlab/CurrentToMassFlow.m`."""
function current_to_massflow(I)
    c = pump_coefficients()
    disc = c.B^2 - 4c.A * (c.DeltaP - c.C * I)
    Q = max((-c.B + sqrt(disc)) / (2c.A), 0.0)
    max(c.rho * Q, 0.0)
end

"""Inverse pump map ṁ [kg/s] → I [A] (exact on ṁ > 0)."""
function massflow_to_current(m)
    c = pump_coefficients()
    Q = m / c.rho
    (c.DeltaP + c.B * Q + c.A * Q^2) / c.C
end

"""One `ode1be` step: backward Euler with a single Newton iteration from x_n."""
function step_ode1be(x, leak, acs, h, p::Params)
    f(z) = rhs(z, leak, acs, p)
    J = ForwardDiff.jacobian(f, x)
    x .+ (I - h * J) \ (h .* f(x))
end

"""One fully converged backward Euler step (reference for the monolithic NLP)."""
function step_be(x, leak, acs, h, p::Params; tol = 1e-10, maxit = 50)
    z = copy(x)
    for _ in 1:maxit
        r = z .- x .- h .* rhs(z, leak, acs, p)
        J = I - h * ForwardDiff.jacobian(w -> rhs(w, leak, acs, p), z)
        dz = J \ r
        z .-= dz
        maximum(abs.(dz) ./ max.(abs.(z), 1.0)) < tol && break
    end
    z
end

"""
    simulate(x0, leak, acs, h, p; scheme = :ode1be, input_at = :next)

Integrate over the grid t₀, t₀+h, …; `leak` and `acs` are sampled on that grid
(length K+1). `input_at = :next` uses inputs at t_{n+1} inside step n (implicit
convention); `:current` uses t_n. Returns the (K+1)×17 state matrix.
"""
function simulate(x0, leak, acs, h, p::Params; scheme = :ode1be, input_at = :next)
    K = length(leak) - 1
    X = zeros(K + 1, NX); X[1, :] = x0
    step = scheme === :ode1be ? step_ode1be : step_be
    for n in 1:K
        j = input_at === :next ? n + 1 : n
        X[n+1, :] = step(X[n, :], leak[j], acs[j], h, p)
    end
    X
end

end # module
