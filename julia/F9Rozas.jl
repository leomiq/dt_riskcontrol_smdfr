"""
    F9Rozas

F9: the one-scenario result of Homem-de-Mello, Valencia, Lagos & Lagos (2026), *Forecasting
outside the Box* (Math. Program., doi:10.1007/s10107-026-02359-9; proposed by H. Rozas, plan V6),
checked on the closed-loop MPC queries of the SMDFR with the paper's own methodology.

Paper (2)–(5): min_{z∈Z} f(z) + E[Q(z, ξ)], Q(z, ξ) = min{qᵀy : W y = h − T z, y ≥ 0}, ξ = (h, T),
W and q fixed (FRFC). Theorem 1: under Assumption 1 (Z ≠ ∅, ri Dom f ≠ ∅), Assumption 2 (dual set
U = {u : Wᵀu ≤ q} non-empty and bounded) and either (i) T deterministic or (ii) T random with
Lemma 1 (n < m; or n = m with 0 ∈ int U; or n = m, 0 an extreme point of U and −TᵀU inside one
orthant w.p. 1), the scenario ξ* = (T*z*, T*) reproduces z*; unique if Assumption 3 holds.

Here z = the flow plan u ∈ ℝ²⁶ of the N2 QP (box [0, 40] kg/s, hard ramp 5 kg/s, R = 0.1),
f = the N2 cost (tracking + λ_pct hinge² + R Σ Δu²), and three second stages:
  P1  ρ Σ_k (ŷ_k + d + ΔT(D) − 1220)₊, total ACS failure starting now, duration D (DET mixture,
      μ×1e4): h random, T = −G fixed → case (i). The loop's `:expected` risk layer.
  P2  partial failures inside the horizon (start j, severity η, duration D): the plant gets η u,
      so T is random → case (ii). As implemented (`:scen_exp`, one row per scenario, m = 1 < n),
      and per step (m = K = n, Lemma 1 (iii) needs every column of G to keep its sign).
  P3  leak-forecast error δL ~ N(0, Σ_L) (Kalman, F5b.0), tracking in L1 (|·| as in the
      newsvendor): h random, T fixed → case (i). The case plan §7.1 said the theorem covers.

Experiments (per captured query; the driver is f9_rozas.jl):
  E0  the captured QP reproduces the loop decision;
  E1  assumptions: λ_min of the Hessian (Asm. 3), dual set (Asm. 2), sign pattern of G (Lemma 1);
  E2  Theorem 1 constructively: ξ* built from the SAA solution must return the SAA plan;
  E3  bilevel (38)–(39) over PHYSICAL one-scenario families (a duration, or (j, η, D)), solved
      derivative-free (grid + golden section / Nelder–Mead), against SAA; CPU time vs N (Table 3);
  E4  contextual case (§5–6): forecast Φ(θ, x) of the scenario from the query features, trained
      with Algorithm 1 + Nelder–Mead on realised costs, evaluated out of sample (phase B);
  E5  closed loop with the learned one-scenario MPC (phase C).
"""
module F9Rozas

include(joinpath(@__DIR__, "F6Tasks.jl"))
using .F6Tasks, JuMP, Clarabel, LinearAlgebra, Statistics, Random, Printf, Serialization
const T6 = F6Tasks; const CL = T6.CL; const F4 = T6.F4; const R4 = CL.R4

const K = 26; const h = 0.5; const ULIM = 1220.0; const DCAP = 400.0; const MUACC = 1e4
const RHOS = (100.0, 1000.0)          # weights of the expected excess (F6 used 10/100/1000)
const RM = 0.1                         # F8 move suppression (paper formulation)
const MMAX = 40.0
const DGRID = collect(-20.0:0.5:130.0) # one-scenario margin δ̂ = ΔT̂ [K] (P1), common grid
const NS_TRAIN = 8                     # realised durations per query (contextual training data)
const ETAS = T6.ETAS; const PETA = T6.PVECS.base; const JSET = [1, 6, 11, 16, 21, 26]
const OKS = (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)

dist() = T6.dist(MUACC)
qD(n) = [min(D, DCAP) for D in T6.qgrid(n, MUACC)]
dTs(curve, Ds) = [CL.dT_at(curve, D) for D in Ds]
maxabs(a, b) = maximum(abs.(a .- b))

# ---------------------------------------------------------------- features (context x)
const FNAMES = (:margin, :d, :m_prev, :leak, :dT10, :dT46, :dT150, :dT400)
"""Context of a query: margin to 1220 K of the bias-corrected PCT, bias, applied flow, leak and
the heat-up without ACS after 10/46/150/400 s (the same function runs inside the loop)."""
features(curve, c, yhat_now) = (margin = ULIM - (yhat_now + c.d), d = c.d, m_prev = c.m_prev,
    leak = c.leak0, dT10 = CL.dT_at(curve, 10.0), dT46 = CL.dT_at(curve, 46.0),
    dT150 = CL.dT_at(curve, 150.0), dT400 = CL.dT_at(curve, 400.0))
fvec(f) = Float64[getfield(f, n) for n in FNAMES]

"""Learned one-scenario policy (E4/E5): margin δ̂ = Φ(θ, x) [K] for the loop's `:scenario` mode.
kinds: `:delta_lin` (δ̂ affine in the standardised features), `:D_lin` (duration D̂ affine,
mapped through the query's heat-up curve), `:D_const` (fixed D̂, no features), `:delta_const`."""
struct Policy
    kind::Symbol; θ::Vector{Float64}; μ::Vector{Float64}; σ::Vector{Float64}
end
function scen_margin(p::Policy, curve, x::Vector{Float64})
    z = (x .- p.μ) ./ p.σ
    p.kind == :delta_const && return p.θ[1]
    p.kind == :D_const && return CL.dT_at(curve, max(p.θ[1], 0.0))
    p.kind == :delta_quad && return p.θ[1] + dot(p.θ[2:length(z)+1], z) + dot(p.θ[length(z)+2:end], z .^ 2)
    v = p.θ[1] + dot(p.θ[2:end], z)
    p.kind == :delta_lin && return v
    p.kind == :D_lin && return CL.dT_at(curve, clamp(v, 0.0, 2DCAP))
    error("policy $(p.kind)")
end
(p::Policy)(curve, c, yhat_now) = scen_margin(p, curve, fvec(features(curve, c, yhat_now)))

# ---------------------------------------------------------------- the captured query
"""Rebuild the N2 QP of a loop query exactly as `ClosedLoop.decide(:N2)` (same linearisation),
plus the heat-up curves and, if `with_GL`, the sensitivity ∂ŷ/∂L of the prediction to the leak."""
function query_record(c, L, act, prev; with_GL = false)
    mid = (act.m_lo + min(act.m_hi, 50.0)) / 2
    m_inc = prev === nothing ? fill(mid, K) : clamp.(vcat(prev[2:end], prev[end]), act.m_lo, act.m_hi)
    X = CL.rom_traj(c, m_inc, L)
    pts = [(X[k+1, :], L[k], m_inc[k]) for k in 1:K]
    C = CL.condense(c, L, X[1, :], pts)
    raw = CL.heatup_curve(c, 0.0, DCAP)                       # cummax ŷ without ACS
    Dl = [min(D, get(c, :Dcap, DCAP)) for D in T6.qgrid(40, MUACC)]
    # the loop's ŷ_now: ΔT = raw − ŷ_now matched to the loop's exp_dT (exact, constant offset)
    off = dTs(raw, Dl) .- c.exp_dT
    yhat_now = mean(off)
    curve0 = raw .- yhat_now
    curves = [e == 0 ? curve0 : CL.heatup_curve(c, yhat_now, DCAP; acs = e * MMAX) for e in ETAS]
    GL = nothing
    if with_GL                                              # condense is affine in L (fixed pts)
        GL = zeros(K, K)
        for j in 1:K
            Lp = copy(L); Lp[j] += 1.0
            GL[:, j] = CL.condense(c, Lp, X[1, :], pts).y0 .- C.y0
        end
    end
    (; y0 = C.y0, G = C.G, yref = c.y_ref, d = c.d, m_prev = c.m_prev, lo = act.m_lo, hi = act.m_hi,
     dm = act.dm, Rm = c.Rmove, curve0, curves, GL, x = fvec(features(curve0, c, yhat_now)),
     yhat_now, off_spread = maximum(off) - minimum(off), rho_loop = c.rho, Dl)
end

# ---------------------------------------------------------------- f(u), Z and the QPs
ypred(q, u) = q.y0 .+ q.G * u
function fval(q, u; quad = true)
    y = ypred(q, u); mm = vcat(q.m_prev, u)
    v = q.Rm * sum(abs2, diff(mm))
    quad && (v += sum(abs2, y .- q.yref) + F4.lam_pct * sum(abs2, max.(y .- F4.U_y, 0.0)))
    v
end
"""P1 objective f(u) + ρ·mean_δ Σ_k (ŷ_k + d + δ − 1220)₊ over the margins δ = ΔT(D)."""
function J1(q, u, δs, ρ)
    y = ypred(q, u); e = 0.0
    for δ in δs, k in 1:K
        v = y[k] + q.d + δ - ULIM; v > 0 && (e += v)
    end
    fval(q, u) + ρ * e / length(δs)
end
"""Generic second stage: f(u) + Σ_i w_i (a_i·u + b_i)₊."""
J_rows(q, u, A, b, w) = fval(q, u) + sum(w .* max.(A * u .+ b, 0.0))
"""P3 objective R Σ Δu² + mean_s Σ_k |ŷ_k + o_ks − y_ref| (offsets o = G_L δL)."""
function J3(q, u, O)
    y = ypred(q, u)
    fval(q, u; quad = false) + sum(abs.(y .+ O .- q.yref)) / size(O, 2)
end

"""N2 model: u in the box, hard ramp from the applied flow, y = y0 + G u; returns the cost f."""
function base(q; quad = true)
    md = Model(Clarabel.Optimizer); set_silent(md)
    # tight tolerances: the objectives reach 3e4 and the differences checked here are 1e-3
    for (a, v) in (("tol_gap_abs", 1e-10), ("tol_gap_rel", 1e-11), ("tol_feas", 1e-10), ("max_iter", 400))
        set_attribute(md, a, v)
    end
    @variable(md, q.lo <= u[1:K] <= q.hi); @variable(md, y[1:K])
    @constraint(md, [k = 1:K], y[k] - sum(q.G[k, j] * u[j] for j in 1:k) == q.y0[k])
    mm = [q.m_prev; u]
    @constraint(md, [k = 1:K], mm[k+1] - mm[k] <= q.dm)
    @constraint(md, [k = 1:K], mm[k] - mm[k+1] <= q.dm)
    obj = @expression(md, q.Rm * sum((mm[k+1] - mm[k])^2 for k in 1:K))
    if quad
        r = @variable(md, [1:K], lower_bound = 0)
        @constraint(md, [k = 1:K], y[k] - r[k] <= F4.U_y)
        obj = @expression(md, obj + sum((y[k] - q.yref)^2 for k in 1:K) + F4.lam_pct * sum(r .^ 2))
    end
    md, u, y, obj
end
okst(md) = termination_status(md) in OKS

"""P1 with offsets O (K×S: margin of scenario s at step k) and probabilities p. Duals λ of the
rows (∈ [0, ρ p_s]); Σ_s λ_ks is the k-th component of the expected subgradient multiplier."""
function solve_p1(q, O::AbstractMatrix, p, ρ; duals = false)
    md, u, y, obj = base(q); S = size(O, 2)
    z = @variable(md, [1:K, 1:S], lower_bound = 0)
    con = @constraint(md, [k = 1:K, s = 1:S], z[k, s] - y[k] >= q.d + O[k, s] - ULIM)
    @objective(md, Min, obj + ρ * sum(p[s] * z[k, s] for k in 1:K, s in 1:S))
    t = @elapsed optimize!(md)
    λ = duals && has_duals(md) ? [dual(con[k, s]) for k in 1:K, s in 1:S] : nothing
    (u = value.(u), ok = okst(md), t, λ)
end
scal(δs) = repeat(reshape(collect(δs), 1, :), K, 1)        # scalar margins → K×S offsets

"""Parametric P1 with ONE scalar margin δ (the one-scenario problem): built once, re-solved."""
mutable struct P1Param; md; u; con; q; end
function P1Param(q, ρ)
    md, u, y, obj = base(q)
    z = @variable(md, [1:K], lower_bound = 0)
    con = @constraint(md, [k = 1:K], z[k] - y[k] >= 0.0)
    @objective(md, Min, obj + ρ * sum(z))
    P1Param(md, u, con, q)
end
function solve_at!(P::P1Param, δ)
    for k in 1:K; set_normalized_rhs(P.con[k], P.q.d + δ - ULIM); end
    optimize!(P.md)
    (u = value.(P.u), ok = okst(P.md))
end

"""min f(u) + Σ_i w_i (A_i u + b_i)₊ (dense rows; P2). Duals λ_i ∈ [0, w_i]."""
function solve_rows(q, A, b, w; duals = false)
    md, u, y, obj = base(q); nr = length(b)
    z = @variable(md, [1:nr], lower_bound = 0)
    con = @constraint(md, [i = 1:nr], z[i] - sum(A[i, j] * u[j] for j in 1:K if A[i, j] != 0) >= b[i])
    @objective(md, Min, obj + sum(w[i] * z[i] for i in 1:nr))
    t = @elapsed optimize!(md)
    λ = duals && has_duals(md) ? [dual(con[i]) for i in 1:nr] : nothing
    (u = value.(u), ok = okst(md), t, λ)
end

"""P3: R Σ Δu² + Σ_s p_s Σ_k |y_k + O_ks − y_ref| (L1 tracking, no quadratic part)."""
function solve_p3(q, O, p)
    md, u, y, obj = base(q; quad = false); S = size(O, 2)
    z = @variable(md, [1:K, 1:S], lower_bound = 0)
    @constraint(md, [k = 1:K, s = 1:S], z[k, s] - y[k] >= O[k, s] - q.yref)
    @constraint(md, [k = 1:K, s = 1:S], z[k, s] + y[k] >= q.yref - O[k, s])
    @objective(md, Min, obj + sum(p[s] * z[k, s] for k in 1:K, s in 1:S))
    t = @elapsed optimize!(md)
    (u = value.(u), ok = okst(md), t)
end

# ---------------------------------------------------------------- P2 scenarios
"""Heat-up curve for a severity η ∈ [0, 1] (linear interpolation between the computed ones)."""
function curve_eta(q, η)
    i = findfirst(==(η), ETAS); i === nothing || return q.curves[i]
    e = sort(collect(ETAS)); η = clamp(η, e[1], e[end])
    k = searchsortedlast(e, η); k == length(e) && (k -= 1)
    w = (η - e[k]) / (e[k+1] - e[k])
    ca = q.curves[findfirst(==(e[k]), ETAS)]; cb = q.curves[findfirst(==(e[k+1]), ETAS)]
    (1 - w) .* ca .+ w .* cb
end
"""Rows of one partial-failure scenario (start j, severity η, duration D). `perstep = false`:
one row at the end of the failure (or of the horizon), as `ClosedLoop` (:scen_exp); true: one
row per step k (m = K): before j the nominal excess, from j on the delivered flow η u until the
failure ends, the remaining heat-up added at that step."""
function p2_scen(q, j, η, D, curve; perstep = false)
    ks = min(K, j + ceil(Int, D / h)); Drem = max(D - (ks - j + 1) * h, 0.0)
    rem = Drem > 0 ? CL.dT_at(curve, Drem) : 0.0
    if !perstep
        e = [i < j ? 1.0 : η for i in 1:K]; e[ks+1:end] .= 0
        return reshape(q.G[ks, :] .* e, 1, K), [q.y0[ks] + q.d + rem - ULIM]
    end
    e = [i < j || i > ks ? 1.0 : η for i in 1:K]
    A = q.G .* e'
    b = [q.y0[k] + q.d + (k == ks ? rem : 0.0) - ULIM for k in 1:K]
    A, b
end
function p2_set(q; perstep = false, Dq = qD(10))
    A = zeros(0, K); b = Float64[]; w = Float64[]
    for j in JSET, (ie, η) in enumerate(ETAS), D in Dq
        a, bb = p2_scen(q, j, η, D, q.curves[ie]; perstep)
        A = vcat(A, a); append!(b, bb)
        append!(w, fill(PETA[ie] / (length(JSET) * length(Dq)), length(bb)))
    end
    A, b, w
end

# ---------------------------------------------------------------- P3 forecast error
const _SL = Ref{Any}(nothing)
"""Leak-forecast covariance Σ_L over the horizon: steady-state local-linear-trend Kalman filter
with the calibrated parameters for σ = 0.1 kg/s measurement noise (F5b.0)."""
function sigmaL()
    _SL[] === nothing || return _SL[]
    f = CL.LeakKF(0.0; CL.KF_PARAMS[0.1]..., thr = Inf)
    for _ in 1:20000; CL.kf_update!(f, 0.0); end
    _SL[] = Symmetric(CL.kf_forecast(f, K).Σ)
end

# ---------------------------------------------------------------- golden section, Nelder–Mead
function golden(f, a, b; it = 14)
    φ = (sqrt(5) - 1) / 2; c = b - φ * (b - a); d = a + φ * (b - a); fc = f(c); fd = f(d)
    for _ in 1:it
        if fc < fd; b, d, fd = d, c, fc; c = b - φ * (b - a); fc = f(c)
        else;       a, c, fc = c, d, fd; d = a + φ * (b - a); fd = f(d); end
    end
    fc < fd ? (c, fc) : (d, fd)
end
"""Nelder–Mead (derivative-free, as Algorithm 1 of the paper with Dias-Garcia et al.); stops when
the best value improves by less than `tol` over `stall` iterations, `restarts` restarts."""
function nelder_mead(f, x0; step = 1.0, tol = 1e-7, maxit = 5000, restarts = 2)
    x = copy(x0); fx = f(x)
    for r in 0:restarts
        n = length(x); S = [copy(x)]
        for i in 1:n; v = copy(x); v[i] += step isa Number ? step : step[i]; push!(S, v); end
        F = map(f, S)
        for it in 1:maxit
            p = sortperm(F); S = S[p]; F = F[p]
            (F[end] - F[1]) <= tol * (abs(F[1]) + 1e-12) && break
            xc = sum(S[1:n]) / n; xr = xc .+ (xc .- S[end]); fr = f(xr)
            if fr < F[1]
                xe = xc .+ 2 .* (xc .- S[end]); fe = f(xe)
                fe < fr ? (S[end] = xe; F[end] = fe) : (S[end] = xr; F[end] = fr)
            elseif fr < F[n]
                S[end] = xr; F[end] = fr
            else
                xk = fr < F[end] ? xc .+ 0.5 .* (xr .- xc) : xc .+ 0.5 .* (S[end] .- xc); fk = f(xk)
                if fk < min(fr, F[end]); S[end] = xk; F[end] = fk
                else
                    for i in 2:n+1; S[i] = S[1] .+ 0.5 .* (S[i] .- S[1]); F[i] = f(S[i]); end
                end
            end
        end
        i = argmin(F)
        F[i] < fx - tol * (abs(fx) + 1e-12) || (F[i] < fx && (x = S[i]; fx = F[i]); break)
        x = S[i]; fx = F[i]
    end
    x, fx
end

# ---------------------------------------------------------------- one query: E0–E3
"""All per-query experiments. `full`: also P2 and P3 (costlier); `scaling`: CPU time vs N."""
function analyse(seed, kq, c, L, act, prev, dec; full = false, scaling = false)
    q = query_record(c, L, act, prev; with_GL = full)
    row = Dict{Symbol,Any}(:seed => seed, :kq => kq, :full => full, :off_spread => q.off_spread)
    for (n, v) in zip(FNAMES, q.x); row[Symbol("x_", n)] = v; end
    grids = Dict{Float64,Any}()
    D200 = qD(200); D1000 = qD(1000); δ200 = dTs(q.curve0, D200); δ1000 = dTs(q.curve0, D1000)
    rng = MersenneTwister(1_000_003 * seed + kq)
    Dn = [min(CL.failure_sample(rng, dist()), DCAP) for _ in 1:NS_TRAIN]; δn = dTs(q.curve0, Dn)
    # E0: the captured QP reproduces the loop decision (40 stratified quantiles, ρ of the loop)
    s0 = solve_p1(q, scal(dTs(q.curve0, q.Dl)), fill(1 / 40, 40), q.rho_loop)
    row[:E0_gap] = maxabs(clamp.(s0.u, q.lo, q.hi), dec.m)
    # E1: Assumption 3 (strict convexity) and the sign pattern of G (Lemma 1 (iii))
    Dm = Matrix(1.0I, K, K) .- diagm(-1 => ones(K - 1))
    row[:lmin_GtG] = eigmin(Symmetric(q.G' * q.G))
    row[:lmin_H] = eigmin(Symmetric(q.G' * q.G + q.Rm * Dm' * Dm))
    lt = [q.G[k, i] for k in 1:K for i in 1:k]
    row[:G_maxpos] = maximum(lt); row[:G_minneg] = minimum(lt)
    row[:G_colsign_ok] = all(j -> all(q.G[k, j] <= 1e-12 for k in j:K) || all(q.G[k, j] >= -1e-12 for k in j:K), 1:K)
    row[:dT0] = q.curve0[1]; row[:dT400] = CL.dT_at(q.curve0, DCAP)
    for ρ in RHOS
        tag = ρ == 100 ? "r100" : "r1000"
        # SAA (200 stratified quantiles) = π_S of the paper for this query
        sS = solve_p1(q, scal(δ200), fill(1 / 200, 200), ρ; duals = true)
        uS = sS.u; JS_in = J1(q, uS, δ200, ρ); JS = J1(q, uS, δ1000, ρ)
        g = sS.λ === nothing ? fill(NaN, K) : vec(sum(sS.λ, dims = 2))   # Σ_s λ_ks ∈ [0, ρ]
        row[Symbol("$(tag)_ok")] = sS.ok; row[Symbol("$(tag)_JS")] = JS; row[Symbol("$(tag)_uS1")] = uS[1]
        row[Symbol("$(tag)_nact")] = count(>(1e-6ρ), g); row[Symbol("$(tag)_gmax")] = maximum(g) / ρ
        Jnom = J1(q, solve_p1(q, scal([-1e3]), [1.0], ρ).u, δ1000, ρ)          # risk ignored
        row[Symbol("$(tag)_Jnom")] = Jnom
        # E2 (case i): ξ* = (T z*, T): margin δ*_k = 1220 − d − ŷ_k(u_S), one per step
        δst = ULIM .- q.d .- ypred(q, uS)
        s2 = solve_p1(q, reshape(δst, K, 1), [1.0], ρ)
        row[Symbol("$(tag)_T1_gap")] = maxabs(s2.u, uS)
        row[Symbol("$(tag)_T1_dJ")] = J1(q, s2.u, δ1000, ρ) - JS
        row[Symbol("$(tag)_T1_spread")] = maximum(δst) - minimum(δst)
        row[Symbol("$(tag)_T1_min")] = minimum(δst); row[Symbol("$(tag)_T1_max")] = maximum(δst)
        # E3: bilevel (38)–(39) over the PHYSICAL family "one total failure of duration D̂"
        # (scalar margin δ̂), derivative-free: grid + golden section on the SAA objective
        P = P1Param(q, ρ)
        tg = @elapsed U = reduce(hcat, [solve_at!(P, δ).u for δ in DGRID])
        Jin = [J1(q, U[:, i], δ200, ρ) for i in eachindex(DGRID)]
        Jex = [J1(q, U[:, i], δ1000, ρ) for i in eachindex(DGRID)]
        Rn = [J1(q, U[:, i], [δ], ρ) for δ in δn, i in eachindex(DGRID)]
        ib = argmin(Jin)
        lo = DGRID[max(ib - 1, 1)]; hi = DGRID[min(ib + 1, length(DGRID))]
        δb, Jb = golden(δ -> J1(q, solve_at!(P, δ).u, δ200, ρ), lo, hi)
        Jb > Jin[ib] && (δb = DGRID[ib]; Jb = Jin[ib])
        ub = solve_at!(P, δb).u
        row[Symbol("$(tag)_BL_delta")] = δb
        row[Symbol("$(tag)_BL_D")] = δb < q.curve0[1] ? -1.0 : δb > CL.dT_at(q.curve0, DCAP) ? Inf :
            (D = range(0, DCAP, length = 4001); D[findfirst(x -> CL.dT_at(q.curve0, x) >= δb, D)])
        row[Symbol("$(tag)_BL_gap")] = maxabs(ub, uS); row[Symbol("$(tag)_BL_gap1")] = abs(ub[1] - uS[1])
        row[Symbol("$(tag)_BL_dJin")] = Jb - JS_in; row[Symbol("$(tag)_BL_dJ")] = J1(q, ub, δ1000, ρ) - JS
        for (nm, D) in (("mean", sum(dist().w ./ dist().mu)), ("median", T6.CL.failure_quantile(0.5, dist())),
                        ("q98", T6.CL.failure_quantile(0.98, dist())))
            un = solve_at!(P, CL.dT_at(q.curve0, min(D, DCAP))).u
            row[Symbol("$(tag)_$(nm)_dJ")] = J1(q, un, δ1000, ρ) - JS
            row[Symbol("$(tag)_$(nm)_gap")] = maxabs(un, uS)
        end
        row[Symbol("$(tag)_tgrid")] = tg
        grids[ρ] = (u1 = Float32.(U[1, :]), Jex = Float32.(Jex), R = Float32.(Rn), JS = JS, uS1 = uS[1], Jnom = Jnom)
        # CPU time vs N (paper Table 3): SAA as one QP vs the one-scenario search
        if scaling && ρ == 1000.0
            for N in (10, 40, 100, 400, 1000, 4000)
                r2 = MersenneTwister(seed + N)
                δN = dTs(q.curve0, [min(CL.failure_sample(r2, dist()), DCAP) for _ in 1:N])
                tS = @elapsed sN = solve_p1(q, scal(δN), fill(1 / N, N), ρ)
                tE = @elapsed JN = [J1(q, U[:, i], δN, ρ) for i in eachindex(DGRID)]
                i2 = argmin(JN)
                row[Symbol("cpu_saa_N$N")] = tS; row[Symbol("cpu_bl_N$N")] = tg + tE
                row[Symbol("cpu_gap_N$N")] = maxabs(U[:, i2], sN.u)
            end
        end
    end
    full || return (row = row, grids = grids, x = q.x, curve = Float32.(q.curve0))
    # ---------------- P2: partial failures inside the horizon (T random), ρ as the loop's 1000
    for ρ in RHOS
        tag = ρ == 100 ? "p2r100" : "p2r1000"
        A, b, w = p2_set(q); wρ = ρ .* w
        sS = solve_rows(q, A, b, wρ; duals = true); uS = sS.u; JS = J_rows(q, uS, A, b, wρ)
        row[Symbol("$(tag)_ok")] = sS.ok
        row[Symbol("$(tag)_VSS")] = J_rows(q, solve_p1(q, scal([-1e3]), [1.0], ρ).u, A, b, wρ) - JS
        # Lemma 1 does not cover the implemented form (m = 1 < n = 26); the scenario only has to
        # reproduce the subgradient at u_S: one row a* = g/ρ with kink at u_S (g = Σ λ_i a_i)
        gsub = sS.λ === nothing ? zeros(K) : vec(sum(sS.λ .* A, dims = 1))
        a1 = gsub ./ ρ
        s1 = solve_rows(q, reshape(a1, 1, K), [-dot(a1, uS)], [ρ])
        row[Symbol("$(tag)_row1_gap")] = maxabs(s1.u, uS)
        row[Symbol("$(tag)_row1_dJ")] = J_rows(q, s1.u, A, b, wρ) - JS
        row[Symbol("$(tag)_gnorm")] = norm(gsub)
        # per-step form (m = K = n): Lemma 1 (iii) ⇒ T* = (R/ρ)·I, h* = T* u*
        Ap, bp, wp = p2_set(q; perstep = true); wpρ = ρ .* wp
        sP = solve_rows(q, Ap, bp, wpρ); uP = sP.u; JP = J_rows(q, uP, Ap, bp, wpρ)
        # E[∂Q] ⊂ simplex conv{0, −R e_i} with R = ρ Σ_s p_s Σ_ki |a_ki|; the lemma's M* = −(R/ρ) I,
        # so the scenario is the row set (c (u*_i − u_i))₊, c = R/ρ, weight ρ (a floor on each flow)
        cst = sum(wp[i] * sum(abs, Ap[i, :]) for i in eachindex(bp))
        sL = solve_rows(q, -cst .* Matrix(1.0I, K, K), cst .* uP, fill(ρ, K))
        row[Symbol("$(tag)_lemma_gap")] = maxabs(sL.u, uP)
        row[Symbol("$(tag)_lemma_dJ")] = J_rows(q, sL.u, Ap, bp, wpρ) - JP
        row[Symbol("$(tag)_lemma_c")] = cst
        # E3: bilevel over PHYSICAL single scenarios (ĵ, η̂, D̂), on the implemented SAA objective
        Dg = vcat([min(T6.CL.failure_quantile(p, dist()), DCAP) for p in 0.05:0.1:0.95], DCAP)
        best = (J = Inf, j = 0, η = NaN, D = NaN, u = uS)
        for j in JSET, (ie, η) in enumerate(ETAS), D in Dg
            a, bb = p2_scen(q, j, η, D, q.curves[ie])
            s = solve_rows(q, a, bb, [ρ]); Jv = J_rows(q, s.u, A, b, wρ)
            Jv < best.J && (best = (J = Jv, j, η, D, u = s.u))
        end
        jb = best.j
        function fphys(v)
            η = clamp(v[1], 0.0, 1.0); D = clamp(v[2], 0.0, DCAP)
            a, bb = p2_scen(q, jb, η, D, curve_eta(q, η))
            J_rows(q, solve_rows(q, a, bb, [ρ]).u, A, b, wρ)
        end
        v, Jv = nelder_mead(fphys, [best.η, best.D]; step = [0.1, 20.0], tol = 1e-9, maxit = 60, restarts = 0)
        if Jv < best.J
            η = clamp(v[1], 0, 1); D = clamp(v[2], 0, DCAP); a, bb = p2_scen(q, best.j, η, D, curve_eta(q, η))
            best = (J = Jv, j = best.j, η, D, u = solve_rows(q, a, bb, [ρ]).u)
        end
        row[Symbol("$(tag)_phys_dJ")] = best.J - JS; row[Symbol("$(tag)_phys_gap")] = maxabs(best.u, uS)
        row[Symbol("$(tag)_phys_j")] = best.j; row[Symbol("$(tag)_phys_eta")] = best.η; row[Symbol("$(tag)_phys_D")] = best.D
    end
    # ---------------- P3: leak-forecast error, L1 tracking (h random, T fixed)
    Σ = sigmaL(); Lc = cholesky(Symmetric(Matrix(Σ) + 1e-14I)).L; σL = sqrt.(diag(Σ))
    r3 = MersenneTwister(77 * seed + kq)
    O200 = q.GL * (Lc * randn(r3, K, 200)); O2000 = q.GL * (Lc * randn(r3, K, 2000))
    s3 = solve_p3(q, O200, fill(1 / 200, 200)); u3 = s3.u; J3S = J3(q, u3, O2000)
    uCE = solve_p3(q, zeros(K, 1), [1.0]).u
    row[:p3_ok] = s3.ok; row[:p3_CE_gap] = maxabs(uCE, u3); row[:p3_CE_dJ] = J3(q, uCE, O2000) - J3S
    ost = q.yref .- ypred(q, u3)                          # ξ*: every kink at the SAA plan
    uT = solve_p3(q, reshape(ost, K, 1), [1.0]).u
    row[:p3_T1_gap] = maxabs(uT, u3); row[:p3_T1_dJ] = J3(q, uT, O2000) - J3S
    δLst = q.GL \ ost
    row[:p3_T1_zmax] = maximum(abs.(δLst) ./ max.(σL, 1e-12))   # is ξ* a plausible leak error?
    row[:p3_sigma_y] = maximum(sqrt.(max.(diag(q.GL * Matrix(Σ) * q.GL'), 0)))
    sg = -3.0:0.1:3.0
    Js = [J3(q, solve_p3(q, reshape(q.GL * (s .* σL), K, 1), [1.0]).u, O200) for s in sg]
    sb = sg[argmin(Js)]; ub3 = solve_p3(q, reshape(q.GL * (sb .* σL), K, 1), [1.0]).u
    row[:p3_BL_s] = sb; row[:p3_BL_gap] = maxabs(ub3, u3); row[:p3_BL_dJ] = J3(q, ub3, O2000) - J3S
    row[:p3_J] = J3S
    (row = row, grids = grids, x = q.x, curve = Float32.(q.curve0))
end

# ---------------------------------------------------------------- phase A: capture + analyse
"""Closed loop of seed `seed` (DET-accelerated family, 40 kg/s, N2 + expected excess ρ = 1000,
R = 0.1, DET failures) with every `every`-th query analysed; P2/P3 on every 4th analysed query.
Saves results/f9_rozas/A_<seed>.jls (resumable)."""
function phaseA(seed, outdir; every = 6, ρloop = 1000.0)
    file = joinpath(outdir, "A_$(seed).jls"); isfile(file) && return :done
    tk = T6.mk("f9", :fail, :detacc, seed, :cap40, "N2_exp_r$(Int(ρloop))_R0.1", :N2;
               spec = (:exp, ρloop, MUACC), mu_design = MUACC, R = RM)
    recs = Any[]; errs = String[]; kq = Ref(0); na = Ref(0)
    T6.CAPTURE[] = (t, c, L, act, prev, dec) -> begin
        kq[] += 1
        (kq[] - 1) % every == 0 || return
        na[] += 1
        try
            push!(recs, analyse(seed, kq[], c, L, act, prev, dec; full = (na[] - 1) % 4 == 0, scaling = na[] == 1))
        catch err
            push!(errs, "kq=$(kq[]): " * sprint(showerror, err))
        end
    end
    tl = @elapsed try
        T6.run_task(tk)
    finally
        T6.CAPTURE[] = nothing
    end
    serialize(file * ".tmp", (recs = recs, errs = errs, wall = tl)); mv(file * ".tmp", file; force = true)
    (seed = seed, n = length(recs), errs = length(errs), wall = tl)
end

# ---------------------------------------------------------------- phase C: closed loop
"""Failure-free closed loop + Rao–Blackwell P(reach 1220 K | failure) (as F6 :rb, μ×1e4)."""
function run_rb(seed, case, riskkw)
    bb = T6.bench(:detacc, seed, :cap40); B = bb.B
    B === nothing && return nothing
    act = T6.actuator(B, :cap40)
    tr = @elapsed lg = CL.run_loop(B, :N2; n = T6.NQ_DET, act, bias = true, keep_states = true, Rmove = RM, riskkw...)
    idx = unique(vcat(1:20:length(lg.t_q), length(lg.t_q)))
    τ = [CL.time_to_limit(B.pH, lg.states[i][1], lg.states[i][2], lg.t_q[i]) for i in idx]
    tq = lg.t_q[idx] .- lg.t_q[1]
    p = T6.rb_pair(tq, k -> CL.failure_surv(τ[k], dist()))
    M = CL.metrics(lg, B.y_ref)
    (; seed, case, p_cond = p.cond, p_unc = p.unc, mean_dev = mean(lg.pct .- B.y_ref), mass = M.mass_plant,
     ise = M.ise, iae = M.iae, pct_max = M.pct_max, tv = sum(abs.(diff(vcat(0.0, lg.m)))),
     t_p50 = M.t_p50, t_p99 = M.t_p99, fails = M.fails, n = M.n, wall = tr)
end

end # module
