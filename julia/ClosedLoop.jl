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
using .HRMQueries, .F4Lib, MAT, JuMP, Clarabel, ForwardDiff, LinearAlgebra, Statistics

export Bench, setup, run_loop, ACT_ARCH, ACT_PHYS

const ROOT = joinpath(@__DIR__, "..")
const h = 0.5
const K = 26
const DELAY = 10.0
const R4 = F4Lib.SMDFRROM

const ACT_ARCH = (name = :arch, m_lo = NaN, m_hi = NaN, ramp = :soft)   # filled in setup
const ACT_PHYS = (name = :phys, m_lo = 0.0, m_hi = 5543.0, ramp = :hard)

"""Everything a closed-loop run needs: plant parameters, the plant history before the
first query, the ROM parameters and the archived references."""
struct Bench
    pH; pR; pR4; y_ref::Float64; act_arch; cs::Matrix{Float64}
    t_start::Float64
    hist::Dict{Float64,Any}     # HRM (x, logic, pct) on the grid before t_start
end

function setup(; t_start = 1160.0)
    arc = matread(joinpath(ROOT, "results", "f1b_hrm_replica", "archive_1180.mat"))
    cs = arc["control_series"]
    pH = hrm_params(A = arc["sampled_A"], B = arc["sampled_B"], mdot_c_drift = arc["sampled_mdot_c_drift"])
    ref = matread(joinpath(ROOT, "results", "f1_rom_replica", "rom_reference_1180.mat"))
    ub = ref["u_bounds"]
    act_arch = (name = :arch, m_lo = current_to_massflow(minimum(ub[:, 1])),
                m_hi = current_to_massflow(maximum(ub[:, 2])), ramp = :soft)
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
          cs, t_start, hist)
end

cost_step(y, y_ref) = (y - y_ref)^2

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
    if act.ramp == :soft
        s = @variable(md, [1:K-1], lower_bound = 0)
        @constraint(md, [k = 1:K-1], u[k+1] - u[k] <= F4Lib.dm_lim + s[k])
        @constraint(md, [k = 1:K-1], u[k] - u[k+1] <= F4Lib.dm_lim + s[k])
        obj = @expression(md, obj + F4Lib.lam_m * sum(s .^ 2))
    else
        mm = [c.m_prev; u]
        @constraint(md, [k = 1:K], mm[k+1] - mm[k] <= F4Lib.dm_lim)
        @constraint(md, [k = 1:K], mm[k] - mm[k+1] <= F4Lib.dm_lim)
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

# ---------------------------------------------------------------- the loop
"""
    run_loop(B, ctrl; n = 61, act = B.act_arch, bias = false, lj = nothing)

Closed loop from t_start over `n` queries. Returns a NamedTuple of per-query logs.
`lj`: archived (t, current) pairs used by `ctrl = :LJ`.
"""
function run_loop(B::Bench, ctrl::Symbol; n = 61, act = B.act_arch, bias = false)
    pH = B.pH; hist = copy(B.hist)
    last = findfirst(==(B.t_start - h), B.cs[:, 1])
    T = B.cs[2:last, 1]; U = B.cs[2:last, 2:3]          # applied inport series (grows)
    m_applied = 0.0; prev = nothing
    log = (t_q = Float64[], pct = Float64[], m = Float64[], u = Float64[], yhat_now = Float64[],
           d = Float64[], time = Float64[], status = String[], m_plan = Vector{Float64}[])
    for k in 1:n
        t_q = B.t_start + (k - 1) * h; τ = t_q - h
        x, Lg, pct_now = hist[τ]
        xH10, L10, _ = hist[t_q-DELAY]
        ext = ExtInput(T, U)
        w = warm_rom(pH, B.pR, xH10, L10.xI, ext, t_q; delay = DELAY, h)
        L = [SMDFRHRM.leak_flow(t_q + (j - 1) * h, pH) for j in 1:K]        # perfect forecast
        yhat_now = R4.output(w.x_init, w.leak_warm[end], w.acs_warm[end], B.pR4)
        d = pct_now - yhat_now
        c = (; p = B.pR4, h, K, S = max.(abs.(w.x_init), 1.0), x_init = w.x_init,
             leak0 = w.leak_warm[end], y_ref = bias ? B.y_ref - d : B.y_ref, m_prev = m_applied)
        if ctrl in (:N1, :N2, :N3)
            dec = decide(ctrl, c, L, act, prev)
            m1 = dec.m[1]; prev = dec.m
        else
            m1 = ctrl == :zero ? 0.0 : ctrl == :max ? act.m_hi :
                 ctrl == :LJ ? current_to_massflow(B.cs[findfirst(==(t_q), B.cs[:, 1]), 2]) : error("ctrl")
            dec = (m = fill(m1, K), status = "-", time = 0.0)
        end
        u1 = m1 > 0 ? F4Lib.massflow_to_current(m1) : 56.0                  # 56 A: pump off
        T = vcat(T, t_q); U = vcat(U, [u1 1.0])
        s = SMDFRHRM.simulate(pH, ExtInput(T, U), t_q; t0 = τ, x0 = x, logic = Lg)
        @assert abs(s.acs[end] - m1) <= 1e-6 * max(1.0, m1) "plant did not receive the decided flow"
        hist[t_q] = (s.x, s.logic, s.pct[end])
        m_applied = m1
        push!(log.t_q, t_q); push!(log.pct, s.pct[end]); push!(log.m, m1); push!(log.u, u1)
        push!(log.yhat_now, yhat_now); push!(log.d, d); push!(log.time, dec.time)
        push!(log.status, dec.status); push!(log.m_plan, dec.m)
    end
    log
end

end # module
