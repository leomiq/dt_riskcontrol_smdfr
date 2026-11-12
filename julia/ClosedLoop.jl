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

export Bench, setup, run_loop, metrics, ACT_ARCH, ACT_PHYS

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
    t_act::Float64              # ACS activation by the DET logic (NaN if unknown)
end

"""
    setup(; t_start = 1160.0, pH = nothing, lead = 17.0, t_max = 2500.0)

Archived realisation (default) or a new one with HRM parameters `pH` (leak A, B, drift).
For a new realisation the MPC starts `lead` s after the DET logic activates the ACS (the
archived loop: activation 1143 s, first query 1160 s); returns `nothing` if the ACS is never
demanded before `t_max`.
"""
function setup(; t_start = 1160.0, pH = nothing, lead = 17.0, t_max = 2500.0)
    ref = matread(joinpath(ROOT, "results", "f1_rom_replica", "rom_reference_1180.mat"))
    ub = ref["u_bounds"]
    act_arch = (name = :arch, m_lo = current_to_massflow(minimum(ub[:, 1])),
                m_hi = current_to_massflow(maximum(ub[:, 2])), ramp = :soft)
    t_act = NaN
    if pH === nothing
        arc = matread(joinpath(ROOT, "results", "f1b_hrm_replica", "archive_1180.mat"))
        cs = arc["control_series"]
        pH = hrm_params(A = arc["sampled_A"], B = arc["sampled_B"], mdot_c_drift = arc["sampled_mdot_c_drift"])
    else
        s = SMDFRHRM.simulate(pH, ExtInput([0.0, 1.0], zeros(2, 2)), t_max)
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
    run_loop(B, ctrl; n = 61, act = B.act_arch, bias = false, eta = 1.0, fail = nothing,
             fbias = 0.0)

Closed loop from t_start over `n` queries. Returns a NamedTuple of per-query logs.
Plant perturbations the MPC does not know about (F5a.2):
  `eta`   ACS efficiency: the plant receives η·ṁ for a commanded ṁ (degraded ACS);
  `fail`  (t_a, t_b): the ACS delivers nothing for t_a < t ≤ t_b (forced failure);
  `fbias` constant error of the leak forecast [kg/s] (ŷ uses ṁ_leak + fbias).
  `audit(c, L, act, prev, dec)` optional callback after each MPC decision (same query data).
`ctrl = :PI` leaves the plant to the DET's own PI (active_mpc = 0); only `fail` applies to it
(the pump is forced off), `eta` does not. If the plant model breaks down (coolant inventory
exhausted, the HRM algebra fails) the loop stops and `log.breakdown[1]` is set.
When the PCT exceeds the upper edge of the DET monitor band (1220 K) the loop stops and
`log.limit[1]` is set: in Simulink the limintm block stops the simulation there (end state of
the event tree), so the replica's continuation (protection off) has no meaning.
The ROM warm-up uses the COMMANDED currents (what the twin knows), the plant the real ones.
"""
function run_loop(B::Bench, ctrl::Symbol; n = 61, act = B.act_arch, bias = false, eta = 1.0,
                  fail = nothing, fbias = 0.0, audit = nothing)
    pH = B.pH; hist = copy(B.hist)
    last = findfirst(==(B.t_start - h), B.cs[:, 1])
    T = B.cs[2:last, 1]; U = B.cs[2:last, 2:3]          # commanded inport series (grows)
    Up = copy(U)                                        # what the plant actually receives
    m_applied = 0.0; prev = nothing
    log = (t_q = Float64[], pct = Float64[], m = Float64[], m_plant = Float64[], u = Float64[],
           yhat_now = Float64[], d = Float64[], time = Float64[], status = String[],
           m_plan = Vector{Float64}[], breakdown = [false], limit = [false])
    for k in 1:n
        t_q = B.t_start + (k - 1) * h; τ = t_q - h
        x, Lg, pct_now = hist[τ]
        xH10, L10, _ = hist[t_q-DELAY]
        failed = fail !== nothing && fail[1] < t_q <= fail[2]
        if ctrl == :PI                       # the DET's own PI drives the pump (no twin needed)
            T = vcat(T, t_q); U = vcat(U, [0.0 0.0]); Up = vcat(Up, failed ? [56.0 1.0] : [0.0 0.0])
            s = try
                SMDFRHRM.simulate(pH, ExtInput(T, Up), t_q; t0 = τ, x0 = x, logic = Lg)
            catch err
                err isa DomainError || rethrow(); log.breakdown[1] = true; break
            end
            hist[t_q] = (s.x, s.logic, s.pct[end])
            push!(log.t_q, t_q); push!(log.pct, s.pct[end]); push!(log.m, s.acs[end])
            push!(log.m_plant, s.acs[end]); push!(log.u, s.current[end]); push!(log.yhat_now, NaN)
            push!(log.d, NaN); push!(log.time, 0.0); push!(log.status, "-"); push!(log.m_plan, [s.acs[end]])
            s.pct[end] > pH.monitor_band[2] && (log.limit[1] = true; break)
            continue
        end
        ext = ExtInput(T, U)
        w = try
            warm_rom(pH, B.pR, xH10, L10.xI, ext, t_q; delay = DELAY, h)
        catch err
            err isa DomainError || rethrow(); log.breakdown[1] = true; break
        end
        L = [SMDFRHRM.leak_flow(t_q + (j - 1) * h, pH) + fbias for j in 1:K]   # forecast
        yhat_now = R4.output(w.x_init, w.leak_warm[end], w.acs_warm[end], B.pR4)
        d = pct_now - yhat_now
        c = (; p = B.pR4, h, K, S = max.(abs.(w.x_init), 1.0), x_init = w.x_init,
             leak0 = w.leak_warm[end], y_ref = bias ? B.y_ref - d : B.y_ref, m_prev = m_applied)
        if ctrl in (:N1, :N2, :N3)
            dec = decide(ctrl, c, L, act, prev)
            audit === nothing || audit(c, L, act, prev, dec)
            m1 = dec.m[1]; prev = dec.m
        else
            m1 = ctrl == :zero ? 0.0 : ctrl == :max ? act.m_hi :
                 ctrl == :LJ ? current_to_massflow(B.cs[findfirst(==(t_q), B.cs[:, 1]), 2]) : error("ctrl")
            dec = (m = fill(m1, K), status = "-", time = 0.0)
        end
        u1 = m1 > 0 ? F4Lib.massflow_to_current(m1) : 56.0                  # 56 A: pump off
        mp = failed ? 0.0 : eta * m1                                         # real flow
        up = mp > 0 ? F4Lib.massflow_to_current(mp) : 56.0
        T = vcat(T, t_q); U = vcat(U, [u1 1.0]); Up = vcat(Up, [up 1.0])
        s = try
            SMDFRHRM.simulate(pH, ExtInput(T, Up), t_q; t0 = τ, x0 = x, logic = Lg)
        catch err
            err isa DomainError || rethrow(); log.breakdown[1] = true; break
        end
        # the DET logic passes the MPC flow only while the protection is active (e.g. not yet
        # in the first steps after the activation, scheduler delay): log what the plant got
        if s.act[end] > 0
            @assert abs(s.acs[end] - mp) <= 1e-6 * max(1.0, mp) "plant did not receive the intended flow"
        else
            mp = s.acs[end]
        end
        hist[t_q] = (s.x, s.logic, s.pct[end])
        m_applied = m1
        push!(log.t_q, t_q); push!(log.pct, s.pct[end]); push!(log.m, m1); push!(log.m_plant, mp)
        push!(log.u, u1)
        push!(log.yhat_now, yhat_now); push!(log.d, d); push!(log.time, dec.time)
        push!(log.status, dec.status); push!(log.m_plan, dec.m)
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

end # module
