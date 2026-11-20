"""
    F7Cert

Global-optimality certificate of closed-loop MPC queries (plan §10, certification point).

The problem certified is exactly the one N3 solves in the loop (`F4Lib.solve_det`):

    P(c, L, act):  min_u  J(u) = Σ_k (ŷ_k(u) − y_ref)² + λ_pct Σ_k [ŷ_k(u) − U_y]₊²  (+ soft ramp)
                   s.t.   act.m_lo ≤ u_k ≤ act.m_hi,  |u_k − u_{k−1}| ≤ 5 kg/s (u_0 = ṁ_prev, hard ramp)

with ŷ(u) from the converged implicit-Euler ROM (the NLP's model; states eliminated).

  * Upper bounds (UB): J of any feasible u. Candidates: the loop's own N3 decision, N1, N2,
    an Ipopt multistart and, when a licence is found, Knitro (multistart + its four
    algorithms). A local solver can only give UBs: it never proves that nothing better exists.
  * Lower bounds (LB): (i) LB = 0, valid because J ≥ 0; it already certifies the queries whose
    best plan has J ≈ 0. (ii) Spatial branch and bound over the control box. On each node
    (a sub-box B) the residual of an affine model of ŷ around the box centre is enclosed
    RIGOROUSLY with interval arithmetic (Taylor model with interval remainder, Proposition
    "Encierro del resto no lineal" of the status report, f3c_rigorous.jl) and the convex
    relaxation of Proposition "Cota inferior por relajación del residuo" is solved as a QP
    (Clarabel). min over open nodes of the node LBs is a global LB of P; the gap UB − LB
    closes as the boxes shrink (the residual falls as width²).

Certificate: J(u_method) − J* ≤ UB_method − LB, valid for any method's decision.
"""
module F7Cert

using LinearAlgebra, Printf, Random, ForwardDiff, IntervalArithmetic, JuMP, Ipopt, Clarabel

const DM = 5.0          # ramp limit per step [kg/s]
const U_Y = 1244.0; const LAM_PCT = 10.0; const LAM_M = 1e-4

# ------------------------------------------------------------------------------------------
# the problem
# ------------------------------------------------------------------------------------------
"""Problem data of one query. `S` is the ROM module (F4Lib.SMDFRROM), `F4` the F4Lib module."""
function problem(c, L, act; S, F4)
    K = c.K
    lo = fill(float(act.m_lo), K); hi = fill(float(act.m_hi), K)
    hard = act.ramp == :hard
    ok = hard ? propagate!(lo, hi, c.m_prev) : true
    ok || error("empty feasible box (ramp from m_prev outside the bounds)")
    (; c, L, act, K, S, F4, lo, hi, hard, m_prev = c.m_prev, y_ref = c.y_ref, h = c.h, p = c.p)
end

"""Ramp propagation on a box: tightens lo/hi so every point is reachable from ṁ_prev with
|Δu| ≤ DM (forward) and consistent with its successors (backward). Returns false if empty."""
function propagate!(lo, hi, m_prev)
    K = length(lo); tol = 1e-9
    a, b = m_prev, m_prev
    for k in 1:K
        lo[k] = max(lo[k], a - DM); hi[k] = min(hi[k], b + DM)
        lo[k] > hi[k] + tol && return false
        hi[k] = max(hi[k], lo[k]); a, b = lo[k], hi[k]
    end
    for k in K-1:-1:1
        lo[k] = max(lo[k], lo[k+1] - DM); hi[k] = min(hi[k], hi[k+1] + DM)
        lo[k] > hi[k] + tol && return false
        hi[k] = max(hi[k], lo[k])
    end
    true
end

"""ŷ(u) with the converged implicit-Euler ROM (as the NLP)."""
function predict(P, u)
    S = P.S; x = P.c.x_init; y = zeros(P.K)
    for k in 1:P.K
        x = S.step_be(x, P.L[k], u[k], P.h, P.p)
        y[k] = S.output(x, P.L[k], u[k], P.p)
    end
    y
end

hinge2(z) = z > 0 ? z^2 : 0.0
function cost(P, u, y)
    J = sum((y .- P.y_ref) .^ 2) + LAM_PCT * sum(hinge2.(y .- U_Y))
    P.hard || (J += LAM_M * sum(hinge2(abs(u[k+1] - u[k]) - DM) for k in 1:P.K-1))
    J
end
J(P, u) = cost(P, u, predict(P, u))

"""
Exact repair of a solver's point into the feasible set: interior-point solvers return points
that violate the ramp by ~1e-7 (their feasibility tolerance), and only an exactly feasible point
gives a valid UB. Forward pass: u_k ← clamp(u_k, u_{k−1} − DM, u_{k−1} + DM) ∩ box.
"""
function repair(P, u)
    v = clamp.(float.(u), P.act.m_lo, P.act.m_hi)
    P.hard || return v
    prev = P.m_prev
    for k in eachindex(v)
        v[k] = clamp(v[k], max(P.act.m_lo, prev - DM), min(P.act.m_hi, prev + DM)); prev = v[k]
    end
    v
end

"""Max violation of the feasible set (box + hard ramp)."""
function infeas(P, u)
    v = max(maximum(P.act.m_lo .- u), maximum(u .- P.act.m_hi), 0.0)
    if P.hard
        mm = [P.m_prev; u]
        v = max(v, maximum(abs.(diff(mm))) - DM)
    end
    max(v, 0.0)
end

# ------------------------------------------------------------------------------------------
# local solvers (upper bounds)
# ------------------------------------------------------------------------------------------
const KN = Ref{Any}(nothing)      # the KNITRO module, if it loads and finds a licence

"""Try to load KNITRO.jl and create a context (needs a licence). Returns true if usable."""
function knitro_init()
    try
        @eval Main using KNITRO
        mod = getfield(Main, :KNITRO)
        md = Model(() -> mod.Optimizer()); set_silent(md)
        @variable(md, x >= 1); @objective(md, Min, (x - 2)^2)
        optimize!(md)
        termination_status(md) in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL) || return false
        KN[] = mod
        return true
    catch e
        @info "Knitro not usable" exception = sprint(showerror, e)[1:min(end, 300)]
        return false
    end
end

"""One local solve of the full-space NLP (same model as N3) from `start`."""
function local_solve(P, start; solver = :ipopt, attrs = Dict{String,Any}(), tlimit = 120.0)
    c, act, F4 = P.c, P.act, P.F4
    g = clamp.(start, act.m_lo, act.m_hi)
    md = if solver == :ipopt
        m_ = Model(Ipopt.Optimizer)
        set_optimizer_attribute(m_, "tol", 1e-10); set_optimizer_attribute(m_, "max_iter", 3000)
        set_optimizer_attribute(m_, "max_wall_time", tlimit)
        m_
    else
        m_ = Model(() -> KN[].Optimizer())
        set_optimizer_attribute(m_, "outlev", 0)
        set_optimizer_attribute(m_, "feastol", 1e-10); set_optimizer_attribute(m_, "opttol", 1e-10)
        set_optimizer_attribute(m_, "maxtime", tlimit)
        m_
    end
    set_silent(md)
    for (k, v) in attrs; set_optimizer_attribute(md, k, v); end
    @variable(md, act.m_lo <= m[k = 1:c.K] <= act.m_hi, start = g[k])
    y = F4.add_scenario!(md, c, m, P.L, g)
    @objective(md, Min, F4.scenario_cost!(md, c, m, y, act.ramp))
    t = @elapsed optimize!(md)
    st = string(termination_status(md))
    u0 = has_values(md) ? value.(m) : g
    u = repair(P, u0)                       # exactly feasible; inf = how far the solver's point was
    (u = u, J = J(P, u), status = st, time = t, inf = infeas(P, u0), moved = maximum(abs.(u .- u0)))
end

"""Random feasible start: a random walk with |Δu| ≤ DM inside the propagated box."""
function random_start(P, rng)
    K = P.K; u = zeros(K)
    if P.hard
        prev = P.m_prev
        for k in 1:K
            a = max(P.lo[k], prev - DM); b = min(P.hi[k], prev + DM)
            u[k] = a + (b - a) * rand(rng); prev = u[k]
        end
    else
        for k in 1:K; u[k] = P.lo[k] + (P.hi[k] - P.lo[k]) * rand(rng); end
    end
    u
end

"""Ipopt multistart: given starts + `n_rand` random feasible walks. Returns all solves."""
function ipopt_multistart(P, starts; n_rand = 20, seed = 1, tlimit = 60.0)
    rng = MersenneTwister(seed)
    S = vcat(starts, [random_start(P, rng) for _ in 1:n_rand])
    [(; local_solve(P, s; tlimit)..., start = i) for (i, s) in enumerate(S)]
end

"""Knitro: its own multistart (ms_enable) and each algorithm from the warm start."""
function knitro_runs(P, warm; ms = 50, tlimit = 300.0)
    KN[] === nothing && return NamedTuple[]
    out = NamedTuple[]
    for alg in 1:4      # 1 IP-direct, 2 IP-CG, 3 active set, 4 SQP
        r = try local_solve(P, warm; solver = :knitro, attrs = Dict("algorithm" => alg), tlimit)
        catch e; (u = warm, J = Inf, status = "ERROR " * sprint(showerror, e)[1:min(end, 80)], time = 0.0, inf = Inf, moved = 0.0) end
        push!(out, (; r..., run = "alg$alg"))
    end
    r = try local_solve(P, warm; solver = :knitro, tlimit,
                        attrs = Dict("ms_enable" => 1, "ms_maxsolves" => ms, "ms_seed" => 1, "ms_maxbndrange" => 100.0))
    catch e; (u = warm, J = Inf, status = "ERROR " * sprint(showerror, e)[1:min(end, 80)], time = 0.0, inf = Inf, moved = 0.0) end
    push!(out, (; r..., run = "ms$ms"))
    out
end

# ------------------------------------------------------------------------------------------
# rigorous enclosure of the residual on a box (generalises f3c_rigorous.jl, part 2)
# ------------------------------------------------------------------------------------------
I_(x) = interval(x)

"""
Affine model a(u) = ŷ(uc) + G (u − uc) of the PCT predictions on the box [lo, hi] and a
rigorous bound ε_k ≥ |ŷ_k(u) − a_k(u)| for every u in the box. Returns nothing if the
remainder enclosure cannot be verified (the box is then split further).
"""
function enclose(P, lo, hi)
    S, p, h, K, c, L = P.S, P.p, P.h, P.K, P.c, P.L
    NX = S.NX; iTw, iM = S.IDX.Tw, S.IDX.dM
    uc = (lo .+ hi) ./ 2
    Sc = max.(abs.(c.x_init), 1.0)
    Xs = zeros(K + 1, NX); Xs[1, :] = c.x_init
    for k in 1:K; Xs[k+1, :] = S.step_be(Xs[k, :], L[k], uc[k], h, p); end
    Jf(z, d) = ForwardDiff.jacobian(w -> S.rhs(w, d, 0.0, p), z)
    dU = [interval(lo[j] - uc[j], hi[j] - uc[j]) for j in 1:K]
    Bs = zeros(NX); Bs[iM] = 1 / Sc[iM]
    Sinv = Diagonal(1 ./ Sc); Sd = Diagonal(Sc)
    pv = p.poly
    p4d(x) = 4pv[1] * x^3 + 3pv[2] * x^2 + 2pv[3] * x + pv[4]
    p4dd(x) = 12pv[1] * x^2 + 6pv[2] * x + 2pv[3]
    T = zeros(NX, K); r = fill(I_(0.0), NX)
    eps = zeros(K); G = zeros(K, K); ystar = zeros(K)
    for k in 1:K
        d = L[k]; xs = Xs[k+1, :]; xp = Xs[k, :]
        As = Sinv * Jf(xs, d) * Sd
        Ms = I - h * As; C = inv(Ms)
        ek = zeros(K); ek[k] = 1.0
        Tn = C * (T .+ h .* Bs * ek')
        Dk = I_.(Ms) * I_.(Tn) .- (I_.(T) .+ h .* I_.(Bs) * I_.(ek)')
        res = (I_.(xs) .- I_.(xp) .- h .* S.rhs(I_.(xs), d, I_(uc[k]), p)) ./ Sc
        E = Matrix{Float64}(I, NX, NX) .- I_.(C) * I_.(Ms)
        lin = I_.(Tn) * dU
        fixed = I_.(C) * (r .- res .- Dk * dU)
        R = fixed .+ interval(-1e-14, 1e-14)
        accepted = false
        for _ in 1:60
            δ = lin .+ R
            Z = I_.(xs) .+ Sc .* δ
            ΔA = (Sinv * I_.(Jf(Z, d)) * Sd) .- I_.(As)
            Fr = fixed .+ I_.(C) * (h .* (ΔA * δ)) .+ E * R
            if all(issubset_interval.(Fr, R))
                R = Fr; accepted = true; break
            end
            any(x -> !isfinite(sup(x)) || !isfinite(inf(x)), Fr) && return nothing
            R = [interval(min(inf(R[j]), inf(Fr[j])) - 0.1 * (sup(Fr[j]) - inf(Fr[j])) - 1e-15,
                          max(sup(R[j]), sup(Fr[j])) + 0.1 * (sup(Fr[j]) - inf(Fr[j])) + 1e-15) for j in 1:NX]
        end
        accepted || return nothing
        r = R; T = Tn
        dd = d + 5550
        e_k = Sc[iTw] * r[iTw] + 0.5 * p4dd(I_(dd) + interval(lo[k], hi[k])) * pown(dU[k], 2)
        eps[k] = mag(e_k)
        ystar[k] = S.output(xs, d, uc[k], p)
        G[k, :] = Sc[iTw] .* T[iTw, :]; G[k, k] += p4d(dd + uc[k])
    end
    all(isfinite, eps) || return nothing
    (; uc, ystar, G, eps)
end

"""Convex relaxation (QP) on the box: a valid LB of J over {u in box, feasible}."""
function relax(P, lo, hi, A)
    K = P.K
    md = Model(Clarabel.Optimizer); set_silent(md)
    set_optimizer_attribute(md, "tol_gap_abs", 1e-12); set_optimizer_attribute(md, "tol_gap_rel", 1e-10)
    set_optimizer_attribute(md, "tol_feas", 1e-10)
    @variable(md, lo[i] <= u[i = 1:K] <= hi[i])
    @variable(md, t[1:K] >= 0); @variable(md, rr[1:K] >= 0)
    @expression(md, aa[k = 1:K], A.ystar[k] + sum(A.G[k, j] * (u[j] - A.uc[j]) for j in 1:k))
    @constraint(md, [k = 1:K], t[k] >= aa[k] - P.y_ref - A.eps[k])
    @constraint(md, [k = 1:K], t[k] >= -(aa[k] - P.y_ref) - A.eps[k])
    @constraint(md, [k = 1:K], rr[k] >= aa[k] - A.eps[k] - U_Y)
    obj = @expression(md, sum(t .^ 2) + LAM_PCT * sum(rr .^ 2))
    if P.hard
        mm = [P.m_prev; u]
        @constraint(md, [k = 1:K], mm[k+1] - mm[k] <= DM)
        @constraint(md, [k = 1:K], mm[k] - mm[k+1] <= DM)
    else
        @variable(md, s[1:K-1] >= 0)
        @constraint(md, [k = 1:K-1], u[k+1] - u[k] <= DM + s[k])
        @constraint(md, [k = 1:K-1], u[k] - u[k+1] <= DM + s[k])
        obj = @expression(md, obj + LAM_M * sum(s .^ 2))
    end
    @objective(md, Min, obj)
    optimize!(md)
    st = termination_status(md)
    st in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE) && return (lb = Inf, u = nothing, status = string(st))
    st in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL) || return (lb = 0.0, u = nothing, status = string(st))
    lb = max(0.0, min(dual_objective_value(md), objective_value(md)))
    (lb = lb, u = value.(u), status = string(st))
end

# ------------------------------------------------------------------------------------------
# spatial branch and bound
# ------------------------------------------------------------------------------------------
"""
Best-first spatial B&B. Returns the certified LB, the final UB (and its u), node counts and
the history (time, LB, UB). Stops when UB − LB ≤ max(tol_abs, tol_rel·UB), on the node limit
or on the time limit (the LB returned is valid in every case).
"""
function bnb(P, UB, ubest; tol_abs = 1e-6, tol_rel = 1e-4, tlimit = 900.0, maxnodes = 20000,
             polish_every = 50, verbose = false)
    t0 = time()
    hist = Tuple{Float64,Float64,Float64}[]
    gapok(lb, ub) = ub - lb <= max(tol_abs, tol_rel * ub)
    if gapok(0.0, UB)          # LB = 0 already certifies (J ≥ 0)
        return (LB = 0.0, UB = UB, u = ubest, nodes = 0, open = 0, failed = 0, status = "trivial (J≥0)",
                hist = [(0.0, 0.0, UB)], root_lb = 0.0, time = 0.0)
    end
    nodes = Tuple{Float64,Vector{Float64},Vector{Float64}}[]   # (lb, lo, hi)
    failed = 0; nexp = 0
    function evalnode(lo, hi)
        A = enclose(P, lo, hi)
        A === nothing && return (lb = 0.0, u = nothing, ok = false)       # parent bound kept by caller
        R = relax(P, lo, hi, A)
        (lb = R.lb, u = R.u, ok = true)
    end
    function try_ub!(u)
        u === nothing && return
        u = repair(P, u)
        Ju = J(P, u)
        if Ju < UB && infeas(P, u) <= 1e-9
            UB = Ju; ubest = copy(u)
        end
    end
    r0 = evalnode(P.lo, P.hi)
    root_lb = r0.ok ? r0.lb : 0.0
    r0.ok || (failed += 1)
    try_ub!(r0.u)
    push!(nodes, (root_lb, copy(P.lo), copy(P.hi)))
    push!(hist, (time() - t0, root_lb, UB))
    status = "limit"
    while !isempty(nodes)
        i = argmin(first.(nodes))
        LBg = nodes[i][1]
        if gapok(LBg, UB)
            status = "certified"; break
        end
        if time() - t0 > tlimit || nexp >= maxnodes
            status = time() - t0 > tlimit ? "time limit" : "node limit"; break
        end
        lbp, lo, hi = nodes[i]; deleteat!(nodes, i); nexp += 1
        # branching variable: widest box side weighted by its curvature reach (later steps
        # accumulate more residual): w_j = (hi_j − lo_j) · (K − j + 1)
        w = [(hi[j] - lo[j]) * (P.K - j + 1) for j in 1:P.K]
        j = argmax(w)
        w[j] <= 1e-9 && continue           # degenerate box: nothing left to split
        mid = (lo[j] + hi[j]) / 2
        for side in 1:2
            l2, h2 = copy(lo), copy(hi)
            side == 1 ? (h2[j] = mid) : (l2[j] = mid)
            P.hard && (propagate!(l2, h2, P.m_prev) || continue)
            r = evalnode(l2, h2)
            r.ok || (failed += 1)
            lb = r.ok ? max(r.lb, lbp) : lbp     # a child's bound can only improve the parent's
            try_ub!(r.u)
            lb < Inf && push!(nodes, (lb, l2, h2))
        end
        if polish_every > 0 && nexp % polish_every == 0 && !isempty(nodes)
            k = argmin(first.(nodes))
            s = local_solve(P, (nodes[k][2] .+ nodes[k][3]) ./ 2; tlimit = 30.0)
            s.inf <= 1e-5 && try_ub!(s.u)
        end
        if nexp % 25 == 0
            lbg = isempty(nodes) ? UB : minimum(first.(nodes))
            push!(hist, (time() - t0, lbg, UB))
            verbose && @printf("  nodes %5d  open %5d  LB %.10g  UB %.10g  gap %.3e\n", nexp, length(nodes), lbg, UB, UB - lbg)
        end
    end
    LB = isempty(nodes) ? UB : min(minimum(first.(nodes)), UB)
    isempty(nodes) && (status = "certified (tree exhausted)")
    push!(hist, (time() - t0, LB, UB))
    (LB = LB, UB = UB, u = ubest, nodes = nexp, open = length(nodes), failed, status, hist, root_lb, time = time() - t0)
end

end # module
