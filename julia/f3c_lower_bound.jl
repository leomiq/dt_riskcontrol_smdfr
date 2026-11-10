# F3c, step 4a (docs/monolithic_formulation_plan.md §8.5): SEMI-RIGOROUS lower bound on the
# global optimum of the nonlinear P0 (N3), hence an optimality gap UB − LB.
#
# 4a.1 Affine reference: a(u) = ŷ(u*) + G (u − u*), G = ∂ŷ/∂u at the NLP solution u* (exact
#      AD through the nonlinear simulation). Any affine a gives a valid bound; this one is tight
#      near u*.
# 4a.2 Residual bound: e(u) = ŷ(u) − a(u); ε_k ≥ max_{u∈box} |e_k(u)| estimated from random
#      samples, extreme patterns and adversarial maximization of ±e_k (Ipopt, exact gradient),
#      times a safety factor. SEMI-rigorous: sampling can miss the maximum (step 4b replaces
#      this by interval arithmetic).
# 4a.3 Convex relaxation: |e_k| ≤ ε_k ⇒ (ŷ_k − y_ref)² ≥ (|a_k − y_ref| − ε_k)_+² and
#      [ŷ_k − U_y]_+² ≥ [a_k − ε_k − U_y]_+²; the ramp penalty depends on u only (exact).
#      L(u) ≤ J(u) on the whole box and L is convex: LB = min L is a QP (Clarabel, dual bound).
# 4a.4 Gap (UB − LB)/UB, for several safety factors.
# 4a.5 Diagnostic: ε and the bound on sub-boxes of half-width ρ around u*, to estimate how much
#      box subdivision (branch and bound) a certified gap would need.
#
# Usage, from the repository root:  julia --project=julia julia/f3c_lower_bound.jl [t_q] [n_samples]
include(joinpath(@__DIR__, "F4Lib.jl"))
using .F4Lib
include(joinpath(@__DIR__, "ReducedP0.jl"))
using JuMP, Ipopt, Clarabel, ForwardDiff, LinearAlgebra, Printf, Random, Statistics

TQ = length(ARGS) >= 1 ? ARGS[1] : "1180"
NS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3000
c = load_query(TQ)
P = reduced_p0(c)
K, m_lo, m_hi, y_ref, U_y = P.K, P.m_lo, P.m_hi, P.y_ref, P.U_y
outdir = joinpath(@__DIR__, "..", "results", "f3c_optimality"); mkpath(outdir)
log_io = open(joinpath(outdir, "lower_bound_$TQ.txt"), "w")
say(s...) = (println(s...); println(log_io, s...); flush(log_io))
fmt(x) = @sprintf("%.3e", x)

say("=== F3c step 4a, t_q = $TQ s, box ṁ ∈ [$(round(m_lo; digits = 4)), $(round(m_hi; digits = 4))]^$K kg/s")

# ---- UB: the monolithic NLP solution, costed on the reduced BE model
nlp = solve_det(c, c.used_leak); nlp = solve_det(c, c.used_leak)
ustar = clamp.(nlp.m, m_lo, m_hi)
UB = P.J(ustar)
say(@sprintf("UB = J(u*) = %.10f (monolithic NLP, Ipopt, %s)", UB, nlp.status))

# ---- 4a.1 affine reference
ystar = P.predict(ustar)
G = ForwardDiff.jacobian(u -> P.predict(u), ustar)
a(u) = ystar .+ G * (u .- ustar)
e(u) = P.predict(u) .- a(u)
say(@sprintf("\n4a.1 G = ∂ŷ/∂u at u*: lower-triangular %s, |G| ∈ [%s, %s] K/(kg/s)",
             istril(G, 0) || maximum(abs.(triu(G, 1))) < 1e-12 ? "yes" : "no",
             fmt(minimum(abs.(G[G .!= 0]))), fmt(maximum(abs.(G)))))

# ---- 4a.2 residual bound over the whole box
say("\n4a.2 max_{u∈box} |e_k(u)|, e = ŷ − a")
rng = MersenneTwister(4)
emax = zeros(K); argbest = [[copy(ustar), copy(ustar)] for _ in 1:K]   # (for +e, for −e)
best_pos = fill(-Inf, K); best_neg = fill(-Inf, K)
function record!(u)
    ev = e(u)
    for k in 1:K
        ev[k] > best_pos[k] && (best_pos[k] = ev[k]; argbest[k][1] = copy(u))
        -ev[k] > best_neg[k] && (best_neg[k] = -ev[k]; argbest[k][2] = copy(u))
    end
end
structured = [[fill(l, K) for l in range(m_lo, m_hi; length = 11)];
              [[k <= s ? m_hi : m_lo for k in 1:K] for s in 1:K-1];
              [[k <= s ? m_lo : m_hi for k in 1:K] for s in 1:K-1];
              [[isodd(k) ? m_hi : m_lo for k in 1:K], [isodd(k) ? m_lo : m_hi for k in 1:K]]]
ts = @elapsed begin
    foreach(record!, structured)
    for _ in 1:NS; record!(m_lo .+ (m_hi - m_lo) .* rand(rng, K)); end
end
eps_sample = max.(best_pos, best_neg)
say(@sprintf("  sampling: %d random + %d extreme patterns in %.1f s; max_k ε_k = %s K (step %d)",
             NS, length(structured), ts, fmt(maximum(eps_sample)), argmax(eps_sample)))

"""Maximize s·e_k(u) over the box (adversarial), exact gradient, L-BFGS Hessian."""
function adversary(k, s, start)
    f(u) = s * e(u)[k]
    md = Model(Ipopt.Optimizer); set_silent(md)
    set_optimizer_attribute(md, "hessian_approximation", "limited-memory")
    set_optimizer_attribute(md, "tol", 1e-9); set_optimizer_attribute(md, "max_iter", 300)
    fop(v...) = f(collect(v))
    ∇f!(gr, v...) = (gr .= ForwardDiff.gradient(f, collect(v)); nothing)
    @operator(md, opf, K, fop, ∇f!)
    @variable(md, m_lo <= u[i = 1:K] <= m_hi, start = start[i])
    @objective(md, Max, opf(u...))
    optimize!(md)
    uu = clamp.(value.(u), m_lo, m_hi)
    (val = f(uu), u = uu)
end
eps_adv = zeros(K)
ta = @elapsed for k in 1:K
    vp = adversary(k, +1, argbest[k][1]); vn = adversary(k, -1, argbest[k][2])
    eps_adv[k] = max(vp.val, vn.val, 0.0)
end
eps_hat = max.(eps_sample, eps_adv)
say(@sprintf("  adversarial (2 Ipopt solves per k) in %.1f s; max_k ε_k = %s K; adversary raised ε above sampling in %d of %d steps (max ratio %.2f)",
             ta, fmt(maximum(eps_hat)), count(eps_adv .> eps_sample .* (1 + 1e-6)), K, maximum(eps_hat ./ max.(eps_sample, 1e-15))))
say("  ε_k [K] by step: " * join(fmt.(eps_hat), " "))

# ---- 4a.3 convex relaxation (QP) on a box
"""LB = min L(u) over lo ≤ u ≤ hi, with residual bound ε (vector). Returns primal, dual, u."""
function relaxation(ε; lo = fill(m_lo, K), hi = fill(m_hi, K))
    md = Model(Clarabel.Optimizer); set_silent(md)
    set_optimizer_attribute(md, "tol_gap_abs", 1e-12); set_optimizer_attribute(md, "tol_gap_rel", 1e-12)
    @variable(md, lo[i] <= u[i = 1:K] <= hi[i])
    @variable(md, t[1:K] >= 0); @variable(md, s[1:K-1] >= 0); @variable(md, r[1:K] >= 0)
    @expression(md, aa[k = 1:K], ystar[k] + sum(G[k, j] * (u[j] - ustar[j]) for j in 1:k))
    @constraint(md, [k = 1:K], t[k] >= aa[k] - y_ref - ε[k])
    @constraint(md, [k = 1:K], t[k] >= -(aa[k] - y_ref) - ε[k])
    @constraint(md, [k = 1:K], r[k] >= aa[k] - ε[k] - U_y)
    @constraint(md, [k = 1:K-1], u[k+1] - u[k] <= P.dm_lim + s[k])
    @constraint(md, [k = 1:K-1], u[k] - u[k+1] <= P.dm_lim + s[k])
    @objective(md, Min, sum(t .^ 2) + P.lam_m * sum(s .^ 2) + P.lam_pct * sum(r .^ 2))
    optimize!(md)
    (primal = objective_value(md), dual = dual_objective_value(md), u = value.(u),
     status = termination_status(md))
end

# ---- 4a.4 gap
say("\n4a.3–4a.4 convex relaxation (Clarabel) and gap; LB = dual bound of the relaxation")
say("  safety  LB              gap = UB − LB    gap / UB")
rows = Any[]
for sf in (0.0, 1.0, 1.5, 2.0)
    R = relaxation(sf .* eps_hat)
    LB = min(R.dual, R.primal)
    say(@sprintf("  %4.1f    %.10f    %.4e       %7.3f %%%s", sf, LB, UB - LB, 100 * (UB - LB) / UB,
                 sf == 0 ? "   (no residual bound: the linear model alone, NOT a valid bound)" : ""))
    push!(rows, (sf, LB, UB - LB, (UB - LB) / UB))
end

# ---- 4a.5 diagnostic on sub-boxes around u*
say("\n4a.5 sub-boxes around u* (half-width ρ·(ṁ_max − ṁ_min)); ε by sampling only, safety 1.5")
say("  ρ       max ε [K]    LB(sub-box)      gap/UB (sub-box)")
diag = Any[]
for ρ in (1.0, 0.5, 0.25, 0.1, 0.05, 0.02)
    lo = max.(ustar .- ρ * (m_hi - m_lo), m_lo); hi = min.(ustar .+ ρ * (m_hi - m_lo), m_hi)
    em = zeros(K)
    for _ in 1:max(300, NS ÷ 5)
        em .= max.(em, abs.(e(lo .+ (hi .- lo) .* rand(rng, K))))
    end
    for v in ([lo], [hi]); em .= max.(em, abs.(e(v[1]))); end
    R = relaxation(1.5 .* em; lo = lo, hi = hi)
    LB = min(R.dual, R.primal)
    say(@sprintf("  %-6.2f  %s    %.10f     %7.3f %%", ρ, fmt(maximum(em)), LB, 100 * (UB - LB) / UB))
    push!(diag, (ρ, maximum(em), LB, (UB - LB) / UB))
end
say("  (a sub-box bound is only valid on that sub-box: certifying the whole box needs covering it,")
say("   i.e. branch and bound; the ρ at which the gap becomes small sets the subdivision depth)")

open(joinpath(outdir, "lower_bound_eps_$TQ.csv"), "w") do io
    println(io, "k,eps_sample,eps_adversarial,eps_hat")
    for k in 1:K; println(io, join((k, eps_sample[k], eps_adv[k], eps_hat[k]), ",")); end
end
open(joinpath(outdir, "lower_bound_gap_$TQ.csv"), "w") do io
    println(io, "safety_factor,LB,gap_abs,gap_rel"); for r in rows; println(io, join(r, ",")); end
end
open(joinpath(outdir, "lower_bound_subbox_$TQ.csv"), "w") do io
    println(io, "rho,max_eps,LB,gap_rel"); for r in diag; println(io, join(r, ",")); end
end
close(log_io)
println("\nwritten to ", outdir)
