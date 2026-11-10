# F3c, steps 1-3 (docs/monolithic_formulation_plan.md §8.5): optimality evidence for the
# NONLINEAR monolithic P0 (N3), independent of the linearized QP.
#
# Reduced space: the states are a function of the controls through the simulation, so
#   J(u) = Σ_k (ŷ_k(u) − y_ref)² + λ_m Σ_k [|u_{k+1}−u_k| − Δṁ]_+² + λ_pct Σ_k [ŷ_k(u) − U_y]_+²,
#   u ∈ [ṁ_min, ṁ_max]^K,
# with ŷ(u) from the replica under one of two discretizations:
#   :be      converged implicit Euler (the model of the monolithic NLP, F2)
#   :ode1be  one Newton iteration per step (exactly what Simulink computes)
# J is written generically, so ForwardDiff differentiates through the simulation
# (gradients and exact Hessians).
#
# Step 1, local certificate: at the NLP solution, first-order conditions with the bound
#   multipliers (μ_k = −∂J/∂u_k ≥ 0 at the upper bound, ≥ 0 at the lower bound with the
#   opposite sign), strict complementarity, and second-order sufficiency: the Hessian
#   restricted to the free variables is positive definite.
# Step 2, P0 on ode1be: reduced-space NLP (Ipopt with exact gradient and Hessian through the
#   ode1be simulation), compared with the BE optimum.
# Step 3, global evidence: (a) multistart of the monolithic NLP from random, structured and
#   LJ/GA starts; (b) differential evolution on J(u), derivative-free; (c) eigenvalues of
#   ∇²J(u) sampled over the control box.
#
# Usage, from the repository root:  julia --project=julia julia/f3c_optimality.jl [t_q] [n_starts] [de_gens]
include(joinpath(@__DIR__, "F4Lib.jl"))
using .F4Lib, .F4Lib.SMDFRROM, JuMP, Ipopt, ForwardDiff, LinearAlgebra, Printf, Random, Statistics

TQ = length(ARGS) >= 1 ? ARGS[1] : "1180"
N_STARTS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 500
DE_GENS = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1500
c = load_query(TQ)
p, h, K, x0 = c.p, c.h, c.K, c.x_init
leak = c.used_leak
m_lo, m_hi = c.m_lo_arch, c.m_hi_arch
y_ref = c.y_ref; U_y = 1244.0; lam_m = 1e-4; dm_lim = 5.0; lam_pct = 10.0
outdir = joinpath(@__DIR__, "..", "results", "f3c_optimality"); mkpath(outdir)
log_io = open(joinpath(outdir, "summary_$TQ.txt"), "w")
say(s...) = (println(s...); println(log_io, s...); flush(log_io))
fmt(x) = @sprintf("%.3e", x)

# ---------------------------------------------------------------- reduced-space cost
primal(v) = v isa ForwardDiff.Dual ? primal(ForwardDiff.value(v)) : v

"""Converged BE step, differentiable: Newton on primal values, then dual Newton refinements."""
function step_be_ad(x, d, u)
    z = SMDFRROM.step_be(primal.(x), d, primal(u), h, p)      # Float64, converged
    T = promote_type(eltype(x), typeof(u))
    zT = T.(z)
    for _ in 1:3                                              # exact 1st and 2nd derivatives
        r = zT .- x .- h .* rhs(zT, d, u, p)
        Jz = I - h * ForwardDiff.jacobian(w -> rhs(w, d, u, p), zT)
        zT = zT .- Jz \ r
    end
    zT
end
step_ode1be_ad(x, d, u) = SMDFRROM.step_ode1be(x, d, u, h, p)

function predict(u::AbstractVector{T}, scheme) where {T}
    step = scheme === :be ? step_be_ad : step_ode1be_ad
    x = T.(x0); y = Vector{T}(undef, K)
    for k in 1:K
        x = step(x, leak[k], u[k])
        y[k] = SMDFRROM.output(x, leak[k], u[k], p)
    end
    y
end
hinge2(z) = z > 0 ? z^2 : zero(z)
function Jred(u, scheme)
    y = predict(u, scheme)
    sum((y .- y_ref) .^ 2) + lam_m * sum(hinge2(abs(u[k+1] - u[k]) - dm_lim) for k in 1:K-1) +
        lam_pct * sum(hinge2.(y .- U_y))
end
Jbe(u) = Jred(u, :be); Jode(u) = Jred(u, :ode1be)

"""
Hessian of J. :ad = forward-over-forward (exact; three nested dual levels through the
implicit step, expensive to compile); :fd = central differences of the exact AD gradient,
H[:, i] = (∇J(u + ε e_i) − ∇J(u − ε e_i)) / 2ε, error O(ε²) plus gradient round-off / ε.
"""
const HESS = get(ENV, "F3C_HESS", "fd") == "ad" ? :ad : :fd
function hess(J, u; ε = 1e-4)
    HESS === :ad && return ForwardDiff.hessian(J, u)
    n = length(u); H = zeros(n, n)
    for i in 1:n
        e = zeros(n); e[i] = ε
        H[:, i] = (ForwardDiff.gradient(J, u .+ e) .- ForwardDiff.gradient(J, u .- e)) ./ (2ε)
    end
    (H + H') / 2
end

# ---------------------------------------------------------------- step 1: local certificate
"""First- and second-order conditions of the box-constrained reduced problem at u."""
function local_certificate(J, u; name, tol_act = 1e-6)
    g = ForwardDiff.gradient(J, u)
    H = hess(J, u); H = (H + H') / 2
    at_hi = u .>= m_hi - tol_act; at_lo = u .<= m_lo + tol_act
    free = .!(at_hi .| at_lo)
    mu_hi = -g[at_hi]; mu_lo = g[at_lo]                       # multipliers of active bounds
    gfree = isempty(g[free]) ? 0.0 : maximum(abs.(g[free]))
    mu_min = isempty([mu_hi; mu_lo]) ? Inf : minimum([mu_hi; mu_lo])
    Hf = H[free, free]
    ev = isempty(Hf) ? Float64[] : eigvals(Symmetric(Hf))
    lmin = isempty(ev) ? Inf : minimum(ev); lmax = isempty(ev) ? NaN : maximum(ev)
    say(@sprintf("  [%s] J = %.10f | active: %d upper, %d lower, %d free", name, J(u), count(at_hi), count(at_lo), count(free)))
    say(@sprintf("  [%s] stationarity on free vars max|∂J/∂u| = %s; min bound multiplier = %s (must be > 0)",
                 name, fmt(gfree), fmt(mu_min)))
    say(@sprintf("  [%s] reduced Hessian on free vars: λ_min = %s, λ_max = %s, cond = %s",
                 name, fmt(lmin), fmt(lmax), fmt(lmax / lmin)))
    ok = (mu_min > 0) && (lmin > 0)
    # first-order decrease that the gradient residual could still allow: |g_free|² / (2 λ_min)
    slack = isempty(ev) ? 0.0 : gfree^2 * count(free) / (2 * lmin)
    say(@sprintf("  [%s] strict local minimum: %s; quadratic-model bound on remaining local improvement ≈ %s",
                 name, ok ? "YES (KKT + strict complementarity + SOSC)" : "NO", fmt(slack)))
    (; g, H, free, mu_min, gfree, lmin, lmax, ok, slack)
end

say("=== F3c, t_q = $TQ s; ṁ ∈ [$(round(m_lo; digits = 4)), $(round(m_hi; digits = 4))] kg/s, K = $K")
nlp = solve_det(c, leak); nlp = solve_det(c, leak)            # monolithic NLP (BE), as F2
u_be = clamp.(nlp.m, m_lo, m_hi)
say(@sprintf("Monolithic NLP (BE, Ipopt): %s, J_BE = %.10f, %.3f s", nlp.status, nlp.obj, nlp.time))
say(@sprintf("Consistency: reduced-space J_BE(u_NLP) = %.10f (must match the NLP objective)", Jbe(u_be)))

say("\n--- Step 1: local optimality certificate of the NLP solution (BE model)")
t1 = @elapsed cert_be = local_certificate(Jbe, u_be; name = "BE")
say(@sprintf("  (%.1f s)", t1))

# ---------------------------------------------------------------- step 2: P0 on ode1be
say("\n--- Step 2: P0 on the exact Simulink discretization (ode1be), reduced-space NLP")
function solve_reduced(J; start)
    md = Model(Ipopt.Optimizer); set_silent(md)
    # exact gradient through the simulation; the Hessian is not needed to find the optimum,
    # Ipopt approximates it with L-BFGS (the Hessian is used only to certify, in step 1)
    set_optimizer_attribute(md, "tol", 1e-10); set_optimizer_attribute(md, "max_iter", 1000)
    set_optimizer_attribute(md, "hessian_approximation", "limited-memory")
    f(v...) = J(collect(v))
    ∇f!(gr, v...) = (gr .= ForwardDiff.gradient(J, collect(v)); nothing)
    @operator(md, opJ, K, f, ∇f!)
    @variable(md, m_lo <= u[k = 1:K] <= m_hi, start = start[k])
    @objective(md, Min, opJ(u...))
    t = @elapsed optimize!(md)
    (u = clamp.(value.(u), m_lo, m_hi), obj = objective_value(md), status = termination_status(md),
     time = t, iters = barrier_iterations(md))
end
r_ode = solve_reduced(Jode; start = u_be)
u_ode = r_ode.u
say(@sprintf("  ode1be NLP (start = BE optimum): %s, %d it, %.2f s, J_ode1be = %.10f", r_ode.status, r_ode.iters, r_ode.time, r_ode.obj))
r_ode2 = solve_reduced(Jode; start = fill((m_lo + m_hi) / 2, K))
say(@sprintf("  ode1be NLP (start = mid-range):  %s, %d it, %.2f s, J_ode1be = %.10f, max|Δu| to previous = %s",
             r_ode2.status, r_ode2.iters, r_ode2.time, r_ode2.obj, fmt(maximum(abs.(r_ode2.u .- u_ode)))))
say(@sprintf("  J_ode1be(u_BE*) = %.10f vs J_ode1be(u_ode1be*) = %.10f → BE optimum is %.2e (%.4f %%) above the ode1be optimum",
             Jode(u_be), r_ode.obj, Jode(u_be) - r_ode.obj, 100 * (Jode(u_be) - r_ode.obj) / r_ode.obj))
say(@sprintf("  max_k |u_BE − u_ode1be| = %s kg/s; Δu₀ = %s kg/s", fmt(maximum(abs.(u_be .- u_ode))), fmt(u_be[1] - u_ode[1])))
cert_ode = local_certificate(Jode, u_ode; name = "ode1be")

# ---------------------------------------------------------------- step 3a: multistart
say("\n--- Step 3a: multistart of the monolithic NLP (BE), $(N_STARTS) random + structured + LJ/GA starts")
rng = MersenneTwister(2026)
starts = Pair{String,Vector{Float64}}[]
for i in 1:N_STARTS; push!(starts, "uniform" => m_lo .+ (m_hi - m_lo) .* rand(rng, K)); end
for i in 1:50                                                 # smooth random walks
    w = cumsum(randn(rng, K)) .* (m_hi - m_lo) / 6 .+ (m_lo + m_hi) / 2
    push!(starts, "walk" => clamp.(w, m_lo, m_hi))
end
for lev in range(m_lo, m_hi; length = 11); push!(starts, "constant" => fill(lev, K)); end
for s in 1:K-1; push!(starts, "bang-bang" => [k <= s ? m_hi : m_lo for k in 1:K]); end
push!(starts, "alternating" => [isodd(k) ? m_hi : m_lo for k in 1:K])
e2 = joinpath(@__DIR__, "..", "results", "e2_optimizers_$TQ", "e2_actions.csv")
if isfile(e2)
    lines = readlines(e2); hdr = split(lines[1], ",")
    vals = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    for (j, name) in enumerate(hdr)
        name in ("t", "NLP") && continue
        push!(starts, name => clamp.(current_to_massflow.([v[j] for v in vals]), m_lo, m_hi))
    end
end
res = Vector{Any}(undef, length(starts))
t3a = @elapsed for (i, (kind, g)) in enumerate(starts)
    s = solve_det(c, leak; guess = g)
    res[i] = (kind = kind, status = string(s.status), J = Jbe(clamp.(s.m, m_lo, m_hi)), u = clamp.(s.m, m_lo, m_hi), time = s.time)
end
ok = [r for r in res if r.status in ("LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED")]
Js = [r.J for r in ok]
# cluster distinct local optima by control distance
clusters = Vector{Any}()
for r in sort(ok; by = r -> r.J)
    j = findfirst(cl -> maximum(abs.(cl.u .- r.u)) < 1e-3, clusters)
    j === nothing ? push!(clusters, (u = r.u, J = r.J, n = 1, kinds = Set([r.kind]))) :
                    (clusters[j] = (u = clusters[j].u, J = clusters[j].J, n = clusters[j].n + 1, kinds = push!(clusters[j].kinds, r.kind)))
end
say(@sprintf("  %d starts in %.1f s; converged %d; failed %d", length(starts), t3a, length(ok), length(starts) - length(ok)))
say(@sprintf("  best J = %.10f, worst converged J = %.10f, NLP J = %.10f", minimum(Js), maximum(Js), Jbe(u_be)))
say(@sprintf("  distinct local optima (max|Δu| > 1e-3 kg/s): %d", length(clusters)))
for (i, cl) in enumerate(clusters[1:min(end, 8)])
    say(@sprintf("    #%d: J = %.10f, found %d times (%s), max|u − u_NLP| = %s", i, cl.J, cl.n,
                 join(sort(collect(cl.kinds)), ","), fmt(maximum(abs.(cl.u .- u_be)))))
end
improve = minimum(Js) < Jbe(u_be) - 1e-8
say("  any start better than the NLP by > 1e-8: ", improve ? "YES" : "no")
open(joinpath(outdir, "multistart_$TQ.csv"), "w") do io
    println(io, "start_kind,status,J_BE,max_du_vs_NLP,time_s")
    for r in res; println(io, join((r.kind, r.status, r.J, maximum(abs.(r.u .- u_be)), r.time), ",")); end
end

# ---------------------------------------------------------------- step 3b: differential evolution
say("\n--- Step 3b: differential evolution on J_BE(u) (derivative-free, DE/rand/1/bin)")
function de(J; NP = 60, G = DE_GENS, F = 0.7, CR = 0.9, rng = MersenneTwister(7))
    P = [m_lo .+ (m_hi - m_lo) .* rand(rng, K) for _ in 1:NP]
    fP = map(J, P); hist = Float64[]
    for gen in 1:G
        for i in 1:NP
            a, b, d = P[rand(rng, setdiff(1:NP, i), 3)]
            v = clamp.(a .+ F .* (b .- d), m_lo, m_hi)
            jr = rand(rng, 1:K)
            trial = [rand(rng) < CR || k == jr ? v[k] : P[i][k] for k in 1:K]
            ft = J(trial)
            ft <= fP[i] && (P[i] = trial; fP[i] = ft)
        end
        push!(hist, minimum(fP))
    end
    ib = argmin(fP)
    (u = P[ib], J = fP[ib], hist = hist, evals = NP * (G + 1))
end
Jbe_f(u) = Jbe(u)                                             # Float64 path
t3b = @elapsed rde = de(Jbe_f)
say(@sprintf("  %d evaluations in %.1f s: best J = %.10f (NLP: %.10f; DE − NLP = %+.3e, %+.3f %%)",
             rde.evals, t3b, rde.J, Jbe(u_be), rde.J - Jbe(u_be), 100 * (rde.J - Jbe(u_be)) / Jbe(u_be)))
pol = solve_det(c, leak; guess = rde.u)
say(@sprintf("  NLP polished from the DE best: J = %.10f, max|u − u_NLP| = %s", Jbe(clamp.(pol.m, m_lo, m_hi)),
             fmt(maximum(abs.(clamp.(pol.m, m_lo, m_hi) .- u_be)))))
open(joinpath(outdir, "de_history_$TQ.csv"), "w") do io
    println(io, "generation,best_J"); for (g, v) in enumerate(rde.hist); println(io, g, ",", v); end
end

# ---------------------------------------------------------------- step 3c: Hessian sampling
NH = parse(Int, get(ENV, "F3C_NH", "300"))
say("\n--- Step 3c: eigenvalues of ∇²J(u) at $NH points of the control box (ode1be model)")
rngh = MersenneTwister(11)
pts = [[m_lo .+ (m_hi - m_lo) .* rand(rngh, K) for _ in 1:NH-20];
       [clamp.(u_be .+ 0.5 .* (m_hi - m_lo) .* (rand(rngh, K) .- 0.5) .* s, m_lo, m_hi) for s in range(0.05, 1; length = 20)]]
lmins = Float64[]; lmaxs = Float64[]; nneg = 0
t3c = @elapsed for u in pts
    Hs = hess(Jode, u); ev = eigvals(Symmetric((Hs + Hs') / 2))
    push!(lmins, minimum(ev)); push!(lmaxs, maximum(ev))
end
nneg = count(<(0), lmins)
say(@sprintf("  %.1f s. λ_min over samples: min %s, median %s, max %s; λ_max median %s",
             t3c, fmt(minimum(lmins)), fmt(median(lmins)), fmt(maximum(lmins)), fmt(median(lmaxs))))
say(@sprintf("  points with a negative eigenvalue: %d of %d (most negative %s, relative to λ_max %s)",
             nneg, NH, fmt(minimum(lmins)), fmt(minimum(lmins) / median(lmaxs))))
open(joinpath(outdir, "hessian_samples_$TQ.csv"), "w") do io
    println(io, "dist_to_NLP_maxnorm,lambda_min,lambda_max")
    for (u, a, b) in zip(pts, lmins, lmaxs); println(io, maximum(abs.(u .- u_be)), ",", a, ",", b); end
end

# ---------------------------------------------------------------- outputs
open(joinpath(outdir, "solutions_$TQ.csv"), "w") do io
    println(io, "t,u_NLP_BE,u_NLP_ode1be,u_DE_best,grad_BE")
    for k in 1:K
        println(io, join((c.tX[c.hor[k]], u_be[k], u_ode[k], rde.u[k], cert_be.g[k]), ","))
    end
end
close(log_io)
println("\nwritten to ", outdir)
