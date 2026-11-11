# Optimality audit of the closed-loop MPC solves (N1, N2, N3).
#
# Along an N3+bias closed loop, at every query the same problem data (c, L) is solved by
#   N3      as in the loop (Ipopt, tol 1e-10, warm start from the shifted previous plan);
#   N3pol   Ipopt again from the N3 solution with the "acceptable" exit disabled
#           (acceptable_iter = 0), so it only stops at tol = 1e-10 or max_iter;
#   N3cold  Ipopt from the default cold start (another local optimum?);
#   N1, N2  the condensed QPs (Clarabel).
# All plans are scored with the same function: the P0 cost on the ode1be ROM (F4Lib.evaluate).
# The gaps are relative to the best of the N3 solutions.
#
# Usage: julia --heap-size-hint=1500M --project=julia julia/f5a_optimality_audit.jl [n=121]
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop, JuMP, Ipopt, Random, Printf, Statistics
const F4 = ClosedLoop.F4Lib

NQ = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 121
OUT = joinpath(@__DIR__, "..", "results", "f5a_closed_loop"); mkpath(OUT)

function solve_n3(c, L, act, guess; acceptable = true)
    g = clamp.(guess, act.m_lo, act.m_hi)
    md = Model(Ipopt.Optimizer); set_silent(md)
    set_optimizer_attribute(md, "tol", 1e-10); set_optimizer_attribute(md, "max_iter", 3000)
    acceptable || set_optimizer_attribute(md, "acceptable_iter", 0)
    @variable(md, act.m_lo <= m[k = 1:c.K] <= act.m_hi, start = g[k])
    y = F4.add_scenario!(md, c, m, L, g)
    @objective(md, Min, F4.scenario_cost!(md, c, m, y, act.ramp))
    optimize!(md)
    (m = value.(m), status = string(termination_status(md)), iters = barrier_iterations(md))
end

cold_guess(c, act) = act.ramp == :hard ? c.m_prev .+ 5.0 .* (1:c.K) : fill(act.m_hi, c.K)

function realisation_bench(seed)
    rng = MersenneTwister(seed)
    A = 0.2 + 0.2rand(rng); Bp = 9 + 6rand(rng); d = 8 + 2rand(rng)
    eta = 0.8 + 0.2rand(rng); fs = 5 + 15rand(rng); fl = 5 + 10rand(rng); fb = -0.25 + 0.3rand(rng)
    B = setup(; pH = ClosedLoop.HRMQueries.hrm_params(; A, B = Bp, mdot_c_drift = d))
    pert = isodd(seed) ? (eta = 1.0, fail = nothing, fbias = 0.0) :
           (eta = eta, fail = (B.t_start + fs, B.t_start + fs + fl), fbias = fb)
    B, pert
end

cases = Any[("archived", setup(), (eta = 1.0, fail = nothing, fbias = 0.0))]
for s in (1, 2)
    B, pert = realisation_bench(s); push!(cases, ("seed$s", B, pert))
end

rows = NamedTuple[]
for (name, B, pert) in cases, act in (B.act_arch, ACT_PHYS)
    k = Ref(0)
    function audit(c, L, act, prev, dec)
        k[] += 1
        J(m) = F4.evaluate(c, m, L).J
        pol = solve_n3(c, L, act, dec.m; acceptable = false)
        cold = solve_n3(c, L, act, cold_guess(c, act))
        n1 = ClosedLoop.decide(:N1, c, L, act, prev)
        n2 = ClosedLoop.decide(:N2, c, L, act, prev)
        Js = (N3 = J(dec.m), N3pol = J(pol.m), N3cold = J(cold.m), N1 = J(n1.m), N2 = J(n2.m))
        Jbest = min(Js.N3, Js.N3pol, Js.N3cold)
        push!(rows, (; case = name, act = act.name, k = k[], st3 = dec.status, st_pol = pol.status,
                     it_pol = pol.iters, st_cold = cold.status, st1 = n1.status, st2 = n2.status,
                     Jbest, gap3 = Js.N3 / Jbest - 1, gap_pol = Js.N3pol / Jbest - 1,
                     gap_cold = Js.N3cold / Jbest - 1, gap1 = Js.N1 / Jbest - 1, gap2 = Js.N2 / Jbest - 1,
                     dm0_pol = abs(dec.m[1] - pol.m[1]), dm0_cold = abs(dec.m[1] - cold.m[1]),
                     dm0_1 = abs(dec.m[1] - n1.m[1]), dm0_2 = abs(dec.m[1] - n2.m[1]),
                     dm_pol = maximum(abs.(dec.m .- pol.m))))
    end
    t = @elapsed run_loop(B, :N3; n = NQ, act, bias = true, pert..., audit)
    @printf("%-9s %-5s %d queries audited in %.0f s\n", name, act.name, k[], t)
end

open(joinpath(OUT, "optimality_audit.csv"), "w") do io
    println(io, join(keys(rows[1]), ","))
    for r in rows
        println(io, join((x isa AbstractFloat ? @sprintf("%.6g", x) : string(x) for x in values(r)), ","))
    end
end

cnt(v) = join(["$s×$(count(==(s), v))" for s in unique(v)], " ")
q(v, p) = quantile(v, p)
println("\ncase      act   statuses N3 (loop) | N3 polished | N1 | N2")
for g in unique((r.case, r.act) for r in rows)
    S = filter(r -> (r.case, r.act) == g, rows)
    @printf("%-9s %-5s %s | %s | %s | %s\n", g..., cnt([r.st3 for r in S]), cnt([r.st_pol for r in S]),
            cnt([r.st1 for r in S]), cnt([r.st2 for r in S]))
end
println("\nrelative cost gap vs best N3 (ode1be P0 cost): median / max     |Δṁ₀| vs N3 loop [kg/s]: max")
for g in unique((r.case, r.act) for r in rows)
    S = filter(r -> (r.case, r.act) == g, rows)
    f(s) = (v = [getfield(r, s) for r in S]; @sprintf("%.1e / %.1e", q(v, 0.5), maximum(v)))
    @printf("%-9s %-5s N3 %s  pol %s  cold %s  N1 %s  N2 %s | pol %.1e cold %.1e N1 %.1e N2 %.1e\n",
            g..., f(:gap3), f(:gap_pol), f(:gap_cold), f(:gap1), f(:gap2),
            maximum(r.dm0_pol for r in S), maximum(r.dm0_cold for r in S),
            maximum(r.dm0_1 for r in S), maximum(r.dm0_2 for r in S))
end
println("written results/f5a_closed_loop/optimality_audit.csv")
