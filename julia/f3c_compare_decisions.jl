# F3c: compare the DECISIONS of every method on the same yardstick.
#
# Each method delivers a control sequence u (26 ACS flows). Every u is evaluated with the
# nonlinear ROM in both discretizations (converged implicit Euler = the NLP model; ode1be =
# Simulink) and placed against
#   - the best cost found by any method (J_best),
#   - the certified lower bound LB of the nonlinear P0 (step 4b, midpoint center), so that
#     J(u) − J* ≤ J(u) − LB is a GUARANTEED bound on the excess of u over the global optimum,
#   - the semi-rigorous LB of step 4a (ε by sampling, safety factor 1).
# Controls come from the stored results of F2/E2, F3b and F3c.
#
# Usage, from the repository root:  julia --project=julia julia/f3c_compare_decisions.jl
include(joinpath(@__DIR__, "F4Lib.jl"))
using .F4Lib
include(joinpath(@__DIR__, "ReducedP0.jl"))
using Printf

root = joinpath(@__DIR__, "..")
function readcsv(path)
    lines = readlines(path); hdr = String.(split(lines[1], ","))
    rows = [parse.(Float64, split(l, ",")) for l in lines[2:end] if !isempty(l)]
    Dict(h => [r[j] for r in rows] for (j, h) in enumerate(hdr))
end

out = open(joinpath(root, "results", "f3c_optimality", "decisions_compared.csv"), "w")
println(out, "t_q,method,J_BE,J_ode1be,excess_vs_best_pct,guaranteed_excess_vs_global_pct,u0_A,ymax_K,lb_certified,lb_semirigorous")
for TQ in ("1180", "1161")
    c = load_query(TQ); P = reduced_p0(c)
    m_lo, m_hi = P.m_lo, P.m_hi
    dec = Pair{String,Vector{Float64}}[]
    tr = readcsv(joinpath(root, "results", "f3_certificates", "trajectories_$TQ.csv"))
    so = readcsv(joinpath(root, "results", "f3c_optimality", "solutions_$TQ.csv"))
    push!(dec, "N3 NLP (Ipopt, Euler implícito)" => so["u_NLP_BE"])
    push!(dec, "N3 NLP sobre ode1be" => so["u_NLP_ode1be"])
    push!(dec, "N2 QP-LTV iterado (SQP)" => tr["N2_ode1be"])
    push!(dec, "N1 QP-LTI (una linealización)" => tr["N1"])
    push!(dec, "Evolución diferencial (90k eval.)" => so["u_DE_best"])
    e2f = joinpath(root, "results", "e2_optimizers_$TQ", "e2_actions.csv")
    if isfile(e2f)
        e2 = readcsv(e2f)
        for (k, name) in (("LJ_seed1", "LJ corregido, semilla 1"), ("LJ_seed2", "LJ corregido, semilla 2"),
                          ("LJ_seed3", "LJ corregido, semilla 3"), ("GA_seed1", "GA, semilla 1"),
                          ("GA_seed2", "GA, semilla 2"), ("LJ_archived", "LJ archivado (lazo original)"))
            push!(dec, name => current_to_massflow.(e2[k]))
        end
    else
        f2 = readcsv(joinpath(root, "results", "f1_rom_replica", "f2_nlp_vs_lj_$TQ.csv"))
        push!(dec, "Acción archivada del lazo" => f2["mdot_LJ"])
    end
    lbc = readcsv(joinpath(root, "results", "f3c_optimality", "rigorous_gap_$(TQ)_mid.csv"))["LB"][1]
    gsr = readcsv(joinpath(root, "results", "f3c_optimality", "lower_bound_gap_$TQ.csv"))
    lbs = gsr["LB"][findfirst(==(1.0), gsr["safety_factor"])]
    vals = [(name, P.J(clamp.(u, m_lo, m_hi)), P.J(clamp.(u, m_lo, m_hi); scheme = :ode1be), u) for (name, u) in dec]
    Jbest = minimum(v[2] for v in vals)
    println("\n=== t_q = $TQ s  (ROM no lineal; LB certificada = $(round(lbc; digits = 6)), LB semi-rigurosa = $(round(lbs; digits = 6)))")
    @printf("%-36s %13s %13s %16s %20s %9s\n", "método", "J (Euler impl.)", "J (ode1be)", "vs mejor [%]", "exceso garantizado [%]", "u₀ [A]")
    for (name, Jb, Jo, u) in sort(vals; by = v -> v[2])
        ex = 100 * (Jb - Jbest) / Jbest
        gex = 100 * (Jb - lbc) / lbc                       # J(u) − J* ≤ J(u) − LB, relative to LB ≤ J*
        y = P.predict(clamp.(u, m_lo, m_hi))
        @printf("%-36s %13.7f %13.7f %16.4f %20.2f %9.4f\n", name, Jb, Jo, ex, gex, massflow_to_current(max(u[1], 1e-9)))
        println(out, join((TQ, "\"" * name * "\"", Jb, Jo, ex, gex, massflow_to_current(max(u[1], 1e-9)), maximum(y), lbc, lbs), ","))
    end
end
close(out)
