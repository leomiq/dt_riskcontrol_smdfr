# F4 / actuator: how do the F2–F4 conclusions depend on the ACS limits?
#
# The archived bound 56.75 A has no physical basis (docs/monolithic_formulation_plan.md
# §8.3): the ECM of CurrentToMassFlow.m is a main-circuit pump (0.515 m duct, 15 kPa,
# nominal flow at ≈ 2 m/s), and the ACS make-up flow uses the foot of its curve. So the
# limits are swept directly in flow units:
#   ṁ_max ∈ {12.86 (archived), 25, 43.55 (56.90 A), 63.99 (57.00 A), 125, 5543 (2 m/s in
#            the ECM duct)} kg/s,  ṁ_min = 0 (a make-up pump does not reverse)
#   ramp ∈ {soft (P0, λ_m = 1e-4), hard 5 kg/s per step from the applied flow}
# For each case: deterministic NLP with the p50 leak; first action under p05 / p95 /
# two-stage SAA (59 members); open-loop out-of-sample cost on held-out members; and the
# model-validity checks (primary flow ṁ_c/ṁ_nom, extrapolation of p₄ beyond its data).
#
# Usage, from the repository root:  julia --project=julia julia/f4_actuator_sweep.jl [t_q ...]
include(joinpath(@__DIR__, "F4Lib.jl"))
using .F4Lib, Printf, Statistics, Random

TQS = isempty(ARGS) ? ["1161", "1180"] : ARGS
M_MAX = [12.8648, 25.0, 43.55, 63.99, 125.0, 5543.0]
N_OOS = 300
out = joinpath(@__DIR__, "..", "results", "f4_risk", "actuator_sweep.csv")
io = open(out, "w")
println(io, "t_q,ramp,m_max,status,J_p50,m0_p50,m0_p05,m0_p95,m0_SAA,binding,y_end,y_max,mdotc_end_over_nom," *
            "p4_arg_max,p4_extrap_K,oos_J_p50,oos_J_p05,oos_J_p95,oos_J_SAA,t_det,t_saa")

for TQ in TQS
    c = load_query(TQ)
    p = c.p
    pv = p.poly
    p4(x) = foldl((a, ci) -> a * x + ci, pv; init = 0.0)
    Random.seed!(1)
    perm = randperm(c.NENS); ins = perm[1:59]; held = perm[60:59+N_OOS]
    a0 = c.x_init[17]
    @printf("\n=== t_q = %s s: inventory %.0f kg below nominal (%.1f %%), ṁ_c = %.3f ṁ_nom; ACS flow applied at t_q−Δt: %.3f kg/s; leak %.3f kg/s\n",
            TQ, -a0, -100a0 / p.M_c, (p.M_c + a0) / p.M_c, c.m_prev, c.leak_now)
    solve_det(c, quantile_leak(c, 50))                    # compile
    println("ramp  ṁ_max    J(p50)     ṁ₀: p50     p05      p95      SAA    | binding  ŷ_end    ṁ_c/nom  p₄ extrap [K] | oos J: p50  p05  p95  SAA")
    for ramp in (:soft, :hard), mh in M_MAX
        kw = (m_lo = 0.0, m_hi = mh, ramp = ramp)
        d50 = solve_det(c, quantile_leak(c, 50); kw...)
        d05 = solve_det(c, quantile_leak(c, 5); kw..., guess = d50.m)
        d95 = solve_det(c, quantile_leak(c, 95); kw..., guess = d50.m)
        sa = solve_saa(c, [scen(c, i) for i in ins]; kw..., guess = d50.m)
        ev = evaluate(c, d50.m, quantile_leak(c, 50))
        m0 = d50.m[1]
        binding = abs(m0 - mh) < 1e-3 * max(1, mh) ? "ṁ_max" :
                  (ramp == :hard && abs(m0 - (c.m_prev + 5.0)) < 1e-3) ? "ramp" : "none"
        mdotc_end = (p.M_c + ev.X[end, 17]) / p.M_c
        arg_max = maximum(quantile_leak(c, 50) .+ d50.m .+ 5550)
        extrap = p4(arg_max) - p4(min(arg_max, p4_fit_range[2]))   # K added by p₄ beyond its data
        oosJ(m) = mean(evaluate(c, m, scen(c, i)).J for i in held)
        o = [oosJ(d50.m), oosJ(d05.m), oosJ(d95.m), oosJ(sa.m[1])]
        @printf("%-5s %7.2f  %9.4f  %8.3f %8.3f %8.3f %8.3f | %-6s %8.3f  %6.4f   %+8.3f     | %.4f %.4f %.4f %.4f\n",
                ramp, mh, d50.obj, m0, d05.m[1], d95.m[1], sa.m0, binding, ev.y[end], mdotc_end, extrap, o...)
        println(io, join(Any[TQ, ramp, mh, d50.status, d50.obj, m0, d05.m[1], d95.m[1], sa.m0, binding,
                             ev.y[end], maximum(ev.y), mdotc_end, arg_max, extrap, o..., d50.time, sa.time], ","))
        flush(io)
    end
end
close(io)
println("\nwritten ", out)
