# F5a.1: deterministic closed loop on the HRM replica (no risk layers, perfect leak forecast).
#
# Controllers: N1 (LTI QP), N2 (real-time-iteration QP), N3 (NLP), each nominal and with the
# offset-free correction; baselines: archived LJ currents (1160–1179.5 s only), no ACS,
# constant maximum flow. Actuator configurations: archived current box with the soft ramp
# of P0, and the physical one (0–5543 kg/s, hard ramp 5 kg/s per step). See
# julia/ClosedLoop.jl for the timing and docs/monolithic_formulation_plan.md §8.6 for the
# metrics fixed in advance.
#
# Usage, from the repository root:
#   julia --project=julia julia/f5a_closed_loop.jl [n_queries=121]   (§8.9 used 61)
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop, Printf, Statistics

const OUT = joinpath(@__DIR__, "..", "results", "f5a_closed_loop"); mkpath(OUT)
n = isempty(ARGS) ? 121 : parse(Int, ARGS[1])
B = setup()
y_ref = B.y_ref
acts = (B.act_arch, ACT_PHYS)

# compile every controller once (excluded from the timing statistics)
for a in acts, ctrl in (:N1, :N2, :N3); run_loop(B, ctrl; n = 2, act = a); end

runs = Any[]
for a in acts
    cases = [(:zero, false), (:N1, false), (:N1, true), (:N2, false), (:N2, true), (:N3, false), (:N3, true)]
    a.name == :arch && append!(cases, [(:LJ, false), (:max, false)])
    for (ctrl, bias) in cases
        nn = ctrl == :LJ ? 40 : n
        t = @elapsed lg = run_loop(B, ctrl; n = nn, act = a, bias)
        push!(runs, (; act = a.name, ctrl, bias, log = lg, wall = t))
        @printf("%-4s %-4s bias=%-5s  %5.1f s  PCT max %.3f K\n", a.name, ctrl, bias, t, maximum(lg.pct))
    end
end

"""Metrics on the queries with t_q ≤ t_end (PCT logged at t_q, after the applied step)."""
function metrics(lg, t_end)
    i = findall(<=(t_end), lg.t_q)
    pct = lg.pct[i]; m = lg.m[i]; e = pct .- y_ref
    tm = lg.time[i]
    (; n = length(i), pct_max = maximum(pct), t_above = 0.5 * count(>(y_ref), pct),
     ise = sum(e .^ 2), iae = 0.5 * sum(abs.(e)), mean_abs_err_last10 = mean(abs.(e[max(1, end - 19):end])),
     mass = 0.5 * sum(m), ramp_max = length(m) > 1 ? maximum(abs.(diff(vcat(0.0, m)))) : 0.0,
     t_p50 = median(tm), t_p99 = quantile(tm, 0.99), t_max = maximum(tm),
     fails = count(s -> !(s in ("-", "LOCALLY_SOLVED", "OPTIMAL", "ALMOST_OPTIMAL", "ALMOST_LOCALLY_SOLVED")), lg.status[i]),
     d_mean = mean(lg.d[i]))
end

open(joinpath(OUT, "summary.csv"), "w") do io
    println(io, "act,ctrl,bias,window_end,n,pct_max,t_above_yref_s,ise,iae,mean_abs_err_last10s,mass_kg,ramp_max,solve_p50_s,solve_p99_s,solve_max_s,solver_fails,mean_nowcast_bias")
    for r in runs, t_end in (1179.5, B.t_start + (n - 1) * 0.5)
        r.ctrl == :LJ && t_end > 1180 && continue
        M = metrics(r.log, t_end)
        @printf(io, "%s,%s,%s,%.1f,%d,%.6f,%.1f,%.6f,%.6f,%.6f,%.4f,%.4f,%.5f,%.5f,%.5f,%d,%.5f\n",
                r.act, r.ctrl, r.bias, t_end, M.n, M.pct_max, M.t_above, M.ise, M.iae, M.mean_abs_err_last10,
                M.mass, M.ramp_max, M.t_p50, M.t_p99, M.t_max, M.fails, M.d_mean)
    end
end
open(joinpath(OUT, "trajectories.csv"), "w") do io
    println(io, "act,ctrl,bias,t_q,pct,m,u,yhat_now,d,solve_s,status")
    for r in runs, i in eachindex(r.log.t_q)
        L = r.log
        @printf(io, "%s,%s,%s,%.1f,%.9f,%.9f,%.9f,%.9f,%.9f,%.6f,%s\n", r.act, r.ctrl, r.bias, L.t_q[i],
                L.pct[i], L.m[i], L.u[i], L.yhat_now[i], L.d[i], L.time[i], L.status[i])
    end
end

println("\nWindow 1160–1179.5 s (41 queries, comparable with the archived LJ loop):")
println("act  ctrl bias   PCT max   t>yref[s]   ISE        IAE      mass[kg]  solve p50/p99 [s]  fails")
for r in runs
    M = metrics(r.log, 1179.5)
    @printf("%-4s %-4s %-5s %9.3f %8.1f %10.3f %9.3f %9.2f   %.4f / %.4f   %d\n", r.act, r.ctrl, r.bias,
            M.pct_max, M.t_above, M.ise, M.iae, M.mass, M.t_p50, M.t_p99, M.fails)
end
println("\nFull window 1160–$(B.t_start + (n - 1) * 0.5) s ($(n) queries):")
for r in runs
    r.ctrl == :LJ && continue
    M = metrics(r.log, B.t_start + (n - 1) * 0.5)
    @printf("%-4s %-4s %-5s %9.3f %8.1f %10.3f %9.3f %9.2f   %.4f / %.4f   %d   |e| last 10 s %.4f K\n", r.act, r.ctrl, r.bias,
            M.pct_max, M.t_above, M.ise, M.iae, M.mass, M.t_p50, M.t_p99, M.fails, M.mean_abs_err_last10)
end
println("written results/f5a_closed_loop/{summary,trajectories}.csv")
