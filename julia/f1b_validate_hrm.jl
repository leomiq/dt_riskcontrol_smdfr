# F1b.3: validate the Julia HRM replica against the archived Simulink HRM run (no MATLAB).
#
# Reproduces what `matlab/DET_scheduler_v1_perc_online_LJ.m` did for the archived loop
# (`models/LJ_test_sampled_vel_ACS_2_loop_1180.mat`, extracted by
# `python/f1b_extract_archive.py`):
#   query 1 (t = 1160 s): one HRM run 0 → 1160 s with zero inputs; the whole trajectory
#                          of PCT, leak and Te_c(12) was logged;
#   query k ≥ 2:          a NEW run from t = 0 to t_k with the control series known up to
#                          t_{k-1} (the inport extrapolates linearly to t_k); only the END
#                          values were logged.
#
# Usage, from the repository root:
#   julia --project=julia julia/f1b_validate_hrm.jl
include(joinpath(@__DIR__, "SMDFRHRM.jl"))
using .SMDFRHRM, MAT, Printf, Statistics

const OUT = joinpath(@__DIR__, "..", "results", "f1b_hrm_replica")
arc = matread(joinpath(OUT, "archive_1180.mat"))
cs = arc["control_series"]                 # [t u_ACS active_mpc]; row 1 = initial zero row
tw = arc["tw_max_V2P_trajectory"]; X1 = vec(arc["X_1_trajectory"]); X2 = vec(arc["X_2_trajectory"])
p = hrm_params(A = arc["sampled_A"], B = arc["sampled_B"], mdot_c_drift = arc["sampled_mdot_c_drift"])
@printf("Leak parameters: A = %.6f, B = %g, mdot_c_drift = %g kg/s, fault at %g s\n",
        p.A, p.B, p.mdot_c_drift, p.fault_time)

ext_upto(r) = ExtInput(cs[2:r, 1], cs[2:r, 2:3])

# ---- query 1: full trajectory 0 → 1160 s
t1 = @elapsed sim1 = simulate(p, ext_upto(2321), 1160.0)
ref_t = tw[2:2322, 1]
@assert sim1.t ≈ ref_t
e_pct = sim1.pct .- tw[2:2322, 2]
e_leak = sim1.leak .- X1[2:2322]
e_tec = sim1.Te_c12 .- X2[2:2322]
@printf("\nQuery 1 (run 0-1160 s, %d steps, %.2f s):\n", length(ref_t) - 1, t1)
for (name, e) in (("PCT max_i T_w", e_pct), ("leak", e_leak), ("Te_c(12)", e_tec))
    @printf("  %-14s max |err| = %.3e   (t of max: %.1f s)\n", name, maximum(abs.(e)),
            ref_t[argmax(abs.(e))])
end
w = ref_t .>= 1100
@printf("  PCT, window 1100-1160 s: max |err| = %.3e K\n", maximum(abs.(e_pct[w])))
ia = findfirst(>(0), sim1.act)
@printf("  ACS activated by the DET logic at t = %.1f s (PCT > 1180 K at %.1f s); max ACS flow before 1160 s = %.3g kg/s\n",
        sim1.t[ia], sim1.t[findfirst(>(1180), sim1.pct)], maximum(sim1.acs))

# ---- queries 2..41: end values of each re-simulation
rows = NamedTuple[]
push!(rows, (; t = 1160.0, pct_sim = sim1.pct[end], pct_arc = tw[2322, 2], leak_sim = sim1.leak[end],
           leak_arc = X1[2322], tec_sim = sim1.Te_c12[end], tec_arc = X2[2322], u_end = 0.0))
t2 = @elapsed for k in 2:41
    tk = 1160 + 0.5 * (k - 1)
    ext = ext_upto(2320 + k)
    s = simulate(p, ext, tk)
    push!(rows, (; t = tk, pct_sim = s.pct[end], pct_arc = tw[2321+k, 2], leak_sim = s.leak[end],
               leak_arc = X1[2321+k], tec_sim = s.Te_c12[end], tec_arc = X2[2321+k],
               u_end = value(ext, tk, 1)))
end
@printf("\nQueries 2-41 (40 re-simulations from t = 0, %.1f s):\n", t2)
ep = [r.pct_sim - r.pct_arc for r in rows]
el = [r.leak_sim - r.leak_arc for r in rows]
ec = [r.tec_sim - r.tec_arc for r in rows]
@printf("  PCT end value      max |err| = %.3e K\n", maximum(abs.(ep)))
@printf("  leak end value     max |err| = %.3e kg/s\n", maximum(abs.(el)))
@printf("  Te_c(12) end value max |err| = %.3e K\n", maximum(abs.(ec)))
r2 = rows[2]
@printf("  t = 1160.5 s (extrapolated current %.2f A): PCT sim %.4f K, archive %.4f K\n",
        r2.u_end, r2.pct_sim, r2.pct_arc)

# ---- "true" trajectory: the same inputs, all known (no extrapolation artefact)
cs_full = vcat(cs[2:end, :], [1180.0 cs[end, 2] 1.0])
true_run = simulate(p, ExtInput(cs_full[:, 1], cs_full[:, 2:3]), 1180.0)
it = findall(>=(1160), true_run.t)
dart = [true_run.pct[findfirst(==(r.t), true_run.t)] - r.pct_arc for r in rows]
@printf("\nArchived PCT at the queries vs PCT with all inputs known: max |diff| = %.3f K (at %.1f s), median %.3f K\n",
        maximum(abs.(dart)), rows[argmax(abs.(dart))].t, median(abs.(dart)))

open(joinpath(OUT, "validation_1180_queries.csv"), "w") do io
    println(io, "t,pct_replica,pct_archive,leak_replica,leak_archive,Tec12_replica,Tec12_archive,u_at_t_extrapolated,pct_all_inputs_known")
    for (r, d) in zip(rows, dart)
        @printf(io, "%.1f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.6f,%.9f\n", r.t, r.pct_sim, r.pct_arc,
                r.leak_sim, r.leak_arc, r.tec_sim, r.tec_arc, r.u_end, r.pct_arc + d)
    end
end
open(joinpath(OUT, "validation_1180_run1.csv"), "w") do io
    println(io, "t,pct_replica,pct_archive,leak_replica,leak_archive,Tec12_replica,Tec12_archive,act")
    for i in eachindex(sim1.t)
        @printf(io, "%.1f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%d\n", sim1.t[i], sim1.pct[i], tw[i+1, 2],
                sim1.leak[i], X1[i+1], sim1.Te_c12[i], X2[i+1], sim1.act[i])
    end
end
println("\nwritten results/f1b_hrm_replica/validation_1180_{run1,queries}.csv")
