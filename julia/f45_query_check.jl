# F4.5, step 0: check that queries built from the Julia HRM replica (julia/HRMQueries.jl)
# reproduce the Simulink queries exported in F1 (results/f1_rom_replica/rom_reference_*.mat).
#
# Compared at t_q − 10 s: the HRM state mapped onto the ROM (N, Cᵢ, Cᵉᵢ, node-mean
# temperatures, ΔM) and the PI integrator — a check of the HRM replica's INTERNAL states,
# beyond the PCT validated in F1b. Compared at t_q − Δt: the ROM state that initialises the
# NLP.
#
# Usage, from the repository root:
#   julia --project=julia julia/f45_query_check.jl
include(joinpath(@__DIR__, "HRMQueries.jl"))
using .HRMQueries, MAT, Printf

const ROOT = joinpath(@__DIR__, "..")
arc = matread(joinpath(ROOT, "results", "f1b_hrm_replica", "archive_1180.mat"))
cs = arc["control_series"]
pH = hrm_params(A = arc["sampled_A"], B = arc["sampled_B"], mdot_c_drift = arc["sampled_mdot_c_drift"])
pR = rom_params()
J = SMDFRROM.IDX
names = ["N"; ["C$i" for i in 1:6]; ["Ce$i" for i in 1:6]; "T_f"; "T_w"; "T_c"; "ΔM"]

for TQ in (1161.0, 1180.0)
    r = matread(joinpath(ROOT, "results", "f1_rom_replica", "rom_reference_$(Int(TQ)).mat"))
    # control series known at t_q: rows up to t_q − Δt (row 1 of cs is the initial zero row)
    last = findfirst(==(TQ - 0.5), cs[:, 1])
    ext = ExtInput(cs[2:last, 1], cs[2:last, 2:3])
    q = build_query(pH, pR, ext, TQ)

    # Simulink values
    x0 = r["x0"]; blocks = vec(x0["block"]); vals = vec(x0["value"])
    get0(n) = vals[findfirst(b -> endswith(b, "/" * n), blocks)]
    xS0 = vcat(get0("Int_N"), vec(get0("Int-C_i")), vec(get0("Int-Ce_i")), get0("Int-T_f"),
               get0("Int-T_w"), get0("Int-T_c"), get0("Int-mdot_c"))
    xI_S = get0("Integrator")
    xs = r["xs"]; blk(n) = findfirst(b -> endswith(b, "/" * n), vec(xs["block"]))
    col(n) = xs["data"][blk(n)]
    tX = vec(xs["time"][blk("Int_N")])
    Xs = hcat(col("Int_N"), col("Int-C_i"), col("Int-Ce_i"), col("Int-T_f"), col("Int-T_w"),
              col("Int-T_c"), col("Int-mdot_c"))
    xS_init = Xs[findfirst(≈(TQ - 0.5), tX), :]

    @printf("\n=== t_q = %.1f s ===\n", TQ)
    @printf("HRM → ROM state at t_q − 10 s (replica vs Simulink):\n")
    for (k, n) in enumerate(names)
        rel = abs(q.xR0[k] - xS0[k]) / max(abs(xS0[k]), 1.0)
        abs_ = abs(q.xR0[k] - xS0[k])
        (k in (J.Tf, J.Tw, J.Tc, J.dM) || k == J.N) &&
            @printf("  %-4s Simulink %-14.8g replica %-14.8g |diff| %.2e (rel %.1e)\n", n, xS0[k], q.xR0[k], abs_, rel)
    end
    relC = maximum(abs.(q.xR0[2:13] .- xS0[2:13]) ./ abs.(xS0[2:13]))
    @printf("  C, Ce: max relative diff %.2e\n", relC)
    @printf("  PI integrator: Simulink %.10f, replica %.10f\n", xI_S, q.xI0)
    @printf("ROM state at t_q − Δt (NLP initial state):\n")
    for k in (J.Tf, J.Tw, J.Tc, J.dM)
        @printf("  %-4s Simulink %-14.8g replica %-14.8g |diff| %.2e\n", names[k], xS_init[k], q.x_init[k], abs(q.x_init[k] - xS_init[k]))
    end
    @printf("  N    relative diff %.2e\n", abs(q.x_init[1] - xS_init[1]) / xS_init[1])
    # ROM PCT estimate ŷ at t_q − Δt with the same leak/ACS used by Simulink
    ip = r["inport_data"]; k0 = findfirst(≈(TQ - 0.5), ip[:, 1])
    yS = SMDFRROM.output(xS_init, ip[k0, 2], SMDFRROM.current_to_massflow(ip[k0, 4]), pR)
    yR = SMDFRROM.output(q.x_init, q.leak_warm[end], q.acs_warm[end], pR)
    @printf("  ŷ (ROM PCT estimate) at t_q − Δt: Simulink %.6f K, replica %.6f K, diff %.2e K\n", yS, yR, yR - yS)
end
