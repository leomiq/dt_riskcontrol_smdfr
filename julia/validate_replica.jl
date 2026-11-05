# F1: validate the algebraic ROM replica against the Simulink ROM.
# Usage, from the repository root:
#   julia --project=julia julia/validate_replica.jl [results/f1_rom_replica/rom_reference_1180.mat]
include(joinpath(@__DIR__, "SMDFRROM.jl"))
using .SMDFRROM, MAT, Printf

file = isempty(ARGS) ? joinpath(@__DIR__, "..", "results", "f1_rom_replica", "rom_reference_1180.mat") : ARGS[1]
r = matread(file)
p = load_params(r["par"])
h = Float64(r["t_sam"]); t_q = Float64(r["t_q"])

# --- logged Simulink states, by block name -> replica order
xs = r["xs"]
blk(name) = findfirst(b -> endswith(b, "/" * name), vec(xs["block"]))
col(name) = xs["data"][blk(name)]
tX = vec(xs["time"][blk("Int_N")])
Xs = hcat(col("Int_N"), col("Int-C_i"), col("Int-Ce_i"), col("Int-T_f"),
          col("Int-T_w"), col("Int-T_c"), col("Int-mdot_c"))
@assert size(Xs, 2) == NX

# --- inputs on the grid
ip = r["inport_data"]
@assert ip[:, 1] ≈ tX "input grid and state grid differ"
leak = ip[:, 2]
u = ip[:, 4]
net = r["mdot_c_net_signal"]; @assert net[:, 1] ≈ tX
acs_sim = net[:, 2] .- leak                     # ACS flow actually applied in Simulink
acs_map = current_to_massflow.(u)
@printf("ACS flow: Simulink vs pump map, max |diff| = %.3e kg/s (range %.2f-%.2f kg/s)\n",
        maximum(abs.(acs_sim .- acs_map)), extrema(acs_sim)...)

# --- algebraic checks at the logged states
ypred = r["pred_tw_max_signal"][:, 2]
yrep0 = [SMDFRROM.output(Xs[k, :], leak[k], acs_sim[k], p) for k in eachindex(tX)]
@printf("Output map on Simulink states: max |ŷ_rep - ŷ_sim| = %.3e K\n", maximum(abs.(yrep0 .- ypred)))
rho_s = r["rho_tot_signal"][:, 2]
rho_r = [SMDFRROM.algebraic(Xs[k, :], p).rho for k in eachindex(tX)]
@printf("Reactivity on Simulink states: max |Δρ| = %.3e (ρ ≈ %.3e)\n", maximum(abs.(rho_r .- rho_s)), rho_s[1])

# --- free simulation from x0
x0 = Xs[1, :]
names = ["N", "C1", "C2", "C3", "C4", "C5", "C6", "Ce1", "Ce2", "Ce3", "Ce4", "Ce5", "Ce6",
         "T_f", "T_w", "T_c", "ΔM"]
hor = tX .> t_q - h                               # MPC horizon samples used by the cost
best = nothing
println("\nscheme   inputs   max|ŷ err| horizon   max|ŷ err| all   worst rel. state err")
for scheme in (:ode1be, :be), input_at in (:next, :current)
    X = simulate(x0, leak, acs_sim, h, p; scheme, input_at)
    y = [SMDFRROM.output(X[k, :], leak[k], acs_sim[k], p) for k in eachindex(tX)]
    e = abs.(y .- ypred)
    rel = vec(maximum(abs.(X .- Xs) ./ max.(abs.(Xs), 1.0); dims = 1))
    i = argmax(rel)
    @printf("%-8s %-8s %12.3e K       %10.3e K   %.2e (%s)\n", scheme, input_at,
            maximum(e[hor]), maximum(e), rel[i], names[i])
    if best === nothing || maximum(e) < best[3]
        global best = (scheme, input_at, maximum(e), X, y)
    end
end

scheme, input_at, emax, X, y = best
println("\nBest: $scheme / $input_at. Per-state max abs error:")
for (j, n) in enumerate(names)
    @printf("  %-4s %.3e  (scale %.3e)\n", n, maximum(abs.(X[:, j] .- Xs[:, j])), maximum(abs.(Xs[:, j])))
end
println("\n  t [s]    ŷ Simulink   ŷ replica   diff [K]")
for k in findall(hor)[1:5:end]
    @printf("  %7.1f  %10.4f  %10.4f  %+.2e\n", tX[k], ypred[k], y[k], y[k] - ypred[k])
end

open(joinpath(dirname(file), "replica_vs_simulink_$(round(Int, t_q)).csv"), "w") do io
    println(io, "t,y_simulink,y_replica,leak,acs,u")
    for k in eachindex(tX)
        @printf(io, "%.1f,%.8f,%.8f,%.8f,%.8f,%.8f\n", tX[k], ypred[k], y[k], leak[k], acs_sim[k], u[k])
    end
end
