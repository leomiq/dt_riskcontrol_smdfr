# F5a.2: Monte Carlo of deterministic closed loops on the HRM replica (no risk layers).
#
# Each realisation r (seeded, common random numbers for all controllers) draws:
#   leak (Park et al., 2021): A ~ U(0.2, 0.4), B ~ U(9, 15), drift d ~ U(8, 10) kg/s — the
#     sampling ranges of the model workspace (A_l/A_s, B_l/B_s, unif_dist_drift);
#   family :nominal   — only the leak varies; perfect leak forecast;
#   family :stressed  — plus perturbations the MPC does not know about:
#     ACS efficiency η ~ U(0.8, 1.0), an ACS failure starting U(5, 20) s after the first
#     query and lasting U(5, 15) s, and a constant leak-forecast error ~ U(−0.25, 0.05) kg/s
#     (the BiLSTM over-predicts the leak by ≈ 0.12 kg/s, plan §8.2).
# The MPC starts 17 s after the ACS activation (as in the archived loop) and runs `n`
# queries. Realisations whose ACS is never demanded before 2500 s are skipped (reported).
#
# Controllers: no ACS, constant maximum (archived box), and N1+bias, N2, N2+bias, N3+bias
# with the archived and the physical actuator configurations.
#
# Usage, from the repository root (W worker processes; 0 = run in this process):
#   julia --project=julia julia/f5a2_monte_carlo.jl [R=6] [W=0] [n=121] [first_seed=1]
using Distributed, Printf, Statistics, Random

R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 6
W = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 0
NQ = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 121
SEED0 = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 1
const OUT = joinpath(@__DIR__, "..", "results", "f5a2_monte_carlo"); mkpath(OUT)

if W > 0
    addprocs(W; exeflags = ["--project=$(joinpath(@__DIR__))", "--heap-size-hint=1200M"])
end
@everywhere include(joinpath($(@__DIR__), "ClosedLoop.jl"))
@everywhere using .ClosedLoop, Random

@everywhere function realisation(seed, n)
    rng = MersenneTwister(seed)
    family = isodd(seed) ? :nominal : :stressed
    A = 0.2 + 0.2rand(rng); Bp = 9 + 6rand(rng); d = 8 + 2rand(rng)
    eta = 0.8 + 0.2rand(rng); fs = 5 + 15rand(rng); fl = 5 + 10rand(rng); fb = -0.25 + 0.3rand(rng)
    pH = ClosedLoop.HRMQueries.hrm_params(; A, B = Bp, mdot_c_drift = d)
    t0 = time()
    B = setup(; pH)
    info = (; seed, family, A, B = Bp, drift = d)
    B === nothing && return (; info..., skipped = true, rows = NamedTuple[])
    pert = family == :nominal ? (eta = 1.0, fail = nothing, fbias = 0.0) :
           (eta = eta, fail = (B.t_start + fs, B.t_start + fs + fl), fbias = fb)
    cases = [(:zero, ACT_PHYS, false), (:max, B.act_arch, false)]
    for a in (B.act_arch, ACT_PHYS)
        append!(cases, [(:N1, a, true), (:N2, a, false), (:N2, a, true), (:N3, a, true)])
    end
    rows = NamedTuple[]
    for (ctrl, a, bias) in cases
        lg = run_loop(B, ctrl; n, act = a, bias, pert...)
        M = metrics(lg, B.y_ref)
        push!(rows, (; info..., t_act = B.t_act, t_start = B.t_start,
                     eta = pert.eta, fail_start = pert.fail === nothing ? NaN : pert.fail[1],
                     fail_len = pert.fail === nothing ? 0.0 : pert.fail[2] - pert.fail[1],
                     fbias = pert.fbias, act = ctrl == :zero ? :none : a.name, ctrl, bias, M...))
    end
    (; info..., skipped = false, rows, wall = time() - t0)
end

# compile on every process before timing (one short archived loop per controller)
@everywhere begin
    Bc = setup()
    for c in (:N1, :N2, :N3); run_loop(Bc, c; n = 2); run_loop(Bc, c; n = 2, act = ACT_PHYS); end
end

seeds = SEED0:(SEED0 + R - 1)
t_all = @elapsed res = W > 0 ? pmap(s -> realisation(s, NQ), seeds) : map(s -> realisation(s, NQ), seeds)
rows = reduce(vcat, [r.rows for r in res]; init = NamedTuple[])
skipped = [r.seed for r in res if r.skipped]
@printf("%d realisations (%d skipped: no ACS demand) in %.1f s with %d workers\n", R, length(skipped), t_all, W)

# ---- CSV of every run (appends if the file exists, so batches can be combined)
file = joinpath(OUT, "runs.csv")
newfile = !isfile(file)
open(file, "a") do io
    newfile && println(io, join(keys(rows[1]), ","))
    for r in rows
        println(io, join((x isa AbstractFloat ? @sprintf("%.8g", x) : string(x) for x in values(r)), ","))
    end
end

# ---- paired summary against N3+bias with the same actuator (common random numbers)
key(r) = (r.act, r.ctrl, r.bias)
println("\nfamily    actuator ctrl  bias   n   ISE mean   ΔISE vs N3+bias (mean ± 95% CI)   PCT max   t>yref[s]  |e| last10   p99 solve  fails")
for fam in (:nominal, :stressed)
    F = filter(r -> r.family == fam, rows)
    isempty(F) && continue
    for k in unique(key.(F))
        S = filter(r -> key(r) == k, F)
        refact = k[1] == :none ? :phys : (k[1] == :arch && k[2] == :max ? :arch : k[1])
        Δ = Float64[]
        for r in S
            o = filter(q -> q.seed == r.seed && q.act == refact && q.ctrl == :N3 && q.bias, F)
            isempty(o) || push!(Δ, r.ise - o[1].ise)
        end
        ci = length(Δ) > 1 ? 1.96 * std(Δ) / sqrt(length(Δ)) : NaN
        @printf("%-9s %-8s %-5s %-5s %3d %10.2f   %+9.3f ± %-8.3f            %8.3f %8.1f %10.4f %10.3f %5d\n",
                fam, k[1], k[2], k[3], length(S), mean(r.ise for r in S), mean(Δ), ci,
                mean(r.pct_max for r in S), mean(r.t_above for r in S), mean(r.err_last10 for r in S),
                maximum(r.t_p99 for r in S), sum(r.fails for r in S))
    end
end
println("appended to results/f5a2_monte_carlo/runs.csv")
