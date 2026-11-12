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
# mode = extreme (F5a.3, plan §8.11): family :extreme for every seed — severe leaks
#   d ~ U(15, 25) kg/s, A ~ U(0.1, 0.2), B ~ U(10, 15) (the leak at activation, 10–17 kg/s,
#   reaches or exceeds the archived pump box of 12.86 kg/s), the MPC starting U(0, 17) s after
#   the activation (acute phase), η ~ U(0.6, 0.9), an ACS failure of U(30, 90) s starting
#   U(5, 20) s after the first query, and a leak-forecast error ~ U(−0.5, 0.1) kg/s. Adds the
#   DET's own PI as a baseline (`ctrl = :PI`, actuator `det`) and drops N1 and N3 (equal to N2
#   in F5a.1/F5a.2; the paired reference is N2+bias). Per-query PCT and flows of every case
#   go to runs_extreme_traj.csv. If the plant
#   breaks down (coolant inventory exhausted) the run stops and is flagged.
#   Results go to runs_extreme.csv.
#
# Usage, from the repository root (W worker processes; 0 = run in this process):
#   julia --project=julia julia/f5a2_monte_carlo.jl [R=6] [W=0] [n=121] [first_seed=1] [mode=standard|extreme]
using Distributed, Printf, Statistics, Random

R = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 6
W = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 0
NQ = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 121
SEED0 = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 1
MODE = length(ARGS) >= 5 ? Symbol(ARGS[5]) : :standard
MODE in (:standard, :extreme) || error("mode must be standard or extreme")
const OUT = joinpath(@__DIR__, "..", "results", "f5a2_monte_carlo"); mkpath(OUT)

if W > 0
    addprocs(W; exeflags = ["--project=$(joinpath(@__DIR__))", "--heap-size-hint=1200M"])
end
@everywhere include(joinpath($(@__DIR__), "ClosedLoop.jl"))
@everywhere using .ClosedLoop, Random

@everywhere function realisation(seed, n, mode)
    rng = MersenneTwister(seed)
    if mode == :extreme
        family = :extreme
        A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng)
        eta = 0.6 + 0.3rand(rng); fs = 5 + 15rand(rng); fl = 30 + 60rand(rng); fb = -0.5 + 0.6rand(rng)
        lead = 17rand(rng)
    else
        family = isodd(seed) ? :nominal : :stressed
        A = 0.2 + 0.2rand(rng); Bp = 9 + 6rand(rng); d = 8 + 2rand(rng)
        eta = 0.8 + 0.2rand(rng); fs = 5 + 15rand(rng); fl = 5 + 10rand(rng); fb = -0.25 + 0.3rand(rng)
        lead = 17.0
    end
    pH = ClosedLoop.HRMQueries.hrm_params(; A, B = Bp, mdot_c_drift = d)
    t0 = time()
    B = setup(; pH, lead)
    info = (; seed, family, A, B = Bp, drift = d)
    B === nothing && return (; info..., skipped = true, rows = NamedTuple[])
    pert = family == :nominal ? (eta = 1.0, fail = nothing, fbias = 0.0) :
           (eta = eta, fail = (B.t_start + fs, B.t_start + fs + fl), fbias = fb)
    cases = [(:zero, ACT_PHYS, false), (:max, B.act_arch, false)]
    mode == :extreme && push!(cases, (:PI, ACT_PHYS, false))
    for a in (B.act_arch, ACT_PHYS)
        mode == :extreme || push!(cases, (:N1, a, true))
        append!(cases, [(:N2, a, false), (:N2, a, true)])
        mode == :extreme || push!(cases, (:N3, a, true))
    end
    rows = NamedTuple[]; traj = NamedTuple[]
    for (ctrl, a, bias) in cases
        lg = run_loop(B, ctrl; n, act = a, bias, pert...)
        M = metrics(lg, B.y_ref)
        actname = ctrl == :zero ? :none : ctrl == :PI ? :det : a.name
        mode == :extreme && for j in eachindex(lg.t_q)
            push!(traj, (; seed, act = actname, ctrl, bias, t = lg.t_q[j] - B.t_act, pct = lg.pct[j],
                         m = lg.m[j], m_plant = lg.m_plant[j]))
        end
        push!(rows, (; info..., lead = B.t_start - B.t_act, t_act = B.t_act, t_start = B.t_start,
                     eta = pert.eta, fail_start = pert.fail === nothing ? NaN : pert.fail[1],
                     fail_len = pert.fail === nothing ? 0.0 : pert.fail[2] - pert.fail[1],
                     fbias = pert.fbias,
                     act = ctrl == :zero ? :none : ctrl == :PI ? :det : a.name, ctrl, bias, M...))
    end
    (; info..., skipped = false, rows, traj, wall = time() - t0)
end

# compile on every process before timing (one short archived loop per controller)
@everywhere begin
    Bc = setup()
    for c in (:N1, :N2, :N3); run_loop(Bc, c; n = 2); run_loop(Bc, c; n = 2, act = ACT_PHYS); end
end

seeds = SEED0:(SEED0 + R - 1)

# ---- CSV of every run (appends, so batches can be combined). Each realisation is written as
# soon as it finishes, so a killed batch keeps the realisations already done.
file = joinpath(OUT, MODE == :extreme ? "runs_extreme.csv" : "runs.csv")
function write_rows(rows)
    isempty(rows) && return
    global file
    header = join(keys(rows[1]), ",")
    if isfile(file) && readline(file) != header       # columns changed since the file was started
        file = replace(file, ".csv" => "_v2.csv")
        @warn "columns differ from the existing file; writing to $file"
    end
    newfile = !isfile(file)
    open(file, "a") do io
        newfile && println(io, header)
        for r in rows
            println(io, join((x isa AbstractFloat ? @sprintf("%.8g", x) : string(x) for x in values(r)), ","))
        end
    end
end
function write_traj(r)
    (r.skipped || isempty(r.traj)) && return
    tf = joinpath(OUT, "runs_extreme_traj.csv"); nf = !isfile(tf)
    open(tf, "a") do io
        nf && println(io, join(keys(r.traj[1]), ","))
        for q in r.traj
            println(io, join((x isa AbstractFloat ? @sprintf("%.8g", x) : string(x) for x in values(q)), ","))
        end
    end
end
run_one(s) = (r = realisation(s, NQ, MODE); write_rows(r.rows); write_traj(r);
              @printf("seed %d done (%s)\n", s, r.skipped ? "skipped" : @sprintf("%.0f s", r.wall)); flush(stdout); r)
t_all = @elapsed res = W > 0 ? pmap(s -> realisation(s, NQ, MODE), seeds) : map(run_one, seeds)
W > 0 && foreach(r -> (write_rows(r.rows); write_traj(r)), res)
rows = reduce(vcat, [r.rows for r in res]; init = NamedTuple[])
skipped = [r.seed for r in res if r.skipped]
@printf("%d realisations (%d skipped: no ACS demand) in %.1f s with %d workers\n", R, length(skipped), t_all, W)

const REFCTRL = MODE == :extreme ? :N2 : :N3
# ---- paired summary against N3+bias (N2+bias in extreme mode) with the same actuator (common random numbers)
key(r) = (r.act, r.ctrl, r.bias)
println("\nfamily    actuator ctrl  bias   n   ISE mean   ΔISE vs ref+bias (mean ± 95% CI)   PCT max   t>yref[s]  |e| last10   p99 solve  fails  #limit(1220K) t_limit[s] breakdowns  (runs stopped at the limit are shorter: ISE not comparable)")
for fam in (:nominal, :stressed, :extreme)
    F = filter(r -> r.family == fam, rows)
    isempty(F) && continue
    for k in unique(key.(F))
        S = filter(r -> key(r) == k, F)
        refact = k[1] in (:none, :det) ? :phys : (k[1] == :arch && k[2] == :max ? :arch : k[1])
        Δ = Float64[]
        for r in S
            o = filter(q -> q.seed == r.seed && q.act == refact && q.ctrl == REFCTRL && q.bias, F)
            isempty(o) || push!(Δ, r.ise - o[1].ise)
        end
        ci = length(Δ) > 1 ? 1.96 * std(Δ) / sqrt(length(Δ)) : NaN
        tl = filter(!isnan, [r.t_limit for r in S])
        @printf("%-9s %-8s %-5s %-5s %3d %10.2f   %+9.3f ± %-8.3f            %8.3f %8.1f %10.4f %10.3f %5d %9d %9.1f %5d\n",
                fam, k[1], k[2], k[3], length(S), mean(r.ise for r in S), mean(Δ), ci,
                mean(r.pct_max for r in S), mean(r.t_above for r in S), mean(r.err_last10 for r in S),
                maximum(r.t_p99 for r in S), sum(r.fails for r in S), count(r.limit for r in S),
                isempty(tl) ? NaN : mean(tl), count(r.breakdown for r in S))
    end
end
println("appended to ", relpath(file, joinpath(@__DIR__, "..")))
