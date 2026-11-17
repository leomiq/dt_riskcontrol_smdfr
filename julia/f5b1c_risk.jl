# F5b.1c (plan §8.12 step 3.ii, §8.14): risk layer on N2 against ACS failures drawn as in the
# DET (the "DET-accelerated" failure family).
#
# Failure model (DET, `matlab/mex_sources/msfcn_indirectps_hrm_c.c` and
# `matlab/simulation_core_dualfluid_FI3_sim_hrm_extract_MCDET.m`): the ACS has 5 components in
# series (any failed component fails the system; `failed_states`); each one is a two-state
# Markov chain up ⇄ down with failure rate λᵢ and repair rate μᵢ (exponential times;
# `models/MCDET_conferencepaper_database.mat`). With the real rates (Σλ ≈ 1.6e-8 s⁻¹, repairs
# of days) a failure inside a 240 s transient has probability ≈ 4e-6, so the rates are
# accelerated as in the commented importance-sampling option of the model workspace:
# λ×1e5, μ×1e4 → Σλ = 1.6e-3 s⁻¹ (MTBF 625 s), P(≥1 failure in 240 s) = 0.32, failure
# length D = mixture of exponentials (weights λᵢ/Σλ, means 46, 46, 25, 6.5 s and 2.9 h):
# median 29 s, mean 64 s, P(D > 130 s) = 0.057; q₀.₉₀ = 103 s, q₀.₉₅ = 136 s,
# q₀.₉₈ = 182 s, q₀.₉₉ = 219 s.
# The chains start all-up at the ACS activation t_act. Realisations are CONDITIONED on at least
# one system failure starting in [t_start, t_start + 240 − 20] s (rejection sampling, as the
# DET branches on the failure event); P(condition) is reported for the unconditional rate.
# Leak: the extreme family (A, B, d, lead as f5a2_monte_carlo.jl extreme); η = 1, no forecast
# bias (the DET has neither). Plant capacity m_max (default 40 kg/s), ramp 5 kg/s per step.
#
# Controllers: N2+bias; N2+bias+risk bound (`risk = (tau = q₁₋α(D), lam)`, α ∈ {0.10, 0.05,
# 0.02, 0.01}); the DET PI. Rows: results/f5b1c_risk/runs.csv (resumable, `budget` seconds).
#
# Usage, from the repository root:
#   julia --project=julia --pkgimages=no --heap-size-hint=1200M julia/f5b1c_risk.jl
#         [first_seed=201] [last_seed=260] [budget_s=480] [m_max=40] [ramp=5]
using Printf, Random, Statistics, MAT
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop

arg(i, d) = length(ARGS) >= i ? ARGS[i] : d
S0 = parse(Int, arg(1, "201")); S1 = parse(Int, arg(2, "260"))
BUDGET = parse(Float64, arg(3, "480")); MMAX = parse(Float64, arg(4, "40")); RAMP = parse(Float64, arg(5, "5"))
const NQ = 481; const WIN = 240.0
const OUT = joinpath(@__DIR__, "..", "results", "f5b1c_risk"); mkpath(OUT)
const FILE = joinpath(OUT, "runs.csv")

db = matread(joinpath(@__DIR__, "..", "models", "MCDET_conferencepaper_database.mat"))
const LAM = vec(db["Lambda_t"]) .* 1e5          # accelerated failure rates [1/s]
const MU = vec(db["Mu_t"]) .* 1e4               # accelerated repair rates [1/s]
const W = LAM ./ sum(LAM)
surv(t) = sum(W .* exp.(-MU .* t))               # P(D > t) of one system failure
function quantile_D(α)                           # q₁₋α(D): surv(q) = α (bisection)
    lo, hi = 0.0, 1e5
    for _ in 1:200; m = (lo + hi) / 2; surv(m) > α ? (lo = m) : (hi = m); end
    (lo + hi) / 2
end
const ALPHAS = (0.10, 0.05, 0.02, 0.01)
const TAU = Dict(a => quantile_D(a) for a in ALPHAS)

"""System down windows of the 5-component series ACS over [t0, t1]: each component alternates
exponential up (λᵢ) and down (μᵢ) times starting up at t0; the system is down when any is."""
function ctmc_windows(rng, t0, t1)
    iv = Tuple{Float64,Float64}[]
    for i in eachindex(LAM)
        t = t0
        while true
            t += randexp(rng) / LAM[i]; t > t1 && break
            d = randexp(rng) / MU[i]; push!(iv, (t, t + d)); t += d
        end
    end
    sort!(iv); out = Tuple{Float64,Float64}[]
    for (a, b) in iv                              # union of the component down intervals
        if !isempty(out) && a <= out[end][2]
            out[end] = (out[end][1], max(out[end][2], b))
        else
            push!(out, (a, b))
        end
    end
    out
end

function draw(seed)
    rng = MersenneTwister(seed)
    A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng); lead = 17rand(rng)
    (; A, Bp, d, lead)
end

"""Failure windows conditioned on ≥ 1 failure starting in [t_start, t_start + WIN − 20]."""
function failures(seed, t_act, t_start)
    rng = MersenneTwister(50_000 + seed); tries = 0
    while true
        tries += 1
        wins = ctmc_windows(rng, t_act, t_start + WIN)
        any(t_start <= a <= t_start + WIN - 20 for (a, _) in wins) && return (wins, tries)
    end
end

cases = vcat([(name = "N2", ctrl = :N2, kw = NamedTuple()), (name = "PI", ctrl = :PI, kw = NamedTuple())],
             [(name = "N2_risk_a$(a)", ctrl = :N2, kw = (risk = (tau = TAU[a], lam = 1e3),)) for a in ALPHAS])

fmt(x) = x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
done = Set{Tuple{Int,String}}()
if isfile(FILE)
    L = readlines(FILE); hdr = split(L[1], ","); is = findfirst(==("seed"), hdr); ic = findfirst(==("case"), hdr)
    for l in L[2:end]; f = split(l, ","); push!(done, (parse(Int, f[is]), f[ic])); end
end
@printf("q(1-α) of D: %s\n", join([@sprintf("α=%.2f: %.1f s", a, TAU[a]) for a in ALPHAS], ", "))

Bc = setup(); run_loop(Bc, :N2; n = 2, act = act_cap(30, 5), bias = true, risk = (tau = 20.0, lam = 1e4),
                       fail = [(1e9, 1e9 + 1)])
t0 = time(); nrun = 0
for seed in S0:S1
    todo = [c for c in cases if !((seed, c.name) in done)]
    isempty(todo) && continue
    r = draw(seed)
    pH = ClosedLoop.HRMQueries.hrm_params(; A = r.A, B = r.Bp, mdot_c_drift = r.d, acs_max = MMAX, acs_ramp = RAMP)
    B = setup(; pH, lead = r.lead)
    B === nothing && continue
    wins, tries = failures(seed, B.t_act, B.t_start)
    inwin = [(a, b) for (a, b) in wins if b > B.t_start && a < B.t_start + WIN]
    for c in todo
        if time() - t0 > BUDGET
            @printf("budget reached after %d runs (%.0f s); run again to continue\n", nrun, time() - t0); exit(0)
        end
        tr = @elapsed lg = run_loop(B, c.ctrl; n = NQ, act = act_cap(MMAX, RAMP), bias = c.ctrl == :N2,
                                    fail = wins, c.kw...)
        M = metrics(lg, B.y_ref)
        ub = filter(!isnan, lg.ubound)
        up = .!isnan.(lg.ubound)
        row = (; seed, case = c.name, m_max = MMAX, ramp = RAMP, A = r.A, B = r.Bp, drift = r.d,
               lead = B.t_start - B.t_act, tries, n_fail = length(inwin),
               first_fail = isempty(inwin) ? NaN : inwin[1][1] - B.t_start,
               max_fail_len = isempty(inwin) ? 0.0 : maximum(b - a for (a, b) in inwin),
               down_time = sum((min(b, B.t_start + WIN) - max(a, B.t_start) for (a, b) in inwin); init = 0.0),
               M..., margin = 1220 - M.pct_max, mean_dev = mean(lg.pct .- B.y_ref),
               ubound_mean = isempty(ub) ? NaN : mean(ub), ubound_min = isempty(ub) ? NaN : minimum(ub),
               frac_above_ubound = isempty(ub) ? NaN : count(lg.pct[up] .> lg.ubound[up]) / count(up), wall = tr)
        nf = !isfile(FILE)
        open(FILE, "a") do io
            nf && println(io, join(keys(row), ","))
            println(io, join(fmt.(values(row)), ","))
        end
        global nrun += 1
        @printf("seed %d %-16s fails %d (max %5.1f s) limit %-5s PCTmax %7.2f ISE %9.1f mass %6.0f U_α %7.1f (%.0f s)\n",
                seed, c.name, row.n_fail, row.max_fail_len, M.limit, M.pct_max, M.ise, M.mass_plant,
                row.ubound_mean, tr); flush(stdout)
    end
end
@printf("complete: %d runs in %.0f s\n", nrun, time() - t0)
