# F5b.1c-RB (plan §8.16): Rao–Blackwell estimator of the probability that an ACS failure
# takes the PCT above 1220 K, per controller, without counting crossings.
#
# For the FIRST system failure (Poisson rate λ = Σλᵢ of the DET-accelerated family, §8.14),
# starting at t along the failure-free trajectory x_t of a controller and lasting D:
#     P(cross) = ∫₀ᵀ λ e^{−λt} S(τ(x_t)) dt,     S(τ) = P(D > τ) = Σᵢ wᵢ e^{−μᵢ τ},
# where τ(x_t) is the time the HRM (the plant) takes to exceed 1220 K with the ACS off from x_t
# (`time_to_limit`). The duration D is integrated out analytically instead of sampled, so the
# estimator has far less variance than counting crossings; it is exact for the first failure
# and ignores later failures (second order). Also reported, conditioned on ≥ 1 failure starting
# in [0, T_c] (T_c = 220 s, the conditioning of f5b1c_risk.jl):
#     P(cross | failure in [0, T_c]) ≈ ∫₀^{T_c} λ e^{−λt} S(τ(x_t)) dt / (1 − e^{−λ T_c}).
# Validation: the DET PI gave 7/60 crossings in f5b1c_risk.jl.
#
# Realisations: the leak draws of f5b1c_risk.jl (seeds 201…), 40 kg/s, ramp 5, η = 1, NO
# failures in the run itself. Controllers: N2+bias, N2+bias+risk bound α ∈ {0.10, 0.05, 0.02,
# 0.01}, the thesis PID and the DET PI. τ is evaluated every `every` queries.
# Rows: results/f5b1c_rb_risk/runs.csv; τ series: tau.csv (resumable, `budget` seconds).
#
# Usage, from the repository root:
#   julia --project=julia --pkgimages=no --heap-size-hint=1200M julia/f5b1c_rb_risk.jl
#         [first_seed=201] [last_seed=220] [budget_s=100000] [every=20]
using Printf, Random, Statistics, MAT
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop

arg(i, d) = length(ARGS) >= i ? ARGS[i] : d
S0 = parse(Int, arg(1, "201")); S1 = parse(Int, arg(2, "220"))
BUDGET = parse(Float64, arg(3, "100000")); EVERY = parse(Int, arg(4, "20"))
const NQ = 481; const TC = 220.0; const MMAX = 40.0; const RAMP = 5.0
const OUT = joinpath(@__DIR__, "..", "results", "f5b1c_rb_risk"); mkpath(OUT)
const FILE = joinpath(OUT, "runs.csv"); const TAUF = joinpath(OUT, "tau.csv")

db = matread(joinpath(@__DIR__, "..", "models", "MCDET_conferencepaper_database.mat"))
const LAM = vec(db["Lambda_t"]) .* 1e5; const MU = vec(db["Mu_t"]) .* 1e4
const LS = sum(LAM); const W = LAM ./ LS
surv(t) = sum(W .* exp.(-MU .* t))
function quantile_D(α)
    lo, hi = 0.0, 1e5
    for _ in 1:200; m = (lo + hi) / 2; surv(m) > α ? (lo = m) : (hi = m); end
    (lo + hi) / 2
end
const ALPHAS = (0.10, 0.05, 0.02, 0.01)
cases = vcat([(name = "N2", ctrl = :N2, kw = NamedTuple())],
             [(name = "N2_risk_a$(a)", ctrl = :N2, kw = (risk = (tau = quantile_D(a), lam = 1e3),)) for a in ALPHAS],
             [(name = "PID", ctrl = :PID, kw = NamedTuple()), (name = "PI", ctrl = :PI, kw = NamedTuple())])

fmt(x) = x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
function append(file, row)
    nf = !isfile(file)
    open(file, "a") do io
        nf && println(io, join(keys(row), ","))
        println(io, join(fmt.(values(row)), ","))
    end
end
done = Set{Tuple{Int,String}}()
if isfile(FILE)
    L = readlines(FILE); hdr = split(L[1], ","); is = findfirst(==("seed"), hdr); ic = findfirst(==("case"), hdr)
    for l in L[2:end]; f = split(l, ","); push!(done, (parse(Int, f[is]), f[ic])); end
end

"""∫ λ e^{−λt} S(τ_t) dt by the trapezoid rule on the evaluation grid (t from the first query)."""
function hazard_integral(t, τ, Tmax)
    g = [LS * exp(-LS * ti) * surv(τi) for (ti, τi) in zip(t, τ)]
    k = findall(<=(Tmax + 1e-9), t)
    tk, gk = t[k], g[k]
    sum((tk[2:end] .- tk[1:end-1]) .* (gk[2:end] .+ gk[1:end-1]) ./ 2)
end

Bc = setup(); run_loop(Bc, :N2; n = 2, bias = true, keep_states = true)
time_to_limit(Bc.pH, Bc.hist[Bc.t_start-0.5][1], Bc.hist[Bc.t_start-0.5][2], Bc.t_start - 0.5; t_cap = 1.0)
t0 = time(); nrun = 0
for seed in S0:S1
    todo = [c for c in cases if !((seed, c.name) in done)]
    isempty(todo) && continue
    rng = MersenneTwister(seed)                                     # f5b1c draw
    A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng); lead = 17rand(rng)
    pH = ClosedLoop.HRMQueries.hrm_params(; A, B = Bp, mdot_c_drift = d, acs_max = MMAX, acs_ramp = RAMP)
    B = setup(; pH, lead); B === nothing && continue
    for c in todo
        if time() - t0 > BUDGET
            @printf("budget reached after %d runs (%.0f s); run again to continue\n", nrun, time() - t0); exit(0)
        end
        tr = @elapsed lg = run_loop(B, c.ctrl; n = NQ, act = act_cap(MMAX, RAMP), bias = c.ctrl == :N2,
                                    keep_states = true, c.kw...)
        M = metrics(lg, B.y_ref)
        idx = unique(vcat(1:EVERY:length(lg.t_q), length(lg.t_q)))
        tt = @elapsed τ = [time_to_limit(pH, lg.states[i][1], lg.states[i][2], lg.t_q[i]) for i in idx]
        t = lg.t_q[idx] .- lg.t_q[1]
        p_unc = hazard_integral(t, τ, t[end])
        p_cond = hazard_integral(t, τ, TC) / (1 - exp(-LS * TC))
        append(FILE, (; seed, case = c.name, A, B = Bp, drift = d, lead = B.t_start - B.t_act, n = M.n,
                      limit = M.limit, p_cross = p_unc, p_cross_cond = p_cond, tau_mean = mean(τ),
                      tau_min = minimum(τ), pct_max = M.pct_max, ise = M.ise, mass = M.mass_plant,
                      mean_dev = mean(lg.pct .- B.y_ref), wall_loop = tr, wall_tau = tt))
        open(TAUF, "a") do io
            isfile(TAUF) && filesize(TAUF) == 0 && println(io, "seed,case,t,tau,pct")
            for (j, i) in enumerate(idx)
                println(io, join(fmt.((seed, c.name, t[j], τ[j], lg.pct[i])), ","))
            end
        end
        global nrun += 1
        @printf("seed %d %-16s P(cross|fail) %.4f  P(cross) %.5f  τ mean %6.1f min %6.1f  PCTmax %7.2f (%.0f + %.0f s)\n",
                seed, c.name, p_cond, p_unc, mean(τ), minimum(τ), M.pct_max, tr, tt); flush(stdout)
    end
end
@printf("complete: %d runs in %.0f s\n", nrun, time() - t0)
