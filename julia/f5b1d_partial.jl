# F5b.1d (plan §8.18): PARTIAL ACS failures inside the planning horizon, with a severity η.
#
# Failure model (assumption to confirm with the team; the DET has no partial failures): the
# timing and lengths of the DET-accelerated family (§8.14: 5 components in series, Markov
# chains λ×1e5, μ×1e4), and each failure takes out 1–4 of 4 parallel ACS trains, so the ACS
# delivers η·ṁ with η ∈ {0.75, 0.5, 0.25, 0} and probabilities {0.4, 0.3, 0.2, 0.1}.
# During a failure N2 pushes its maximum, so the plant receives η·ṁ_max.
#
# Evaluation: Rao–Blackwell (§8.16) generalised to severities, on failure-free trajectories:
#     P(cross | failure in [0, T_c]) = ∫ λ e^{−λt} Σ_η p_η S(τ_η(x_t)) dt / (1 − e^{−λ T_c}),
# τ_η = HRM time to 1220 K with the ACS delivering η·ṁ_max (`time_to_limit(…; flow)`); a
# failure that does not reach 1220 K within 400 s counts as no crossing.
#
# Controllers (40 kg/s, ramp 5): N2+bias; the total-failure bound of §8.14 (α = 0.02, assumes
# η = 0); the η-aware quantile bound (`:var_eta`, α ∈ {0.05, 0.01}); and the two formulations
# that are NOT a bound because η multiplies the decision inside the horizon: scenario CVaR_β
# (`:scen_cvar`, β ∈ {0.90, 0.98}) and scenario expected excess (`:scen_exp`, ρ ∈ {100, 1000}),
# on scenarios (start j ∈ {1, 6, 11, 16, 21, 26}, η, D on 10 quantiles).
# Rows: results/f5b1d_partial/runs.csv (resumable, `budget` seconds).
#
# Usage, from the repository root:
#   julia --project=julia --pkgimages=no --heap-size-hint=1200M julia/f5b1d_partial.jl
#         [first_seed=201] [last_seed=208] [budget_s=360] [every=20] [m_max=40] [cases=a,b,…]
using Printf, Random, Statistics, MAT
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop
const CL = ClosedLoop

arg(i, d) = length(ARGS) >= i ? ARGS[i] : d
S0 = parse(Int, arg(1, "201")); S1 = parse(Int, arg(2, "208"))
BUDGET = parse(Float64, arg(3, "360")); EVERY = parse(Int, arg(4, "20"))
MMAX = parse(Float64, arg(5, "40"))
ONLY = length(ARGS) >= 6 ? split(ARGS[6], ",") : nothing           # optional subset of case names
const NQ = 481; const TC = 220.0; const RAMP = 5.0; const TCAP = 400.0
const OUT = joinpath(@__DIR__, "..", "results", "f5b1d_partial"); mkpath(OUT)
const FILE = joinpath(OUT, MMAX == 40 ? "runs.csv" : "runs_m$(Int(MMAX)).csv")

db = matread(joinpath(@__DIR__, "..", "models", "MCDET_conferencepaper_database.mat"))
const LAM = vec(db["Lambda_t"]) .* 1e5; const MU = vec(db["Mu_t"]) .* 1e4
const LS = sum(LAM); const DIST = (w = LAM ./ LS, mu = MU)
const ETAS = [0.75, 0.5, 0.25, 0.0]; const PETA = [0.4, 0.3, 0.2, 0.1]
const DQ10 = [CL.failure_quantile((i - 0.5) / 10, DIST) for i in 1:10]
const J = [1, 6, 11, 16, 21, 26]
base = (etas = ETAS, peta = PETA, mmax = MMAX)
cases = [(name = "N2", kw = NamedTuple()),
         (name = "N2_total_a0.02", kw = (risk = (tau = CL.failure_quantile(0.98, DIST), lam = 1e3),)),
         (name = "N2_vareta_a0.05", kw = (risk = (mode = :var_eta, alpha = 0.05, dist = DIST, lam = 1e3, base...),)),
         (name = "N2_vareta_a0.01", kw = (risk = (mode = :var_eta, alpha = 0.01, dist = DIST, lam = 1e3, base...),)),
         (name = "N2_vareta_a0.002", kw = (risk = (mode = :var_eta, alpha = 0.002, dist = DIST, lam = 1e3, base...),)),
         (name = "N2_scencvar_b0.9", kw = (risk = (mode = :scen_cvar, beta = 0.90, J, Dq = DQ10, lam = 1e3, base...),)),
         (name = "N2_scencvar_b0.98", kw = (risk = (mode = :scen_cvar, beta = 0.98, J, Dq = DQ10, lam = 1e3, base...),)),
         (name = "N2_scenexp_r100", kw = (risk = (mode = :scen_exp, rho = 100.0, J, Dq = DQ10, base...),)),
         (name = "N2_scenexp_r1000", kw = (risk = (mode = :scen_exp, rho = 1000.0, J, Dq = DQ10, base...),))]

fmt(x) = x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
done = Set{Tuple{Int,String}}()
if isfile(FILE)
    L = readlines(FILE); hdr = split(L[1], ","); is = findfirst(==("seed"), hdr); ic = findfirst(==("case"), hdr)
    for l in L[2:end]; f = split(l, ","); push!(done, (parse(Int, f[is]), f[ic])); end
end

function rb(t, taus)          # taus[η index][time index]
    g = [LS * exp(-LS * ti) * sum(pe * (τ[k] >= TCAP ? 0.0 : CL.failure_surv(τ[k], DIST)) for (pe, τ) in zip(PETA, taus))
         for (k, ti) in enumerate(t)]
    k = findall(<=(TC + 1e-9), t); tk, gk = t[k], g[k]
    sum((tk[2:end] .- tk[1:end-1]) .* (gk[2:end] .+ gk[1:end-1]) ./ 2) / (1 - exp(-LS * TC))
end

Bc = setup(); run_loop(Bc, :N2; n = 2, bias = true, keep_states = true,
                       risk = (mode = :scen_cvar, beta = 0.9, J, Dq = DQ10, lam = 1e3, base...))
t0 = time(); nrun = 0
for seed in S0:S1
    todo = [c for c in cases if !((seed, c.name) in done) && (ONLY === nothing || c.name in ONLY)]
    isempty(todo) && continue
    rng = MersenneTwister(seed)                                     # f5b1c draw
    A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng); lead = 17rand(rng)
    pH = CL.HRMQueries.hrm_params(; A, B = Bp, mdot_c_drift = d, acs_max = MMAX, acs_ramp = RAMP)
    B = setup(; pH, lead); B === nothing && continue
    for c in todo
        if time() - t0 > BUDGET
            @printf("budget reached after %d runs (%.0f s); run again to continue\n", nrun, time() - t0); exit(0)
        end
        tr = @elapsed lg = run_loop(B, :N2; n = NQ, act = act_cap(MMAX, RAMP), bias = true, keep_states = true, c.kw...)
        idx = unique(vcat(1:EVERY:length(lg.t_q), length(lg.t_q)))
        tt = @elapsed taus = [[time_to_limit(pH, lg.states[i][1], lg.states[i][2], lg.t_q[i]; t_cap = TCAP, flow = e * MMAX)
                               for i in idx] for e in ETAS]
        t = lg.t_q[idx] .- lg.t_q[1]
        M = metrics(lg, B.y_ref)
        row = (; seed, case = c.name, n = M.n, p_cross_cond = rb(t, taus),
               p_total_only = rb(t, [i == 4 ? taus[4] : fill(TCAP, length(t)) for i in 1:4]) / PETA[4],
               tau0_mean = mean(taus[4]), pct_max = M.pct_max, ise = M.ise, mass = M.mass_plant,
               mean_dev = mean(lg.pct .- B.y_ref), t_p99 = M.t_p99, fails = M.fails, wall_loop = tr, wall_tau = tt)
        nf = !isfile(FILE)
        open(FILE, "a") do io
            nf && println(io, join(keys(row), ","))
            println(io, join(fmt.(values(row)), ","))
        end
        global nrun += 1
        @printf("seed %d %-18s P(cross|fail) %.4f  K below y_ref %6.2f  mass %5.0f  p99 solve %.2f s fails %d (%.0f + %.0f s)\n",
                seed, c.name, row.p_cross_cond, -row.mean_dev, M.mass_plant, M.t_p99, M.fails, tr, tt); flush(stdout)
    end
end
@printf("complete: %d runs in %.0f s\n", nrun, time() - t0)
