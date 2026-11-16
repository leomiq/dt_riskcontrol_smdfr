# F5b.1b (plan §8.12, step 3.i): deterministic margins, no probabilistic view.
#
# Extreme family (same seeds and draw as f5b1a_capacity_sweep.jl) on a plant with finite ACS
# capacity (m_max, ramp) from the F5b.1a map; η and the ACS failure act on every controller.
# Cases (common random numbers):
#   N2        N2+bias (reference)
#   N2_DSM    N2+bias + λ_s Σ[DSM_k − τ]₊², DSM_k = (ŷ_k − y_ref)/(U − y_ref),
#             U ∈ {1244, 1220} K, τ ∈ {0.00238 (thesis), 0.1, 0.25}, λ_s ∈ {1e3, 1e4, 1e5}
#   PI_1144   the DET PI (setpoint 1144 K, i.e. 35 K below y_ref = 1179.3 K)
#   PI_yref   the same PI with setpoint y_ref (same target as N2)
# Optional twin settings of F5b.0 (`twin = measured` → measured ACS flow in the warm-up;
# `twin = kalman` → also noisy leak σ = 0.1 kg/s and the Kalman forecaster).
# Rows: results/f5b1b_dsm/runs.csv (resumable, `budget` seconds).
#
# Usage, from the repository root:
#   julia --project=julia --pkgimages=no --heap-size-hint=1200M julia/f5b1b_dsm.jl
#         [m_max=40] [ramp=5] [first_seed=101] [last_seed=114] [budget_s=480] [twin=commanded]
using Printf, Random, Statistics
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop

arg(i, d) = length(ARGS) >= i ? ARGS[i] : d
MMAX = parse(Float64, arg(1, "40")); RAMP = parse(Float64, arg(2, "5"))
S0 = parse(Int, arg(3, "101")); S1 = parse(Int, arg(4, "114"))
BUDGET = parse(Float64, arg(5, "480")); TWIN = arg(6, "commanded"); NQ = 481
const OUT = joinpath(@__DIR__, "..", "results", "f5b1b_dsm"); mkpath(OUT)
const FILE = joinpath(OUT, "runs.csv")
twin_kw = TWIN == "commanded" ? NamedTuple() : TWIN == "measured" ? (measured_acs = true,) :
          TWIN == "kalman" ? (measured_acs = true, leak_noise = 0.1, forecast = :kalman) : error("twin")

function draw(seed)
    rng = MersenneTwister(seed)
    A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng)
    eta = 0.6 + 0.3rand(rng); fs = 5 + 15rand(rng); fl = 30 + 60rand(rng); fb = -0.5 + 0.6rand(rng)
    lead = 17rand(rng)
    (; A, Bp, d, eta, fs, fl, fb, lead)
end

cases = [(name = "N2", ctrl = :N2, kw = NamedTuple()),
         (name = "PI_1144", ctrl = :PI, kw = NamedTuple()),
         (name = "PI_yref", ctrl = :PI, kw = (pi_setpoint = :yref,))]
for U in (1244.0, 1220.0), thr in (0.00238, 0.1, 0.25), lam in (1e3, 1e4, 1e5)
    push!(cases, (name = "N2_DSM_U$(Int(U))_t$(thr)_l$(Int(lam))", ctrl = :N2,
                  kw = (dsm = (lam = lam, thr = thr, U = U),)))
end

fmt(x) = x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
done = Set{Tuple{Int,String,String}}()
if isfile(FILE)
    L = readlines(FILE); hdr = split(L[1], ",")
    is, ic, it = (findfirst(==(k), hdr) for k in ("seed", "case", "twin"))
    for l in L[2:end]; f = split(l, ","); push!(done, (parse(Int, f[is]), f[ic], f[it])); end
end

Bc = setup(); run_loop(Bc, :N2; n = 2, act = act_cap(30, 5), bias = true, dsm = (lam = 1e4, thr = 0.1, U = 1220.0))
t0 = time(); nrun = 0
for seed in S0:S1
    todo = [c for c in cases if !((seed, c.name, TWIN) in done)]
    isempty(todo) && continue
    r = draw(seed)
    pH = ClosedLoop.HRMQueries.hrm_params(; A = r.A, B = r.Bp, mdot_c_drift = r.d, acs_max = MMAX, acs_ramp = RAMP)
    B = setup(; pH, lead = r.lead)
    B === nothing && continue
    pert = (eta = r.eta, fail = (B.t_start + r.fs, B.t_start + r.fs + r.fl))
    for c in todo
        if time() - t0 > BUDGET
            @printf("budget reached after %d runs (%.0f s); run again to continue\n", nrun, time() - t0); exit(0)
        end
        kw = haskey(c.kw, :pi_setpoint) ? (pi_setpoint = B.y_ref,) : c.kw
        tw = c.ctrl == :N2 ? twin_kw : NamedTuple()
        tr = @elapsed lg = run_loop(B, c.ctrl; n = NQ, act = act_cap(MMAX, RAMP), bias = c.ctrl == :N2,
                                    eta_pi = true, noise_seed = seed, pert..., kw..., tw...)
        M = metrics(lg, B.y_ref)
        e = lg.pct .- B.y_ref
        row = (; seed, case = c.name, twin = TWIN, m_max = MMAX, ramp = RAMP, A = r.A, B = r.Bp, drift = r.d,
               eta = r.eta, fail_len = r.fl, M..., margin = 1220 - M.pct_max,
               mean_dev = mean(e), t_above_1200 = 0.5 * count(>(1200), lg.pct), wall = tr)
        nf = !isfile(FILE)
        open(FILE, "a") do io
            nf && println(io, join(keys(row), ","))
            println(io, join(fmt.(values(row)), ","))
        end
        global nrun += 1
        @printf("seed %d %-28s limit %-5s PCTmax %7.2f ISE %9.1f mass %6.0f (%.0f s)\n",
                seed, c.name, M.limit, M.pct_max, M.ise, M.mass_plant, tr); flush(stdout)
    end
end
@printf("complete: %d runs in %.0f s\n", nrun, time() - t0)
