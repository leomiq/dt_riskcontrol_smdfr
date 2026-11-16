# F5b.0 (plan §8.12, step 2): what the twin does not know, in closed loop.
#
# Extreme family of F5a.3 (same seeds and draw as f5b1a_capacity_sweep.jl), on a plant with
# finite ACS capacity (m_max, ramp) chosen from the F5b.1a map, with the η and ACS-failure
# perturbations of the family. Controller: N2+bias. Cases (common random numbers):
#   A  measured ACS flow in the ROM warm-up: commanded (as F5a) vs measured; oracle leak;
#   B  leak forecast without Park: oracle, persistence, Kalman (with / without jump
#      intervention), with measurement noise σ ∈ {0.1, 0.3} kg/s, without and with jumps of
#      the break (1–2 per run at U(0, 240) s after the first query, Δ ~ −U(1, 5) kg/s);
#      measured ACS flow on;
#   PI limited to the same capacity, without and with jumps (reference).
# Rows are appended to results/f5b0_closed_loop/runs.csv (resumable, `budget` seconds).
#
# Usage, from the repository root:
#   julia --project=julia --pkgimages=no --heap-size-hint=1200M julia/f5b0_closed_loop.jl
#         [m_max=40] [ramp=5] [first_seed=101] [last_seed=114] [budget_s=480] [n=481]
using Printf, Random, Statistics
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop

arg(i, d) = length(ARGS) >= i ? ARGS[i] : d
MMAX = parse(Float64, arg(1, "40")); RAMP = parse(Float64, arg(2, "5"))
S0 = parse(Int, arg(3, "101")); S1 = parse(Int, arg(4, "114"))
BUDGET = parse(Float64, arg(5, "480")); NQ = parse(Int, arg(6, "481"))
const OUT = joinpath(@__DIR__, "..", "results", "f5b0_closed_loop"); mkpath(OUT)
const FILE = joinpath(OUT, "runs.csv")

function draw(seed)                       # = f5b1a_capacity_sweep.jl / f5a2_monte_carlo.jl extreme
    rng = MersenneTwister(seed)
    A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng)
    eta = 0.6 + 0.3rand(rng); fs = 5 + 15rand(rng); fl = 30 + 60rand(rng); fb = -0.5 + 0.6rand(rng)
    lead = 17rand(rng)
    rj = MersenneTwister(10_000 + seed)   # jumps of the break: separate stream
    nj = rand(rj, 1:2); jumps = [(240rand(rj), -(1 + 4rand(rj))) for _ in 1:nj]
    (; A, Bp, d, eta, fs, fl, fb, lead, jumps)
end

cases = [(name = "A_commanded", ctrl = :N2, jumps = false, kw = (measured_acs = false,)),
         (name = "A_measured", ctrl = :N2, jumps = false, kw = (measured_acs = true,)),
         (name = "PI", ctrl = :PI, jumps = false, kw = NamedTuple()),
         (name = "PI_jumps", ctrl = :PI, jumps = true, kw = NamedTuple())]
for σ in (0.1, 0.3), j in (false, true), f in (:oracle, :persist, :kalman, :kalman_noint)
    push!(cases, (name = "B_$(f)_s$(σ)$(j ? "_jumps" : "")", ctrl = :N2, jumps = j,
                  kw = (measured_acs = true, leak_noise = σ, forecast = f)))
end

fmt(x) = x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
done = Set{Tuple{Int,String}}()
if isfile(FILE)
    L = readlines(FILE); hdr = split(L[1], ","); is = findfirst(==("seed"), hdr); ic = findfirst(==("case"), hdr)
    for l in L[2:end]; f = split(l, ","); push!(done, (parse(Int, f[is]), f[ic])); end
end

Bc = setup(); run_loop(Bc, :N2; n = 2, act = act_cap(30, 5), bias = true, measured_acs = true,
                       leak_noise = 0.1, forecast = :kalman)
t0 = time(); nrun = 0
for seed in S0:S1, jumps in (false, true)
    todo = [c for c in cases if c.jumps == jumps && !((seed, c.name) in done)]
    isempty(todo) && continue
    r = draw(seed)
    pH = ClosedLoop.HRMQueries.hrm_params(; A = r.A, B = r.Bp, mdot_c_drift = r.d, acs_max = MMAX,
                                           acs_ramp = RAMP, leak_jumps = Tuple{Float64,Float64}[])
    B = setup(; pH, lead = r.lead)
    B === nothing && continue
    if jumps        # jump times relative to the first query: rebuild the bench with them
        pH = ClosedLoop.HRMQueries.hrm_params(; A = r.A, B = r.Bp, mdot_c_drift = r.d, acs_max = MMAX,
                 acs_ramp = RAMP, leak_jumps = [(B.t_start + tj, dj) for (tj, dj) in r.jumps])
        B = setup(; pH, lead = r.lead)
    end
    pert = (eta = r.eta, fail = (B.t_start + r.fs, B.t_start + r.fs + r.fl))
    for c in todo
        if time() - t0 > BUDGET
            @printf("budget reached after %d runs (%.0f s); run again to continue\n", nrun, time() - t0); exit(0)
        end
        tr = @elapsed lg = run_loop(B, c.ctrl; n = NQ, act = act_cap(MMAX, RAMP), bias = c.ctrl == :N2,
                                    noise_seed = seed, eta_pi = true, pert..., c.kw...)
        M = metrics(lg, B.y_ref)
        fe = filter(!isnan, lg.ferr)
        row = (; seed, case = c.name, m_max = MMAX, ramp = RAMP, A = r.A, B = r.Bp, drift = r.d,
               eta = r.eta, fail_len = r.fl, njumps = jumps ? length(r.jumps) : 0,
               jump_sum = jumps ? sum(last, r.jumps) : 0.0, M..., margin = 1220 - M.pct_max,
               ferr_p50 = isempty(fe) ? NaN : median(fe), ferr_p99 = isempty(fe) ? NaN : quantile(fe, 0.99),
               ferr_max = isempty(fe) ? NaN : maximum(fe), wall = tr)
        nf = !isfile(FILE)
        open(FILE, "a") do io
            nf && println(io, join(keys(row), ","))
            println(io, join(fmt.(values(row)), ","))
        end
        global nrun += 1
        @printf("seed %d %-24s limit %-5s PCTmax %7.2f ISE %9.1f ferr p99 %.3f  (%.0f s)\n",
                seed, c.name, M.limit, M.pct_max, M.ise, row.ferr_p99, tr); flush(stdout)
    end
end
@printf("complete: %d runs in %.0f s\n", nrun, time() - t0)
