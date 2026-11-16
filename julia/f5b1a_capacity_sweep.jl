# F5b.1a (plan §8.12, step 1): sweep of the ACS capacity m_max and hard ramp on the extreme
# family of F5a.3, N2+bias against the DET's PI limited to the SAME capacity.
#
# The capacity is a PLANT parameter (`hrm_params(…; acs_max, acs_ramp)`): it limits whatever
# drives the pump — the DET PI before the first MPC query and as a baseline, and the MPC flow.
# The ramp limits increases only (a failure stops the pump at once; the recovery ramps up).
# The MPC knows the same bounds (`act_cap(m_max, dm)`), but not the failure or η (as in F5a.3).
# η applies to the PI too (`eta_pi = true`; in F5a.3 the PI was not degraded).
#
# Realisations: those of `f5a2_monte_carlo.jl … extreme` with the same seeds (same draw).
# Per seed: PI without limit (reference, once), and for every (m_max, ramp): PI capped, N2+bias.
#
# Rows are appended to results/f5b1a_capacity/runs.csv as soon as each run ends; rows already
# there are skipped, so the sweep can be run in short foreground batches (`budget` seconds).
#
# Usage, from the repository root:
#   julia --project=julia --heap-size-hint=1200M julia/f5b1a_capacity_sweep.jl [first_seed=101]
#         [last_seed=114] [budget_s=480] [n=481] [mmax=15,20,25,30,40,50] [ramps=2.5,5,10]
using Printf, Random
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop

arg(i, d) = length(ARGS) >= i ? ARGS[i] : d
S0 = parse(Int, arg(1, "101")); S1 = parse(Int, arg(2, "114"))
BUDGET = parse(Float64, arg(3, "480")); NQ = parse(Int, arg(4, "481"))
MMAX = parse.(Float64, split(arg(5, "15,20,25,30,40,50"), ","))
RAMPS = parse.(Float64, split(arg(6, "2.5,5,10"), ","))
const OUT = joinpath(@__DIR__, "..", "results", "f5b1a_capacity"); mkpath(OUT)
const FILE = joinpath(OUT, "runs.csv"); const TRAJ = joinpath(OUT, "traj.csv")   # gzipped for the repository (traj.csv.gz)
const TOL = 1e-3

"""Extreme-family draw of f5a2_monte_carlo.jl (identical RNG sequence)."""
function draw(seed)
    rng = MersenneTwister(seed)
    A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng)
    eta = 0.6 + 0.3rand(rng); fs = 5 + 15rand(rng); fl = 30 + 60rand(rng); fb = -0.5 + 0.6rand(rng)
    lead = 17rand(rng)
    (; A, Bp, d, eta, fs, fl, fb, lead)
end

"""Fractions of the queries whose first action sits on each bound: magnitude (hi), ramp-up,
ramp-down, zero; `free` = none of them (interior decision)."""
function saturation(m, m_hi, dm)
    Δ = diff(vcat(0.0, m)); n = length(m)
    hi = m .>= m_hi - TOL; lo = m .<= TOL
    up = .!hi .& (Δ .>= dm - TOL); dn = .!lo .& (Δ .<= -dm + TOL)
    free = .!(hi .| lo .| up .| dn)
    (sat_hi = count(hi) / n, sat_up = count(up) / n, sat_dn = count(dn) / n,
     sat_lo = count(lo) / n, free = count(free) / n)
end

fmt(x) = x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
function append(file, row)
    nf = !isfile(file)
    open(file, "a") do io
        nf && println(io, join(keys(row), ","))
        println(io, join(fmt.(values(row)), ","))
    end
end
done = Set{Tuple{Int,Float64,Float64,String}}()
if isfile(FILE)
    lines = readlines(FILE); hdr = split(lines[1], ",")
    is, im, ir, ic = (findfirst(==(c), hdr) for c in ("seed", "m_max", "ramp", "ctrl"))
    for l in lines[2:end]
        f = split(l, ","); push!(done, (parse(Int, f[is]), parse(Float64, f[im]), parse(Float64, f[ir]), f[ic]))
    end
end

# compile
Bc = setup(); run_loop(Bc, :N2; n = 2, act = act_cap(30, 5), bias = true)

t0 = time(); nrun = 0
configs = vcat([(Inf, Inf)], [(m, r) for m in MMAX for r in RAMPS])
for seed in S0:S1, (mmax, ramp) in configs
    ctrls = isinf(mmax) ? ["PI"] : ["PI", "N2"]
    all(c -> (seed, mmax, ramp, c) in done, ctrls) && continue
    if time() - t0 > BUDGET
        @printf("budget reached after %d runs (%.0f s); run again to continue\n", nrun, time() - t0); exit(0)
    end
    r = draw(seed)
    pH = ClosedLoop.HRMQueries.hrm_params(; A = r.A, B = r.Bp, mdot_c_drift = r.d,
                                           acs_max = mmax, acs_ramp = ramp)
    B = setup(; pH, lead = r.lead)
    B === nothing && (@printf("seed %d: no ACS demand, skipped\n", seed); continue)
    pert = (eta = r.eta, fail = (B.t_start + r.fs, B.t_start + r.fs + r.fl), fbias = r.fb)
    act = isinf(mmax) ? ACT_PHYS : act_cap(mmax, ramp)
    for c in ctrls
        (seed, mmax, ramp, c) in done && continue
        tr = @elapsed lg = run_loop(B, Symbol(c); n = NQ, act, bias = c == "N2", eta_pi = true, pert...)
        M = metrics(lg, B.y_ref)
        sat = saturation(lg.m, act.m_hi, isinf(ramp) ? Inf : ramp)
        row = (; seed, m_max = mmax, ramp, ctrl = c, A = r.A, B = r.Bp, drift = r.d,
               lead = B.t_start - B.t_act, t_act = B.t_act, eta = r.eta, fail_start = r.fs,
               fail_len = r.fl, fbias = r.fb, M..., margin = 1220 - M.pct_max, sat..., wall = tr)
        append(FILE, row)
        nf = !isfile(TRAJ)
        open(TRAJ, "a") do io
            nf && println(io, "seed,m_max,ramp,ctrl,t,pct,m,m_plant")
            for j in eachindex(lg.t_q)
                println(io, join(fmt.((seed, mmax, ramp, c, lg.t_q[j] - B.t_act, lg.pct[j], lg.m[j], lg.m_plant[j])), ","))
            end
        end
        global nrun += 1
        @printf("seed %d m_max %5.1f ramp %4.1f %-2s: limit %-5s PCTmax %7.2f free %.2f  (%.0f s)\n",
                seed, mmax, ramp, c, M.limit, M.pct_max, sat.free, tr); flush(stdout)
    end
end
@printf("sweep complete: %d runs in %.0f s\n", nrun, time() - t0)
