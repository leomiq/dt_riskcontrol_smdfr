# F5b.1b0 (plan §8.12 step 3.i.0, §8.15): sweep of the DSM penalty weight λ_dsm as in the
# thesis (metaheuristic loop: `matlab/fitness_cost_dsm.m`, results/results_penalty_saturation
# and results/results_retuned_lambda1e4), N2+DSM against the DET PI, with the thesis indices
# (`thesis_metrics`: IAE, reach / settling time in a 0.30 K band, steady-state error, minimum
# overshoot, minimum safety margin, peak PCT, control total variation, …).
#
# Penalty (thesis form): λ_dsm Σ_k [DSM_k − thr]₊², DSM_k = (ŷ_k − y_ref)/(U − y_ref),
# thr = 0.00237737, U ∈ {1244 (thesis), 1220 (DET end state)} K,
# λ_dsm ∈ {0, 1e2, 1e3, 1e4, 1e5, 1e6} (1e4 was the retuned thesis value, 1e6 over-penalised).
#
# Sets (common random numbers inside each set):
#   archived  the archived realisation (the thesis grid: first query 1160 s), n = 61 (thesis)
#             and n = 241; actuator `arch` (pump box 56.69–56.75 A → ≤ 12.86 kg/s) and `wide`
#             (≤ 56.90 A → 43.55 kg/s, the thesis "widened" bound), soft ramp as in the thesis;
#             the plant is capped at the same capacity.
#   nominal   nominal family of f5a2_monte_carlo.jl (odd seeds 1…39), no perturbations, n = 121,
#             the same two actuators.
#   detacc    DET-accelerated failure family of f5b1c_risk.jl (seeds 201…230), 40 kg/s, hard
#             ramp 5 kg/s per step, n = 481; λ_dsm ∈ {0, 1e3, 1e4, 1e5, 1e6}.
# PID: the THESIS PID (`matlab/pid_controller.m`: incremental, K_p = −0.5, K_i = −0.4 A/K,
# setpoint y_ref, [56.68, 100] A, 60 A/s; `ctrl = :PID`), on the plant capped at the same
# capacity (`PID`) and without cap (`PID_free`, as in the thesis). PI: the DET's internal PI
# (setpoint 1144 K, `ctrl = :PI`), same two variants (`PI`, `PI_free`). detacc: capped only.
# Rows: results/f5b1b0_dsm_lambda/runs.csv; trajectories of the archived set: traj_archived.csv.
#
# Usage, from the repository root:
#   julia --project=julia --pkgimages=no --heap-size-hint=1200M julia/f5b1b0_dsm_lambda.jl
#         [set=archived|nominal|detacc|all] [budget_s=100000]
using Printf, Random, Statistics, MAT
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop
const CL = ClosedLoop

SET = length(ARGS) >= 1 ? ARGS[1] : "all"
BUDGET = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 100000.0
const OUT = joinpath(@__DIR__, "..", "results", "f5b1b0_dsm_lambda"); mkpath(OUT)
const FILE = joinpath(OUT, "runs.csv"); const TRAJ = joinpath(OUT, "traj_archived.csv")
const THR = 0.00237737
const M_WIDE = CL.F4Lib.SMDFRROM.current_to_massflow(56.90)

fmt(x) = x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
function append(file, row)
    nf = !isfile(file)
    open(file, "a") do io
        nf && println(io, join(keys(row), ","))
        println(io, join(fmt.(values(row)), ","))
    end
end
done = Set{String}()
if isfile(FILE)
    L = readlines(FILE); hdr = split(L[1], ","); ik = findfirst(==("key"), hdr)
    for l in L[2:end]; push!(done, split(l, ",")[ik]); end
end

cases(lams) = vcat([(name = "N2", ctrl = :N2, lam = 0.0, U = NaN)],
                   [(name = "N2_DSM_U$(Int(U))_l$(Int(l))", ctrl = :N2, lam = l, U = U)
                    for U in (1244.0, 1220.0) for l in lams])

# ---------------- DET-accelerated failures (same as f5b1c_risk.jl)
db = matread(joinpath(@__DIR__, "..", "models", "MCDET_conferencepaper_database.mat"))
const LAM = vec(db["Lambda_t"]) .* 1e5; const MU = vec(db["Mu_t"]) .* 1e4
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
    for (a, b) in iv
        !isempty(out) && a <= out[end][2] ? (out[end] = (out[end][1], max(out[end][2], b))) : push!(out, (a, b))
    end
    out
end
function failures(seed, t_act, t_start)
    rng = MersenneTwister(50_000 + seed)
    while true
        w = ctmc_windows(rng, t_act, t_start + 240.0)
        any(t_start <= a <= t_start + 220.0 for (a, _) in w) && return w
    end
end

t0 = time(); nrun = 0
function run_case!(set, seed, B, act, actname, c, n; pert = NamedTuple(), traj = false, ctrl_kw = NamedTuple())
    key = "$(set)|$(seed)|$(actname)|$(c.name)|$(n)"
    key in done && return
    if time() - t0 > BUDGET
        @printf("budget reached after %d runs (%.0f s); run again to continue\n", nrun, time() - t0); exit(0)
    end
    dsm = c.lam > 0 ? (lam = c.lam, thr = THR, U = c.U) : nothing
    tr = @elapsed lg = c.ctrl == :N2 ? run_loop(B, :N2; n, act, bias = true, dsm, pert...) :
                                      run_loop(B, c.ctrl; n, act, pert..., ctrl_kw...)
    M = metrics(lg, B.y_ref); T = thesis_metrics(lg, B.y_ref)
    append(FILE, (; key, set, seed, act = actname, case = c.name, lam = c.lam, U = c.U, n_req = n,
                  n = M.n, limit = M.limit, ise = M.ise, fails = M.fails, T..., wall = tr))
    if traj
        nf = !isfile(TRAJ)
        open(TRAJ, "a") do io
            nf && println(io, "act,case,n_req,t,pct,m,m_plant")
            for j in eachindex(lg.t_q)
                println(io, join(fmt.((actname, c.name, n, lg.t_q[j], lg.pct[j], lg.m[j], lg.m_plant[j])), ","))
            end
        end
    end
    push!(done, key); global nrun += 1
    @printf("%-8s %4d %-5s %-22s n=%3d IAE %8.1f settle %6.1f reach %6.1f ss %6.2f peak %7.2f TV %7.1f mass %6.0f (%.0f s)\n",
            set, seed, actname, c.name, M.n, T.IAE, T.settling_time, T.reach_time, T.steady_state_error,
            T.peak_PCT, T.total_variation, T.coolant_mass, tr); flush(stdout)
end

Bc = setup(); run_loop(Bc, :N2; n = 2, bias = true, dsm = (lam = 1e4, thr = THR, U = 1244.0))

LAMS = (1e2, 1e3, 1e4, 1e5, 1e6)
if SET in ("archived", "all")
    for (actname, mh) in (("arch", NaN), ("wide", M_WIDE))
        cap = isnan(mh) ? nothing : (acs_max = mh, acs_ramp = Inf)
        Ba = setup(; cap = cap === nothing ? (acs_max = Bc.act_arch.m_hi, acs_ramp = Inf) : cap)
        Bf = setup()                                                   # PI without cap
        act = isnan(mh) ? Ba.act_arch : (name = :wide, m_lo = Ba.act_arch.m_lo, m_hi = mh, ramp = :soft)
        for n in (61, 241)
            for c in cases(LAMS); run_case!("archived", 0, Ba, act, actname, c, n; traj = true); end
            run_case!("archived", 0, Ba, act, actname, (name = "PI", ctrl = :PI, lam = 0.0, U = NaN), n; traj = true)
            run_case!("archived", 0, Bf, act, actname, (name = "PI_free", ctrl = :PI, lam = 0.0, U = NaN), n; traj = true)
            run_case!("archived", 0, Ba, act, actname, (name = "PID", ctrl = :PID, lam = 0.0, U = NaN), n; traj = true)
            run_case!("archived", 0, Bf, act, actname, (name = "PID_free", ctrl = :PID, lam = 0.0, U = NaN), n; traj = true)
        end
    end
end
if SET in ("nominal", "all")
    for seed in 1:2:39, (actname, mh) in (("arch", NaN), ("wide", M_WIDE))
        rng = MersenneTwister(seed)                                     # f5a2 standard draw
        A = 0.2 + 0.2rand(rng); Bp = 9 + 6rand(rng); d = 8 + 2rand(rng)
        mhi = isnan(mh) ? Bc.act_arch.m_hi : mh
        pH = CL.HRMQueries.hrm_params(; A, B = Bp, mdot_c_drift = d, acs_max = mhi)
        B = setup(; pH, lead = 17.0); B === nothing && continue
        Bf = setup(; pH = CL.HRMQueries.hrm_params(; A, B = Bp, mdot_c_drift = d), lead = 17.0)
        act = isnan(mh) ? B.act_arch : (name = :wide, m_lo = B.act_arch.m_lo, m_hi = mh, ramp = :soft)
        for c in cases(LAMS); run_case!("nominal", seed, B, act, actname, c, 121); end
        run_case!("nominal", seed, B, act, actname, (name = "PI", ctrl = :PI, lam = 0.0, U = NaN), 121)
        run_case!("nominal", seed, Bf, act, actname, (name = "PI_free", ctrl = :PI, lam = 0.0, U = NaN), 121)
        run_case!("nominal", seed, B, act, actname, (name = "PID", ctrl = :PID, lam = 0.0, U = NaN), 121)
        run_case!("nominal", seed, Bf, act, actname, (name = "PID_free", ctrl = :PID, lam = 0.0, U = NaN), 121)
    end
end
if SET in ("detacc", "all")
    for seed in 201:230
        rng = MersenneTwister(seed)                                     # f5b1c draw
        A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng); lead = 17rand(rng)
        pH = CL.HRMQueries.hrm_params(; A, B = Bp, mdot_c_drift = d, acs_max = 40.0, acs_ramp = 5.0)
        B = setup(; pH, lead); B === nothing && continue
        pert = (fail = failures(seed, B.t_act, B.t_start),)
        act = act_cap(40, 5)
        for c in cases((1e3, 1e4, 1e5, 1e6)); run_case!("detacc", seed, B, act, "cap40", c, 481; pert); end
        run_case!("detacc", seed, B, act, "cap40", (name = "PI", ctrl = :PI, lam = 0.0, U = NaN), 481; pert)
        run_case!("detacc", seed, B, act, "cap40", (name = "PID", ctrl = :PID, lam = 0.0, U = NaN), 481; pert)
    end
end
@printf("complete: %d runs in %.0f s\n", nrun, time() - t0)
