# F7: certificate of global optimality for a subset of closed-loop MPC queries (F7Cert.jl).
#
# For each seed: re-run the F6 closed loop of N3 (same task as the F6 campaign, so the same
# plant, leak and failures), capture the queries that F6 audited and that are "hard" (best
# plan with J > 0.01, i.e. not trivially certified by J ≥ 0), and for each query:
#   1. decisions of the loop (N3), N1 and N2, scored with the NLP's cost J (implicit Euler);
#   2. upper bounds: Ipopt multistart (4 structured + F7_NRAND random feasible starts) and,
#      if a licence is found, Knitro (algorithms 1–4 and its own multistart);
#   3. lower bound: spatial branch and bound with rigorous interval enclosures (F7_TLIMIT s).
#
# Usage, from the repository root:
#   julia --project=julia julia/f7_certify.jl <campaign: nom|det|ext> <seeds, e.g. 1,3,5> [nk]
# Environment: F7_TLIMIT (B&B seconds per query, default 900), F7_NRAND (default 20),
#              F7_KNITRO=0 to skip Knitro, F7_NK (queries per seed, default 3).
# Output: results/f7_cert/<campaign>_<seeds>.csv (one row per query), *_hist.csv (B&B
#         history), *_local.csv (every local solve), queries/*.jls (serialised queries).
include(joinpath(@__DIR__, "F6Tasks.jl"))
using .F6Tasks
include(joinpath(@__DIR__, "F7Cert.jl"))
using .F7Cert, Printf, Serialization, Dates, Logging
global_logger(ConsoleLogger(stderr, Logging.Error))   # interval NaI warnings of failed enclosures
const CL = F6Tasks.CL; const F4 = F6Tasks.F4; const FC = F7Cert

camp = ARGS[1]
seeds = parse.(Int, split(ARGS[2], ","))
NK = parse(Int, get(ENV, "F7_NK", "3"))
TLIM = parse(Float64, get(ENV, "F7_TLIMIT", "900"))
NRAND = parse(Int, get(ENV, "F7_NRAND", "20"))
ROOT = joinpath(@__DIR__, "..")
out = joinpath(ROOT, "results", "f7_cert"); mkpath(joinpath(out, "queries"))
tag = "$(camp)_$(join(seeds, "-"))"
use_kn = get(ENV, "F7_KNITRO", "1") == "1" && FC.knitro_init()
println("F7 $tag  Knitro: ", use_kn ? "yes" : "no", "  B&B limit $(TLIM) s  ", now())

# ---- audited hard queries of F6 (N3 loops)
function hard_ks(key)
    lines = readlines(joinpath(ROOT, "results", "f6_paper", "$(camp)_audit.csv"))
    hd = split(lines[1], ","); ik = findfirst(==("k"), hd); ij = findfirst(==("Jbest"), hd)
    ks = Int[]
    for l in lines[2:end]
        f = split(l, ",")
        f[1] == key && parse(Float64, f[ij]) > 0.01 && push!(ks, parse(Int, f[ik]))
    end
    sort!(ks)
    length(ks) <= NK && return ks
    unique([ks[round(Int, 1 + (length(ks) - 1) * i / (NK - 1))] for i in 0:NK-1])
end

csvrow(io, nt) = println(io, join((x isa AbstractFloat ? @sprintf("%.10g", x) : string(x) for x in values(nt)), ","))
function writerows(path, rows)
    isempty(rows) && return
    open(path, "w") do io
        println(io, join(keys(rows[1]), ","))
        foreach(r -> csvrow(io, r), rows)
    end
end

rows = Any[]; hrows = Any[]; lrows = Any[]
for s in seeds
    T = F6Tasks.campaign_tasks(camp; seeds = [s])
    T = filter(t -> t.case == "N3b" && (camp != "nom" || t.act == :cap40), T)
    length(T) == 1 || error("expected one N3b task for $camp seed $s, got $(length(T))")
    t = merge(T[1], (audit = false,))
    ks = hard_ks(t.key)
    isempty(ks) && (println("  $(t.key): no hard audited query"); continue)
    cap = Dict{Int,Any}(); cnt = Ref(0)
    F6Tasks.CAPTURE[] = (tk, c, L, act, prev, dec) -> begin
        cnt[] += 1
        cnt[] in ks && (cap[cnt[]] = (c = c, L = copy(L), act = act, prev = prev === nothing ? nothing : copy(prev), m = copy(dec.m)))
    end
    tl = @elapsed F6Tasks.run_task(t)
    F6Tasks.CAPTURE[] = nothing
    println("  $(t.key): loop $(round(tl, digits = 1)) s, queries $(sort(collect(keys(cap))))")
    for k in sort(collect(keys(cap)))
        q = cap[k]
        serialize(joinpath(out, "queries", "$(camp)_$(s)_k$(k).jls"), q)
        P = FC.problem(q.c, q.L, q.act; S = F4.SMDFRROM, F4 = F4)
        n1 = CL.decide(:N1, q.c, q.L, q.act, q.prev); n2 = CL.decide(:N2, q.c, q.L, q.act, q.prev)
        # decisions as the solvers returned them violate the ramp by up to ~1e-7: repaired
        # (exactly feasible) points are scored; the raw violation is reported
        ul, u1, u2 = FC.repair(P, q.m), FC.repair(P, n1.m), FC.repair(P, n2.m)
        Jl, J1, J2 = FC.J(P, ul), FC.J(P, u1), FC.J(P, u2)
        viol = max(FC.infeas(P, q.m), FC.infeas(P, n1.m), FC.infeas(P, n2.m))
        cold = q.act.ramp == :hard ? q.c.m_prev .+ 5.0 .* (1:P.K) : fill(min(q.act.m_hi, 50.0), P.K)
        ti = @elapsed ms = FC.ipopt_multistart(P, [q.m, n1.m, n2.m, cold]; n_rand = NRAND, seed = s * 1000 + k)
        okms = [r for r in ms if r.inf <= 1e-5]
        Jip = minimum(r.J for r in okms; init = Inf)
        for r in ms
            push!(lrows, (; key = t.key, k, solver = "ipopt", run = "start$(r.start)", J = r.J, status = r.status,
                          time = r.time, inf = r.inf, moved = r.moved, u1 = r.u[1]))
        end
        kn = use_kn ? FC.knitro_runs(P, q.m) : NamedTuple[]
        for r in kn
            push!(lrows, (; key = t.key, k, solver = "knitro", run = r.run, J = r.J, status = r.status,
                          time = r.time, inf = r.inf, moved = r.moved, u1 = r.u[1]))
        end
        Jkn = minimum((r.J for r in kn if r.inf <= 1e-5); init = Inf)
        cands = [(Jl, ul), (J1, u1), (J2, u2), [(r.J, r.u) for r in ms]..., [(r.J, r.u) for r in kn if isfinite(r.J)]...]
        cands = [(Jc, u) for (Jc, u) in cands if FC.infeas(P, u) <= 1e-9]
        UB0, i0 = findmin(first.(cands)); u0 = cands[i0][2]
        B = FC.bnb(P, UB0, u0; tlimit = TLIM)
        for (tt, lb, ub) in B.hist
            push!(hrows, (; key = t.key, k, t = tt, LB = lb, UB = ub))
        end
        LB, UB = B.LB, B.UB
        rel(x) = UB > 0 ? x / UB : 0.0
        row = (; key = t.key, campaign = camp, seed = s, k, ramp = string(q.act.ramp),
               box_width_mean = sum(P.hi .- P.lo) / P.K, raw_violation = viol,
               J_loopN3 = Jl, J_N1 = J1, J_N2 = J2, J_ipopt_ms = Jip,
               n_ipopt = length(ms), n_ipopt_ok = length(okms),
               n_ipopt_distinct = length(unique(round.(sort([r.J for r in okms]), sigdigits = 7))),
               J_knitro = Jkn, n_knitro = length(kn), UB, LB, root_LB = B.root_lb,
               gap_abs = UB - LB, gap_rel = rel(UB - LB),
               cert_loopN3_abs = Jl - LB, cert_N1_abs = J1 - LB, cert_N2_abs = J2 - LB,
               cert_loopN3_rel = rel(Jl - LB), cert_N1_rel = rel(J1 - LB), cert_N2_rel = rel(J2 - LB),
               rms_gap_K = sqrt(max(UB - LB, 0) / P.K),
               bnb_status = B.status, nodes = B.nodes, open = B.open, failed_enclosures = B.failed,
               t_ipopt = ti, t_bnb = B.time)
        push!(rows, row)
        @printf("    k=%3d  J: N3 %.6g  N1 %.6g  N2 %.6g  ipopt-ms %.6g  knitro %.6g | LB %.6g (root %.6g)  gap %.2e (%.3f %%)  %s, %d nodes, %.0f s\n",
                k, Jl, J1, J2, Jip, Jkn, LB, B.root_lb, UB - LB, 100 * rel(UB - LB), B.status, B.nodes, B.time)
        writerows(joinpath(out, "$tag.csv"), rows)
        writerows(joinpath(out, "$(tag)_hist.csv"), hrows)
        writerows(joinpath(out, "$(tag)_local.csv"), lrows)
    end
end
println("done ", now())
