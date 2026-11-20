# F7, Knitro part: re-solve the certified closed-loop queries with Knitro (needs a licence).
#
# Reads the queries serialised by f7_certify.jl (results/f7_cert/queries/*.jls) and the
# certificates (results/f7_cert/<campaign>_<seed>.csv), runs Knitro's four algorithms and its
# multistart from the loop's N3 decision, and writes results/f7_cert/knitro.csv: Knitro's best J
# against the Ipopt/B&B upper bound and the certified lower bound of the same query.
# A Knitro J below the certified LB would reveal an error (it cannot happen if both are right).
#
# Usage, from the repository root (KNITRO.jl in the default environment, licence found by
# Knitro, e.g. ARTELYS_LICENSE=<dir>):
#   julia --project=julia julia/f7_knitro.jl [pattern]   (pattern filters query files, e.g. ext_10)
include(joinpath(@__DIR__, "F6Tasks.jl"))
using .F6Tasks
include(joinpath(@__DIR__, "F7Cert.jl"))
using .F7Cert, Printf, Serialization
const F4 = F6Tasks.F4; const FC = F7Cert

FC.knitro_init() || error("Knitro is not usable (licence?)")
ROOT = joinpath(@__DIR__, ".."); out = joinpath(ROOT, "results", "f7_cert")
pat = length(ARGS) >= 1 ? ARGS[1] : ""
files = sort(filter(f -> endswith(f, ".jls") && occursin(pat, f), readdir(joinpath(out, "queries"))))

"""Certified (UB, LB) of a query from the f7_certify.jl tables."""
function cert_of(camp, seed, k)
    for f in readdir(out)
        (startswith(f, "$(camp)_") && endswith(f, ".csv") && !occursin("_hist", f) && !occursin("_local", f)) || continue
        L = readlines(joinpath(out, f)); hd = split(L[1], ",")
        ix(n) = findfirst(==(n), hd)
        for l in L[2:end]
            v = split(l, ",")
            v[ix("seed")] == string(seed) && v[ix("k")] == string(k) &&
                return (UB = parse(Float64, v[ix("UB")]), LB = parse(Float64, v[ix("LB")]))
        end
    end
    (UB = NaN, LB = NaN)
end

rows = String["file,run,J_knitro,status,time,raw_violation,UB_cert,LB_cert,excess_over_UB,below_LB"]
for f in files
    m = match(r"^(\w+)_(\d+)_k(\d+)\.jls$", f)
    camp, seed, k = m[1], parse(Int, m[2]), parse(Int, m[3])
    q = deserialize(joinpath(out, "queries", f))
    P = FC.problem(q.c, q.L, q.act; S = F4.SMDFRROM, F4 = F4)
    cert = cert_of(camp, seed, k)
    for r in FC.knitro_runs(P, q.m)
        push!(rows, join((f, r.run, @sprintf("%.10g", r.J), r.status, @sprintf("%.3g", r.time),
                          @sprintf("%.2e", r.inf), @sprintf("%.10g", cert.UB), @sprintf("%.10g", cert.LB),
                          @sprintf("%.3e", r.J - cert.UB), r.J < cert.LB - 1e-9), ","))
        @printf("%-22s %-6s J %.10g (UB %.10g, LB %.10g)  %s\n", f, r.run, r.J, cert.UB, cert.LB, r.status)
    end
    write(joinpath(out, "knitro.csv"), join(rows, "\n") * "\n")
end
