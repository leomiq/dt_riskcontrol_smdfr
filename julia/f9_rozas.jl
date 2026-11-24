# F9: the one-scenario result of Homem-de-Mello et al. (2026) (proposed by H. Rozas) checked on
# the SMDFR MPC with the paper's methodology (F9Rozas.jl has the formulations and experiments).
#
# Usage, from the repository root:
#   julia --project=julia julia/f9_rozas.jl A [workers=20] [first:last seeds] [every=6]
#       phase A: closed loops (N2 + expected excess ρ = 1000, R = 0.1, DET failures, 40 kg/s),
#       every `every`-th query analysed (E0–E3, P1 on all, P2/P3 on every 4th), resumable:
#       results/f9_rozas/A_<seed>.jls, then A_queries.csv (one row per analysed query)
#   julia --project=julia julia/f9_rozas.jl B
#       phase B (master only): contextual one-scenario forecasts Φ(θ, x) (paper §5–6,
#       Algorithm 1 + Nelder–Mead) trained on seeds ≤ 235, tested on 236–250:
#       B_results.csv, B_learning.csv, B_policies.jls
#   julia --project=julia julia/f9_rozas.jl C [workers=20] [first:last seeds]
#       phase C: closed loops with the learned one-scenario MPC vs SAA vs N2, Rao–Blackwell:
#       C_rb.csv (resumable)
# Environment: F9_OUT (output folder), F6_NQ (queries per loop; short smoke tests).
using Distributed, Printf, Dates, Serialization, Random, Statistics

PHASE = ARGS[1]
W = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20
SEEDS = length(ARGS) >= 3 ? (p = parse.(Int, split(ARGS[3], ":")); p[1]:p[end]) : (201:250)
EVERY = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 6
const OUT = get(ENV, "F9_OUT", joinpath(@__DIR__, "..", "results", "f9_rozas")); mkpath(OUT)

PHASE in ("A", "C") && W > 0 && addprocs(W; exeflags = ["--project=$(@__DIR__)", "--heap-size-hint=2000M"],
                                       env = ["OPENBLAS_NUM_THREADS" => "1", "OMP_NUM_THREADS" => "1"])
@everywhere include(joinpath($(@__DIR__), "F9Rozas.jl"))
@everywhere using .F9Rozas, Dates
const F9 = F9Rozas

fmt(x) = x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
function write_rows(file, rows)
    isempty(rows) && return
    cols = unique(reduce(vcat, [collect(keys(r)) for r in rows]))
    open(file, "w") do io
        println(io, join(cols, ","))
        for r in rows; println(io, join([fmt(get(r, c, NaN)) for c in cols], ",")); end
    end
end
ts() = Dates.format(now(), "HH:MM:SS")

if PHASE == "A"
    todo = [s for s in SEEDS if !isfile(joinpath(OUT, "A_$(s).jls"))]
    @printf("[%s] phase A: %d seeds to run (%d workers, every %d)\n", ts(), length(todo), nworkers(), EVERY)
    res = pmap(todo) do s
        r = try
            F9.phaseA(s, OUT; every = EVERY)
        catch err
            (seed = s, error = sprint(showerror, err))
        end
        println("[", Dates.format(now(), "HH:MM:SS"), "] ", r); flush(stdout); r
    end
    rows = Any[]
    for s in SEEDS
        f = joinpath(OUT, "A_$(s).jls"); isfile(f) || continue
        d = deserialize(f)
        isempty(d.errs) || println("seed $s: ", length(d.errs), " query errors, e.g. ", first(d.errs))
        append!(rows, [r.row for r in d.recs])
    end
    write_rows(joinpath(OUT, "A_queries.csv"), [Dict(string(k) => v for (k, v) in r) for r in rows])
    @printf("[%s] A_queries.csv: %d queries\n", ts(), length(rows))

elseif PHASE == "B"
    recs = Any[]
    for s in 201:250
        f = joinpath(OUT, "A_$(s).jls"); isfile(f) || continue
        append!(recs, deserialize(f).recs)
    end
    SPLIT = parse(Int, get(ENV, "F9_SPLIT", "235"))          # last training seed
    tr = [r for r in recs if r.row[:seed] <= SPLIT]; te = [r for r in recs if r.row[:seed] > SPLIT]
    @printf("[%s] phase B: %d train queries, %d test queries\n", ts(), length(tr), length(te))
    Xtr = reduce(hcat, [r.x for r in tr])'; μ = vec(mean(Xtr, dims = 1)); σ = max.(vec(std(Xtr, dims = 1)), 1e-9)
    G = F9.DGRID; δlo, δhi = G[1], G[end]
    itp(v, δ) = (δ = clamp(δ, δlo, δhi); x = (δ - δlo) / (G[2] - G[1]) + 1; i = min(floor(Int, x), length(G) - 1);
                 w = x - i; (1 - w) * v[i] + w * v[i+1])
    curve(r) = Float64.(r.curve)
    dmean = sum(F9.dist().w ./ F9.dist().mu); dmed = F9.CL.failure_quantile(0.5, F9.dist())
    results = Any[]; learning = Any[]; policies = Dict{Any,Any}()
    for ρ in F9.RHOS
        margin(p, r) = F9.scen_margin(p, curve(r), r.x)
        # Algorithm 1: forecast → plan (grid interpolation of the one-scenario QP) → realised cost
        loss(p, R) = mean(mean(itp(view(r.grids[ρ].R, s, :), margin(p, r)) for s in axes(r.grids[ρ].R, 1)) for r in R)
        excost(p, r) = itp(r.grids[ρ].Jex, margin(p, r))
        mk(kind, θ) = F9.Policy(kind, θ, μ, σ)
        nf = length(μ)
        fit(kind, θ0, step, R) = (F9.nelder_mead(θ -> loss(mk(kind, θ), R), θ0; step, tol = 1e-7, maxit = 4000, restarts = 3)[1])
        # variance-reduced variant (not in the paper): the realised cost replaced by its expectation
        # over D given the query (Rao–Blackwell, the distribution of D is known; 1000 quantiles)
        lossE(p, R) = mean(excost(p, r) for r in R)
        fitE(kind, θ0, step, R) = (F9.nelder_mead(θ -> lossE(mk(kind, θ), R), θ0; step, tol = 1e-7, maxit = 4000, restarts = 3)[1])
        # 1-D searches: the loss is not unimodal → grid, then golden section around the best node
        function search1(f, a, b, n)
            xs = range(a, b, length = n); fs = map(f, xs); i = argmin(fs)
            x, fx = F9.golden(f, xs[max(i - 1, 1)], xs[min(i + 1, n)]; it = 30)
            fx < fs[i] ? x : xs[i]
        end
        δc = search1(v -> loss(mk(:delta_const, [v]), tr), δlo, δhi, 151)
        Dc = search1(v -> loss(mk(:D_const, [v]), tr), 0.0, 2F9.DCAP, 161)
        tlin = @elapsed θd = fit(:delta_lin, vcat(δc, zeros(nf)), vcat(5.0, fill(2.0, nf)), tr)
        tlinD = @elapsed θD = fit(:D_lin, vcat(Dc, zeros(nf)), vcat(20.0, fill(10.0, nf)), tr)
        tq = @elapsed θq = fit(:delta_quad, vcat(θd, zeros(nf)), vcat(2.0, fill(1.0, 2nf)), tr)
        δcE = search1(v -> lossE(mk(:delta_const, [v]), tr), δlo, δhi, 151)
        DcE = search1(v -> lossE(mk(:D_const, [v]), tr), 0.0, 2F9.DCAP, 161)
        θdE = fitE(:delta_lin, vcat(δcE, zeros(nf)), vcat(5.0, fill(2.0, nf)), tr)
        θDE = fitE(:D_lin, vcat(DcE, zeros(nf)), vcat(20.0, fill(10.0, nf)), tr)
        pols = [("nominal (risk ignored)", mk(:delta_const, [δlo])),
                ("mean duration E[D] (point forecast)", mk(:D_const, [dmean])),
                ("median duration", mk(:D_const, [dmed])),
                ("constant margin (no features)", mk(:delta_const, [δc])),
                ("constant duration D̂ (no features)", mk(:D_const, [Dc])),
                ("AD linear Φ(θ,x) → margin", mk(:delta_lin, θd)),
                ("AD linear Φ(θ,x) → duration", mk(:D_lin, θD)),
                ("AD quadratic Φ(θ,x) → margin", mk(:delta_quad, θq)),
                ("[E-loss] constant duration D̂", mk(:D_const, [DcE])),
                ("[E-loss] linear Φ(θ,x) → margin", mk(:delta_lin, θdE)),
                ("[E-loss] linear Φ(θ,x) → duration", mk(:D_lin, θDE))]
        policies[ρ] = Dict(name => p for (name, p) in pols)
        JS = [r.grids[ρ].JS for r in te]; Jn = [itp(r.grids[ρ].Jex, δlo) for r in te]
        Jor = [minimum(r.grids[ρ].Jex) for r in te]
        for (name, p) in vcat(pols, [("oracle ξ*_x (best one scenario per query)", nothing)])
            J = p === nothing ? Jor : [excost(p, r) for r in te]
            u1 = p === nothing ? [r.grids[ρ].u1[argmin(r.grids[ρ].Jex)] for r in te] : [itp(r.grids[ρ].u1, margin(p, r)) for r in te]
            Jtr = p === nothing ? NaN : loss(p, tr)
            push!(results, Dict("rho" => ρ, "policy" => replace(name, "," => ";"), "J_test" => mean(J), "regret_vs_SAA" => mean(J .- JS),
                "regret_max" => maximum(J .- JS), "VSS_captured" => (sum(Jn) - sum(J)) / (sum(Jn) - sum(JS)),
                "u1_mae_vs_SAA" => mean(abs.(u1 .- [r.grids[ρ].uS1 for r in te])),
                "loss_train" => Jtr, "theta" => p === nothing ? "" : join(round.(p.θ; sigdigits = 4), " ")))
        end
        push!(results, Dict("rho" => ρ, "policy" => "SAA pi_S (200 quantiles; target)", "J_test" => mean(JS), "regret_vs_SAA" => 0.0,
            "regret_max" => 0.0, "VSS_captured" => 1.0, "u1_mae_vs_SAA" => 0.0, "loss_train" => NaN, "theta" => ""))
        @printf("[%s] ρ = %g: fitted (%.0f s, %.0f s, %.0f s); test VSS = %.4g per query\n", ts(), ρ, tlin, tlinD, tq, mean(Jn .- JS))
        # Theorem 2: the learned forecast converges as the training set grows (test cost)
        for Nt in (25, 100, 300, 1000, length(tr)), rep in 1:(Nt == length(tr) ? 1 : 3)
            Nt > length(tr) && continue
            sub = tr[randperm(MersenneTwister(100rep + Nt), length(tr))[1:Nt]]
            θs = fit(:delta_lin, vcat(δc, zeros(nf)), vcat(5.0, fill(2.0, nf)), sub)
            J = [excost(mk(:delta_lin, θs), r) for r in te]
            push!(learning, Dict("rho" => ρ, "N_train" => Nt, "rep" => rep, "regret_vs_SAA" => mean(J .- JS),
                                 "VSS_captured" => (sum(Jn) - sum(J)) / (sum(Jn) - sum(JS))))
        end
    end
    write_rows(joinpath(OUT, "B_results.csv"), results)
    write_rows(joinpath(OUT, "B_learning.csv"), learning)
    serialize(joinpath(OUT, "B_policies.jls"), (policies = policies, μ = μ, σ = σ, names = F9.FNAMES))
    for r in results; @printf("  ρ=%-5g %-45s regret %10.4g  VSS captured %6.3f\n", r["rho"], r["policy"], r["regret_vs_SAA"], r["VSS_captured"]); end

elseif PHASE == "B2"
    # E5 found a closed-loop instability of the contextual policies at ρ = 1000: the margin
    # depends on features that the decision itself moves within one query (margin, d, applied
    # flow, 10-s heat-up), i.e. the context is endogenous, which the paper's i.i.d. setting
    # ignores. Refit the expected-loss linear forecast on the slow features only.
    recs = Any[]
    for s in 201:250
        f = joinpath(OUT, "A_$(s).jls"); isfile(f) && append!(recs, deserialize(f).recs)
    end
    SPLIT = parse(Int, get(ENV, "F9_SPLIT", "235"))
    tr = [r for r in recs if r.row[:seed] <= SPLIT]; te = [r for r in recs if r.row[:seed] > SPLIT]
    B0 = deserialize(joinpath(OUT, "B_policies.jls")); μ, σ = B0.μ, B0.σ
    G = F9.DGRID; δlo, δhi = G[1], G[end]
    itp(v, δ) = (δ = clamp(δ, δlo, δhi); x = (δ - δlo) / (G[2] - G[1]) + 1; i = min(floor(Int, x), length(G) - 1);
                 w = x - i; (1 - w) * v[i] + w * v[i+1])
    slow = [findfirst(==(n), F9.FNAMES) for n in (:leak, :dT46, :dT150, :dT400)]
    pols = Dict{Any,Any}(); res = Any[]
    for ρ in F9.RHOS
        excost(p, r) = itp(r.grids[ρ].Jex, F9.scen_margin(p, Float64.(r.curve), r.x))
        full(v) = (θ = zeros(1 + length(μ)); θ[1] = v[1]; θ[1 .+ slow] = v[2:end]; F9.Policy(:delta_lin, θ, μ, σ))
        θ0 = B0.policies[ρ]["[E-loss] linear Φ(θ,x) → margin"].θ
        v, _ = F9.nelder_mead(v -> mean(excost(full(v), r) for r in tr), vcat(θ0[1], θ0[1 .+ slow]);
                              step = vcat(5.0, fill(2.0, length(slow))), tol = 1e-7, maxit = 4000, restarts = 3)
        p = full(v); pols[ρ] = p
        JS = [r.grids[ρ].JS for r in te]; Jn = [itp(r.grids[ρ].Jex, δlo) for r in te]; J = [excost(p, r) for r in te]
        push!(res, Dict("rho" => ρ, "policy" => "[E-loss] linear Φ(θ;x_slow) → margin", "regret_vs_SAA" => mean(J .- JS),
                        "VSS_captured" => (sum(Jn) - sum(J)) / (sum(Jn) - sum(JS)), "theta" => join(round.(p.θ; sigdigits = 4), " ")))
        @printf("ρ = %g: slow-feature policy captures %.3f of the VSS (test)\n", ρ, res[end]["VSS_captured"])
    end
    write_rows(joinpath(OUT, "B2_results.csv"), res)
    serialize(joinpath(OUT, "B2_policies.jls"), pols)

elseif PHASE == "C"
    P = deserialize(joinpath(OUT, "B_policies.jls")).policies
    P2 = isfile(joinpath(OUT, "B2_policies.jls")) ? deserialize(joinpath(OUT, "B2_policies.jls")) : nothing
    file = joinpath(OUT, "C_rb.csv")
    done = Set{String}()
    if isfile(file)
        L = readlines(file); ic = findfirst(==("case"), split(L[1], ",")); is = findfirst(==("seed"), split(L[1], ","))
        for l in L[2:end]; v = split(l, ","); push!(done, "$(v[is])|$(v[ic])"); end
    end
    tasks = Any[]
    for s in SEEDS
        push!(tasks, (s, "N2", NamedTuple()))
        for ρ in F9.RHOS
            r = Int(ρ)
            push!(tasks, (s, "SAA_r$r", (risk = (mode = :expected, rho = ρ, Dq = F9.T6.qgrid(40, F9.MUACC)),)))
            # E[D] as the scenario is not run: its margin never activates (= N2, phase B)
            for (tag, name) in (("scen_ADlin", "AD linear Φ(θ,x) → margin"),
                                ("scen_ADlinE", "[E-loss] linear Φ(θ,x) → margin"),
                                ("scen_ADlinDE", "[E-loss] linear Φ(θ,x) → duration"),
                                ("scen_Dconst", "constant duration D̂ (no features)"))
                push!(tasks, (s, "$(tag)_r$r", (risk = (mode = :scenario, fn = P[ρ][name], rho = ρ),)))
            end
            P2 === nothing || push!(tasks, (s, "scen_ADslowE_r$r", (risk = (mode = :scenario, fn = P2[ρ], rho = ρ),)))
        end
    end
    tasks = filter(t -> !("$(t[1])|$(t[2])" in done), tasks)
    @printf("[%s] phase C: %d loops to run (%d workers)\n", ts(), length(tasks), nworkers())
    lk = ReentrantLock()
    pmap(tasks) do t
        r = try
            F9.run_rb(t[1], t[2], t[3])
        catch err
            (seed = t[1], case = t[2], error = replace(sprint(showerror, err), "," => ";", "\n" => " "))
        end
        r === nothing && return nothing
        println("[", Dates.format(now(), "HH:MM:SS"), "] ", t[1], " ", t[2]); flush(stdout)
        r
    end |> res -> begin
        rows = [Dict(string(k) => v for (k, v) in pairs(r)) for r in res if r !== nothing && !haskey(r, :error)]
        errs = [r for r in res if r !== nothing && haskey(r, :error)]
        isempty(errs) || foreach(println, errs)
        old = isfile(file) ? readlines(file) : String[]
        tmp = file * ".new"; write_rows(tmp, rows)
        if !isempty(old) && isfile(tmp)          # append to the existing table (same columns)
            new = readlines(tmp); open(file, "a") do io; foreach(l -> println(io, l), new[2:end]); end; rm(tmp)
        elseif isfile(tmp)
            mv(tmp, file; force = true)
        end
    end
end
