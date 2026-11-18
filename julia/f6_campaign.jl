# F6: paper campaign (plan §10, points 2 and 3) — parallel, resumable driver.
#
# Campaigns (task definitions in F6Tasks.jl):
#   sens  point 2: failure acceleration μ×{3e3, 1e4, 3e4} × ACS capacity {30,…,50} kg/s,
#         Rao–Blackwell P(cross | failure) of N2, PI, PID and the risk bounds (VaR α, CVaR 0.95)
#   form  point 3: every risk formulation of §8.17 on 50 realisations (40 kg/s, μ×1e4)
#   eta   point 2: severity probabilities of partial failures (§8.18), 4 p_η vectors
#   det   point 3: closed loop with DET-accelerated failures: N1/N2/N3, risk layers, DSM,
#         DET PI, thesis PID; thesis control metrics + per-query optimality audit
#   nom   point 3: nominal family (thesis setting), archived box and 40 kg/s; metrics + audit
#   ext   point 3: extreme family of F5a.3; metrics + audit
#
# One task = one closed-loop run. Workers pull tasks (most expensive first); the master
# appends every result to results/f6_paper/<campaign>*.csv as soon as it arrives, so the run
# can be stopped and restarted at any time (finished keys are skipped). Failed tasks are
# logged to <campaign>_errors.log and retried on the next start.
#
# Usage, from the repository root:
#   julia --project=julia julia/f6_campaign.jl <campaign>[,<campaign>…] [workers=30] [first:last seeds]
# e.g.  julia --project=julia julia/f6_campaign.jl sens,form,eta,det,nom,ext 30
#       julia --project=julia julia/f6_campaign.jl det 4 201:202        (smoke test)
using Distributed, Printf, Dates

CAMPS = split(ARGS[1], ",")
W = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 30
SEEDS = length(ARGS) >= 3 ? (p = parse.(Int, split(ARGS[3], ":")); p[1]:p[end]) : nothing
const OUT = get(ENV, "F6_OUT", joinpath(@__DIR__, "..", "results", "f6_paper")); mkpath(OUT)   # F6_OUT: smoke tests

W > 0 && addprocs(W; exeflags = ["--project=$(@__DIR__)", "--heap-size-hint=1500M"],
                  env = ["OPENBLAS_NUM_THREADS" => "1", "OMP_NUM_THREADS" => "1"])
@everywhere include(joinpath($(@__DIR__), "F6Tasks.jl"))
@everywhere using .F6Tasks

# ---- tasks still to do
fmt(x) = x isa AbstractFloat ? @sprintf("%.8g", x) : string(x)
function done_keys(file)
    s = Set{String}()
    isfile(file) || return s
    L = readlines(file); isempty(L) && return s
    ik = findfirst(==("key"), split(L[1], ","))
    for l in L[2:end]; push!(s, split(l, ",")[ik]); end
    s
end
tasks = Any[]
for c in CAMPS
    T = F6Tasks.campaign_tasks(String(c); seeds = SEEDS)
    done = union(done_keys(joinpath(OUT, "$(c).csv")), done_keys(joinpath(OUT, "$(c)_skipped.csv")))
    todo = filter(t -> !(t.key in done), T)
    @printf("%-5s %5d tasks, %5d done, %5d to run\n", c, length(T), length(T) - length(todo), length(todo))
    append!(tasks, todo)
end
sort!(tasks; by = t -> -t.cost)                       # long tasks first (load balance)
isempty(tasks) && (println("nothing to do"); exit(0))

# ---- writer (master only; one lock, whole rows)
const LK = ReentrantLock()
function write_out(out::Dict)
    lock(LK) do
        for (f, rows) in out
            isempty(rows) && continue
            file = joinpath(OUT, f); nf = !isfile(file) || filesize(file) == 0
            open(file, "a") do io
                nf && println(io, join(keys(rows[1]), ","))
                for r in rows; println(io, join(fmt.(values(r)), ",")); end
            end
        end
    end
end
function log_error(t, msg)
    lock(LK) do
        open(joinpath(OUT, "$(t.campaign)_errors.log"), "a") do io
            println(io, "[$(now())] $(t.key)\n$(msg)\n")
        end
    end
end

# ---- run
@everywhere function safe_run(t)
    try
        F6Tasks.run_task(t)
    catch err
        "ERROR: " * sprint(showerror, err, catch_backtrace())[1:min(end, 4000)]
    end
end
println("warming up $(max(W, 1)) process(es)…"); flush(stdout)
tw = @elapsed (W > 0 ? (@sync for p in workers(); @async remotecall_wait(F6Tasks.warmup, p); end) : F6Tasks.warmup())
@printf("warm-up %.0f s; %d tasks on %d workers, started %s\n", tw, length(tasks), max(W, 1), now()); flush(stdout)

t0 = time(); ndone = Ref(0); nerr = Ref(0); next = Ref(0)
function take_task()
    lock(LK) do
        next[] >= length(tasks) && return nothing
        next[] += 1; tasks[next[]]
    end
end
function report(t, res, dt)
    ndone[] += 1
    if res isa AbstractString
        nerr[] += 1; log_error(t, res)
        @printf("[%s] %5d/%d ERROR %s\n", Dates.format(now(), "HH:MM:SS"), ndone[], length(tasks), t.key)
    else
        write_out(res)
        el = time() - t0
        @printf("[%s] %5d/%d %-48s %6.0f s  (ETA %.0f min)\n", Dates.format(now(), "HH:MM:SS"), ndone[],
                length(tasks), t.key, dt, el / ndone[] * (length(tasks) - ndone[]) / 60)
    end
    flush(stdout)
end
if W > 0
    @sync for p in workers()
        @async while (t = take_task()) !== nothing
            dt = @elapsed res = remotecall_fetch(safe_run, p, t)
            report(t, res, dt)
        end
    end
else
    while (t = take_task()) !== nothing
        dt = @elapsed res = safe_run(t); report(t, res, dt)
    end
end
@printf("finished %d tasks (%d errors) in %.1f min\n", ndone[], nerr[], (time() - t0) / 60)
