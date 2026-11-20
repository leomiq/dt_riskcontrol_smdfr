"""
    F6Tasks

Task definitions of the paper campaign (plan §10, points 2 and 3), run in parallel by
`f6_campaign.jl` (one closed-loop run per task). Three task kinds:

  :rb    failure-free closed loop + Rao–Blackwell crossing probability (as f5b1c_rb_risk.jl,
         plan §8.16), evaluated for SEVERAL failure accelerations μ× at once: the loop of a
         controller that does not use the failure model does not depend on μ, only the integral
         does. τ series are written so any other rate can be re-weighted afterwards.
  :eta   failure-free loop + generalised Rao–Blackwell over severities η (as f5b1d_partial.jl,
         §8.18), evaluated for several severity probability vectors p_η at once.
  :fail  closed loop with the perturbations of a family (DET-accelerated failures as
         f5b1c_risk.jl; nominal and extreme families of f5a2_monte_carlo.jl), with the control
         metrics of the thesis (`thesis_metrics`) and, for N1/N2/N3 without a risk layer, a
         per-query optimality audit (every `AUDIT_EVERY` queries: the loop decision, N1, N2,
         N3 warm-started and polished, N3 cold, all scored with the same P0 cost on the
         ode1be ROM; gap = J/J_best − 1, as f5a_optimality_audit.jl).

Families (draws identical to the original scripts, common random numbers across controllers):
  :detacc   seeds 201…: A ~ U(.1,.2), B ~ U(10,15), d ~ U(15,25), lead ~ U(0,17) (f5b1c draw);
            plant capacity m_max, ramp 5 kg/s per step.
  :nominal  odd seeds: A ~ U(.2,.4), B ~ U(9,15), d ~ U(8,10), lead 17 s (f5a2 standard).
  :extreme  seeds 101…: f5a2 extreme draw (η, one failure of U(30,90) s, forecast bias).
"""
module F6Tasks

include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop, MAT, JuMP, Ipopt, Random, Statistics, Printf
const CL = ClosedLoop
const F4 = CL.F4Lib
const h = 0.5

export run_task

# ---------------------------------------------------------------- failure model (DET, §8.14)
const DB = matread(joinpath(@__DIR__, "..", "models", "MCDET_conferencepaper_database.mat"))
const LAM0 = vec(DB["Lambda_t"]); const MU0 = vec(DB["Mu_t"])
const LAM_ACC = 1e5                                   # failure-rate acceleration (fixed)
const LAM = LAM0 .* LAM_ACC; const LS = sum(LAM)
dist(muacc) = (w = LAM ./ LS, mu = MU0 .* muacc)
const MUS = (3e3, 1e4, 3e4)                           # repair-rate accelerations evaluated
const TC = 220.0                                      # conditioning window [s] (f5b1c)
const NQ_DET = parse(Int, get(ENV, "F6_NQ", "481"))   # F6_NQ: short smoke tests

const ETAS = [0.75, 0.5, 0.25, 0.0]
const PVECS = (base = [0.4, 0.3, 0.2, 0.1], uniform = [0.25, 0.25, 0.25, 0.25],
               pess = [0.1, 0.2, 0.3, 0.4], twotrain = [0.0, 0.7, 0.0, 0.3])
const TCAP_ETA = 400.0

const _cache = Dict{Any,Any}()
# optional extra callback CAPTURE[](t, c, L, act, prev, dec) after each MPC decision of a :fail task
# (f7_certify.jl uses it to capture closed-loop queries without changing the loop)
const CAPTURE = Ref{Any}(nothing)
cached(f, k) = get!(f, _cache, k)
qgrid(n, muacc) = cached(() -> [CL.failure_quantile((i - 0.5) / n, dist(muacc)) for i in 1:n], (:dq, n, muacc))
const DSM_THR = 0.00237737

"""Keyword arguments of `run_loop` for a risk / DSM specification (see f5b1c_rb_risk.jl,
f5b1d_partial.jl, f5b1b0_dsm_lambda.jl)."""
function riskkw(spec, mmax)
    spec === nothing && return NamedTuple()
    k = spec[1]
    if k == :var
        _, a, mu = spec; return (risk = (tau = CL.failure_quantile(1 - a, dist(mu)), lam = 1e3),)
    elseif k == :cvar
        _, b, mu = spec; return (risk = (mode = :cvar, beta = b, Dq = qgrid(200, mu), lam = 1e3),)
    elseif k == :wilks
        _, mu = spec; return (risk = (mode = :wilks, dist = dist(mu), n = 59, lam = 1e3),)
    elseif k == :exp
        _, r, mu = spec; return (risk = (mode = :expected, rho = r, Dq = qgrid(40, mu)),)
    elseif k == :vareta
        _, a, pn = spec
        return (risk = (mode = :var_eta, alpha = a, dist = dist(1e4), lam = 1e3, etas = ETAS,
                        peta = PVECS[pn], mmax = mmax),)
    elseif k in (:scencvar, :scenexp)
        _, v = spec
        base = (etas = ETAS, peta = PVECS.base, mmax = mmax, J = [1, 6, 11, 16, 21, 26], Dq = qgrid(10, 1e4))
        return k == :scencvar ? (risk = (mode = :scen_cvar, beta = v, lam = 1e3, base...),) :
                                (risk = (mode = :scen_exp, rho = v, base...),)
    elseif k == :dsm
        _, l, U = spec; return (dsm = (lam = l, thr = DSM_THR, U = U),)
    end
    error("risk spec $spec")
end

# ---------------------------------------------------------------- realisations
m_arch_hi() = cached(() -> setup().act_arch.m_hi, :m_arch_hi)

function draw(family, seed)
    rng = MersenneTwister(seed)
    if family == :detacc
        A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng); lead = 17rand(rng)
        return (; A, Bp, d, lead, eta = 1.0, fs = NaN, fl = NaN, fb = 0.0)
    elseif family == :nominal
        A = 0.2 + 0.2rand(rng); Bp = 9 + 6rand(rng); d = 8 + 2rand(rng)
        return (; A, Bp, d, lead = 17.0, eta = 1.0, fs = NaN, fl = NaN, fb = 0.0)
    elseif family == :extreme
        A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng)
        eta = 0.6 + 0.3rand(rng); fs = 5 + 15rand(rng); fl = 30 + 60rand(rng); fb = -0.5 + 0.6rand(rng)
        lead = 17rand(rng)
        return (; A, Bp, d, lead, eta, fs, fl, fb)
    end
    error("family $family")
end

"""Plant capacity of an actuator name: `capNN` = NN kg/s with ramp 5; `arch` = archived
pump box (plant capped at its maximum); `phys`/`none` = unlimited plant."""
function plant_cap(act)
    s = string(act)
    startswith(s, "cap") && return (parse(Float64, s[4:end]), 5.0)
    act == :arch && return (m_arch_hi(), Inf)
    (Inf, Inf)
end

function bench(family, seed, act)
    cached((:bench, family, seed, act)) do
        r = draw(family, seed)
        mmax, ramp = plant_cap(act)
        pH = CL.HRMQueries.hrm_params(; A = r.A, B = r.Bp, mdot_c_drift = r.d, acs_max = mmax, acs_ramp = ramp)
        B = setup(; pH, lead = r.lead)
        (; B, r, mmax)
    end
end

function actuator(B, act)
    s = string(act)
    startswith(s, "cap") && return act_cap(parse(Float64, s[4:end]), 5.0)
    act == :arch && return B.act_arch
    ACT_PHYS
end

# DET-accelerated failure windows conditioned on ≥ 1 failure starting in [t_start, t_start+220]
function ctmc_windows(rng, t0, t1, mu)
    iv = Tuple{Float64,Float64}[]
    for i in eachindex(LAM)
        t = t0
        while true
            t += randexp(rng) / LAM[i]; t > t1 && break
            d = randexp(rng) / mu[i]; push!(iv, (t, t + d)); t += d
        end
    end
    sort!(iv); out = Tuple{Float64,Float64}[]
    for (a, b) in iv
        !isempty(out) && a <= out[end][2] ? (out[end] = (out[end][1], max(out[end][2], b))) : push!(out, (a, b))
    end
    out
end
function det_failures(seed, t_act, t_start; muacc = 1e4)
    rng = MersenneTwister(50_000 + seed); tries = 0
    while true
        tries += 1
        w = ctmc_windows(rng, t_act, t_start + 240.0, MU0 .* muacc)
        any(t_start <= a <= t_start + 220.0 for (a, _) in w) && return (w, tries)
    end
end

# ---------------------------------------------------------------- Rao–Blackwell
trapz_upto(t, g, Tmax) = (k = findall(<=(Tmax + 1e-9), t); tk = t[k]; gk = g[k];
                          sum((tk[2:end] .- tk[1:end-1]) .* (gk[2:end] .+ gk[1:end-1]) ./ 2))
"""P(cross) (unconditional, over the run) and P(cross | failure in [0, TC]) for failure-length
distribution `ds`; `survfun(k)` = Σ over severities of p·S(τ) at evaluation point k."""
function rb_pair(t, survfun)
    g = [LS * exp(-LS * ti) * survfun(k) for (k, ti) in enumerate(t)]
    (unc = trapz_upto(t, g, t[end]), cond = trapz_upto(t, g, TC) / (1 - exp(-LS * TC)))
end

# ---------------------------------------------------------------- optimality audit
function solve_n3(c, L, act, guess; acceptable = true)
    g = clamp.(guess, act.m_lo, act.m_hi)
    md = Model(Ipopt.Optimizer); set_silent(md)
    set_optimizer_attribute(md, "tol", 1e-10); set_optimizer_attribute(md, "max_iter", 3000)
    acceptable || set_optimizer_attribute(md, "acceptable_iter", 0)
    @variable(md, act.m_lo <= m[k = 1:c.K] <= act.m_hi, start = g[k])
    y = F4.add_scenario!(md, c, m, L, g)
    @objective(md, Min, F4.scenario_cost!(md, c, m, y, act.ramp))
    t = @elapsed optimize!(md)
    (m = value.(m), status = string(termination_status(md)), time = t)
end
cold_guess(c, act) = act.ramp == :hard ? c.m_prev .+ 5.0 .* (1:c.K) : fill(min(act.m_hi, 50.0), c.K)

function make_audit(every, rows, key)
    k = Ref(0)
    function audit(c, L, act, prev, dec)
        k[] += 1
        (k[] - 1) % every == 0 || return
        J(m) = F4.evaluate(c, clamp.(m, act.m_lo, act.m_hi), L).J
        pol = solve_n3(c, L, act, dec.m; acceptable = false)
        cold = solve_n3(c, L, act, cold_guess(c, act))
        n1 = CL.decide(:N1, c, L, act, prev); n2 = CL.decide(:N2, c, L, act, prev)
        n3 = CL.decide(:N3, c, L, act, prev)
        Js = (loop = J(dec.m), N1 = J(n1.m), N2 = J(n2.m), N3 = J(n3.m), N3pol = J(pol.m), N3cold = J(cold.m))
        Jb = minimum(values(Js))
        gap(x) = Jb > 0 ? x / Jb - 1 : x - Jb
        push!(rows, (; key, k = k[], Jbest = Jb, gap_loop = gap(Js.loop), gap_N1 = gap(Js.N1),
                     gap_N2 = gap(Js.N2), gap_N3 = gap(Js.N3), gap_N3pol = gap(Js.N3pol),
                     gap_N3cold = gap(Js.N3cold), st_N3 = n3.status, st_pol = pol.status,
                     st_cold = cold.status, dm0_N3pol = abs(dec.m[1] - pol.m[1]),
                     t_N1 = n1.time, t_N2 = n2.time, t_N3 = n3.time))
    end
end

# ---------------------------------------------------------------- tasks
nanrow(names) = NamedTuple{names}(ntuple(_ -> NaN, length(names)))
const MNAMES = (:n, :breakdown, :limit, :t_limit, :pct_max, :pct_final, :t_above, :ise, :iae, :err_last10,
                :mass_cmd, :mass_plant, :ramp_max, :t_p50, :t_p99, :t_max, :fails)
const TNAMES = (:min_safety_margin, :margin_1220, :min_overshoot, :max_abs_dev_pt, :reach_time,
                :settling_time, :steady_state_error, :IAE, :cumulative_error, :peak_PCT, :final_PCT,
                :total_variation, :max_control_rate, :tv_plant, :coolant_mass, :runtime_median)

"""Run one task; returns a Dict file => vector of rows (NamedTuples) to append."""
function run_task(t)
    out = Dict{String,Vector{Any}}()
    put(f, r) = push!(get!(out, f, Any[]), r)
    bb = bench(t.family, t.seed, t.act)
    B, r = bb.B, bb.r
    head = (; key = t.key, campaign = t.campaign, family = t.family, seed = t.seed, case = t.case,
            ctrl = t.ctrl, bias = t.bias, act = t.act, mmax = bb.mmax, mu_design = t.mu_design,
            A = r.A, B = r.Bp, drift = r.d)
    if B === nothing
        put(t.campaign * "_skipped.csv", (; key = t.key, seed = t.seed, reason = "no ACS demand"))
        return out
    end
    act = actuator(B, t.act)
    kw = riskkw(t.spec, bb.mmax)
    lead = B.t_start - B.t_act
    if t.kind == :rb
        tr = @elapsed lg = run_loop(B, t.ctrl; n = NQ_DET, act, bias = t.bias, keep_states = true, kw...)
        idx = unique(vcat(1:20:length(lg.t_q), length(lg.t_q)))
        tt = @elapsed τ = [time_to_limit(B.pH, lg.states[i][1], lg.states[i][2], lg.t_q[i]) for i in idx]
        tq = lg.t_q[idx] .- lg.t_q[1]
        ps = [rb_pair(tq, k -> CL.failure_surv(τ[k], dist(mu))) for mu in MUS]
        M = metrics(lg, B.y_ref); T = thesis_metrics(lg, B.y_ref)
        put(t.campaign * ".csv", (; head..., lead,
            p_cond_mu3e3 = ps[1].cond, p_cond_mu1e4 = ps[2].cond, p_cond_mu3e4 = ps[3].cond,
            p_unc_mu3e3 = ps[1].unc, p_unc_mu1e4 = ps[2].unc, p_unc_mu3e4 = ps[3].unc,
            tau_mean = mean(τ), tau_min = minimum(τ), mean_dev = mean(lg.pct .- B.y_ref),
            ubound_mean = (u = filter(!isnan, lg.ubound); isempty(u) ? NaN : mean(u)),
            M..., T..., wall_loop = tr, wall_tau = tt))
        for (j, i) in enumerate(idx)
            put(t.campaign * "_tau.csv", (; key = t.key, t = tq[j], tau = τ[j], pct = lg.pct[i]))
        end
    elseif t.kind == :eta
        tr = @elapsed lg = run_loop(B, :N2; n = NQ_DET, act, bias = true, keep_states = true, kw...)
        idx = unique(vcat(1:20:length(lg.t_q), length(lg.t_q)))
        tt = @elapsed taus = [[time_to_limit(B.pH, lg.states[i][1], lg.states[i][2], lg.t_q[i];
                                              t_cap = TCAP_ETA, flow = e * bb.mmax) for i in idx] for e in ETAS]
        tq = lg.t_q[idx] .- lg.t_q[1]
        D = dist(1e4)
        S(e, k) = taus[e][k] >= TCAP_ETA ? 0.0 : CL.failure_surv(taus[e][k], D)
        pc = map(PVECS) do p
            rb_pair(tq, k -> sum(p[e] * S(e, k) for e in eachindex(ETAS))).cond
        end
        ptot = rb_pair(tq, k -> S(4, k)).cond
        M = metrics(lg, B.y_ref); T = thesis_metrics(lg, B.y_ref)
        put(t.campaign * ".csv", (; head..., lead, p_total_only = ptot,
            p_cond_base = pc.base, p_cond_uniform = pc.uniform, p_cond_pess = pc.pess, p_cond_twotrain = pc.twotrain,
            tau0_mean = mean(taus[4]), tau25_mean = mean(taus[3]), mean_dev = mean(lg.pct .- B.y_ref),
            M..., T..., wall_loop = tr, wall_tau = tt))
        for (e, η) in enumerate(ETAS), (j, i) in enumerate(idx)
            put(t.campaign * "_tau.csv", (; key = t.key, eta = η, t = tq[j], tau = taus[e][j], pct = lg.pct[i]))
        end
    elseif t.kind == :fail
        pert = NamedTuple(); finfo = (n_fail = 0, first_fail = NaN, max_fail_len = 0.0, down_time = 0.0, tries = 0)
        if t.family == :detacc
            wins, tries = det_failures(t.seed, B.t_act, B.t_start)
            inw = [(a, b) for (a, b) in wins if b > B.t_start && a < B.t_start + 240.0]
            pert = (fail = wins,)
            finfo = (n_fail = length(inw), first_fail = isempty(inw) ? NaN : inw[1][1] - B.t_start,
                     max_fail_len = isempty(inw) ? 0.0 : maximum(b - a for (a, b) in inw),
                     down_time = sum((min(b, B.t_start + 240) - max(a, B.t_start) for (a, b) in inw); init = 0.0),
                     tries)
        elseif t.family == :extreme
            pert = (eta = r.eta, fail = (B.t_start + r.fs, B.t_start + r.fs + r.fl), fbias = r.fb)
            finfo = (n_fail = 1, first_fail = r.fs, max_fail_len = r.fl, down_time = r.fl, tries = 1)
        end
        arows = Any[]
        audit = t.audit ? make_audit(t.audit_every, arows, t.key) : nothing
        if CAPTURE[] !== nothing
            a0 = audit; cap = CAPTURE[]
            audit = (c, L, act, prev, dec) -> (a0 === nothing || a0(c, L, act, prev, dec); cap(t, c, L, act, prev, dec))
        end
        tr = @elapsed lg = run_loop(B, t.ctrl; n = t.n, act, bias = t.bias, audit, pert..., kw...)
        M = metrics(lg, B.y_ref); T = thesis_metrics(lg, B.y_ref)
        au = if isempty(arows)
            (n_audit = 0, gap_loop_med = NaN, gap_loop_max = NaN, gap_N1_max = NaN, gap_N2_max = NaN,
             gap_N3_max = NaN, gap_N3cold_max = NaN)
        else
            g(s) = [getfield(a, s) for a in arows]
            (n_audit = length(arows), gap_loop_med = median(g(:gap_loop)), gap_loop_max = maximum(g(:gap_loop)),
             gap_N1_max = maximum(g(:gap_N1)), gap_N2_max = maximum(g(:gap_N2)), gap_N3_max = maximum(g(:gap_N3)),
             gap_N3cold_max = maximum(g(:gap_N3cold)))
        end
        ub = filter(!isnan, lg.ubound)
        put(t.campaign * ".csv", (; head..., lead, n_req = t.n, eta = r.eta, fbias = r.fb, finfo...,
            M..., T..., ubound_mean = isempty(ub) ? NaN : mean(ub), au..., wall = tr))
        foreach(a -> put(t.campaign * "_audit.csv", a), arows)
        for j in eachindex(lg.t_q)
            put(t.campaign * "_traj.csv", (; key = t.key, t = lg.t_q[j] - lg.t_q[1], pct = lg.pct[j], m = lg.m[j],
                                          m_plant = lg.m_plant[j], ubound = lg.ubound[j], solve = lg.time[j]))
        end
    else
        error("task kind $(t.kind)")
    end
    out
end

# ---------------------------------------------------------------- campaigns
function mk(campaign, kind, family, seed, act, case, ctrl; bias = ctrl in (:N1, :N2, :N3), spec = nothing,
            mu_design = NaN, n = NQ_DET, audit = false, audit_every = 20, cost = 1.0)
    key = "$(campaign)|$(family)|$(seed)|$(act)|$(case)"
    (; key, campaign, kind, family, seed, act, case, ctrl, bias, spec, mu_design, n, audit, audit_every, cost)
end
ecode(x) = (s = @sprintf("%.0e", x); replace(s, "e+0" => "e", "e+" => "e"))   # 3e3, 1e4, …

"""Task list of a campaign; `seeds` overrides the default seeds of its family."""
function campaign_tasks(name; seeds = nothing)
    T = Any[]
    if name == "sens"          # point 2: failure acceleration μ × ACS capacity (Rao–Blackwell)
        S = seeds === nothing ? (201:250) : seeds
        for s in S, m in (30, 35, 40, 45, 50)
            act = Symbol("cap$m")
            push!(T, mk(name, :rb, :detacc, s, act, "N2", :N2; cost = 1))
            push!(T, mk(name, :rb, :detacc, s, act, "PI", :PI; cost = 2))
            push!(T, mk(name, :rb, :detacc, s, act, "PID", :PID; cost = 2))
            for mu in MUS
                for a in (0.05, 0.02, 0.01)
                    push!(T, mk(name, :rb, :detacc, s, act, "N2_var_a$(a)_mu$(ecode(mu))", :N2;
                                spec = (:var, a, mu), mu_design = mu, cost = 1.5))
                end
                push!(T, mk(name, :rb, :detacc, s, act, "N2_cvar_b0.95_mu$(ecode(mu))", :N2;
                            spec = (:cvar, 0.95, mu), mu_design = mu, cost = 2))
            end
        end
    elseif name == "form"      # point 3: all risk formulations, 40 kg/s, μ×1e4 (Rao–Blackwell)
        S = seeds === nothing ? (201:250) : seeds
        for s in S
            push!(T, mk(name, :rb, :detacc, s, :cap40, "N2_var_a0.1", :N2; spec = (:var, 0.1, 1e4), mu_design = 1e4))
            for b in (0.90, 0.98)
                push!(T, mk(name, :rb, :detacc, s, :cap40, "N2_cvar_b$(b)", :N2; spec = (:cvar, b, 1e4), mu_design = 1e4, cost = 2))
            end
            push!(T, mk(name, :rb, :detacc, s, :cap40, "N2_wilks59", :N2; spec = (:wilks, 1e4), mu_design = 1e4, cost = 2))
            for rho in (10.0, 100.0, 1000.0)
                push!(T, mk(name, :rb, :detacc, s, :cap40, "N2_exp_r$(Int(rho))", :N2; spec = (:exp, rho, 1e4), mu_design = 1e4, cost = 3))
            end
        end
    elseif name == "eta"       # point 2: severity probabilities (partial failures, §8.18)
        S = seeds === nothing ? (201:250) : seeds
        for s in S
            push!(T, mk(name, :eta, :detacc, s, :cap40, "N2", :N2; cost = 3))
            push!(T, mk(name, :eta, :detacc, s, :cap40, "N2_total_a0.02", :N2; spec = (:var, 0.02, 1e4), mu_design = 1e4, cost = 3))
            for pn in keys(PVECS)
                push!(T, mk(name, :eta, :detacc, s, :cap40, "N2_vareta_a0.01_$(pn)", :N2; spec = (:vareta, 0.01, pn), mu_design = 1e4, cost = 3))
            end
            push!(T, mk(name, :eta, :detacc, s, :cap40, "N2_scencvar_b0.9", :N2; spec = (:scencvar, 0.9), mu_design = 1e4, cost = 6))
            push!(T, mk(name, :eta, :detacc, s, :cap40, "N2_scenexp_r100", :N2; spec = (:scenexp, 100.0), mu_design = 1e4, cost = 6))
        end
    elseif name == "det"       # point 3: closed loop WITH DET-accelerated failures, thesis metrics + audit
        S = seeds === nothing ? (201:250) : seeds
        for s in S
            for c in (:N1, :N2, :N3)
                push!(T, mk(name, :fail, :detacc, s, :cap40, "$(c)b", c; audit = true, cost = c == :N3 ? 8 : 4))
            end
            push!(T, mk(name, :fail, :detacc, s, :cap40, "N2_var_a0.02", :N2; spec = (:var, 0.02, 1e4), mu_design = 1e4, cost = 2))
            push!(T, mk(name, :fail, :detacc, s, :cap40, "N2_cvar_b0.95", :N2; spec = (:cvar, 0.95, 1e4), mu_design = 1e4, cost = 2))
            push!(T, mk(name, :fail, :detacc, s, :cap40, "N2_wilks59", :N2; spec = (:wilks, 1e4), mu_design = 1e4, cost = 2))
            push!(T, mk(name, :fail, :detacc, s, :cap40, "N2_exp_r100", :N2; spec = (:exp, 100.0, 1e4), mu_design = 1e4, cost = 3))
            push!(T, mk(name, :fail, :detacc, s, :cap40, "N2_DSM_U1244_l1e4", :N2; spec = (:dsm, 1e4, 1244.0), cost = 2))
            push!(T, mk(name, :fail, :detacc, s, :cap40, "PI", :PI; cost = 2))
            push!(T, mk(name, :fail, :detacc, s, :cap40, "PID", :PID; cost = 2))
        end
    elseif name == "nom"       # point 3: nominal family (thesis setting), thesis metrics + audit
        S = seeds === nothing ? (1:2:99) : seeds
        for s in S, act in (:arch, :cap40)
            for c in (:N1, :N2, :N3)
                push!(T, mk(name, :fail, :nominal, s, act, "$(c)b", c; n = min(121, NQ_DET), audit = true, audit_every = 10, cost = c == :N3 ? 3 : 1.5))
            end
            push!(T, mk(name, :fail, :nominal, s, act, "PID", :PID; n = min(121, NQ_DET), cost = 0.5))
            push!(T, mk(name, :fail, :nominal, s, act, "PI", :PI; n = min(121, NQ_DET), cost = 0.5))
        end
    elseif name == "ext"       # point 3: extreme family of F5a.3 (unlimited plant), thesis metrics + audit
        S = seeds === nothing ? (101:150) : seeds
        for s in S
            push!(T, mk(name, :fail, :extreme, s, :arch, "N2b", :N2; cost = 3))
            for c in (:N1, :N2, :N3)
                push!(T, mk(name, :fail, :extreme, s, :phys, "$(c)b", c; audit = true, cost = c == :N3 ? 8 : 4))
            end
            push!(T, mk(name, :fail, :extreme, s, :phys, "N2", :N2; bias = false, cost = 3))
            push!(T, mk(name, :fail, :extreme, s, :phys, "PI", :PI; cost = 2))
            push!(T, mk(name, :fail, :extreme, s, :phys, "PID", :PID; cost = 2))
        end
    else
        error("unknown campaign $name (sens, form, eta, det, nom, ext)")
    end
    T
end

"""Compile everything on this process with two-query runs (excluded from the timings)."""
function warmup()
    Bc = setup()
    for c in (:N1, :N2, :N3); run_loop(Bc, c; n = 2, bias = true); end
    run_loop(Bc, :N2; n = 2, bias = true, keep_states = true, risk = (tau = 20.0, lam = 1e3))
    time_to_limit(Bc.pH, Bc.hist[Bc.t_start-0.5][1], Bc.hist[Bc.t_start-0.5][2], Bc.t_start - 0.5; t_cap = 1.0)
    nothing
end

end # module
