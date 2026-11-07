# Shared code of the F4 scripts: query loading, MC-dropout leak scenarios, the
# monolithic NLP (deterministic and two-stage SAA) with configurable actuator limits,
# and the ode1be re-simulation used to cost every decision.
#
# Actuator limits (see docs/monolithic_formulation_plan.md §8.3):
#   m_lo, m_hi   bounds on the ACS make-up flow ṁ_ACS [kg/s] (default: archived current bounds)
#   ramp = :soft  P0 as in fitness_cost.m: |Δṁ| ≤ 5 kg/s per step penalized with λ_m = 1e-4,
#                 only between horizon steps
#   ramp = :hard  |ṁ_k − ṁ_{k−1}| ≤ 5 kg/s per step as a constraint, k = 0 included
#                 (from the flow applied at t_q − Δt): a rate limiter of 10 kg/s²
module F4Lib
include(joinpath(@__DIR__, "SMDFRROM.jl"))
using .SMDFRROM, MAT, JuMP, Ipopt, Statistics
export load_query, evaluate, solve_det, solve_saa, scen, quantile_leak, current_to_massflow,
       massflow_to_current, p4_fit_range

const U_y = 1244.0; const lam_m = 1e-4; const dm_lim = 5.0; const lam_pct = 10.0   # fitness_cost.m
const p4_fit_range = (1737.6, 5550.5)    # ṁ range of the data behind p₄ (models/fit_poly.mat)

"""Load the ROM query of F1 and the MC-dropout ensemble of python/f4_leak_ensemble.py."""
function load_query(TQ::AbstractString; root = joinpath(@__DIR__, ".."))
    r = matread(joinpath(root, "results", "f1_rom_replica", "rom_reference_$TQ.mat"))
    e = matread(joinpath(root, "results", "f4_risk", "leak_ensemble_$TQ.mat"))
    p = load_params(r["par"])
    h = Float64(r["t_sam"]); t_q = Float64(r["t_q"]); y_ref = Float64(r["target"])
    xs = r["xs"]; blk(n) = findfirst(b -> endswith(b, "/" * n), vec(xs["block"]))
    col(n) = xs["data"][blk(n)]
    tX = vec(xs["time"][blk("Int_N")])
    Xs = hcat(col("Int_N"), col("Int-C_i"), col("Int-Ce_i"), col("Int-T_f"),
              col("Int-T_w"), col("Int-T_c"), col("Int-mdot_c"))
    ip = r["inport_data"]
    k0 = findfirst(≈(t_q - h), tX); hor = k0+1:length(tX); K = length(hor)
    ub = r["u_bounds"]
    leak_now = ip[hor[1], 2]
    @assert isapprox(Float64(e["leak_now"]), leak_now; atol = 1e-9) "ensemble window does not end at t_q"
    Y = e["Y"][:, 1:K-1, 1]
    real = vec(e["real"])
    (; TQ, p, h, t_q, y_ref, tX, hor, K, x_init = Xs[k0, :], leak0 = ip[k0, 2], leak_now,
       m_prev = current_to_massflow(ip[k0, 4]),            # ACS flow applied at t_q − Δt
       m_lo_arch = current_to_massflow(minimum(ub[:, 1])), m_hi_arch = current_to_massflow(maximum(ub[:, 2])),
       S = max.(abs.(Xs[k0, :]), 1.0), Y, NENS = size(Y, 1),
       det_leak = vcat(leak_now, e["y_det"][1:K-1, 1]), used_leak = ip[hor, 2],
       real_leak = length(real) >= K - 1 ? vcat(leak_now, real[1:K-1]) : nothing)
end

scen(c, i) = vcat(c.leak_now, c.Y[i, :])
quantile_leak(c, pc) = vcat(c.leak_now, [quantile(c.Y[:, k], pc / 100) for k in 1:c.K-1])

cost_P0(c, y, m) = sum((y .- c.y_ref) .^ 2) + lam_m * sum(max.(0, abs.(diff(m)) .- dm_lim) .^ 2) +
                   lam_pct * sum(max.(0, y .- U_y) .^ 2)

"""ode1be re-simulation (= Simulink): PCT, P0 cost, states."""
function evaluate(c, m, L)
    X = simulate(c.x_init, vcat(c.leak0, L), vcat(0.0, m), c.h, c.p; scheme = :ode1be, input_at = :next)
    y = [SMDFRROM.output(X[k+1, :], L[k], m[k], c.p) for k in 1:c.K]
    (y = y, J = cost_P0(c, y, m), X = X)
end

"""Implicit-Euler ROM constraints for one scenario; returns (PCT expressions, scaled states)."""
function add_scenario!(md, c, m, L, guess)
    p, h, K, S = c.p, c.h, c.K, c.S
    Xg = simulate(c.x_init, vcat(c.leak0, L), vcat(0.0, guess), h, p; scheme = :be)
    z = @variable(md, [k = 0:K, j = 1:NX], start = Xg[k+1, j] / S[j])
    for j in 1:NX; fix(z[0, j], c.x_init[j] / S[j]; force = true); end
    x(k, j) = S[j] * z[k, j]
    for k in 1:K
        N = x(k, 1); C = [x(k, 1 + i) for i in 1:6]; Ce = [x(k, 7 + i) for i in 1:6]
        Tf, Tw, Tc, dM = x(k, 14), x(k, 15), x(k, 16), x(k, 17)
        M_act  = @expression(md, p.M_c + dM)
        mdot_c = @expression(md, p.mdot_c * M_act / p.M_c)
        h_wc   = @expression(md, p.h_wc * (mdot_c / 5550)^0.8)
        Tin_c  = @expression(md, p.Te_c0 - (p.Q0 * mdot_c / 2775 - (5500 - mdot_c)) / (mdot_c * p.c_pc))
        rho    = @expression(md, p.rho_0 + p.alfa_f * (Tf - p.Tbar_f0) + p.alfa_c * (Tc - p.Tbar_c0) + p.rhoex)
        Pi     = @expression(md, p.P_0 / p.N_0 * N * p.Pi_0)
        f = Vector{Any}(undef, NX)
        f[1] = @expression(md, (rho - p.beta) * N / p.Lambda + sum(p.lambda_i[i] * C[i] for i in 1:6))
        for i in 1:6
            f[1+i] = @expression(md, p.beta_i[i] / p.Lambda * N - (p.lambda_i[i] + p.lambda_c) * C[i] + p.lambda_e * Ce[i])
            f[7+i] = @expression(md, p.lambda_c * C[i] - (p.lambda_i[i] + p.lambda_e) * Ce[i])
        end
        f[14] = @expression(md, (Pi + p.mdot_f * p.c_pf * 2 * (p.Te_f0 - Tf) - p.h_fw * p.A_fw * (Tf - Tw)) / (p.M_f * p.c_pf))
        f[15] = @expression(md, (p.h_fw * p.A_fw * (Tf - Tw) - h_wc * p.A_wc * (Tw - Tc)) / (p.M_w * p.c_pw))
        f[16] = @expression(md, (mdot_c * p.c_pc * 2 * (Tin_c - Tc) + h_wc * p.A_wc * (Tw - Tc)) / (M_act * p.c_pc))
        f[17] = @expression(md, L[k] + m[k])
        for j in 1:NX
            @constraint(md, z[k, j] - z[k-1, j] == h * f[j] / S[j])
        end
    end
    pv = p.poly
    [@expression(md, x(k, 15) + sum(pv[i] * (L[k] + m[k] + 5550)^(length(pv) - i) for i in eachindex(pv)) + 0.4) for k in 1:K]
end

"""P0 cost of one scenario; the rate limit is soft (P0) or hard (from ṁ_prev)."""
function scenario_cost!(md, c, m, y, ramp)
    K = c.K
    rr = @variable(md, [1:K], lower_bound = 0)
    @constraint(md, [k = 1:K], y[k] - U_y <= rr[k])
    track = @expression(md, sum((y[k] - c.y_ref)^2 for k in 1:K) + lam_pct * sum(rr .^ 2))
    if ramp == :soft
        s = @variable(md, [1:K-1], lower_bound = 0)
        @constraint(md, [k = 1:K-1], m[k+1] - m[k] <= dm_lim + s[k])
        @constraint(md, [k = 1:K-1], m[k] - m[k+1] <= dm_lim + s[k])
        return @expression(md, track + lam_m * sum(s .^ 2))
    end
    mm = [c.m_prev; m]
    @constraint(md, [k = 1:K], mm[k+1] - mm[k] <= dm_lim)
    @constraint(md, [k = 1:K], mm[k] - mm[k+1] <= dm_lim)
    track
end

function new_model()
    md = Model(Ipopt.Optimizer); set_silent(md)
    set_optimizer_attribute(md, "tol", 1e-10); set_optimizer_attribute(md, "max_iter", 1000)
    md
end

clampg(c, g, lo, hi) = clamp.(g, lo, hi)

"""Deterministic NLP (F2) with leak L."""
function solve_det(c, L; m_lo = c.m_lo_arch, m_hi = c.m_hi_arch, ramp = :soft, guess = nothing)
    g = clampg(c, guess === nothing ? fill(m_hi, c.K) : guess, m_lo, m_hi)
    ramp == :hard && guess === nothing && (g = clampg(c, c.m_prev .+ 5.0 .* (1:c.K), m_lo, m_hi))
    md = new_model()
    @variable(md, m_lo <= m[k = 1:c.K] <= m_hi, start = g[k])
    y = add_scenario!(md, c, m, L, g)
    @objective(md, Min, scenario_cost!(md, c, m, y, ramp))
    t = @elapsed optimize!(md)
    (m = value.(m), obj = objective_value(md), status = termination_status(md), time = t)
end

"""Two-stage SAA: ṁ₀ shared, ṁ_k (k ≥ 2) per scenario. With m0_fixed, only the recourse."""
function solve_saa(c, Ls; m0_fixed = nothing, m_lo = c.m_lo_arch, m_hi = c.m_hi_arch, ramp = :soft, guess = nothing)
    K = c.K
    g = clampg(c, guess === nothing ? fill(m_hi, K) : guess, m_lo, m_hi)
    md = new_model()
    @variable(md, m_lo <= m0 <= m_hi, start = g[1])
    m0_fixed === nothing || fix(m0, m0_fixed; force = true)
    costs = Any[]; rec = Any[]
    for L in Ls
        mr = @variable(md, [k = 2:K], lower_bound = m_lo, upper_bound = m_hi)
        for k in 2:K; set_start_value(mr[k], g[k]); end
        m = [m0; [mr[k] for k in 2:K]]
        push!(costs, scenario_cost!(md, c, m, add_scenario!(md, c, m, L, g), ramp)); push!(rec, mr)
    end
    @objective(md, Min, sum(costs) / length(Ls))
    t = @elapsed optimize!(md)
    (m0 = value(m0), obj = objective_value(md), status = termination_status(md), time = t,
     m = [vcat(value(m0), [value(r_[k]) for k in 2:K]) for r_ in rec])
end

end
