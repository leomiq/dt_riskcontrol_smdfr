# F2: P0 written as a monolithic NLP (N3 level) and solved with Ipopt at the
# archived query t_q = 1180 s; compared with the archived Luus-Jaakola action.
#
#   min  Σ_k (ŷ_k − y_ref)² + λ_m Σ_k s_k² + λ_pct Σ_k r_k²
#   s.t. x_k = x_{k−1} + Δt f(x_k, ṁ_leak,k, ṁ_ACS,k)         (implicit Euler, k = 1..K)
#        ŷ_k = T_w,k + p₄(ṁ_leak,k + ṁ_ACS,k + 5550) + 0.4
#        |ṁ_ACS,k − ṁ_ACS,k−1| ≤ Δṁ + s_k,  ŷ_k − U_y ≤ r_k,  s, r ≥ 0
#        ṁ(u_min) ≤ ṁ_ACS,k ≤ ṁ(u_max),   x_0 = Simulink state at t_q − Δt
#
# The hinges [·]₊² of fitness_cost.m become slacks (exact). The decision is the ACS
# flow; the current is recovered with the inverse pump map. Every candidate is
# then re-simulated with the ode1be replica (= Simulink, see validate_replica.jl),
# so all costs are compared on the same model.
#
# Usage, from the repository root:  julia --project=julia julia/f2_monolithic_nlp.jl
include(joinpath(@__DIR__, "SMDFRROM.jl"))
using .SMDFRROM, MAT, JuMP, Ipopt, Printf

file = joinpath(@__DIR__, "..", "results", "f1_rom_replica", "rom_reference_1180.mat")
r = matread(file)
p = load_params(r["par"])
h = Float64(r["t_sam"]); t_q = Float64(r["t_q"]); y_ref = Float64(r["target"])
U_y = 1244.0; lam_m = 1e-4; dm_lim = 5.0; lam_pct = 10.0      # fitness_cost.m

xs = r["xs"]; blk(n) = findfirst(b -> endswith(b, "/" * n), vec(xs["block"]))
col(n) = xs["data"][blk(n)]
tX = vec(xs["time"][blk("Int_N")])
Xs = hcat(col("Int_N"), col("Int-C_i"), col("Int-Ce_i"), col("Int-T_f"),
          col("Int-T_w"), col("Int-T_c"), col("Int-mdot_c"))
ip = r["inport_data"]
k0 = findfirst(≈(t_q - h), tX)                 # state at t_q − Δt is the NLP initial state
hor = k0+1:length(tX)                          # 26 samples, t_q … t_q + 12.5 s
K = length(hor)
x_init = Xs[k0, :]
leak = ip[hor, 2]
u_opt = vec(r["u_opt"]); @assert length(u_opt) == K
ub = r["u_bounds"]; u_lo, u_hi = minimum(ub[:, 1]), maximum(ub[:, 2])
m_lo, m_hi = current_to_massflow(u_lo), current_to_massflow(u_hi)
@printf("t_q = %.1f s, K = %d, y_ref = %.4f K, u ∈ [%.4f, %.4f] A → ṁ_ACS ∈ [%.3f, %.3f] kg/s\n",
        t_q, K, y_ref, u_lo, u_hi, m_lo, m_hi)

"""P0 cost of fitness_cost.m for a PCT trajectory y and ACS flows m (both length K)."""
function cost_P0(y, m)
    Jt = sum((y .- y_ref) .^ 2)
    Jm = lam_m * sum(max.(0, abs.(diff(m)) .- dm_lim) .^ 2)
    Jp = lam_pct * sum(max.(0, y .- U_y) .^ 2)
    (J = Jt + Jm + Jp, Jt, Jm, Jp)
end

"""Re-simulate an ACS flow sequence with the ode1be replica and return (y, cost)."""
function evaluate(m)
    X = simulate(x_init, vcat(ip[k0, 2], leak), vcat(0.0, m), h, p; scheme = :ode1be, input_at = :next)
    y = [SMDFRROM.output(X[k+1, :], leak[k], m[k], p) for k in 1:K]
    y, cost_P0(y, m)
end

function solve_nlp(m_guess; print_level = 0)
    S = max.(abs.(x_init), 1.0)                 # state scaling (N ~ 1e13, C ~ 1e16)
    Xg = simulate(x_init, vcat(ip[k0, 2], leak), vcat(0.0, m_guess), h, p; scheme = :be)
    md = Model(Ipopt.Optimizer)
    set_optimizer_attribute(md, "print_level", print_level)
    set_optimizer_attribute(md, "tol", 1e-10)
    set_optimizer_attribute(md, "max_iter", 500)
    @variable(md, m_lo <= m[k = 1:K] <= m_hi, start = m_guess[k])
    @variable(md, z[k = 0:K, j = 1:NX], start = Xg[k+1, j] / S[j])
    @variable(md, s[1:K-1] >= 0)
    @variable(md, rr[1:K] >= 0)
    for j in 1:NX
        fix(z[0, j], x_init[j] / S[j]; force = true)
    end
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
        # Tin_f − Te_f = 2(Te_f0 − T_f);  Tin_c − Te_c = 2(Tin_c − T_c)
        f[14] = @expression(md, (Pi + p.mdot_f * p.c_pf * 2 * (p.Te_f0 - Tf) - p.h_fw * p.A_fw * (Tf - Tw)) / (p.M_f * p.c_pf))
        f[15] = @expression(md, (p.h_fw * p.A_fw * (Tf - Tw) - h_wc * p.A_wc * (Tw - Tc)) / (p.M_w * p.c_pw))
        f[16] = @expression(md, (mdot_c * p.c_pc * 2 * (Tin_c - Tc) + h_wc * p.A_wc * (Tw - Tc)) / (M_act * p.c_pc))
        f[17] = @expression(md, leak[k] + m[k])
        for j in 1:NX
            @constraint(md, z[k, j] - z[k-1, j] == h * f[j] / S[j])
        end
    end
    pv = p.poly
    @expression(md, y[k = 1:K], x(k, 15) +
        sum(pv[i] * (leak[k] + m[k] + 5550)^(length(pv) - i) for i in eachindex(pv)) + 0.4)
    @constraint(md, [k = 1:K-1], m[k+1] - m[k] <= dm_lim + s[k])
    @constraint(md, [k = 1:K-1], m[k] - m[k+1] <= dm_lim + s[k])
    @constraint(md, [k = 1:K], y[k] - U_y <= rr[k])
    @objective(md, Min, sum((y[k] - y_ref)^2 for k in 1:K) + lam_m * sum(s .^ 2) + lam_pct * sum(rr .^ 2))
    t = @elapsed optimize!(md)
    (m = value.(m), obj = objective_value(md), status = termination_status(md),
     iters = barrier_iterations(md), time = t)
end

# --- reference: archived LJ action, Simulink output and replica re-simulation
m_LJ = current_to_massflow.(u_opt)
y_sim = r["pred_tw_max_signal"][hor, 2]
c_sim = cost_P0(y_sim, m_LJ)
y_LJ, c_LJ = evaluate(m_LJ)
@printf("\nLuus–Jaakola (archived u_opt): J = %.6f on Simulink, %.6f on the replica\n", c_sim.J, c_LJ.J)

# --- NLP from several starts (detects local minima)
solve_nlp(m_LJ)                                 # compile
starts = ["LJ" => m_LJ, "u_min" => fill(m_lo, K), "u_max" => fill(m_hi, K), "mid" => fill((m_lo + m_hi) / 2, K)]
println("\nstart   status              iters  time [s]   J_NLP (BE)   J re-sim (ode1be)   J_track     J_mdot")
res = Dict{String,Any}()
for (name, mg) in starts
    sol = solve_nlp(mg)
    yv, c = evaluate(sol.m)
    res[name] = (sol, yv, c)
    @printf("%-7s %-19s %5d  %8.3f   %.6f     %.6f          %.6f   %.2e\n",
            name, string(sol.status), sol.iters, sol.time, sol.obj, c.J, c.Jt, c.Jm)
end

best = argmin(n -> res[n][3].J, collect(keys(res)))
sol, y_best, c_best = res[best]
u_best = massflow_to_current.(sol.m)
@printf("\nBest NLP start: %s.  J_NLP = %.6f vs J_LJ = %.6f (replica)  →  LJ gap = %.2f %%\n",
        best, c_best.J, c_LJ.J, 100 * (c_LJ.J - c_best.J) / c_best.J)
@printf("First action u₀: NLP %.5f A, LJ %.5f A  (Δ = %+.5f A)\n", u_best[1], u_opt[1], u_best[1] - u_opt[1])
@printf("Lower bound from the LJ cost? none: J_NLP is a local optimum on BE; see report.\n")

println("\n   t [s]   u_LJ [A]   u_NLP [A]   ŷ_LJ [K]    ŷ_NLP [K]")
for k in 1:K
    @printf("  %7.1f   %.5f    %.5f    %.4f   %.4f\n", tX[hor[k]], u_opt[k], u_best[k], y_LJ[k], y_best[k])
end

open(joinpath(dirname(file), "f2_nlp_vs_lj_1180.csv"), "w") do io
    println(io, "t,u_LJ,u_NLP,mdot_LJ,mdot_NLP,y_LJ,y_NLP,leak")
    for k in 1:K
        @printf(io, "%.1f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f\n", tX[hor[k]], u_opt[k], u_best[k],
                m_LJ[k], sol.m[k], y_LJ[k], y_best[k], leak[k])
    end
end
