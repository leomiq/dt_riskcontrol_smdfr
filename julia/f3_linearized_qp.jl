# F3: linearized prediction models (N1 LTI, N2 LTV) of the ROM replica, the linear
# MPC they give (a convex QP), and the transfer function ṁ_ACS → ŷ.
#
# N1 (LTI): one linearization at the query state, (x̄, ū, d̄) = (x_init, ṁ_mid, ṁ_leak at t_q−Δt)
#     ẋ ≈ f̄ + A(x−x̄) + B(ṁ_ACS−ū) + E(ṁ_leak−d̄),    ŷ ≈ ȳ + C(x−x̄) + D(ṁ_net−n̄)
# N2 (LTV): the same, linearized at every step along a nominal trajectory
#     (x̄_k from the nonlinear model with ū_k), so the model follows the leak.
#     Iterated (sequential QP) up to convergence, N2 → N3 in the limit.
# Both are discretized with implicit Euler, as the NLP. The QP's action is then
# re-simulated with the nonlinear ode1be replica (= Simulink) and costed with P0.
#
# Usage, from the repository root:  julia --project=julia julia/f3_linearized_qp.jl
include(joinpath(@__DIR__, "SMDFRROM.jl"))
using .SMDFRROM, MAT, JuMP, Ipopt, Printf, LinearAlgebra, ForwardDiff, DelimitedFiles

TQ = isempty(ARGS) ? "1180" : ARGS[1]          # query instant, as in rom_reference_<TQ>.mat
dir = joinpath(@__DIR__, "..", "results", "f1_rom_replica")
r = matread(joinpath(dir, "rom_reference_$TQ.mat"))
p = load_params(r["par"])
h = Float64(r["t_sam"]); t_q = Float64(r["t_q"]); y_ref = Float64(r["target"])
U_y = 1244.0; lam_m = 1e-4; dm_lim = 5.0; lam_pct = 10.0

xs = r["xs"]; blk(n) = findfirst(b -> endswith(b, "/" * n), vec(xs["block"]))
col(n) = xs["data"][blk(n)]
tX = vec(xs["time"][blk("Int_N")])
Xs = hcat(col("Int_N"), col("Int-C_i"), col("Int-Ce_i"), col("Int-T_f"),
          col("Int-T_w"), col("Int-T_c"), col("Int-mdot_c"))
ip = r["inport_data"]
k0 = findfirst(≈(t_q - h), tX); hor = k0+1:length(tX); K = length(hor)
x_init = Xs[k0, :]; leak = ip[hor, 2]; leak0 = ip[k0, 2]
ub = r["u_bounds"]; m_lo = current_to_massflow(minimum(ub[:, 1])); m_hi = current_to_massflow(maximum(ub[:, 2]))
S = max.(abs.(x_init), 1.0)

poly(n) = SMDFRROM.polyval(p.poly, n)
dpoly(n) = ForwardDiff.derivative(poly, n)
cost_P0(y, m) = sum((y .- y_ref) .^ 2) + lam_m * sum(max.(0, abs.(diff(m)) .- dm_lim) .^ 2) +
                lam_pct * sum(max.(0, y .- U_y) .^ 2)
function evaluate(m)        # nonlinear ode1be re-simulation (validated = Simulink)
    X = simulate(x_init, vcat(leak0, leak), vcat(0.0, m), h, p; scheme = :ode1be, input_at = :next)
    y = [SMDFRROM.output(X[k+1, :], leak[k], m[k], p) for k in 1:K]
    y, cost_P0(y, m), X
end

"""Jacobians of f and of the output at (x, leak, acs)."""
function lin(x, d, u)
    A = ForwardDiff.jacobian(z -> rhs(z, d, u, p), x)
    B = ForwardDiff.derivative(v -> rhs(x, d, v, p), u)
    E = ForwardDiff.derivative(v -> rhs(x, v, u, p), d)
    (; f = rhs(x, d, u, p), A, B, E, D = dpoly(d + u + 5550), y = SMDFRROM.output(x, d, u, p))
end

"""
Linear MPC (QP). `pts[k]` = (x̄_k, d̄_k, ū_k) is the linearization point of step k
(the same point for every k gives N1). Implicit Euler on the affine model.
"""
function solve_qp(pts)
    L = [lin(pt...) for pt in pts]
    md = Model(Ipopt.Optimizer); set_silent(md)
    set_optimizer_attribute(md, "tol", 1e-10)
    @variable(md, m_lo <= m[1:K] <= m_hi, start = (m_lo + m_hi) / 2)
    @variable(md, z[0:K, 1:NX])               # scaled state, x = S .* z
    @variable(md, s[1:K-1] >= 0); @variable(md, rr[1:K] >= 0)
    for j in 1:NX; fix(z[0, j], x_init[j] / S[j]; force = true); end
    for k in 1:K
        l = L[k]; xb, db, ub_ = pts[k]
        for j in 1:NX
            fj = l.f[j] + sum(l.A[j, i] * (S[i] * z[k, i] - xb[i]) for i in 1:NX if l.A[j, i] != 0; init = zero(AffExpr)) +
                 l.B[j] * (m[k] - ub_) + l.E[j] * (leak[k] - db)
            @constraint(md, z[k, j] - z[k-1, j] == h * fj / S[j])
        end
    end
    @expression(md, y[k = 1:K], L[k].y + (S[15] * z[k, 15] - pts[k][1][15]) +
                                L[k].D * ((leak[k] + m[k]) - (pts[k][2] + pts[k][3])))
    @constraint(md, [k = 1:K-1], m[k+1] - m[k] <= dm_lim + s[k])
    @constraint(md, [k = 1:K-1], m[k] - m[k+1] <= dm_lim + s[k])
    @constraint(md, [k = 1:K], y[k] - U_y <= rr[k])
    @objective(md, Min, sum((y[k] - y_ref)^2 for k in 1:K) + lam_m * sum(s .^ 2) + lam_pct * sum(rr .^ 2))
    t = @elapsed optimize!(md)
    (m = value.(m), y_pred = value.(y), obj = objective_value(md), status = termination_status(md), time = t)
end

# --- NLP reference (f2_monolithic_nlp.jl output)
csv, hdr = readdlm(joinpath(dir, "f2_nlp_vs_lj_$TQ.csv"), ','; header = true)
m_nlp = Float64.(csv[:, findfirst(==("mdot_NLP"), vec(hdr))])
y_nlp, J_nlp, _ = evaluate(m_nlp)
@printf("NLP (N3) reference: J = %.6f, u₀ ↔ ṁ₀ = %.4f kg/s\n\n", J_nlp, m_nlp[1])

m_mid = (m_lo + m_hi) / 2
rows = Any[]
report(name, sol) = begin
    y, J, _ = evaluate(sol.m)
    pred_err = maximum(abs.(sol.y_pred .- y))
    @printf("%-22s %-10s %7.3f s   J_QP = %.6f   J (nonlinear re-sim) = %.6f   gap vs NLP %+7.2f %%   |ŷ_lin − ŷ| ≤ %.3f K   ṁ₀ = %.3f\n",
            name, string(sol.status), sol.time, sol.obj, J, 100 * (J - J_nlp) / J_nlp, pred_err, sol.m[1])
    push!(rows, (name, sol.time, sol.obj, J, pred_err, sol.m))
    J
end

# --- N1: LTI at the query state
solve_qp(fill((x_init, leak0, m_mid), K))      # compile
sol1 = solve_qp(fill((x_init, leak0, m_mid), K))
report("N1 (LTI)", sol1)

# --- N2: LTV along the nonlinear trajectory of the previous iterate (SQP-like)
m_it = fill(m_mid, K)
for it in 1:6
    _, _, X = evaluate(m_it)
    pts = [(X[k+1, :], leak[k], m_it[k]) for k in 1:K]
    sol = solve_qp(pts)
    J = report("N2 (LTV), iter $it", sol)
    step = maximum(abs.(sol.m .- m_it)); global m_it = sol.m
    step < 1e-6 && break
end

# --- transfer function of N1: ṁ_ACS → ŷ (deviation variables)
l = lin(x_init, leak0, m_mid)
Cv = zeros(NX); Cv[15] = 1.0
G(s) = (Cv' * ((s * I - l.A) \ l.B)) + l.D
ev = eigvals(l.A)
@printf("\nN1 at t_q − Δt: eigenvalues of A (1/s): fastest %.3e, slowest %.3e; unstable: %d\n",
        minimum(real.(ev)), maximum(real.(ev)), count(>(0), real.(ev)))
# ΔM is a pure integrator (dΔM/dt = ṁ_leak + ṁ_ACS): G has a pole at s = 0. Split it:
#   G(s) = D + G_T(s) + G_M(s)/s, with G_T the path with the inventory frozen and
#   G_M(s) = ∂ŷ/∂ΔM through ṁ_c and h_wc (dynamic), G_M(0) its static value.
r16 = 1:NX-1
AT, BT, aM = l.A[r16, r16], l.B[r16], l.A[r16, NX]
GT0 = -(Cv[r16]' * (AT \ BT))                  # inventory frozen
GM0 = -(Cv[r16]' * (AT \ aM))                  # K per kg of inventory, static
@printf("Inventory frozen: static gain ∂ŷ/∂ṁ_ACS = %.5f K/(kg/s)  [direct D = %.5f, through T_w = %.5f]\n",
        l.D + GT0, l.D, GT0)
@printf("Integrating path: ∂ŷ/∂ΔM = %.3e K/kg  →  slope %.3e K/s per kg/s of ACS flow sustained\n", GM0, GM0)
println("  ω [rad/s]    |G(iω)| [K/(kg/s)]   phase [deg]")
for ω in (1e-3, 1e-2, 0.1, 1.0, 10.0)
    g = G(im * ω); @printf("  %8.3f     %.5f            %8.2f\n", ω, abs(g), rad2deg(angle(g)))
end

open(joinpath(dir, "f3_qp_vs_nlp_$TQ.csv"), "w") do io
    println(io, "t,mdot_NLP," * join(["mdot_" * replace(r[1], r"[ (),]" => "") for r in rows], ","))
    for k in 1:K
        println(io, join([@sprintf("%.1f", tX[hor[k]]); @sprintf("%.8f", m_nlp[k]);
                          [@sprintf("%.8f", r[6][k]) for r in rows]], ","))
    end
end
