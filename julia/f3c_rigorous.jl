# F3c, step 4b (docs/monolithic_formulation_plan.md §8.5): RIGOROUS bound on the nonlinear
# residual e_k(u) = ŷ_k(u) − a_k(u) over the whole control box, hence a certified lower
# bound on the global optimum of the nonlinear P0 (converged implicit-Euler model, as the NLP).
#
# First-order Taylor model with interval remainder (scaled states ζ = δ ./ S):
#   δ_k(u) = x_k(u) − x*_k = T_k Δu + r_k,   Δu = u − u*,   r_k ∈ interval box.
# Implicit Euler for the true and nominal trajectories, subtracted, with the mean-value form
# f(z,u) − f(z*,u*) = [A](Z_k)(z − z*) + B (u − u*), [A] = interval Jacobian of f over the
# reachable box Z_k, and B = e_ΔM constant (the control enters only the inventory):
#   M*_k δ_k = δ_{k−1} + h B Δu_k + h ([A] − A*_k) δ_k − res_k,     M*_k = I − h A*_k
# Linear part (point): T_k = C_k (T_{k−1} + h B e_kᵀ), C_k ≈ (M*_k)⁻¹, with its floating-point
# defect D_k = M*_k T_k − (T_{k−1} + h B e_kᵀ) enclosed by intervals. Remainder:
#   r_k ∈ C_k ( r_{k−1} + h ([A] − A*_k) [δ_k] − res_k − D_k [Δu] ) + (I − C_k M*_k) r_k
# r_k appears on both sides: a candidate box R is accepted when the operator maps it into
# itself (existence by Brouwer; the implicit-Euler step is assumed to have a unique solution
# near the nominal one, which Newton's convergence supports). Output:
#   e_k = S_Tw r_k[T_w] + ½ [p₄''](leak_k + U + 5550) Δu_k²     (exact remainder of p₄)
# All interval operations use outward rounding (IntervalArithmetic.jl); res_k and D_k make
# the floating-point point computations rigorous as well.
#
# Parts: 1 = validate the interval Jacobian; 2 = propagate the remainder (ε_k rigorous);
#        3 = certified lower bound and gap.
# Usage, from the repository root:  julia --project=julia julia/f3c_rigorous.jl [t_q] [last_part]
include(joinpath(@__DIR__, "F4Lib.jl"))
using .F4Lib
include(joinpath(@__DIR__, "ReducedP0.jl"))
using IntervalArithmetic, ForwardDiff, LinearAlgebra, Printf, Random, JuMP, Clarabel

TQ = length(ARGS) >= 1 ? ARGS[1] : "1180"
LAST = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
c = load_query(TQ)
P = reduced_p0(c)
S = F4Lib.SMDFRROM
p, h, K, NX = c.p, c.h, c.K, S.NX
leak = c.used_leak
m_lo, m_hi, y_ref, U_y = P.m_lo, P.m_hi, P.y_ref, P.U_y
Sc = c.S; iTw, iM = S.IDX.Tw, S.IDX.dM
outdir = joinpath(@__DIR__, "..", "results", "f3c_optimality"); mkpath(outdir)
log_io = open(joinpath(outdir, "rigorous_$(TQ)_$(get(ENV, "F3C_CENTER", "opt")).txt"), "w")
say(s...) = (println(s...); println(log_io, s...); flush(log_io))
fmt(x) = @sprintf("%.3e", x)
I_(x) = interval(x)
lo(x) = inf(x); hi(x) = sup(x)

say("=== F3c step 4b (rigorous residual bound), t_q = $TQ s, parts 1..$LAST")
nlp = solve_det(c, leak); nlp = solve_det(c, leak)
ustar = clamp.(nlp.m, m_lo, m_hi)
UB = P.J(ustar)
# Expansion point of the Taylor model (any point of the box gives a valid bound; the
# residual grows ~ |u − u_c|², so the box midpoint minimizes the worst case). The UB stays J(u*).
CENTER = get(ENV, "F3C_CENTER", "opt")
uc = CENTER == "mid" ? fill((m_lo + m_hi) / 2, K) : copy(ustar)
Xs = S.simulate(c.x_init, vcat(c.leak0, leak), vcat(0.0, uc), h, p; scheme = :be)   # rows 0..K
say(@sprintf("UB = J(u*) = %.10f; Taylor-model center: %s", UB, CENTER == "mid" ? "box midpoint" : "u* (NLP solution)"))

"""Jacobian of f w.r.t. x (A does not depend on u: the control enters only ΔM)."""
Jf(z, d) = ForwardDiff.jacobian(w -> S.rhs(w, d, 0.0, p), z)

# ------------------------------------------------------------------ part 1
say("\n4b.1 interval Jacobian [A](Z) over sampled reach boxes; it must contain every sampled A(x)")
rng = MersenneTwister(5)
NSAMP = 300
samples = [m_lo .+ (m_hi - m_lo) .* rand(rng, K) for _ in 1:NSAMP]
append!(samples, [fill(m_lo, K), fill(m_hi, K)])
Xsamp = [S.simulate(c.x_init, vcat(c.leak0, leak), vcat(0.0, u), h, p; scheme = :be) for u in samples]
dev = [maximum(abs.(X[k+1, :] .- Xs[k+1, :]) for X in Xsamp) for k in 1:K]   # per state
viol = 0; widths = Float64[]
t1 = @elapsed for k in 1:K
    Z = [interval(Xs[k+1, j] - 1.5dev[k][j], Xs[k+1, j] + 1.5dev[k][j]) for j in 1:NX]
    AI = Jf(Z, leak[k])
    for X in Xsamp[1:50]
        A = Jf(X[k+1, :], leak[k])
        global viol += count(!issubset_interval(I_(A[i, j]), AI[i, j]) for i in 1:NX, j in 1:NX)
    end
    As = Jf(Xs[k+1, :], leak[k])
    # scaled half-width of [A] relative to the scaled A*
    W = [radius(AI[i, j]) * Sc[j] / Sc[i] for i in 1:NX, j in 1:NX]
    push!(widths, opnorm(W, Inf) / opnorm([abs(As[i, j]) * Sc[j] / Sc[i] for i in 1:NX, j in 1:NX], Inf))
end
say(@sprintf("  %d steps, %d sampled Jacobians each: %d entries outside the interval Jacobian (must be 0); %.1f s", K, 50, viol, t1))
say(@sprintf("  relative half-width of [A] (scaled, ∞-norm): step 1 %s, step %d %s, max %s",
             fmt(widths[1]), K, fmt(widths[end]), fmt(maximum(widths))))
say("  max sampled state deviation over the box at step $K, scaled: " * fmt(maximum(dev[K] ./ Sc)))
LAST < 2 && (close(log_io); exit())

# ------------------------------------------------------------------ part 2
say("\n4b.2 propagation of the Taylor model with interval remainder (whole control box)")
dU = [interval(m_lo - uc[j], m_hi - uc[j]) for j in 1:K]
Bs = zeros(NX); Bs[iM] = 1 / Sc[iM]
Sinv = Diagonal(1 ./ Sc); Sd = Diagonal(Sc)
pv = p.poly
p4d(x) = 4pv[1] * x^3 + 3pv[2] * x^2 + 2pv[3] * x + pv[4]
p4dd(x) = 12pv[1] * x^2 + 6pv[2] * x + 2pv[3]
T = zeros(NX, K); r = fill(I_(0.0), NX)
eps_rig = zeros(K); Grec = zeros(K, K); ystar = zeros(K); iters = Int[]; ok_all = true
t2 = @elapsed for k in 1:K
    global T, r
    d = leak[k]; xs = Xs[k+1, :]; xp = Xs[k, :]
    As = Sinv * Jf(xs, d) * Sd
    Ms = I - h * As; C = inv(Ms)
    ek = zeros(K); ek[k] = 1.0
    rhsT = T .+ h .* Bs * ek'
    Tn = C * rhsT
    # rigorous enclosures of the floating-point defects
    Dk = I_.(Ms) * I_.(Tn) .- (I_.(T) .+ h .* I_.(Bs) * I_.(ek)')
    res = (I_.(xs) .- I_.(xp) .- h .* S.rhs(I_.(xs), d, I_(uc[k]), p)) ./ Sc
    E = Matrix{Float64}(I, NX, NX) .- I_.(C) * I_.(Ms)
    lin = I_.(Tn) * dU                                   # exact hull of the linear part
    fixed = I_.(C) * (r .- res .- Dk * dU)               # terms that do not depend on R
    R = fixed .+ interval(-1e-14, 1e-14)
    accepted = false
    for it in 1:60
        δ = lin .+ R
        Z = I_.(xs) .+ Sc .* δ
        ΔA = (Sinv * I_.(Jf(Z, d)) * Sd) .- I_.(As)
        Fr = fixed .+ I_.(C) * (h .* (ΔA * δ)) .+ E * R
        if all(issubset_interval.(Fr, R))
            R = Fr; accepted = true; push!(iters, it); break
        end
        # enlarge the candidate: hull with the image, inflated by 10 %
        R = [interval(min(lo(R[j]), lo(Fr[j])) - 0.1 * (hi(Fr[j]) - lo(Fr[j])) - 1e-15,
                      max(hi(R[j]), hi(Fr[j])) + 0.1 * (hi(Fr[j]) - lo(Fr[j])) + 1e-15) for j in 1:NX]
    end
    if !accepted
        say("  step $k: remainder enclosure NOT verified (no self-mapping box in 60 iterations)")
        global ok_all = false; break
    end
    r = R; T = Tn
    dd = d + 5550
    e_k = Sc[iTw] * r[iTw] + 0.5 * p4dd(I_(dd) + interval(m_lo, m_hi)) * pown(dU[k], 2)
    eps_rig[k] = mag(e_k)
    ystar[k] = S.output(xs, d, uc[k], p)
    Grec[k, :] = Sc[iTw] .* T[iTw, :]; Grec[k, k] += p4d(dd + uc[k])
end
say(@sprintf("  %s in %.1f s; fixed-point iterations per step: %s", ok_all ? "all $K steps verified" : "FAILED", t2,
             isempty(iters) ? "-" : string(minimum(iters), "–", maximum(iters))))
ok_all || (close(log_io); exit())
Gad = ForwardDiff.jacobian(u -> P.predict(u), uc)
say("  consistency: max |G_rec − G_AD| = " * fmt(maximum(abs.(Grec .- Gad))) *
    ", max |ŷ*_rec − ŷ*_AD| = " * fmt(maximum(abs.(ystar .- P.predict(uc)))))
# sampled residual of THIS affine model (same center), over the random + extreme samples
eps_s = zeros(K)
for u in samples
    eps_s .= max.(eps_s, abs.(P.predict(u) .- (ystar .+ Grec * (u .- uc))))
end
say("  k    ε_rigorous [K]   ε_sampled [K]   ratio   (sampled = $(length(samples)) points, same center)")
for k in (1, 2, 5, 10, 15, 20, 26)
    k > K && continue
    say(@sprintf("  %2d   %s        %s            %.2f", k, fmt(eps_rig[k]), fmt(eps_s[k]), eps_rig[k] / eps_s[k]))
end
say(@sprintf("  ε_rigorous ≥ ε_sampled in all steps: %s (it must be, if both are right)",
             all(eps_rig .>= eps_s .* (1 - 1e-9)) ? "yes" : "NO"))
LAST < 3 && (close(log_io); exit())

# ------------------------------------------------------------------ part 3
say("\n4b.3 certified lower bound: convex relaxation with the rigorous ε (Clarabel dual bound)")
md = Model(Clarabel.Optimizer); set_silent(md)
set_optimizer_attribute(md, "tol_gap_abs", 1e-12); set_optimizer_attribute(md, "tol_gap_rel", 1e-12)
@variable(md, m_lo <= u[1:K] <= m_hi)
@variable(md, t[1:K] >= 0); @variable(md, s[1:K-1] >= 0); @variable(md, rr[1:K] >= 0)
@expression(md, aa[k = 1:K], ystar[k] + sum(Grec[k, j] * (u[j] - uc[j]) for j in 1:k))
@constraint(md, [k = 1:K], t[k] >= aa[k] - y_ref - eps_rig[k])
@constraint(md, [k = 1:K], t[k] >= -(aa[k] - y_ref) - eps_rig[k])
@constraint(md, [k = 1:K], rr[k] >= aa[k] - eps_rig[k] - U_y)
@constraint(md, [k = 1:K-1], u[k+1] - u[k] <= P.dm_lim + s[k])
@constraint(md, [k = 1:K-1], u[k] - u[k+1] <= P.dm_lim + s[k])
@objective(md, Min, sum(t .^ 2) + P.lam_m * sum(s .^ 2) + P.lam_pct * sum(rr .^ 2))
optimize!(md)
LB = min(dual_objective_value(md), objective_value(md))
say(@sprintf("  %s: LB = %.10f, UB = %.10f, certified gap = %.4e (%.3f %% of UB)",
             termination_status(md), LB, UB, UB - LB, 100 * (UB - LB) / UB))
open(joinpath(outdir, "rigorous_eps_$(TQ)_$(CENTER).csv"), "w") do io
    println(io, "k,eps_rigorous,eps_sampled")
    for k in 1:K; println(io, k, ",", eps_rig[k], ",", eps_s[k]); end
end
open(joinpath(outdir, "rigorous_gap_$(TQ)_$(CENTER).csv"), "w") do io
    println(io, "UB,LB,gap_abs,gap_rel"); println(io, join((UB, LB, UB - LB, (UB - LB) / UB), ","))
end
close(log_io)
println("\nwritten to ", outdir)
