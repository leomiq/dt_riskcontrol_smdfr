# F3b: the linear MPC (N1, N2) in sparse-monolithic and condensed form, solved with a
# general NLP solver (Ipopt) and a convex conic solver (Clarabel), with dual certificates;
# and a test of why N2 stalls short of the NLP at t_q = 1180 s.
#
# Notation: k = horizon step (1..K), i = SQP iteration. In N2 the linearization point of
# step k at iteration i is the state x̄_k^(i−1) of the nonlinear trajectory simulated with
# the incumbent control of iteration i−1 (all K Jacobians are computed before the QP).
#
# Deviation form (numerically clean: no f̄ − A x̄ cancellation with N ~ 1e13, C ~ 1e16):
#   δ_k = x_k − x̄_k,  Δu_k = u_k − ū_k,  Δd_k = d_k − d̄_k
#   (I − h A_k) δ_k = δ_{k−1} + r_k + h B_k Δu_k + h E_k Δd_k,   r_k = x̄_{k−1} − x̄_k + h f̄_k
#   ŷ_k = ȳ_k + δ_k[T_w] + D_k (Δd_k + Δu_k)
# r_k = 0 when the nominal trajectory is a converged implicit-Euler (BE) solution.
# Sparse form: δ_1..δ_K are variables (monolithic, as the NLP). Condensed form: δ is
# eliminated by the recursion, ŷ = ŷ⁰ + G u with G ∈ R^{K×K} lower triangular, and only
# u (plus the hinge slacks) remain.
#
# Usage, from the repository root:  julia --project=julia julia/f3_certificates_condensed.jl [t_q]
include(joinpath(@__DIR__, "F4Lib.jl"))
using .F4Lib, .F4Lib.SMDFRROM, JuMP, Ipopt, Clarabel, ForwardDiff, LinearAlgebra, Printf, MAT

TQ = isempty(ARGS) ? "1180" : ARGS[1]
c = load_query(TQ)
p, h, K, S, x0 = c.p, c.h, c.K, c.S, c.x_init
leak, leak0 = c.used_leak, c.leak0                    # the forecast used by F2/F3
m_lo, m_hi = c.m_lo_arch, c.m_hi_arch
m_mid = (m_lo + m_hi) / 2
y_ref = c.y_ref; U_y = 1244.0; lam_m = 1e-4; dm_lim = 5.0; lam_pct = 10.0
Sd = Diagonal(S)
Tw = SMDFRROM.IDX.Tw
outdir = joinpath(@__DIR__, "..", "results", "f3_certificates"); mkpath(outdir)

poly(n) = SMDFRROM.polyval(p.poly, n)
cost_P0(y, m) = sum((y .- y_ref) .^ 2) + lam_m * sum(max.(0, abs.(diff(m)) .- dm_lim) .^ 2) +
                lam_pct * sum(max.(0, y .- U_y) .^ 2)
"""Nonlinear re-simulation: scheme :ode1be (= Simulink) or :be (converged, = the NLP model)."""
function resim(m; scheme = :ode1be)
    X = simulate(x0, vcat(leak0, leak), vcat(0.0, m), h, p; scheme = scheme, input_at = :next)
    y = [SMDFRROM.output(X[k+1, :], leak[k], m[k], p) for k in 1:K]
    (y = y, J = cost_P0(y, m), X = X)
end

"""Linearization data of step k at (x̄, d̄, ū)."""
function lin(x, d, u)
    (f = rhs(x, d, u, p), A = ForwardDiff.jacobian(z -> rhs(z, d, u, p), x),
     B = ForwardDiff.derivative(v -> rhs(x, d, v, p), u), E = ForwardDiff.derivative(v -> rhs(x, v, u, p), d),
     D = ForwardDiff.derivative(poly, d + u + 5550), y = SMDFRROM.output(x, d, u, p))
end

"""Linearization points: N1 = the query point for every k; N2 = a nominal trajectory X (rows 1..K+1)."""
points_N1() = (xbar0 = x0, pts = fill((x0, leak0, m_mid), K))
points_N2(m_inc; scheme) = (X = resim(m_inc; scheme = scheme).X;
                            (xbar0 = X[1, :], pts = [(X[k+1, :], leak[k], m_inc[k]) for k in 1:K]))

"""Deviation-form data shared by both QP forms (scaled: δ = S ζ)."""
function qp_data(P)
    L = [lin(pt...) for pt in P.pts]
    xbar = vcat([P.xbar0], [pt[1] for pt in P.pts])       # x̄_0 … x̄_K (vector of states)
    blocks = map(1:K) do k
        l = L[k]; _, db, ub = P.pts[k]
        M = Sd \ (I - h * l.A) * Sd                       # scaled (I − hA_k)
        r = Sd \ (xbar[k] - xbar[k+1] + h * l.f)          # scaled residual
        b = Sd \ (h * l.B); e = Sd \ (h * l.E) * (leak[k] - db)
        (M = M, r = r, b = b, e = e, ub = ub, db = db, D = l.D, ybar = l.y)
    end
    (blocks = blocks, zeta0 = Sd \ (x0 - P.xbar0))
end

function add_cost!(md, u, y)
    s = @variable(md, [1:K-1], lower_bound = 0); r = @variable(md, [1:K], lower_bound = 0)
    @constraint(md, [k = 1:K-1], u[k+1] - u[k] <= dm_lim + s[k])
    @constraint(md, [k = 1:K-1], u[k] - u[k+1] <= dm_lim + s[k])
    @constraint(md, [k = 1:K], y[k] - U_y <= r[k])
    @objective(md, Min, sum((y[k] - y_ref)^2 for k in 1:K) + lam_m * sum(s .^ 2) + lam_pct * sum(r .^ 2))
end

function optimizer(name)
    if name == :ipopt
        md = Model(Ipopt.Optimizer); set_silent(md); set_optimizer_attribute(md, "tol", 1e-10)
    else
        md = Model(Clarabel.Optimizer); set_silent(md)
        set_optimizer_attribute(md, "tol_gap_abs", 1e-12); set_optimizer_attribute(md, "tol_gap_rel", 1e-12)
        set_optimizer_attribute(md, "tol_feas", 1e-12)
    end
    md
end

function finish(md, u, y, t)
    st = termination_status(md)
    primal = objective_value(md)
    dual = try dual_objective_value(md) catch; NaN end
    (u = value.(u), y = value.(y), primal = primal, dual = dual, gap = primal - dual,
     status = st, time = t, nvar = num_variables(md))
end

"""Sparse monolithic QP: δ (scaled ζ) are variables, the linear dynamics are constraints."""
function solve_sparse(Q, solver)
    md = optimizer(solver)
    @variable(md, m_lo <= u[1:K] <= m_hi, start = m_mid)
    @variable(md, ζ[1:K, 1:NX])
    for k in 1:K
        B = Q.blocks[k]; prev = k == 1 ? Q.zeta0 : ζ[k-1, :]
        @constraint(md, B.M * ζ[k, :] .== prev .+ B.r .+ B.b .* (u[k] - B.ub) .+ B.e)
    end
    @expression(md, y[k = 1:K], Q.blocks[k].ybar + S[Tw] * ζ[k, Tw] +
                                Q.blocks[k].D * ((leak[k] - Q.blocks[k].db) + (u[k] - Q.blocks[k].ub)))
    add_cost!(md, u, y)
    t = @elapsed optimize!(md)
    finish(md, u, y, t)
end

"""Condensation: ŷ = y0 + G u, by propagating ζ_k = F_k ζ_{k−1} + (affine in u_k)."""
function condense(Q)
    Z0 = copy(Q.zeta0); Zu = zeros(NX, K)
    y0 = zeros(K); G = zeros(K, K)
    for k in 1:K
        B = Q.blocks[k]; Mf = lu(B.M)
        Z0 = Mf \ (Z0 .+ B.r .- B.b .* B.ub .+ B.e)
        Zu = Mf \ Zu; Zu[:, k] .+= Mf \ B.b
        y0[k] = B.ybar + S[Tw] * Z0[Tw] + B.D * ((leak[k] - B.db) - B.ub)
        G[k, :] = S[Tw] .* Zu[Tw, :]; G[k, k] += B.D
    end
    (y0 = y0, G = G)
end

function solve_condensed(Q, solver; C = condense(Q))
    md = optimizer(solver)
    @variable(md, m_lo <= u[1:K] <= m_hi, start = m_mid)
    @expression(md, y[k = 1:K], C.y0[k] + sum(C.G[k, j] * u[j] for j in 1:k))
    add_cost!(md, u, y)
    t = @elapsed optimize!(md)
    (finish(md, u, y, t)..., G = C.G)
end

# --- NLP reference (N3, converged implicit Euler, Ipopt) and its timing
nlp = solve_det(c, leak)                                  # compile
nlp = solve_det(c, leak)
m_nlp = nlp.m
J_nlp = resim(m_nlp).J
@printf("t_q = %s s. NLP (N3, Ipopt): J_BE = %.8f, J_ode1be = %.8f, %.3f s\n\n", TQ, nlp.obj, J_nlp, nlp.time)

# --- N1: four solver/form combinations, certificates
rows = Any[]
Q1 = qp_data(points_N1())
println("N1 (LTI)   form       solver     status            vars   time [s]   primal          dual            gap        max|u−u_NLP| [kg/s]  J re-sim")
for form in (:sparse, :condensed), solver in (:ipopt, :clarabel)
    f = form == :sparse ? solve_sparse : solve_condensed
    f(Q1, solver)                                         # compile
    s = f(Q1, solver)
    du = maximum(abs.(s.u .- m_nlp)); J = resim(s.u).J
    @printf("           %-10s %-10s %-17s %4d   %8.4f   %.10f  %.10f  %+.1e   %.2e             %.8f\n",
            form, solver, string(s.status), s.nvar, s.time, s.primal, s.dual, s.gap, du, J)
    push!(rows, ("N1", string(form), string(solver), string(s.status), s.nvar, s.time, s.primal, s.dual, s.gap, du, s.u[1] - m_nlp[1], J))
end

# --- N2: SQP iterations, linearizing along ode1be (as F3) or converged BE trajectories
traj = Dict{String,Vector{Float64}}("NLP" => m_nlp)
traj["N1"] = solve_condensed(Q1, :clarabel).u
iters = Any[]
for scheme in (:ode1be, :be)
    println("\nN2 (LTV), nominal trajectories simulated with $scheme (condensed QP, Clarabel):")
    println("  it   max|u_i − u_{i−1}|   max|u_i − u_NLP|   J_BE(u_i)       J_ode1be(u_i)   gap dual    time [s]")
    m_it = fill(m_mid, K)
    for it in 1:8
        s = solve_condensed(qp_data(points_N2(m_it; scheme = scheme)), :clarabel)
        step = maximum(abs.(s.u .- m_it)); du = maximum(abs.(s.u .- m_nlp))
        Jbe, Jo = resim(s.u; scheme = :be).J, resim(s.u).J
        @printf("  %2d   %.3e            %.3e          %.8f      %.8f      %+.1e    %.4f\n", it, step, du, Jbe, Jo, s.gap, s.time)
        push!(iters, (string(scheme), it, step, du, Jbe, Jo, s.gap, s.time))
        m_it = s.u
        step < 1e-8 && break
    end
    traj["N2_" * string(scheme)] = m_it
    if scheme == :be
        Q2 = qp_data(points_N2(m_it; scheme = :be))
        for form in (:sparse, :condensed), solver in (:ipopt, :clarabel)
            f = form == :sparse ? solve_sparse : solve_condensed
            s = f(Q2, solver)
            push!(rows, ("N2_BE", string(form), string(solver), string(s.status), s.nvar, s.time, s.primal, s.dual,
                         s.gap, maximum(abs.(s.u .- m_nlp)), s.u[1] - m_nlp[1], resim(s.u).J))
        end
    end
end
@printf("\nNLP on the same converged-BE model: J_BE = %.8f (the N2-BE fixed point should match it)\n", nlp.obj)

# --- outputs
open(joinpath(outdir, "qp_solvers_$TQ.csv"), "w") do io
    println(io, "model,form,solver,status,nvar,time_s,primal,dual,gap,max_du_vs_NLP,du0_vs_NLP,J_resim_ode1be")
    for r in rows; println(io, join(r, ",")); end
    println(io, join(("N3_NLP", "sparse", "ipopt", "LOCALLY_SOLVED", 0, nlp.time, nlp.obj, NaN, NaN, 0.0, 0.0, J_nlp), ","))
end
open(joinpath(outdir, "n2_iterations_$TQ.csv"), "w") do io
    println(io, "nominal_scheme,iter,step,max_du_vs_NLP,J_BE,J_ode1be,gap,time_s")
    for r in iters; println(io, join(r, ",")); end
end
names = ["NLP", "N1", "N2_ode1be", "N2_be"]
open(joinpath(outdir, "trajectories_$TQ.csv"), "w") do io
    println(io, "t," * join(names, ","))
    for k in 1:K
        println(io, join([@sprintf("%.1f", c.tX[c.hor[k]]); [@sprintf("%.10f", traj[n][k]) for n in names]], ","))
    end
end
# sparsity patterns for the figure: equality-constraint matrix of the sparse N1 QP in
# (ζ_1..ζ_K, u) and the condensed G
B1 = Q1.blocks
I_, J_ = Int[], Int[]
for k in 1:K, (a, b) in Tuple.(findall(!iszero, B1[k].M))
    push!(I_, (k - 1) * NX + a); push!(J_, (k - 1) * NX + b)
end
for k in 2:K, a in 1:NX
    push!(I_, (k - 1) * NX + a); push!(J_, (k - 2) * NX + a)
end
for k in 1:K, a in findall(!iszero, B1[k].b)
    push!(I_, (k - 1) * NX + a); push!(J_, K * NX + k)
end
matwrite(joinpath(outdir, "sparsity_$TQ.mat"), Dict("rows" => I_, "cols" => J_, "nrows" => K * NX,
         "ncols" => K * NX + K, "G" => condense(Q1).G))
println("\nwritten to ", outdir)
