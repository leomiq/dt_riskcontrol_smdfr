# Reduced-space P0 (states eliminated by the simulation), differentiable with ForwardDiff.
# Shared by the F3c scripts. Include after F4Lib.jl and build with `reduced_p0(c)`.
#
#   ŷ(u): K PCT predictions from the replica with controls u (one per horizon step)
#   J(u) = Σ (ŷ_k − y_ref)² + λ_m Σ [|u_{k+1} − u_k| − Δṁ]_+² + λ_pct Σ [ŷ_k − U_y]_+²
# scheme = :be (converged implicit Euler, the monolithic NLP's model) or :ode1be (Simulink).
using ForwardDiff, LinearAlgebra

primal(v) = v isa ForwardDiff.Dual ? primal(ForwardDiff.value(v)) : v
hinge2(z) = z > 0 ? z^2 : zero(z)

function reduced_p0(c; leak = c.used_leak, U_y = 1244.0, lam_m = 1e-4, dm_lim = 5.0, lam_pct = 10.0)
    p, h, K, x0, y_ref = c.p, c.h, c.K, c.x_init, c.y_ref
    S = F4Lib.SMDFRROM
    function step_be_ad(x, d, u)          # Newton on primal values, then dual refinements
        z = S.step_be(primal.(x), d, primal(u), h, p)
        T = promote_type(eltype(x), typeof(u)); zT = T.(z)
        for _ in 1:3
            r = zT .- x .- h .* S.rhs(zT, d, u, p)
            zT = zT .- (I - h * ForwardDiff.jacobian(w -> S.rhs(w, d, u, p), zT)) \ r
        end
        zT
    end
    step_ode1be_ad(x, d, u) = S.step_ode1be(x, d, u, h, p)
    function predict(u::AbstractVector{T}; scheme = :be) where {T}
        step = scheme === :be ? step_be_ad : step_ode1be_ad
        x = T.(x0); y = Vector{T}(undef, K)
        for k in 1:K
            x = step(x, leak[k], u[k]); y[k] = S.output(x, leak[k], u[k], p)
        end
        y
    end
    cost(u, y) = sum((y .- y_ref) .^ 2) + lam_m * sum(hinge2(abs(u[k+1] - u[k]) - dm_lim) for k in 1:K-1) +
                 lam_pct * sum(hinge2.(y .- U_y))
    J(u; scheme = :be) = cost(u, predict(u; scheme = scheme))
    (; predict, J, cost, K, y_ref, U_y, lam_m, dm_lim, lam_pct, m_lo = c.m_lo_arch, m_hi = c.m_hi_arch)
end
