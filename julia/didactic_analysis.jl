# Data for the didactic figures of docs/estado_mpc_monolitico.tex (no new results: the same
# validated replica, analysed so that a newcomer can see the structure of the problem).
#
#   eigenvalues_<t_q>.csv   eigenvalues of A = ∂f/∂x at the query state (stiffness)
#   bode_<t_q>.csv          frequency response G(iω) from ṁ_ACS to ŷ of the linearized ROM
#   pump_curve.csv          ACS pump map I [A] → ṁ [kg/s] (ECM of CurrentToMassFlow.m)
#   authority_<t_q>.csv     nonlinear PCT over the horizon for constant ACS flows
#
# Usage, from the repository root:  julia --project=julia julia/didactic_analysis.jl
include(joinpath(@__DIR__, "F4Lib.jl"))
using .F4Lib, ForwardDiff, LinearAlgebra, Printf
S = F4Lib.SMDFRROM
out = joinpath(@__DIR__, "..", "results", "didactic"); mkpath(out)

open(joinpath(out, "pump_curve.csv"), "w") do io
    println(io, "current_A,mdot_kgs")
    for I in range(56.60, 87.5; length = 400)
        println(io, I, ",", current_to_massflow(I))
    end
end

for TQ in ("1161", "1180")
    c = load_query(TQ); p = c.p; x0 = c.x_init; d0 = c.leak0
    um = (c.m_lo_arch + c.m_hi_arch) / 2
    A = ForwardDiff.jacobian(z -> S.rhs(z, d0, um, p), x0)
    B = ForwardDiff.derivative(v -> S.rhs(x0, d0, v, p), um)
    D = ForwardDiff.derivative(n -> S.polyval(p.poly, n), d0 + um + 5550)
    ev = eigvals(A)
    open(joinpath(out, "eigenvalues_$TQ.csv"), "w") do io
        println(io, "re,im")
        for e in ev; println(io, real(e), ",", imag(e)); end
    end
    Cv = zeros(S.NX); Cv[S.IDX.Tw] = 1.0
    open(joinpath(out, "bode_$TQ.csv"), "w") do io
        println(io, "omega,mag,phase_deg")
        for ω in 10 .^ range(-4, 3; length = 300)
            g = Cv' * ((im * ω * I - A) \ B) + D
            println(io, ω, ",", abs(g), ",", rad2deg(angle(g)))
        end
    end
    # nonlinear response of the PCT to constant ACS flows over the horizon (converged BE)
    levels = range(c.m_lo_arch, c.m_hi_arch; length = 5)
    Ys = Vector{Vector{Float64}}()
    for m in levels
        X = S.simulate(x0, vcat(d0, c.used_leak), vcat(0.0, fill(m, c.K)), c.h, p; scheme = :be)
        push!(Ys, [S.output(X[k+1, :], c.used_leak[k], m, p) for k in 1:c.K])
    end
    open(joinpath(out, "authority_$TQ.csv"), "w") do io
        println(io, "t," * join([@sprintf("m%.3f", m) for m in levels], ","))
        for k in 1:c.K
            println(io, join([c.tX[c.hor[k]]; [Ys[j][k] for j in eachindex(levels)]], ","))
        end
    end
    @printf("t_q=%s: eigenvalues Re ∈ [%.3e, %.3e], %d unstable; |G| at 1e-3 rad/s = %.3e; ΔPCT(end) between u_min and u_max = %.3f K\n",
            TQ, minimum(real.(ev)), maximum(real.(ev)), count(>(1e-9), real.(ev)),
            abs(Cv' * ((im * 1e-3 * I - A) \ B) + D), Ys[1][end] - Ys[end][end])
end
println("written to ", out)
