# F5b.1d (plan §8.18), figure data: time to 1220 K on the HRM when a failure leaves the ACS
# delivering η·ṁ_max, along the N2+bias trajectory of seed 201 (DET-accelerated leak draw),
# for η on a fine grid and ṁ_max ∈ {30, 40} kg/s. Explains why partial failures add almost no
# risk: τ_η grows quickly with η and P(D > τ_η) decays exponentially.
# Output: results/f5b1d_partial/tau_eta.csv.
# Usage: julia --project=julia --pkgimages=no --heap-size-hint=1200M julia/f5b1d_tau_eta.jl
using Printf, Random
include(joinpath(@__DIR__, "ClosedLoop.jl"))
using .ClosedLoop
const CL = ClosedLoop

OUT = joinpath(@__DIR__, "..", "results", "f5b1d_partial"); mkpath(OUT)
seed = 201
rng = MersenneTwister(seed)
A = 0.1 + 0.1rand(rng); Bp = 10 + 5rand(rng); d = 15 + 10rand(rng); lead = 17rand(rng)
open(joinpath(OUT, "tau_eta.csv"), "w") do io
    println(io, "m_max,t,leak,pct,eta,tau")
    for mmax in (30.0, 40.0)
        pH = CL.HRMQueries.hrm_params(; A, B = Bp, mdot_c_drift = d, acs_max = mmax, acs_ramp = 5.0)
        B = setup(; pH, lead)
        lg = run_loop(B, :N2; n = 481, act = act_cap(mmax, 5), bias = true, keep_states = true)
        for i in (1, 61, 121, 241, 361, 481), η in 0:0.125:0.875
            τ = time_to_limit(pH, lg.states[i][1], lg.states[i][2], lg.t_q[i]; t_cap = 600.0, flow = η * mmax)
            leak = -CL.HRMQueries.SMDFRHRM.leak_flow(lg.t_q[i], pH)
            println(io, join([@sprintf("%.6g", x) for x in (mmax, lg.t_q[i] - lg.t_q[1], leak, lg.pct[i], η, τ)], ","))
            @printf("m_max %2.0f t %5.1f leak %5.2f η %.3f τ %6.1f\n", mmax, lg.t_q[i] - lg.t_q[1], leak, η, τ); flush(stdout)
        end
    end
end
