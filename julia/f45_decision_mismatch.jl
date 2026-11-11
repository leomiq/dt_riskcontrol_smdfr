# F4.5, R4 (second part): does the ROM–HRM mismatch change the MPC DECISION?
#
# At a subset of the archived queries, three decisions are compared on the HRM replica
# (the plant), all with a perfect leak forecast:
#   ROM       the monolithic NLP N3 on the ROM (julia/F4Lib.jl), as in F2–F4;
#   ROM+bias  the same NLP with the output corrected by the nowcast error
#             d = PCT_HRM(t_q − Δt) − ŷ_ROM(t_q − Δt), held constant over the horizon
#             (offset-free MPC with a constant output disturbance);
#   HRM       an oracle ("N3-HRM"): the same P0 problem solved directly on the HRM replica
#             as a full-space NLP (50 states per step as variables, converged backward
#             Euler as constraints, exact Hessians, Ipopt), like N3 on the ROM.
# For each, the cost P0 the HRM actually incurs, the first move ṁ₀ and the PCT peak.
#
# Usage, from the repository root:
#   julia --project=julia julia/f45_decision_mismatch.jl
include(joinpath(@__DIR__, "HRMQueries.jl"))
include(joinpath(@__DIR__, "F4Lib.jl"))
using .HRMQueries, .F4Lib, MAT, Printf, Statistics, JuMP, Ipopt, ForwardDiff

const ROOT = joinpath(@__DIR__, "..")
const OUT = joinpath(ROOT, "results", "f45_robustness"); mkpath(OUT)
const h = 0.5; const K = 26; const DELAY = 10.0
const value_ext = SMDFRHRM.value

arc = matread(joinpath(ROOT, "results", "f1b_hrm_replica", "archive_1180.mat"))
cs = arc["control_series"]
pH = hrm_params(A = arc["sampled_A"], B = arc["sampled_B"], mdot_c_drift = arc["sampled_mdot_c_drift"])
ref = matread(joinpath(ROOT, "results", "f1_rom_replica", "rom_reference_1180.mat"))
pR = rom_params(); pR4 = F4Lib.SMDFRROM.load_params(ref["par"])
y_ref = Float64(ref["target"]); ub = ref["u_bounds"]
m_lo_arch = current_to_massflow(minimum(ub[:, 1])); m_hi_arch = current_to_massflow(maximum(ub[:, 2]))
ext_full = ExtInput(cs[2:end, 1], cs[2:end, 2:3])

println("HRM replica run 0 → 1150 s …")
s = SMDFRHRM.simulate(pH, ext_full, 1150.0)
saved = Dict{Float64,Any}(1150.0 => (s.x, s.logic, s.pct[end]))
for t in 1150.0:h:1179.0
    x, L, _ = saved[t]
    s2 = SMDFRHRM.simulate(pH, ext_full, t + h; t0 = t, x0 = x, logic = L)
    saved[t+h] = (s2.x, s2.logic, s2.pct[end])
end

"""HRM PCT over the horizon for flows m (pure function of m, AD-friendly). Same timing as
`SMDFRHRM.simulate`: leak at t_n, ACS flow of step k at t_{n+1} = t_q + (k−1)Δt."""
function hrm_pct(xq, t_q, m)
    x = xq; y = similar(m)
    for k in 1:K
        x = SMDFRHRM.step_ode1be(x, SMDFRHRM.leak_flow(t_q + (k - 2) * h, pH), m[k], h, pH)
        y[k] = SMDFRHRM.pct(x)
        @assert y[k] == x[SMDFRHRM.IDX.Tw[12]] "the hottest wall node is not node 12"
    end
    y
end
cost(y, m) = sum((y .- y_ref) .^ 2) + F4Lib.lam_m * sum(max.(0, abs.(diff(m)) .- F4Lib.dm_lim) .^ 2) +
             F4Lib.lam_pct * sum(max.(0, y .- F4Lib.U_y) .^ 2)

"""HRM right-hand side as JuMP expressions (same equations as `SMDFRHRM.rhs`)."""
function hrm_rhs_expr(md, x, leak, acs, p)
    I_ = SMDFRHRM.IDX; NN = SMDFRHRM.NN
    N = x[I_.N]; C = x[I_.C]; Ce = x[I_.Ce]; Tf = x[I_.Tf]; Tw = x[I_.Tw]; Tc = x[I_.Tc]; dM = x[I_.dM]
    Mtot = sum(p.M_c)
    M_act = [@expression(md, p.M_c[i] * (1 + dM / Mtot)) for i in 1:NN]
    mdot = @expression(md, p.mdot_c[1] * (1 + dM / Mtot))           # identical in every node
    h_wc = [@expression(md, p.h_wc[i] * (mdot / 5550)^0.8) for i in 1:NN]
    cp = [@expression(md, 175.1 - 4.961e-2 * Tc[i] + 1.985e-5 * Tc[i]^2 - 2.099e-9 * Tc[i]^3 - 1.524e6 / Tc[i]^2) for i in 1:NN]
    Tin_c = @expression(md, p.Te_c0[12] - (p.Q0 * mdot / 2775 - (5500 - mdot)) / (mdot * cp[1]))
    Te_c = Vector{Any}(undef, NN); Te_f = Vector{Any}(undef, NN)
    Te_c[1] = @expression(md, 2Tc[1] - Tin_c); Tin_f = p.Te_f0[12]
    Te_f[1] = @expression(md, 2Tf[1] - Tin_f)
    for i in 2:NN
        Te_c[i] = @expression(md, 2Tc[i] - Te_c[i-1]); Te_f[i] = @expression(md, 2Tf[i] - Te_f[i-1])
    end
    rho = @expression(md, p.rho_0 + p.alfa_f * ((Te_f[12] + Tin_f) / 2 - p.Tbar_f0) +
                          p.alfa_c * ((Te_c[12] + Tin_c) / 2 - p.Tbar_c0) + p.rhoex)
    P = @expression(md, p.P_0 / p.N_0 * N)
    f = Vector{Any}(undef, SMDFRHRM.NX)
    f[I_.N] = @expression(md, (rho - p.beta) * N / p.Lambda + sum(p.lambda_i[i] * C[i] for i in 1:6))
    for i in 1:6
        f[I_.C[i]] = @expression(md, p.beta_i[i] / p.Lambda * N - (p.lambda_i[i] + p.lambda_c) * C[i] + p.lambda_e * Ce[i])
        f[I_.Ce[i]] = @expression(md, p.lambda_c * C[i] - (p.lambda_i[i] + p.lambda_e) * Ce[i])
    end
    for i in 1:NN
        Tf_up = i == 1 ? Tin_f : Te_f[i-1]; Tc_up = i == 1 ? Tin_c : Te_c[i-1]
        f[I_.Tf[i]] = @expression(md, (P * p.Pi_0[i] + p.mdot_f[i] * p.c_pf[i] * (Tf_up - Te_f[i]) -
                                       p.h_fw[i] * p.A_fw[i] * (Tf[i] - Tw[i])) / (p.M_f[i] * p.c_pf[i]))
        f[I_.Tw[i]] = @expression(md, (p.h_fw[i] * p.A_fw[i] * (Tf[i] - Tw[i]) -
                                       h_wc[i] * p.A_wc[i] * (Tw[i] - Tc[i])) / (p.M_w[i] * p.c_pw[i]))
        f[I_.Tc[i]] = @expression(md, (mdot * cp[i] * (Tc_up - Te_c[i]) + h_wc[i] * p.A_wc[i] * (Tw[i] - Tc[i])) /
                                      (M_act[i] * cp[i]))
    end
    f[I_.dM] = @expression(md, leak + acs)
    f
end

"""
Oracle N3-HRM: P0 on the HRM as a full-space NLP (states are variables, converged backward
Euler as constraints, exact Hessians), like F4Lib.solve_det on the ROM. The PCT is T_w of
node 12, which is the hottest wall node in every query (checked by the caller).
"""
function solve_hrm(xq, t_q, m_prev, cfg, guess)
    m_lo, m_hi, ramp = cfg
    NX = SMDFRHRM.NX; S = max.(abs.(xq), 1.0); tw12 = SMDFRHRM.IDX.Tw[12]
    g = clamp.(guess, m_lo, m_hi)
    Lh = [SMDFRHRM.leak_flow(t_q + (k - 2) * h, pH) for k in 1:K]   # HRM: leak at t_n
    # initial guess of the states: ode1be trajectory under the guess
    Xg = zeros(K + 1, NX); Xg[1, :] = xq
    for k in 1:K; Xg[k+1, :] = SMDFRHRM.step_ode1be(Xg[k, :], Lh[k], g[k], h, pH); end
    md = Model(Ipopt.Optimizer); set_silent(md)
    set_optimizer_attribute(md, "tol", 1e-10); set_optimizer_attribute(md, "max_iter", 1000)
    @variable(md, m_lo <= m[k = 1:K] <= m_hi, start = g[k])
    z = @variable(md, [k = 0:K, j = 1:NX], start = Xg[k+1, j] / S[j])
    for j in 1:NX; fix(z[0, j], xq[j] / S[j]; force = true); end
    y = Vector{Any}(undef, K)
    for k in 1:K
        xk = [@expression(md, S[j] * z[k, j]) for j in 1:NX]
        f = hrm_rhs_expr(md, xk, Lh[k], m[k], pH)
        for j in 1:NX
            @constraint(md, z[k, j] - z[k-1, j] == h * f[j] / S[j])
        end
        y[k] = xk[tw12]
    end
    rr = @variable(md, [1:K], lower_bound = 0)
    @constraint(md, [k = 1:K], y[k] - F4Lib.U_y <= rr[k])
    obj = @expression(md, sum((y[k] - y_ref)^2 for k in 1:K) + F4Lib.lam_pct * sum(rr .^ 2))
    if ramp == :hard
        mm = [m_prev; m]
        @constraint(md, [k = 1:K], mm[k+1] - mm[k] <= F4Lib.dm_lim)
        @constraint(md, [k = 1:K], mm[k] - mm[k+1] <= F4Lib.dm_lim)
    else
        s = @variable(md, [1:K-1], lower_bound = 0)
        @constraint(md, [k = 1:K-1], m[k+1] - m[k] <= F4Lib.dm_lim + s[k])
        @constraint(md, [k = 1:K-1], m[k] - m[k+1] <= F4Lib.dm_lim + s[k])
        obj = @expression(md, obj + F4Lib.lam_m * sum(s .^ 2))
    end
    @objective(md, Min, obj)
    t = @elapsed optimize!(md)
    (m = clamp.(value.(m), m_lo, m_hi), status = termination_status(md), time = t)
end

rows = NamedTuple[]
for t_q in (1160.0, 1162.5, 1165.0, 1167.5, 1170.0, 1172.5, 1175.0, 1177.5, 1180.0)
    xH10, L10, _ = saved[t_q-DELAY]
    w = warm_rom(pH, pR, xH10, L10.xI, ext_full, t_q; delay = DELAY, h)
    L = [SMDFRHRM.leak_flow(t_q + (k - 1) * h, pH) for k in 1:K]
    m_prev = value_ext(ext_full, t_q - h, 2) > 0 ? current_to_massflow(value_ext(ext_full, t_q - h, 1)) : 0.0
    xq = saved[t_q-h][1]
    d = saved[t_q-h][3] - F4Lib.SMDFRROM.output(w.x_init, w.leak_warm[end], w.acs_warm[end], pR4)
    mk(yr) = (; p = pR4, h, K, S = max.(abs.(w.x_init), 1.0), x_init = w.x_init, leak0 = w.leak_warm[end],
              y_ref = yr, m_prev, m_lo_arch, m_hi_arch)
    for (cfgname, cfg) in ((:arch, (m_lo_arch, m_hi_arch, :soft)), (:phys, (0.0, 5543.0, :hard)))
        kw = (m_lo = cfg[1], m_hi = cfg[2], ramp = cfg[3])
        s_rom = solve_det(mk(y_ref), L; kw...)
        s_bias = solve_det(mk(y_ref - d), L; kw...)        # (ŷ + d − y_ref)² = (ŷ − (y_ref − d))²
        s_hrm = solve_hrm(xq, t_q, m_prev, cfg, s_bias.m)
        J(m) = cost(hrm_pct(xq, t_q, m), m)
        Jo = J(s_hrm.m)
        for (name, sol) in (("ROM", s_rom), ("ROM+bias", s_bias), ("HRM oracle", s_hrm))
            y = hrm_pct(xq, t_q, sol.m)
            push!(rows, (; t_q, cfg = cfgname, method = name, status = string(sol.status), time = sol.time,
                         m0 = sol.m[1], mass = h * sum(sol.m), J_hrm = J(sol.m), excess = J(sol.m) / Jo - 1,
                         pct_max = maximum(y), pct_end = y[end], d))
        end
        r = rows[end-2:end]
        @printf("t_q %.1f %-4s d=%+.3f | ṁ₀: ROM %7.3f, ROM+bias %7.3f, HRM %7.3f | J_hrm: %.4f / %.4f / %.4f (excess %+.1f%% / %+.1f%%) | %s %.1fs\n",
                t_q, cfgname, d, r[1].m0, r[2].m0, r[3].m0, r[1].J_hrm, r[2].J_hrm, r[3].J_hrm,
                100r[1].excess, 100r[2].excess, r[3].status, r[3].time)
    end
end

open(joinpath(OUT, "decision_mismatch.csv"), "w") do io
    ks = keys(rows[1]); println(io, join(ks, ","))
    for r in rows; println(io, join((x isa AbstractFloat ? @sprintf("%.9g", x) : string(x) for x in values(r)), ",")); end
end
println("\nMedian excess cost on the HRM over the HRM oracle:")
for cfg in (:arch, :phys), meth in ("ROM", "ROM+bias")
    R = filter(r -> r.cfg == cfg && r.method == meth, rows)
    O = filter(r -> r.cfg == cfg && r.method == "HRM oracle", rows)
    @printf("  %-4s %-9s median %+.2f%%, max %+.2f%% | max |Δṁ₀| vs oracle %.3f kg/s\n", cfg, meth,
            100median(r.excess for r in R), 100maximum(r.excess for r in R),
            maximum(abs(a.m0 - b.m0) for (a, b) in zip(R, O)))
end
println("written results/f45_robustness/decision_mismatch.csv")
