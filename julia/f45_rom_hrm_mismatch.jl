# F4.5, R4: ROM–HRM mismatch over the 41 archived queries t_q = 1160 … 1180 s.
#
# For every query (built from the Julia HRM replica exactly as the loop builds it):
#   1. solve the monolithic NLP N3 (julia/F4Lib.jl) on the ROM with a PERFECT leak
#      forecast, so that only the model error remains, for two actuator configurations:
#        :arch  archived current box [56.69, 56.75] A (≈ 0–12.9 kg/s), soft ramp (P0);
#        :phys  0 ≤ ṁ ≤ 5543 kg/s with the hard ramp |Δṁ| ≤ 5 kg/s per step (plan §8.3);
#   2. apply the optimal flow sequence to the HRM replica from its state at t_q − Δt;
#   3. compare the HRM PCT with the ROM prediction ŷ over the 26-step horizon.
# Also: the ROM "nowcast" error at t_q − Δt (known when deciding), the bias-corrected
# error, the actuator authority in both models, and the P0 cost the HRM actually incurs.
#
# Usage, from the repository root:
#   julia --project=julia julia/f45_rom_hrm_mismatch.jl
include(joinpath(@__DIR__, "HRMQueries.jl"))
include(joinpath(@__DIR__, "F4Lib.jl"))
using .HRMQueries, .F4Lib, MAT, Printf, Statistics

const ROOT = joinpath(@__DIR__, "..")
const OUT = joinpath(ROOT, "results", "f45_robustness"); mkpath(OUT)
const h = 0.5; const K = 26; const DELAY = 10.0

arc = matread(joinpath(ROOT, "results", "f1b_hrm_replica", "archive_1180.mat"))
cs = arc["control_series"]
pH = hrm_params(A = arc["sampled_A"], B = arc["sampled_B"], mdot_c_drift = arc["sampled_mdot_c_drift"])
ref = matread(joinpath(ROOT, "results", "f1_rom_replica", "rom_reference_1180.mat"))
pR = rom_params()                                   # ROM params for HRMQueries' warm-up
pR4 = F4Lib.SMDFRROM.load_params(ref["par"])        # same values, F4Lib's type
y_ref = Float64(ref["target"])
ub = ref["u_bounds"]
m_lo_arch = current_to_massflow(minimum(ub[:, 1])); m_hi_arch = current_to_massflow(maximum(ub[:, 2]))
ext_full = ExtInput(cs[2:end, 1], cs[2:end, 2:3])  # the applied series (all known)
const value = SMDFRHRM.value
cur_at(t) = value(ext_full, t, 1); act_at(t) = value(ext_full, t, 2)

# ---- true HRM trajectory: save (state, logic) at every 0.5 s from 1150 to 1179.5 s
println("HRM replica run 0 → 1150 s …")
s = SMDFRHRM.simulate(pH, ext_full, 1150.0)
saved = Dict{Float64,Any}(1150.0 => (s.x, s.logic, s.pct[end]))
for t in 1150.0:h:1179.0
    x, L, _ = saved[t]
    s2 = SMDFRHRM.simulate(pH, ext_full, t + h; t0 = t, x0 = x, logic = L)
    saved[t+h] = (s2.x, s2.logic, s2.pct[end])
end

"""HRM PCT over the horizon when the flows `m` (K) are applied from t_q on."""
function hrm_response(t_q, m)
    x, L, _ = saved[t_q-h]
    last = findfirst(==(t_q - h), cs[:, 1])
    tf = [t_q + (k - 1) * h for k in 1:K]
    u = [m[k] > 0 ? F4Lib.massflow_to_current(m[k]) : 56.0 for k in 1:K]   # 56 A: pump off
    ext = ExtInput(vcat(cs[2:last, 1], tf), vcat(cs[2:last, 2:3], hcat(u, ones(K))))
    r = SMDFRHRM.simulate(pH, ext, tf[end]; t0 = t_q - h, x0 = x, logic = L)
    @assert maximum(abs.(r.acs[2:end] .- m)) < 1e-6 * max(1, maximum(m)) "ACS flow not reproduced"
    r.pct[2:end]
end

cost(y, m) = sum((y .- y_ref) .^ 2) + F4Lib.lam_m * sum(max.(0, abs.(diff(m)) .- F4Lib.dm_lim) .^ 2) +
             F4Lib.lam_pct * sum(max.(0, y .- F4Lib.U_y) .^ 2)

rows = NamedTuple[]; steps = NamedTuple[]
for t_q in 1160.0:h:1180.0
    xH10, L10, _ = saved[t_q-DELAY]
    w = warm_rom(pH, pR, xH10, L10.xI, ext_full, t_q; delay = DELAY, h)
    x_init = w.x_init
    L = [SMDFRHRM.leak_flow(t_q + (k - 1) * h, pH) for k in 1:K]      # perfect forecast (ROM inport)
    m_prev = act_at(t_q - h) > 0 ? current_to_massflow(cur_at(t_q - h)) : 0.0
    c = (; p = pR4, h, K, S = max.(abs.(x_init), 1.0), x_init, leak0 = w.leak_warm[end], y_ref,
         m_prev, m_lo_arch, m_hi_arch)
    # nowcast: ROM estimate vs HRM PCT at t_q − Δt (both known at decision time)
    pct_now = saved[t_q-h][3]
    yhat_now = F4Lib.SMDFRROM.output(x_init, w.leak_warm[end], w.acs_warm[end], pR4)
    # authority: PCT at the end of the horizon, full flow vs no flow
    y0R = F4Lib.evaluate(c, zeros(K), L).y; y1R = F4Lib.evaluate(c, fill(m_hi_arch, K), L).y
    y0H = hrm_response(t_q, zeros(K)); y1H = hrm_response(t_q, fill(m_hi_arch, K))
    for (cfg, kw) in ((:arch, (m_lo = m_lo_arch, m_hi = m_hi_arch, ramp = :soft)),
                      (:phys, (m_lo = 0.0, m_hi = 5543.0, ramp = :hard)))
        sol = solve_det(c, L; kw...)
        ev = F4Lib.evaluate(c, sol.m, L)
        yH = hrm_response(t_q, sol.m)
        e = yH .- ev.y
        push!(rows, (; t_q, cfg, status = string(sol.status), solve_s = sol.time, m0 = sol.m[1],
                     mass = h * sum(sol.m), J_rom = ev.J, J_hrm = cost(yH, sol.m),
                     e_now = pct_now - yhat_now, e_max = maximum(abs.(e)), e_end = e[end],
                     e_bc_max = maximum(abs.(e .- (pct_now - yhat_now))),
                     pct_hrm_max = maximum(yH), yhat_max = maximum(ev.y),
                     auth_rom = y0R[end] - y1R[end], auth_hrm = y0H[end] - y1H[end]))
        for k in 1:K
            push!(steps, (; t_q, cfg, k, t = t_q + (k - 1) * h, m = sol.m[k], yhat = ev.y[k], pct_hrm = yH[k]))
        end
    end
    r = rows[end-1]
    @printf("t_q %.1f  nowcast err %+.3f K | arch: m0 %.2f, max|e| %.3f K, bias-corr %.3f K, J_rom %.3f J_hrm %.3f | auth ROM %.3f HRM %.3f K\n",
            t_q, r.e_now, r.m0, r.e_max, r.e_bc_max, r.J_rom, r.J_hrm, r.auth_rom, r.auth_hrm)
end

function writecsv(file, tab)
    open(file, "w") do io
        ks = keys(tab[1]); println(io, join(ks, ","))
        for r in tab; println(io, join((x isa AbstractFloat ? @sprintf("%.9g", x) : string(x) for x in values(r)), ",")); end
    end
end
writecsv(joinpath(OUT, "rom_hrm_mismatch.csv"), rows)
writecsv(joinpath(OUT, "rom_hrm_mismatch_steps.csv"), steps)

println("\nSummary over the 41 queries:")
for cfg in (:arch, :phys)
    R = filter(r -> r.cfg == cfg, rows)
    @printf("  %-5s  nowcast err [%.3f, %.3f] K | max|e| median %.3f, max %.3f K | bias-corrected max %.3f K | J_hrm/J_rom median %.2f | solve median %.3f s | status %s\n",
            cfg, minimum(r.e_now for r in R), maximum(r.e_now for r in R),
            median(r.e_max for r in R), maximum(r.e_max for r in R), maximum(r.e_bc_max for r in R),
            median(r.J_hrm / r.J_rom for r in R), median(r.solve_s for r in R), unique(r.status for r in R))
end
R = filter(r -> r.cfg == :arch, rows)
@printf("  authority (12.9 kg/s for 12.5 s): ROM %.3f K, HRM %.3f K (ratio HRM/ROM median %.3f)\n",
        median(r.auth_rom for r in R), median(r.auth_hrm for r in R), median(r.auth_hrm / r.auth_rom for r in R))
println("written results/f45_robustness/rom_hrm_mismatch{,_steps}.csv")
