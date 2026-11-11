"""
    HRMQueries

Build MPC queries from the Julia HRM replica, without Simulink (F4.5 of
`docs/monolithic_formulation_plan.md`). Reproduces what the loop does at a query time t_q
(`matlab/DET_scheduler_v1_perc_online_LJ.m` and `matlab/DET_2.m`):

1. run the HRM from t = 0 to t_q − `delay` (10 s) with the control series known at t_q;
2. map the HRM state onto the ROM (`Raiz_to_Rom` of `DET_2.m`): N, Cᵢ, Cᵉᵢ and ΔM are
   copied; T_f, T_w and T_c of the 0-D ROM are the MEANS of the 12 HRM nodes;
3. warm the ROM up from t_q − `delay` to t_q − Δt with the HRM leak and the applied ACS
   currents (inputs at t_{n+1}, as in F1). The ROM state at t_q − Δt is the initial state
   of the monolithic NLP (`F4Lib.load_query` uses the same instant).

The ROM parameters are fixed across queries; they are read from the F1 export
(`results/f1_rom_replica/rom_reference_1180.mat`).
"""
module HRMQueries

include(joinpath(@__DIR__, "SMDFRHRM.jl"))
include(joinpath(@__DIR__, "SMDFRROM.jl"))
using .SMDFRHRM, .SMDFRROM, MAT

export rom_params, hrm_to_rom, build_query, warm_rom, ExtInput, hrm_params, SMDFRHRM, SMDFRROM

const ROOT = joinpath(@__DIR__, "..")

"""ROM parameters (identical in every F1 export)."""
rom_params(file = joinpath(ROOT, "results", "f1_rom_replica", "rom_reference_1180.mat")) =
    SMDFRROM.load_params(matread(file)["par"])

"""HRM state (50) → ROM state (17), as `DET_2.m` overwrites the ROM operating point."""
function hrm_to_rom(xH)
    I = SMDFRHRM.IDX; J = SMDFRROM.IDX
    xR = zeros(SMDFRROM.NX)
    xR[J.N] = xH[I.N]; xR[J.C] = xH[I.C]; xR[J.Ce] = xH[I.Ce]
    xR[J.Tf] = sum(xH[I.Tf]) / 12; xR[J.Tw] = sum(xH[I.Tw]) / 12; xR[J.Tc] = sum(xH[I.Tc]) / 12
    xR[J.dM] = xH[I.dM]
    xR
end

"""
    build_query(pH, pR, ext, t_q; delay = 10.0, h = 0.5)

`ext` is the root-inport series (u_ACS [A], active_mpc) known at t_q. Returns the HRM state
and PI integrator at t_q − delay, the ROM state there (`xR0`) and at t_q − Δt (`x_init`),
the warm-up inputs, and the HRM run itself (for the PCT and leak at t_q).
"""
function build_query(pH, pR, ext::ExtInput, t_q; delay = 10.0, h = 0.5)
    t0 = t_q - delay
    sH = SMDFRHRM.simulate(pH, ext, t0; h)
    w = warm_rom(pH, pR, sH.x, sH.logic.xI, ext, t_q; delay, h)
    (; t_q, t0, xH0 = sH.x, xI0 = sH.logic.xI, w..., hrm = sH)
end

"""
    warm_rom(pH, pR, xH0, xI0, ext, t_q; delay = 10.0, h = 0.5)

ROM state at t_q − Δt from the HRM state `xH0` at t_q − delay (PI integrator `xI0`):
projection `hrm_to_rom` and a warm-up with the HRM leak and the applied ACS currents
(ROM inports at t_{n+1}).
"""
function warm_rom(pH, pR, xH0, xI0, ext::ExtInput, t_q; delay = 10.0, h = 0.5)
    t0 = t_q - delay
    xR0 = hrm_to_rom(xH0)
    tw = t0:h:(t_q - h)
    leak = [SMDFRHRM.leak_flow(t, pH) for t in tw]
    cur = [value(ext, t, 1) for t in tw]
    act = [value(ext, t, 2) for t in tw]
    # With active_mpc = 0 the ROM's own PI drives the pump; below the pump threshold
    # (≈ 56.687 A) its flow is 0, which is the only case supported here.
    acs = [a > 0 ? SMDFRROM.current_to_massflow(u) : 0.0 for (u, a) in zip(cur, act)]
    if any(act .<= 0)
        @assert xI0 < 50 "PI current near the pump threshold during the warm-up: not supported"
    end
    X = SMDFRROM.simulate(xR0, leak, acs, h, pR; scheme = :ode1be, input_at = :next)
    (; xR0, x_init = X[end, :], X_warm = X, t_warm = collect(tw), leak_warm = leak, acs_warm = acs)
end

end # module
