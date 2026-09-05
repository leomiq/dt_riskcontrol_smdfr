"""Compare the penalised closed-loop realization against the unpenalised pilot.

Produces the two-curve figure: PCT versus time to t=1190, with and without the
DSM safety-margin penalty, plus the applied current / coolant mass flow and the
DSM trace showing where the margin activates.

Baseline : phase2c_fullloop_pilot_LJ_seed1.mat        (lambda = 0)
Penalised: phase2c_fullloop_penalised_*_LJ_seed1.mat  (dsm_thr 0.00238, lambda 1e6)

Same solver (LJ), same seed, same accident stream, same bounds -- the cost
function is the only difference.
"""
import os
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat

# Repository-relative by default; override with the environment variables when
# the run artefacts live elsewhere (e.g. a scratch directory on a cluster).
REPO = Path(__file__).resolve().parent.parent
P = Path(os.environ.get("P2C_RUN_DIR", REPO / "results"))
OUT = Path(os.environ.get("P2C_OUT_DIR", REPO / "results"))

TARGET = 1179.2712      # CalculateTs(100, 543)
U_Y = 1244.0
SPAN = U_Y - TARGET
DSM_THR = 0.00237737


def load_run(path):
    d = loadmat(path, squeeze_me=True, struct_as_record=False)
    return {
        "t": np.asarray(d["time_snapshots"], float),
        "pct": np.asarray(d["PCT"], float),
        "u": np.asarray(d["u_V2P_trajectory"], float),
        "rt": np.asarray(d.get("optimizer_runtimes", []), float),
    }


def current_to_massflow(u):
    """Mirror of CurrentToMassFlow's operating range, for reporting only.

    The analytical model maps [56.69, 56.75] A onto roughly [0.58, 12.86] kg/s
    strongly nonlinearly; a linear interpolation over the measured endpoints is
    enough for the figure's secondary axis. Measured anchors from
    test_control_authority.m.
    """
    anchors_u = np.array([56.69, 56.696, 56.702, 56.708, 56.714, 56.72,
                          56.726, 56.732, 56.738, 56.744, 56.75])
    anchors_m = np.array([0.5798, 1.8086, 3.0373, 4.2660, 5.4946, 6.7231,
                          7.9516, 9.1800, 10.4083, 11.6366, 12.8648])
    return np.interp(u, anchors_u, anchors_m)


def main():
    base_f = P / "phase2c_fullloop_pilot_LJ_seed1.mat"
    cands = sorted(P.glob("phase2c_fullloop_penalised*_LJ_seed1.mat"))
    if not base_f.exists():
        sys.exit(f"baseline missing: {base_f}")
    if not cands:
        sys.exit("penalised realization not found yet -- still running?")
    pen_f = cands[-1]

    base = load_run(base_f)
    pen = load_run(pen_f)
    print(f"baseline : {base_f.name}")
    print(f"penalised: {pen_f.name}\n")

    for name, r in (("baseline", base), ("penalised", pen)):
        dsm = (r["pct"] - TARGET) / SPAN
        m = current_to_massflow(r["u"])
        print(f"=== {name} ===")
        print(f"  queries      : {len(r['t'])}  t [{r['t'][0]:.1f} .. {r['t'][-1]:.1f}]")
        print(f"  PCT          : first {r['pct'][0]:.4f}  last {r['pct'][-1]:.4f}  "
              f"max {r['pct'].max():.4f}  min {r['pct'].min():.4f}")
        print(f"  DSM          : max {dsm.max():.6f}  last {dsm[-1]:.6f}")
        print(f"  samples over threshold ({DSM_THR:.5f}): {int((dsm > DSM_THR).sum())} of {len(dsm)}")
        print(f"  u            : mean {r['u'].mean():.5f}  range [{r['u'].min():.5f}, {r['u'].max():.5f}]")
        print(f"  mass flow    : mean {m.mean():.4f}  total {m.sum():.2f}")
        print(f"  |PCT-target| : mean {np.abs(r['pct']-TARGET).mean():.4f} K")
        if r["rt"].size:
            print(f"  optimizer    : {r['rt'].sum()/3600:.2f} h")
        print()

    # --- the comparison that matters ---
    n = min(len(base["t"]), len(pen["t"]))
    dm = current_to_massflow(pen["u"][:n]) - current_to_massflow(base["u"][:n])
    dpct = pen["pct"][:n] - base["pct"][:n]
    dsm_b = (base["pct"][:n] - TARGET) / SPAN
    dsm_p = (pen["pct"][:n] - TARGET) / SPAN

    print("=== EFFECT OF THE PENALTY ===")
    print(f"extra coolant injected (sum of per-query difference): {dm.sum():+.2f}")
    print(f"mean extra coolant per query                        : {dm.mean():+.4f}")
    print(f"queries with MORE coolant than baseline             : {int((dm > 0).sum())} of {n}")
    print(f"mean PCT difference                                 : {dpct.mean():+.4f} K")
    print(f"peak DSM  baseline {dsm_b.max():.6f} -> penalised {dsm_p.max():.6f} "
          f"({dsm_p.max()-dsm_b.max():+.6f})")
    print(f"mean DSM  baseline {dsm_b.mean():.6f} -> penalised {dsm_p.mean():.6f}")
    print(f"tracking |PCT-target| baseline {np.abs(base['pct'][:n]-TARGET).mean():.4f} K "
          f"-> penalised {np.abs(pen['pct'][:n]-TARGET).mean():.4f} K")

    print("\n=== per-query detail (first 20) ===")
    print(f"{'t':>8} {'PCT_base':>11} {'PCT_pen':>11} {'dPCT':>9} "
          f"{'DSM_base':>10} {'DSM_pen':>10} {'u_base':>10} {'u_pen':>10} {'dmdot':>9}")
    for i in range(min(20, n)):
        print(f"{base['t'][i]:8.1f} {base['pct'][i]:11.4f} {pen['pct'][i]:11.4f} "
              f"{dpct[i]:+9.4f} {dsm_b[i]:10.6f} {dsm_p[i]:10.6f} "
              f"{base['u'][i]:10.5f} {pen['u'][i]:10.5f} {dm[i]:+9.4f}")

    np.savez(OUT / "penalised_vs_baseline.npz",
             t=base["t"][:n], pct_base=base["pct"][:n], pct_pen=pen["pct"][:n],
             u_base=base["u"][:n], u_pen=pen["u"][:n],
             dsm_base=dsm_b, dsm_pen=dsm_p, dsm_thr=DSM_THR,
             target=TARGET, U_y=U_Y)
    print(f"\nsaved {OUT / 'penalised_vs_baseline.npz'}")


if __name__ == "__main__":
    main()
