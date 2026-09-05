"""Two separate comparisons, kept apart on purpose.

FIGURE 1 -- effect of the DSM safety-margin penalty (the experiment just run)
    LJ lambda=0  vs  LJ lambda=1e6, same solver, same seed, same accident.
    The cost function is the only variable, so any difference in the PCT
    trajectory is attributable to the penalty.

FIGURE 2 -- solver comparison (from the archived KPIs, no new runs needed)
    LJ / GA / PID on the same scenario. Note LJ and GA share min_safety_margin
    and peak_PCT to the 15th digit, because with a single-node DET both
    reconstruct the same root trajectory; they differ in control smoothness.
    PID is the non-optimised reference and cannot carry a J_psm term at all
    (its Kp/Ki/Kd law does not minimise a cost).

Mixing the two axes would confound "effect of the constraint" with "effect of
the search algorithm", which is why they stay in separate figures.
"""
import csv
import os
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat

# Repository-relative by default; override with the environment variables when
# the run artefacts live elsewhere (e.g. a scratch directory on a cluster).
REPO = Path(__file__).resolve().parent.parent
WORK = Path(os.environ.get("P2C_RUN_DIR", REPO / "results"))
EVID = Path(os.environ.get("P2C_OUT_DIR", REPO / "results"))

TARGET = 1179.2712
U_Y = 1244.0
SPAN = U_Y - TARGET
DSM_THR = 0.00237737

# Measured anchors from test_control_authority.m / sweep stage 1
ANCHOR_U = np.array([56.69, 56.696, 56.702, 56.708, 56.714, 56.72,
                     56.726, 56.732, 56.738, 56.744, 56.75])
ANCHOR_M = np.array([0.5798, 1.8086, 3.0373, 4.2660, 5.4946, 6.7231,
                     7.9516, 9.1800, 10.4083, 11.6366, 12.8648])


def massflow(u):
    return np.interp(u, ANCHOR_U, ANCHOR_M)


def load_run(path):
    d = loadmat(path, squeeze_me=True, struct_as_record=False)
    return {
        "t": np.asarray(d["time_snapshots"], float),
        "pct": np.asarray(d["PCT"], float),
        "u": np.asarray(d["u_V2P_trajectory"], float),
        "rt": np.asarray(d.get("optimizer_runtimes", []), float),
    }


def load_checkpoint(path):
    """A checkpoint holds u and timings but not the closed-loop PCT.

    run_phase2c_full_loop writes checkpoints with '-v7.3' (HDF5), which scipy's
    loadmat cannot read, so use h5py and fall back to loadmat for older files.
    """
    try:
        import h5py
        with h5py.File(path, "r") as f:
            n = int(np.asarray(f["snapshot_id_done"]).ravel()[0])
            t = np.asarray(f["time_snapshots"]).ravel().astype(float)
            u = np.asarray(f["u_V2P_trajectory"]).ravel().astype(float)
            rt = np.asarray(f["optimizer_runtimes"]).ravel().astype(float)
    except (ImportError, OSError):
        d = loadmat(path, squeeze_me=True, struct_as_record=False)
        n = int(d["snapshot_id_done"])
        t = np.asarray(d["time_snapshots"], float).ravel()
        u = np.asarray(d["u_V2P_trajectory"], float).ravel()
        rt = np.asarray(d["optimizer_runtimes"], float).ravel()
    return {"n_done": n, "t": t[:n], "u": u[:n], "rt": rt[:n]}


def figure1():
    print("=" * 74)
    print("FIGURE 1 -- DSM penalty effect (LJ lambda=0 vs LJ lambda=1e6)")
    print("=" * 74)

    base_f = WORK / "phase2c_fullloop_pilot_LJ_seed1.mat"
    if not base_f.exists():
        print(f"  baseline missing: {base_f}")
        return
    base = load_run(base_f)

    # run_phase2c_full_loop names the final file "<saveTag>_<solver>_seed<n>.mat",
    # so the penalised run lands as penalised_thr..._lam..._LJ_seed1.mat
    pen_files = [p for p in sorted(WORK.glob("penalised_*_LJ_seed1.mat"))
                 if "checkpoint" not in p.name]
    ckpts = sorted(WORK.glob("penalised_*_LJ_seed1_checkpoint.mat"))

    if pen_files:
        pen = load_run(pen_files[-1])
        complete = True
        print(f"  penalised run: {pen_files[-1].name} (complete)")
    elif ckpts:
        pen = load_checkpoint(ckpts[-1])
        complete = False
        print(f"  penalised run: {ckpts[-1].name} "
              f"(checkpoint, {pen['n_done']} queries -- control only, no closed-loop PCT)")
    else:
        print("  no penalised run or checkpoint found")
        return

    n = min(len(base["t"]), len(pen["t"]))
    mb, mp = massflow(base["u"][:n]), massflow(pen["u"][:n])

    print(f"\n  queries compared: {n}  t [{base['t'][0]:.1f} .. {base['t'][n-1]:.1f}]")
    print(f"\n  {'metric':<34}{'baseline':>14}{'penalised':>14}{'delta':>13}")
    print("  " + "-" * 73)
    rows = [
        ("queries at u_max (56.75)", int((base["u"][:n] > 56.7499).sum()),
         int((pen["u"][:n] > 56.7499).sum()), "count"),
        ("queries at u_min (56.69)", int((base["u"][:n] < 56.6901).sum()),
         int((pen["u"][:n] < 56.6901).sum()), "count"),
    ]
    for name, b, p, _ in rows:
        print(f"  {name:<34}{b:>14d}{p:>14d}{p-b:>+13d}")
    for name, b, p in [
        ("mean u [A]", base["u"][:n].mean(), pen["u"][:n].mean()),
        ("min u [A]", base["u"][:n].min(), pen["u"][:n].min()),
        ("mean coolant mass flow", mb.mean(), mp.mean()),
        ("min coolant mass flow", mb.min(), mp.min()),
        ("total coolant (sum)", mb.sum(), mp.sum()),
        ("control total variation [A]", np.abs(np.diff(base["u"][:n])).sum(),
         np.abs(np.diff(pen["u"][:n])).sum()),
    ]:
        print(f"  {name:<34}{b:>14.4f}{p:>14.4f}{p-b:>+13.4f}")

    if complete:
        dsm_b = (base["pct"][:n] - TARGET) / SPAN
        dsm_p = (pen["pct"][:n] - TARGET) / SPAN
        print()
        for name, b, p in [
            ("peak PCT [K]", base["pct"][:n].max(), pen["pct"][:n].max()),
            ("min PCT [K]", base["pct"][:n].min(), pen["pct"][:n].min()),
            ("max DSM", dsm_b.max(), dsm_p.max()),
            ("mean DSM", dsm_b.mean(), dsm_p.mean()),
            ("mean |PCT-target| [K]", np.abs(base["pct"][:n]-TARGET).mean(),
             np.abs(pen["pct"][:n]-TARGET).mean()),
            ("min safety margin [K]", (U_Y - base["pct"][:n]).min(),
             (U_Y - pen["pct"][:n]).min()),
        ]:
            print(f"  {name:<34}{b:>14.4f}{p:>14.4f}{p-b:>+13.4f}")

        np.savez(EVID / "figure1_penalty_effect.npz",
                 t=base["t"][:n], pct_base=base["pct"][:n], pct_pen=pen["pct"][:n],
                 u_base=base["u"][:n], u_pen=pen["u"][:n],
                 mdot_base=mb, mdot_pen=mp,
                 dsm_base=dsm_b, dsm_pen=dsm_p,
                 dsm_thr=DSM_THR, target=TARGET, U_y=U_Y)
        print(f"\n  saved {EVID / 'figure1_penalty_effect.npz'}")
    else:
        print("\n  (closed-loop PCT unavailable until the run finishes and writes"
              "\n   phase2c_fullloop_penalised*_LJ_seed1.mat)")

    print("\n  per-query control (first 24)")
    print(f"  {'t':>8}{'u_base':>11}{'u_pen':>11}{'du':>10}"
          f"{'mdot_base':>12}{'mdot_pen':>11}{'dmdot':>10}")
    for i in range(min(24, n)):
        print(f"  {base['t'][i]:8.1f}{base['u'][i]:11.5f}{pen['u'][i]:11.5f}"
              f"{pen['u'][i]-base['u'][i]:+10.5f}{mb[i]:12.4f}{mp[i]:11.4f}"
              f"{mp[i]-mb[i]:+10.4f}")


def figure2():
    print("\n" + "=" * 74)
    print("FIGURE 2 -- solver comparison (archived KPIs, no new runs)")
    print("=" * 74)

    f = EVID / "phase2c_archived_trajectory_kpis.csv"
    if not f.exists():
        print(f"  missing: {f}")
        return
    with f.open() as fh:
        rows = list(csv.DictReader(fh))

    keys = ["min_safety_margin_K", "peak_PCT_K", "settling_time_s", "IAE_Ks",
            "steady_state_error_K", "total_variation_A", "max_control_rate_A_s"]
    print(f"\n  {'metric':<26}" + "".join(f"{r['solver']:>16}" for r in rows))
    print("  " + "-" * (26 + 16*len(rows)))
    for k in keys:
        line = f"  {k:<26}"
        for r in rows:
            v = float(r[k])
            line += f"{v:>16.6f}" if abs(v) < 1e6 else f"{'inf':>16}"
        print(line)

    lj = next(r for r in rows if r["solver"] == "LJ")
    ga = next(r for r in rows if r["solver"] == "GA")
    print("\n  Note: LJ and GA agree to the 15th digit on min_safety_margin "
          f"({lj['min_safety_margin_K']})")
    print(f"        and peak_PCT ({lj['peak_PCT_K']}), because with a single-node")
    print("        DET both reconstruct the same root trajectory. They differ in")
    print(f"        control smoothness: total variation {float(lj['total_variation_A']):.3f} (LJ)"
          f" vs {float(ga['total_variation_A']):.3f} (GA).")
    pid = next(r for r in rows if r["solver"] == "PID")
    print(f"        PID pays {float(pid['total_variation_A']):.2f} A of total variation and "
          f"{float(pid['max_control_rate_A_s']):.1f} A/s peak rate --")
    print("        far harsher actuation, and it cannot carry a J_psm term at all.")


if __name__ == "__main__":
    figure1()
    figure2()
