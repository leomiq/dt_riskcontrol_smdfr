"""Extract the archived per-instant control results for the Chapter 7 series.

Reads the 61 MAT workspaces of each archived series (LJ test 11_10, GA test
11_10, PID root series in Version5) and produces:
  * phase2c_archived_trajectory_per_instant.csv  (t_q, u_applied, PCT, f_opt,
    optimizer runtime, passes/generations)
  * phase2c_archived_trajectory_kpis.csv          (closed-loop KPIs computed
    over the 61-instant grid, replicating calculate_performance_metrics.m)
"""
import csv
import os
from pathlib import Path

import numpy as np
from scipy.io import loadmat

# The archived per-solver sweep directories ("LJ test 11_10", "GA test 11_10")
# are raw run output and are not redistributed with this repository; point
# ARCHIVE_DIR at them to regenerate the archived-trajectory CSVs.
REPO = Path(__file__).resolve().parent.parent
BASE = Path(os.environ.get("ARCHIVE_DIR", REPO / "results" / "archive"))
OUT = Path(os.environ.get("P2C_OUT_DIR", REPO / "results"))

SERIES = {
    "LJ": (BASE / "LJ test 11_10", "LJ_test_sampled_vel_ACS_2_loop_{}.mat", 0.3),
    "GA": (BASE / "GA test 11_10", "GA_test_sampled_vel_ACS_2_loop_{}.mat", 0.3),
    "PID": (BASE, "PID_test_sampled_vel_ACS_2_loop_{}.mat", 0.2),
}
PCT_MAX = 1244.0


def kpis(time, pct, control, pct_ref, tolerance):
    time = np.asarray(time, float)
    pct = np.asarray(pct, float)
    control = np.asarray(control, float)
    error = pct_ref - pct
    m = {}
    m["min_safety_margin"] = float((PCT_MAX - pct).min())
    m["safety_violation"] = bool(m["min_safety_margin"] < 0)
    m["min_overshoot"] = float((pct - pct_ref).min())
    m["percent_overshoot"] = m["min_overshoot"] / pct_ref * 100.0
    settled = np.where(np.abs(error) <= tolerance)[0]
    outside = np.where(np.abs(error) > tolerance)[0]
    if settled.size == 0:
        m["settling_time"] = float("inf")
    elif outside.size == 0:
        m["settling_time"] = float(time[settled[0]])
    else:
        last = outside[-1]
        m["settling_time"] = float(time[last + 1]) if last + 1 < len(time) else float("inf")
    m["IAE"] = float(np.trapz(np.abs(error), time))
    m["cumulative_error"] = float(np.trapz(error, time))
    m["steady_state_error"] = float(error[int(round(0.8 * len(error))):].mean())
    m["peak_PCT"] = float(pct.max())
    d = np.diff(control)
    m["total_variation"] = float(np.abs(d).sum())
    m["max_control_rate"] = float(np.abs(d).max() / np.mean(np.diff(time)))
    return m


def main():
    per_instant = []
    kpi_rows = []
    for name, (folder, fmt, tol) in SERIES.items():
        times = [1160.0 + 0.5 * i for i in range(61)]
        last = loadmat(folder / fmt.format("1190"), squeeze_me=True, struct_as_record=False)
        t_snap = np.asarray(last["time_snapshots"], float).ravel()
        u_traj = np.asarray(last["u_V2P_trajectory"], float).ravel()
        tw = np.asarray(last["tw_max_V2P_trajectory"], float)
        ctrl = np.asarray(last["control_series"], float)
        target = float(np.asarray(last["target"]).ravel()[0])

        pct_at_q = np.array([tw[np.argmin(np.abs(tw[:, 0] - tq)), 1] for tq in t_snap])
        u_at_q = np.array([ctrl[np.argmin(np.abs(ctrl[:, 0] - tq)), 1] for tq in t_snap])

        for i, tq in enumerate(times):
            fn = fmt.format(f"{tq:g}")
            d = loadmat(folder / fn, variable_names=["f_opt", "t_vector", "u_opt"],
                        squeeze_me=True, struct_as_record=False)
            f_opt = float(np.asarray(d["f_opt"]).ravel()[0]) if "f_opt" in d else ""
            tv = np.asarray(d["t_vector"], float).ravel() if "t_vector" in d else np.array([])
            rt = float(tv.sum()) if tv.size else ""
            npg = int(np.count_nonzero(tv)) if tv.size else ""
            per_instant.append([name, f"{tq:g}", u_at_q[i], pct_at_q[i], f_opt, rt, npg])

        m = kpis(t_snap, pct_at_q, u_at_q, target, tol)
        kpi_rows.append([name, m["min_safety_margin"], m["safety_violation"],
                         m["min_overshoot"], m["percent_overshoot"], m["settling_time"],
                         m["IAE"], m["cumulative_error"], m["steady_state_error"],
                         m["peak_PCT"], m["total_variation"], m["max_control_rate"]])
        print(f"=== {name} (tolerance {tol}) ===")
        for k, v in m.items():
            print(f"  {k:20s} = {v}")

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "phase2c_archived_trajectory_per_instant.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["solver", "t_q_s", "u_applied_A", "PCT_K", "f_opt", "runtime_s", "passes_generations"])
        w.writerows(per_instant)
    with (OUT / "phase2c_archived_trajectory_kpis.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["solver", "min_safety_margin_K", "safety_violation", "min_overshoot_K",
                    "percent_overshoot", "settling_time_s", "IAE_Ks", "cumulative_error_Ks",
                    "steady_state_error_K", "peak_PCT_K", "total_variation_A", "max_control_rate_A_s"])
        w.writerows(kpi_rows)
    print(f"\nEscrito en {OUT}")


if __name__ == "__main__":
    main()
