"""F5b.1b0 summary (results/f5b1b0_dsm_lambda/runs.csv, julia/f5b1b0_dsm_lambda.jl): thesis
indices of N2+DSM over λ_dsm against the thesis PID and the DET PI.

Per set / actuator / horizon: mean over realisations of IAE, reach and settling time (0.30 K
band; Inf counted separately), steady-state error, minimum overshoot, peak PCT, margin to
1220 K, control total variation [kg/s], coolant mass, realisations reaching 1220 K; and, for
N2+DSM, the largest paired change of the trajectory vs N2 (max over realisations of the
max |ΔPCT| is not logged: the paired Δ IAE and Δ TV are reported instead).
Usage: .venv/Scripts/python.exe python/f5b1b0_dsm_lambda_summary.py
"""
import csv
import os
from collections import defaultdict
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "f5b1b0_dsm_lambda", "runs.csv"))))
f = lambda r, k: float(r[k])
groups = defaultdict(list)
for r in rows:
    groups[(r["set"], r["act"], r["n_req"])].append(r)

ORDER = ["N2"] + [f"N2_DSM_U{U}_l{l}" for U in (1244, 1220) for l in (100, 1000, 10000, 100000, 1000000)] + \
        ["PID", "PID_free", "PI", "PI_free"]
for g in sorted(groups):
    G = groups[g]
    by = defaultdict(dict)
    for r in G:
        by[r["case"]][r["seed"]] = r
    ref = by.get("N2", {})
    print(f"\n=== set {g[0]}, actuator {g[1]}, n = {g[2]} ({len(ref)} realisations) ===")
    print(f"{'case':24s} {'#1220':>5s} {'IAE':>8s} {'dIAE vs N2':>11s} {'reach':>7s} {'settle':>7s} {'#nosettle':>9s} "
          f"{'ss err':>7s} {'min ovs':>8s} {'peak':>8s} {'marg1220':>8s} {'TVplant':>8s} {'dTV':>7s} {'mass':>7s}")
    for case in [c for c in ORDER if c in by]:
        S = by[case]; ss = sorted(S)
        v = lambda k: np.array([f(S[s], k) for s in ss])
        fin = lambda x: x[np.isfinite(x)]
        reach, settle = v("reach_time"), v("settling_time")
        dI = np.mean([f(S[s], "IAE") - f(ref[s], "IAE") for s in ss if s in ref]) if ref else np.nan
        dT = np.mean([f(S[s], "tv_plant") - f(ref[s], "tv_plant") for s in ss if s in ref]) if ref else np.nan
        print(f"{case:24s} {sum(S[s]['limit'] == 'true' for s in ss):5d} {v('IAE').mean():8.1f} {dI:+11.2f} "
              f"{(fin(reach).mean() if fin(reach).size else np.inf):7.1f} {(fin(settle).mean() if fin(settle).size else np.inf):7.1f} "
              f"{np.sum(~np.isfinite(settle)):9d} {v('steady_state_error').mean():7.2f} {v('min_overshoot').mean():8.2f} "
              f"{v('peak_PCT').mean():8.2f} {v('margin_1220').mean():8.2f} {v('tv_plant').mean():8.1f} {dT:+7.2f} "
              f"{v('coolant_mass').mean():7.0f}")
