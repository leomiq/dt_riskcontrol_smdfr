"""F5b.0 closed-loop summary (results/f5b0_closed_loop/runs.csv, julia/f5b0_closed_loop.jl).

Paired differences (common random numbers) per case:
  block A: measured vs commanded ACS flow in the warm-up (reference A_commanded);
  block B: each forecaster against the oracle with the same noise level and jumps.
Usage: .venv/Scripts/python.exe python/f5b0_closed_loop_summary.py
"""
import csv
import os
from collections import defaultdict
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "f5b0_closed_loop", "runs.csv"))))
by = defaultdict(dict)
for r in rows:
    by[r["case"]][r["seed"]] = r


def ref_of(case):
    if case.startswith("A_") or case.startswith("PI"):
        return "A_commanded" if case.startswith("A_") else None
    parts = case.split("_")                      # B_<forecast>_s<σ>[_jumps]
    tail = case[case.index("_s"):]
    return "B_oracle" + tail


def ci(x):
    x = np.asarray(x)
    return (x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else (np.nan, np.nan)


print(f"{'case':28s} {'n':>3s} {'#1220':>5s} {'PCTmax':>8s} {'ISE':>9s} {'mass':>7s} {'ferr p99':>9s} "
      f"{'ΔPCTmax vs ref':>16s} {'ΔISE vs ref (%)':>18s}")
for case in sorted(by, key=lambda c: (c[0], c)):
    S = by[case]; seeds = sorted(S)
    f = lambda k: np.array([float(S[s][k]) for s in seeds])
    lim = sum(S[s]["limit"] == "true" for s in seeds)
    ref = ref_of(case); dp = di = ""
    if ref and ref in by and ref != case:
        R = by[ref]; com = [s for s in seeds if s in R]
        a, b = ci([float(S[s]["pct_max"]) - float(R[s]["pct_max"]) for s in com])
        c, d = ci([100 * (float(S[s]["ise"]) / float(R[s]["ise"]) - 1) for s in com])
        dp = f"{a:+7.2f}±{b:.2f}"; di = f"{c:+7.1f}±{d:.1f}"
    fe = f("ferr_p99")
    print(f"{case:28s} {len(seeds):3d} {lim:5d} {f('pct_max').mean():8.2f} {f('ise').mean():9.0f} "
          f"{f('mass_plant').mean():7.0f} {np.nanmean(fe):9.3f} {dp:>16s} {di:>18s}")
