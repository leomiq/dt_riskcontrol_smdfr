"""F5b.1b summary (results/f5b1b_dsm/runs.csv, julia/f5b1b_dsm.jl): every case paired against
N2+bias on the same realisation (common random numbers).
Columns: PCT max (mean, worst), paired ΔPCTmax and ΔISE (mean ± 95 % CI), injected mass,
mean deviation from y_ref (negative = pre-cooling) and time above 1200 K.
Usage: .venv/Scripts/python.exe python/f5b1b_dsm_summary.py [twin=commanded]
"""
import csv
import os
import sys
from collections import defaultdict
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
TWIN = sys.argv[1] if len(sys.argv) > 1 else "commanded"
rows = [r for r in csv.DictReader(open(os.path.join(ROOT, "results", "f5b1b_dsm", "runs.csv"))) if r["twin"] == TWIN]
by = defaultdict(dict)
for r in rows:
    by[r["case"]][r["seed"]] = r
ref = by["N2"]


def ci(x):
    x = np.asarray(x, float)
    return x.mean(), (1.96 * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan)


print(f"{'case':30s} {'n':>3s} {'#1220':>5s} {'PCTmax':>8s} {'worst':>8s} {'dPCTmax vs N2':>15s} "
      f"{'dISE vs N2 %':>14s} {'mass':>6s} {'mean dev':>8s} {'t>1200':>7s}")
for case in sorted(by, key=lambda c: (not c.startswith("N2"), c)):
    S = by[case]; seeds = sorted(S)
    f = lambda k: np.array([float(S[s][k]) for s in seeds])
    com = [s for s in seeds if s in ref]
    a, b = ci([float(S[s]["pct_max"]) - float(ref[s]["pct_max"]) for s in com])
    c, d = ci([100 * (float(S[s]["ise"]) / float(ref[s]["ise"]) - 1) for s in com])
    print(f"{case:30s} {len(seeds):3d} {sum(S[s]['limit'] == 'true' for s in seeds):5d} {f('pct_max').mean():8.2f} "
          f"{f('pct_max').max():8.2f} {a:+8.3f}±{b:5.3f} {c:+7.2f}±{d:5.2f} {f('mass_plant').mean():6.0f} "
          f"{f('mean_dev').mean():8.2f} {f('t_above_1200').mean():7.1f}")
