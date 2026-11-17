"""F5b.1c summary (results/f5b1c_risk/runs.csv, julia/f5b1c_risk.jl): risk bound on N2 against
DET-accelerated ACS failures. Every case paired against N2+bias on the same realisation.

Columns: realisations reaching 1220 K (with a Wilson 95 % interval), PCT max (mean / worst),
paired Δ PCT max and Δ ISE, injected mass, mean deviation from y_ref (negative = pre-cooling),
mean U_α. Also the crossing rate split by the length of the longest failure in the window, and
P(condition) = 1 / mean(tries) to turn conditional rates into rates per transient.
Usage: .venv/Scripts/python.exe python/f5b1c_risk_summary.py
"""
import csv
import os
from collections import defaultdict
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "f5b1c_risk", "runs.csv"))))
by = defaultdict(dict)
for r in rows:
    by[r["case"]][r["seed"]] = r
ref = by["N2"]
seeds_all = sorted(set.intersection(*[set(v) for v in by.values()]))   # complete realisations only


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def ci(x):
    x = np.asarray(x, float)
    return x.mean(), (1.96 * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan)


tries = np.array([float(ref[s]["tries"]) for s in seeds_all])
print(f"{len(seeds_all)} complete realisations; P(>=1 failure starting in the window) ~ {1 / tries.mean():.3f}")
mf = np.array([float(ref[s]["max_fail_len"]) for s in seeds_all])
print(f"longest failure in window: median {np.median(mf):.1f} s, P(>130 s) {np.mean(mf > 130):.3f}\n")
order = ["N2", "N2_risk_a0.1", "N2_risk_a0.05", "N2_risk_a0.02", "N2_risk_a0.01", "PI"]
print(f"{'case':16s} {'n':>3s} {'#1220':>5s} {'P(cross|fail) [95% CI]':>24s} {'PCTmax':>8s} {'worst':>8s} "
      f"{'dPCTmax vs N2':>15s} {'dISE vs N2':>16s} {'mass':>6s} {'mean dev':>8s} {'U_a':>7s}")
for case in [c for c in order if c in by]:
    S = by[case]; ss = seeds_all
    f = lambda k: np.array([float(S[s][k]) for s in ss])
    k = sum(S[s]["limit"] == "true" for s in ss); lo, hi = wilson(k, len(ss))
    a, b = ci([float(S[s]["pct_max"]) - float(ref[s]["pct_max"]) for s in ss])
    c, d = ci([float(S[s]["ise"]) - float(ref[s]["ise"]) for s in ss])
    print(f"{case:16s} {len(ss):3d} {k:5d} {k / len(ss):8.3f} [{lo:.3f}, {hi:.3f}] {f('pct_max').mean():8.2f} "
          f"{f('pct_max').max():8.2f} {a:+8.3f}±{b:5.3f} {c:+9.0f}±{d:6.0f} {f('mass_plant').mean():6.0f} "
          f"{f('mean_dev').mean():8.2f} {np.nanmean(f('ubound_mean')):7.1f}")

print("\ncrossings by longest failure in the window:")
bins = [(0, 60), (60, 100), (100, 140), (140, 1e9)]
print("case             " + " ".join(f"{f'{a:.0f}-{b:.0f}s' if b < 1e8 else f'>{a:.0f}s':>12s}" for a, b in bins))
for case in [c for c in order if c in by]:
    S = by[case]; cells = []
    for a, b in bins:
        ss = [s for s in seeds_all if a <= float(ref[s]["max_fail_len"]) < b]
        k = sum(S[s]["limit"] == "true" for s in ss)
        cells.append(f"{k:>4d}/{len(ss):<4d}   ")
    print(f"{case:16s} " + " ".join(cells))
