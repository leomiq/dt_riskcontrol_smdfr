"""F5b.1a: regime map (m_max, ramp) -> {always saturated, free decision with risk, no risk}
from results/f5b1a_capacity/runs.csv (julia/f5b1a_capacity_sweep.jl).

Per configuration and controller: realisations reaching the DET end state (PCT > 1220 K),
mean / max PCT, minimum margin to 1220 K, mean fractions of first actions on each bound, and
the injected mass. Paired difference N2 - PI (same capacity) in PCT max.

Regime rule (per configuration, on N2+bias):
  saturated  : free fraction < 0.10 (the first action is on a bound almost always)
  no risk    : no realisation reaches 1220 K and the minimum margin > 5 K
  free+risk  : otherwise (interior decisions AND realisations at or near the limit)
Usage: .venv/Scripts/python.exe python/f5b1a_capacity_map.py
"""
import csv
import os
from collections import defaultdict
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
F = os.path.join(ROOT, "results", "f5b1a_capacity", "runs.csv")
rows = list(csv.DictReader(open(F)))
num = lambda r, k: float(r[k])
by = defaultdict(list)
for r in rows:
    by[(num(r, "m_max"), num(r, "ramp"), r["ctrl"])].append(r)

out = []
print(f"{'m_max':>5} {'ramp':>4} {'ctrl':4} {'n':>3} {'#1220':>5} {'PCTmax':>7} {'worst':>7} {'minmarg':>7} "
      f"{'hi':>5} {'up':>5} {'dn':>5} {'lo':>5} {'free':>5} {'mass':>7} {'dPCT(N2-PI)':>12} regime")
keys = sorted(by, key=lambda k: (k[0], k[1], k[2]))
for k in keys:
    S = by[k]
    pm = np.array([num(r, "pct_max") for r in S])
    lim = sum(r["limit"] == "true" for r in S)
    sat = {c: np.mean([num(r, f"sat_{c}") for r in S]) for c in ("hi", "up", "dn", "lo")}
    free = np.mean([num(r, "free") for r in S])
    mass = np.mean([num(r, "mass_plant") for r in S])
    d = ""
    regime = ""
    if k[2] == "N2":
        P = {r["seed"]: num(r, "pct_max") for r in by.get((k[0], k[1], "PI"), [])}
        dd = [num(r, "pct_max") - P[r["seed"]] for r in S if r["seed"] in P]
        if dd:
            d = f"{np.mean(dd):+6.2f}±{1.96 * np.std(dd, ddof=1) / np.sqrt(len(dd)):.2f}"
        regime = ("saturated" if free < 0.10 else
                  "no risk" if lim == 0 and (1220 - pm.max()) > 5 else "free+risk")
    print(f"{k[0]:5.0f} {k[1]:4.1f} {k[2]:4s} {len(S):3d} {lim:5d} {pm.mean():7.2f} {pm.max():7.2f} {1220 - pm.max():7.2f} "
          f"{sat['hi']:5.2f} {sat['up']:5.2f} {sat['dn']:5.2f} {sat['lo']:5.2f} {free:5.2f} {mass:7.0f} {d:>12s} {regime}")
    out.append(dict(m_max=k[0], ramp=k[1], ctrl=k[2], n=len(S), n_limit=lim, pct_max_mean=pm.mean(),
                    pct_max_worst=pm.max(), margin_min=1220 - pm.max(), **{f"sat_{c}": v for c, v in sat.items()},
                    free=free, mass=mass, dpct_vs_pi=d, regime=regime))
with open(os.path.join(ROOT, "results", "f5b1a_capacity", "map.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader()
    for o in out:
        w.writerow({a: (f"{b:.6g}" if isinstance(b, float) else b) for a, b in o.items()})
