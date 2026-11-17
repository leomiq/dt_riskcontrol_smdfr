"""F5b.1c-RB summary (results/f5b1c_rb_risk/runs.csv, julia/f5b1c_rb_risk.jl): Rao–Blackwell
probability that an ACS failure (DET-accelerated) takes the PCT above 1220 K, per controller,
with paired comparisons against N2+bias (same realisation) and 95 % intervals.

P(cross | failure in the window) and P(cross) per 240 s transient (accelerated rates); with the
real DET rates every failure lasts days, so P(cross) ≈ P(failure) ≈ 3.8e-6 for any controller.
Usage: .venv/Scripts/python.exe python/f5b1c_rb_summary.py
"""
import csv
import os
from collections import defaultdict
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "f5b1c_rb_risk", "runs.csv"))))
by = defaultdict(dict)
for r in rows:
    by[r["case"]][r["seed"]] = r
ref = by["N2"]
seeds = sorted(set.intersection(*[set(v) for v in by.values()]))
f = lambda r, k: float(r[k])


def ci(x):
    x = np.asarray(x, float)
    return x.mean(), 1.96 * x.std(ddof=1) / np.sqrt(len(x))


order = ["N2", "N2_risk_a0.1", "N2_risk_a0.05", "N2_risk_a0.02", "N2_risk_a0.01", "PID", "PI"]
print(f"{len(seeds)} realisations with every case\n")
print(f"{'case':16s} {'P(cross|fail)':>16s} {'ratio vs N2 (paired)':>22s} {'P(cross)/240s':>14s} "
      f"{'tau mean':>8s} {'tau min':>8s} {'ISE/ISE(N2)':>11s} {'mass':>6s} {'mean dev':>8s}")
out = []
for case in [c for c in order if c in by]:
    S = by[case]
    pc = np.array([f(S[s], "p_cross_cond") for s in seeds]); p0 = np.array([f(ref[s], "p_cross_cond") for s in seeds])
    m, h = ci(pc); lr = np.log(pc / p0); rm, rh = ci(lr)
    pu = np.mean([f(S[s], "p_cross") for s in seeds])
    ise = np.mean([f(S[s], "ise") for s in seeds]) / np.mean([f(ref[s], "ise") for s in seeds])
    print(f"{case:16s} {m:8.4f} ± {h:.4f} {np.exp(rm):8.3f} [{np.exp(rm - rh):.3f}, {np.exp(rm + rh):.3f}] "
          f"{pu:14.5f} {np.mean([f(S[s], 'tau_mean') for s in seeds]):8.1f} {np.mean([f(S[s], 'tau_min') for s in seeds]):8.1f} "
          f"{ise:11.2f} {np.mean([f(S[s], 'mass') for s in seeds]):6.0f} {np.mean([f(S[s], 'mean_dev') for s in seeds]):8.2f}")
    out.append(dict(case=case, p_cond=m, p_cond_ci=h, ratio=np.exp(rm), ratio_lo=np.exp(rm - rh), ratio_hi=np.exp(rm + rh),
                    p_uncond=pu, ise_ratio=ise, mass=np.mean([f(S[s], "mass") for s in seeds]),
                    mean_dev=np.mean([f(S[s], "mean_dev") for s in seeds])))
with open(os.path.join(ROOT, "results", "f5b1c_rb_risk", "summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader()
    for o in out:
        w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in o.items()})
