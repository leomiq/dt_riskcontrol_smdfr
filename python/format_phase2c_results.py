"""Format the Phase 2C fixed-state repeated-seed results into the thesis
summary table and statistics block (P-M4, PAT-C-065, PAT-C-066).

Usage:
    python format_phase2c_results.py <path_to_phase2c_fixed_state_results.mat>

The script prints:
  * per-solver mean / SD / median / IQR / range for f_opt, u_apply, minDSM_K,
    runtime_s, passes/generations and evaluations;
  * a ready-to-paste LaTeX paragraph for the Chapter 7 placeholder.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat


def stats(a):
    a = np.asarray(a, dtype=float).ravel()
    q1, q3 = np.percentile(a, [25, 75])
    return dict(
        mean=a.mean(),
        sd=a.std(ddof=1) if a.size > 1 else 0.0,
        median=np.median(a),
        iqr=q3 - q1,
        mn=a.min(),
        mx=a.max(),
    )


def fmt(s, dec=3):
    return (f"{s['mean']:.{dec}f} \u00b1 {s['sd']:.{dec}f} "
            f"(med {s['median']:.{dec}f}, IQR {s['iqr']:.{dec}f}, "
            f"range [{s['mn']:.{dec}f}, {s['mx']:.{dec}f}])")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    f = Path(sys.argv[1])
    d = loadmat(f, squeeze_me=True, struct_as_record=False)
    lj = d["LJ"]
    ga = d["GA"]
    cfg = d["config"]

    rows = []
    for name, runs in (("LJ", lj), ("GA", ga)):
        n = len(runs)
        print(f"=== {name} (n={n}) ===")
        for key, label, dec in (
            ("f_opt", "f_opt", 4),
            ("u_apply", "u_apply (A)", 6),
            ("minDSM_K", "minDSM (K)", 3),
            ("runtime_s", "runtime (s)", 2),
        ):
            s = stats([getattr(r, key) for r in runs])
            print(f"  {label:14s}: {fmt(s, dec)}")
        key = "passes" if name == "LJ" else "generations"
        s = stats([getattr(r, key) for r in runs])
        print(f"  {key:14s}: {fmt(s, 1)}")
        s = stats([getattr(r, "evaluations") for r in runs])
        print(f"  evaluations : {fmt(s, 1)}")
        rows.append((name, n))

    print("\n--- Config ---")
    print(f"query_time = {cfg.query_time_s} s; seeds = {cfg.seeds}; "
          f"LJ(max_passes={cfg.LJ.max_passes}, M={cfg.LJ.M}, R={cfg.LJ.R}, "
          f"r_init={cfg.LJ.r_init}, gamma={cfg.LJ.gamma}); "
          f"GA(generations={cfg.GA.generations}, pop_size={cfg.GA.pop_size}, "
          f"crossover={cfg.GA.crossover_frac}, mutation={cfg.GA.mutation_rate}); "
          f"u in [{cfg.u_min}, {cfg.u_max}] A; target={cfg.target} K")
    print(f"Y_denorm_50 mapping: {cfg.Y_denorm_50_mapping}")

    print("\n--- LaTeX paragraph draft ---")
    ljf = stats([getattr(r, "f_opt") for r in lj])
    ljd = stats([getattr(r, "minDSM_K") for r in lj])
    ljr = stats([getattr(r, "runtime_s") for r in lj])
    gaf = stats([getattr(r, "f_opt") for r in ga])
    gad = stats([getattr(r, "minDSM_K") for r in ga])
    gar = stats([getattr(r, "runtime_s") for r in ga])
    print(
        "At the fixed query state tq = 1180 s with 20 recorded seeds, LJ "
        f"achieved f* = {ljf['median']:.4f} (median; IQR {ljf['iqr']:.4f}), "
        f"first action u* = {np.median([getattr(r,'u_apply') for r in lj]):.4f} A, "
        f"minimum predicted DSM = {ljd['median']:.2f} K (IQR {ljd['iqr']:.2f} K) "
        f"and optimizer runtime = {ljr['median']:.1f} s (IQR {ljr['iqr']:.1f} s). "
        "GA under an equivalent evaluation budget (420 evaluations) gave "
        f"f* = {gaf['median']:.4f} (IQR {gaf['iqr']:.4f}), "
        f"minimum predicted DSM = {gad['median']:.2f} K (IQR {gad['iqr']:.2f} K) "
        f"and runtime = {gar['median']:.1f} s (IQR {gar['iqr']:.1f} s)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
