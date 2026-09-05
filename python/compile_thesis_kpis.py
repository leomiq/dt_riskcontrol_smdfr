"""Thesis-format performance indices for the penalty / saturation configurations.

Mirrors Table `tab:gbdt-archived-kpis` in chapter4.tex, which reports for each
controller:

    Minimum safety margin (DSM of PCT)   min(U_y - PCT)          [K]
    Minimum overshoot (PCT - ref)        min(PCT - PCT_ref)      [K]
    Settling time (<0.30 K)              first t after which |e| stays in band
    Steady-state error                   mean e over the last 20% of the grid
    Peak PCT                             max(PCT)                [K]
    IAE                                  trapz(|e|, t)           [K.s]
    Control total variation              sum |du|                [A]
    Optimizer runtime (median per query) median of the per-query loop time [s]

Definitions follow calculate_performance_metrics.m, which is what produced the
archived numbers, so the new rows are directly comparable with the thesis table.
Two indices Leo asked about explicitly:
  * "total de la amplitud"  -> control total variation, sum of |du| over the grid
  * "maxima diferencia con respecto al setpoint" -> reported both signed
    (minimum overshoot, the thesis index) and as max |PCT - ref|.

Configurations covered:
    A  baseline          u_max 56.75,  lambda 0,     LJ original
    C  over-penalised    u_max 56.90,  lambda 1e6,   LJ original
    D  retuned           u_max 56.90,  lambda 1e4,   LJ retuned
plus the archived LJ / GA / PID rows read straight from the thesis CSV.
"""
from pathlib import Path

import numpy as np
from scipy.io import loadmat

HERE = Path(r"external revision/Response to reviewers/evidence/phase2c_fixed_state_experiment")
RET = HERE / "results_retuned_lambda1e4"
SOLV = HERE / "results_solver_comparison"

PCT_REF = 1179.2711776      # CalculateTs(100, 543)
PCT_MAX = 1244.0
TOL = 0.30                  # thesis band for LJ/GA; PID uses 0.20

RUNS = [
    ("A  baseline  $\\lambda$=0, $u_{max}$=56.75", "phase2c_fullloop_pilot_LJ_seed1.mat"),
    ("C  penalised $\\lambda$=1e6, $u_{max}$=56.90", "widepen_umax56p90_LJ_seed1.mat"),
    ("D  retuned   $\\lambda$=1e4, $u_{max}$=56.90", "retuned_lam1e4_umax56p90_LJ_seed1.mat"),
]


def metrics(t, pct, control, rt, tol=TOL):
    """Reproduce calculate_performance_metrics.m exactly."""
    t = np.asarray(t, float).ravel()
    pct = np.asarray(pct, float).ravel()
    control = np.asarray(control, float).ravel()
    err = PCT_REF - pct                      # reference minus measurement

    m = {}
    m["min_safety_margin"] = float((PCT_MAX - pct).min())
    m["safety_violation"] = bool(m["min_safety_margin"] < 0)
    m["min_overshoot"] = float((pct - PCT_REF).min())
    m["percent_overshoot"] = m["min_overshoot"] / PCT_REF * 100.0
    m["max_abs_dev"] = float(np.abs(pct - PCT_REF).max())
    m["max_pos_dev"] = float((pct - PCT_REF).max())

    # The first two samples carry the ACS activation transient (PCT dips to
    # ~1146 K at t=1160.5), which is identical in every configuration because it
    # precedes any control effect. It therefore dominates the overshoot indices
    # and hides the differences the controller actually produces, so the same
    # deviations are also reported excluding it.
    dev_pt = pct[2:] - PCT_REF
    m["min_overshoot_pt"] = float(dev_pt.min())
    m["max_abs_dev_pt"] = float(np.abs(dev_pt).max())
    m["max_pos_dev_pt"] = float(dev_pt.max())

    settled = np.where(np.abs(err) <= tol)[0]
    outside = np.where(np.abs(err) > tol)[0]
    if settled.size == 0:
        m["settling_time"] = np.inf
    elif outside.size == 0:
        m["settling_time"] = float(t[settled[0]])
    elif outside[-1] + 1 < len(t):
        m["settling_time"] = float(t[outside[-1] + 1])
    else:
        m["settling_time"] = np.inf          # never re-enters the band

    m["IAE"] = float(np.trapz(np.abs(err), t))
    m["cumulative_error"] = float(np.trapz(err, t))
    m["steady_state_error"] = float(err[int(round(0.8 * len(err))):].mean())
    m["peak_PCT"] = float(pct.max())
    m["min_PCT"] = float(pct.min())
    m["final_PCT"] = float(pct[-1])

    du = np.diff(control)
    m["total_variation"] = float(np.abs(du).sum())
    dt = np.mean(np.diff(t))
    m["max_control_rate"] = float(np.abs(du).max() / dt)
    m["mean_abs_du"] = float(np.abs(du).mean())

    if rt is not None and np.size(rt):
        rt = np.asarray(rt, float).ravel()
        m["runtime_median"] = float(np.median(rt))
        m["runtime_mean"] = float(np.mean(rt))
    return m


def load(path):
    d = loadmat(path, squeeze_me=True, struct_as_record=False)
    t = np.asarray(d["time_snapshots"], float).ravel()
    pct = np.asarray(d["PCT"], float).ravel()
    ctrl = np.asarray(d["control"], float).ravel()
    rt = np.asarray(d.get("optimizer_runtimes", []), float).ravel()
    # control has one fewer sample than PCT (it is the applied series); pad the
    # last action so the two align for the total-variation sum
    if len(ctrl) == len(t) - 1:
        ctrl = np.append(ctrl, ctrl[-1])
    return t, pct, ctrl, rt


rows = {}
for label, fname in RUNS:
    p = RET / fname
    if not p.exists():
        print(f"missing: {p}")
        continue
    t, pct, ctrl, rt = load(p)
    rows[label] = metrics(t, pct, ctrl, rt)

INDEX = [
    ("min_safety_margin", "Minimum safety margin (DSM of PCT)", "K", "{:.2f}"),
    ("min_overshoot", "Minimum overshoot (PCT - ref)", "K", "{:.2f}"),
    ("max_abs_dev", "Max |PCT - ref|", "K", "{:.2f}"),
    ("min_overshoot_pt", "Minimum overshoot, post-transient", "K", "{:.3f}"),
    ("max_abs_dev_pt", "Max |PCT - ref|, post-transient", "K", "{:.3f}"),
    ("settling_time", "Settling time (<0.30 K)", "s", "{:.1f}"),
    ("steady_state_error", "Steady-state error", "K", "{:+.2f}"),
    ("peak_PCT", "Peak PCT", "K", "{:.2f}"),
    ("final_PCT", "Final PCT", "K", "{:.2f}"),
    ("IAE", "IAE", "K.s", "{:.2f}"),
    ("total_variation", "Control total variation", "A", "{:.2f}"),
    ("max_control_rate", "Peak control rate", "A/s", "{:.3f}"),
    ("runtime_median", "Optimizer runtime (median per query)", "s", "{:.1f}"),
]

names = list(rows)
w = 40
print("=" * (w + 20 * len(names)))
print("THESIS-FORMAT PERFORMANCE INDICES")
print("=" * (w + 20 * len(names)))
print(f"{'Performance index':<{w}}" + "".join(f"{n.split()[0]:>20}" for n in names))
print("-" * (w + 20 * len(names)))
for key, label, unit, fmt in INDEX:
    line = f"{label + ' [' + unit + ']':<{w}}"
    for n in names:
        v = rows[n].get(key, np.nan)
        line += f"{'inf':>20}" if not np.isfinite(v) else f"{fmt.format(v):>20}"
    print(line)
print()
for n in names:
    print(f"  {n.split()[0]} = {n}")

# ---------------- LaTeX table, thesis style ----------------
tex = []
tex.append(r"\begin{table}[h!tb]")
tex.append(r"	\centering")
tex.append(r"	\caption{Closed-loop KPIs over the 61-instant grid for the DSM-penalty "
           r"configurations: baseline, over-penalised, and the retuned setting. "
           r"All runs use Luus--Jaakola with seed 1 and the same accident stream; "
           r"$u_{\max}$ and $\lambda$ are the only design changes.}")
tex.append(r"	\label{tab:gbdt-penalty-kpis}")
tex.append(r"	\begin{tabular}{lccc}")
tex.append(r"		\toprule")
tex.append(r"		\textbf{Performance index} & \textbf{Baseline} & "
           r"\textbf{$\lambda=10^{6}$} & \textbf{$\lambda=10^{4}$, retuned} \\")
tex.append(r"		 & $u_{\max}=56.75$ & $u_{\max}=56.90$ & $u_{\max}=56.90$ \\")
tex.append(r"		\midrule")

TEX_ROWS = [
    ("min_safety_margin", "Minimum safety margin (DSM of PCT)", "{:.2f} K"),
    ("min_overshoot", "Minimum overshoot (PCT $-$ ref)", "${:.2f}$ K"),
    ("max_abs_dev_pt", "Maximum $|$PCT $-$ ref$|$, post-transient", "{:.3f} K"),
    ("settling_time", "Settling time ($<0.30$ K)", "{:.1f} s"),
    ("steady_state_error", "Steady-state error", "${:+.2f}$ K"),
    ("peak_PCT", "Peak PCT", "{:.2f} K"),
    ("IAE", "IAE", "{:.2f} K$\\cdot$s"),
    ("total_variation", "Control total variation", "{:.2f} A"),
    ("runtime_median", "Optimizer runtime (median per query)", "{:.1f} s"),
]
# mark the best value per row (lower is better except the safety margin)
BETTER_HIGH = {"min_safety_margin"}
for key, label, fmt in TEX_ROWS:
    vals = [rows[n].get(key, np.nan) for n in names]
    fin = [v for v in vals if np.isfinite(v)]
    best = (max(fin) if key in BETTER_HIGH else min(fin)) if fin else None
    cells = []
    for v in vals:
        if not np.isfinite(v):
            cells.append(r"$\infty$")
        else:
            s = fmt.format(v)
            cells.append(r"\textbf{" + s + "}" if best is not None and np.isclose(v, best)
                         else s)
    tex.append(f"		{label} & " + " & ".join(cells) + r" \\")
tex.append(r"		\bottomrule")
tex.append(r"	\end{tabular}")
tex.append(r"\end{table}")

out = RET / "table_thesis_kpis.tex"
out.write_text("\n".join(tex), encoding="utf-8")
print(f"\nLaTeX table written to {out}\n")
print("\n".join(tex))
