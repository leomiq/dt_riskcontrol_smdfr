"""Figure 2 -- solver comparison (LJ / GA / PID) from archived data only.

Deliberately kept separate from the DSM-penalty figures: this axis answers
"which optimizer solves the same problem better", while the penalty figures
answer "what does the safety constraint do to the trajectory". Mixing them
would confound the effect of the search algorithm with the effect of the
constraint.

  figures/solver_kpi_bars.{pdf,png}     KPI bars, grouped by what they measure
  figures/solver_actuation.{pdf,png}    control effort: TV and peak rate (log)
  figures/solver_comparison.{pdf,png}   composed 2x2 overview

Two facts the figures must make legible:
  * LJ and GA agree to the 15th digit on min_safety_margin and peak_PCT --
    with a single-node DET both reconstruct the same root trajectory, so the
    optimizer choice cannot change the thermal outcome here.
  * PID is the non-optimised reference. It tracks slightly better in IAE but
    pays two orders of magnitude more actuation, and it cannot carry a J_psm
    term at all (a Kp/Ki/Kd law does not minimise a cost).
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 10.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "lines.linewidth": 1.5,
    "figure.dpi": 150,
})

HERE = Path(r"external revision/Response to reviewers/evidence/phase2c_fixed_state_experiment")
OUT = HERE / "results_solver_comparison"     # archived LJ / GA / PID study
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

COLORS = {"LJ": "#1f77b4", "GA": "#d62728", "PID": "#2ca02c"}
ORDER = ["LJ", "GA", "PID"]

rows = list(csv.DictReader((HERE / "phase2c_archived_trajectory_kpis.csv").open(encoding="utf-8")))
K = {r["solver"]: r for r in rows}


def val(solver, key):
    v = float(K[solver][key])
    return v


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)
    print("OK", name)


def _bars(ax, key, ylabel, fmt="{:.4f}", ylim=None, logy=False, note=None):
    vals = [val(s, key) for s in ORDER]
    finite = [v for v in vals if np.isfinite(v)]
    x = np.arange(len(ORDER))
    bars = ax.bar(x, [v if np.isfinite(v) else 0 for v in vals],
                  color=[COLORS[s] for s in ORDER], alpha=0.75, width=0.6,
                  edgecolor="0.25", linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(ORDER)
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y", alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    if logy:
        ax.set_yscale("log")
    if ylim:
        ax.set_ylim(*ylim)
    hi = max(finite) if finite else 1.0
    for xi, v, b in zip(x, vals, bars):
        if not np.isfinite(v):
            ax.text(xi, 0.02 * hi, "inf", ha="center", va="bottom",
                    fontsize=9, color=COLORS[ORDER[xi]], fontweight="bold")
            continue
        off = 1.25 if logy else 1.02
        ax.annotate(fmt.format(v), xy=(xi, v),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8.5, color="0.2")
    if note:
        ax.set_title(note, fontsize=9, color="0.35")


# ---- 1. thermal-safety KPIs: identical for LJ and GA ----
fig, axs = plt.subplots(1, 3, figsize=(9.8, 3.5))
_bars(axs[0], "min_safety_margin_K", "Min safety margin (K)", "{:.4f}",
      ylim=(62.0, 63.2), note="thermal safety")
_bars(axs[1], "peak_PCT_K", "Peak PCT (K)", "{:.4f}",
      ylim=(1180.8, 1181.8), note="peak temperature")
_bars(axs[2], "IAE_Ks", r"IAE (K$\cdot$s)", "{:.2f}",
      ylim=(0, 62), note="tracking error")
axs[0].annotate("LJ = GA to the\n15th digit", xy=(0.5, 62.786805), xytext=(1.6, 63.02),
                ha="center", fontsize=8, color="0.3",
                arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8,
                                shrinkA=2, shrinkB=2))
fig.tight_layout()
save(fig, "solver_kpi_bars")

# ---- 2. actuation effort ----
fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.6))
_bars(axs[0], "total_variation_A", "Control total variation (A)", "{:.3f}",
      logy=True, note="cumulative actuation")
_bars(axs[1], "max_control_rate_A_s", r"Peak control rate (A s$^{-1}$)", "{:.3f}",
      logy=True, note="peak actuation rate")
for ax in axs:
    ax.set_ylim(0.02, 3e3)
fig.tight_layout()
save(fig, "solver_actuation")

# ---- composed 2x2 ----
plt.rcParams.update({"font.size": 10.5, "axes.labelsize": 11})
fig, axes = plt.subplots(2, 2, figsize=(7.6, 6.6))
_bars(axes[0, 0], "min_safety_margin_K", "Min safety margin (K)", "{:.4f}",
      ylim=(62.0, 63.2), note="thermal safety")
axes[0, 0].annotate("LJ = GA to the 15th digit",
                    xy=(0.5, 62.786805), xytext=(1.55, 63.05),
                    ha="center", fontsize=7.5, color="0.3",
                    arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8,
                                    shrinkA=2, shrinkB=2))
_bars(axes[0, 1], "IAE_Ks", r"IAE (K$\cdot$s)", "{:.2f}", ylim=(0, 62),
      note="tracking error")
_bars(axes[1, 0], "total_variation_A", "Total variation (A)", "{:.3f}",
      logy=True, note="cumulative actuation")
axes[1, 0].set_ylim(0.2, 3e3)
_bars(axes[1, 1], "max_control_rate_A_s", r"Peak rate (A s$^{-1}$)", "{:.3f}",
      logy=True, note="peak actuation rate")
axes[1, 1].set_ylim(0.02, 3e3)
for ax, lab in zip(axes.ravel(), ["(a)", "(b)", "(c)", "(d)"]):
    ax.text(0.02, 0.97, lab, transform=ax.transAxes, va="top", ha="left",
            fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "solver_comparison.pdf")
fig.savefig(FIG / "solver_comparison.png")
plt.close(fig)
print("OK solver_comparison")

# ---- table for the text ----
keys = [
    ("min_safety_margin_K", "Min safety margin (K)", "{:.6f}"),
    ("peak_PCT_K", "Peak PCT (K)", "{:.6f}"),
    ("percent_overshoot", "Overshoot (%)", "{:+.4f}"),
    ("settling_time_s", "Settling time (s)", "{:.1f}"),
    ("IAE_Ks", "IAE (K.s)", "{:.4f}"),
    ("steady_state_error_K", "Steady-state error (K)", "{:+.6f}"),
    ("total_variation_A", "Control total variation (A)", "{:.6f}"),
    ("max_control_rate_A_s", "Peak control rate (A/s)", "{:.4f}"),
]
print("\n" + "=" * 78)
print("SOLVER COMPARISON (archived realization)")
print("=" * 78)
print(f"{'metric':<30}" + "".join(f"{s:>16}" for s in ORDER))
print("-" * 78)
for key, label, fmt in keys:
    line = f"{label:<30}"
    for s in ORDER:
        v = val(s, key)
        line += f"{fmt.format(v):>16}" if np.isfinite(v) else f"{'inf':>16}"
    print(line)
print("-" * 78)
print(f"LJ/GA share min_safety_margin and peak_PCT exactly "
      f"({val('LJ','min_safety_margin_K'):.15f} K)")
print(f"PID actuation vs LJ: total variation x{val('PID','total_variation_A')/val('LJ','total_variation_A'):.1f}, "
      f"peak rate x{val('PID','max_control_rate_A_s')/val('LJ','max_control_rate_A_s'):.0f}")
print(f"\nfigures in {FIG}")
