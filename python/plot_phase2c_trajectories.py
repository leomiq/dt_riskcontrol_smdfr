"""Publication-quality Phase 2C comparison figures (vector PDF + PNG previews).

Per-indicator figures and a composed 2x3 overview (thesis):
  figures/pct_vs_tq.{pdf,png}      PCT(t_q)          LJ / GA / PID
  figures/dsm_vs_tq.{pdf,png}      DSM(t_q)          LJ / GA / PID
  figures/u_vs_tq.{pdf,png}        applied control   LJ / GA / PID
  figures/fopt_vs_tq.{pdf,png}     objective f*(t_q) LJ / GA
  figures/runtime_vs_tq.{pdf,png}  optimizer runtime LJ / GA
  figures/variance_1180.{pdf,png}  20-seed boxplots at t_q = 1180 s
  figures/phase2c_comparison.pdf   composed 2x3 overview (thesis)
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

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
MARK = {"LJ": "o", "GA": "s", "PID": "^"}
ME = 6
TARGET = 1179.2711776
XTICKS = list(range(1160, 1191, 10))
XLIM = (1158, 1192)

rows = list(csv.DictReader((HERE / "phase2c_archived_trajectory_per_instant.csv").open(encoding="utf-8")))
data = {s: {} for s in ["LJ", "GA", "PID"]}
for r in rows:
    d = data[r["solver"]]
    d.setdefault("tq", []).append(float(r["t_q_s"]))
    d.setdefault("u", []).append(float(r["u_applied_A"]))
    d.setdefault("pct", []).append(float(r["PCT_K"]))
    d.setdefault("f", []).append(float(r["f_opt"]) if r["f_opt"] else np.nan)
    d.setdefault("rt", []).append(float(r["runtime_s"]) if r["runtime_s"] else np.nan)

mat = loadmat(HERE / "phase2c_fixed_state_results.mat", squeeze_me=True, struct_as_record=False)
VAR = {
    "f_opt": {k: np.array([getattr(r, "f_opt") for r in mat[k]], float) for k in ["LJ", "GA"]},
    "u_apply": {k: np.array([getattr(r, "u_apply") for r in mat[k]], float) for k in ["LJ", "GA"]},
    "minDSM_K": {k: np.array([getattr(r, "minDSM_K") for r in mat[k]], float) for k in ["LJ", "GA"]},
    "runtime_s": {k: np.array([getattr(r, "runtime_s") for r in mat[k]], float) for k in ["LJ", "GA"]},
}

def _base_ax(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(*XLIM)
    ax.set_xticks(XTICKS)
    ax.grid(True, alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)

def _line(ax, solvers, key, ylog=False, ylim=None):
    for s in solvers:
        d = data[s]
        ax.plot(d["tq"], d[key], color=COLORS[s], marker=MARK[s], markersize=4,
                linewidth=1.5, markevery=ME, label=s)
    if ylog:
        ax.set_yscale("log")
    if ylim:
        ax.set_ylim(*ylim)

def _target(ax):
    ax.axhline(TARGET, color="0.45", ls="--", lw=1.1,
               label=r"setpoint $T_p^*=1179.27$ K")

def _plot_pid_zoomed(ax, ylo=56.55, yhi=56.85):
    """PID inside the LJ/GA-comparison scale; off-scale samples are marked at
    the axis edges with their values (30 A at t=1160; up to 57.26 A)."""
    d = data["PID"]
    tq = np.array(d["tq"], float)
    u = np.array(d["u"], float)
    vis = (u >= ylo) & (u <= yhi)
    ax.plot(tq[vis], u[vis], color=COLORS["PID"], marker=MARK["PID"], markersize=4,
            linewidth=1.5, markevery=ME, label="PID")
    below = tq[u < ylo]
    above = tq[u > yhi]
    if below.size:
        ax.scatter(below, np.full_like(below, ylo), marker="v", s=45, color=COLORS["PID"],
                   clip_on=False, zorder=6, edgecolor="white", linewidth=0.5)
        ax.annotate("PID 30 A", xy=(below[0], ylo), xytext=(below[0] + 1.2, ylo + 0.045),
                    arrowprops=dict(arrowstyle="->", color=COLORS["PID"], lw=1.0),
                    fontsize=7.5, color=COLORS["PID"])
    if above.size:
        ax.scatter(above, np.full_like(above, yhi), marker="^", s=45, color=COLORS["PID"],
                   clip_on=False, zorder=6, edgecolor="white", linewidth=0.5)
        ax.annotate("PID up to 57.26 A", xy=(above[-1], yhi), xytext=(above[-1] - 5.2, yhi - 0.055),
                    arrowprops=dict(arrowstyle="->", color=COLORS["PID"], lw=1.0),
                    fontsize=7.5, color=COLORS["PID"])

def _acs_band(ax):
    ax.axvspan(1160, 1161.5, color="0.85", alpha=0.5, lw=0)

def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)
    print("OK", name)

# ---- 1. PCT vs t_q ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_base_ax(ax, r"Query time $t_q$ (s)", "PCT (K)")
_acs_band(ax)
_line(ax, ["LJ", "GA", "PID"], "pct", ylim=(1135, 1185))
_target(ax)
ax.legend(frameon=False, loc="lower left", ncol=2)
save(fig, "pct_vs_tq")

# ---- 2. DSM vs t_q ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_base_ax(ax, r"Query time $t_q$ (s)", "DSM (K)")
_acs_band(ax)
for s in ["LJ", "GA", "PID"]:
    d = data[s]
    ax.plot(d["tq"], 1244 - np.array(d["pct"]), color=COLORS[s], marker=MARK[s],
            markersize=4, linewidth=1.5, markevery=ME, label=s)
ax.set_ylim(45, 110)
ax.annotate("ACS activation transient", xy=(1160.5, 97.5), xytext=(1164, 100),
            arrowprops=dict(arrowstyle="->", color="0.3", lw=1.0), fontsize=8.5, color="0.25")
ax.legend(frameon=False, loc="lower right")
save(fig, "dsm_vs_tq")

# ---- 3. u vs t_q (comparison scale; PID off-scale values annotated) ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_base_ax(ax, r"Query time $t_q$ (s)", "Applied ACS current (A)")
_acs_band(ax)
for s in ["LJ", "GA"]:
    d = data[s]
    ax.plot(d["tq"], d["u"], color=COLORS[s], marker=MARK[s], markersize=4,
            linewidth=1.5, markevery=ME, label=s)
_plot_pid_zoomed(ax)
ax.axhline(56.75, color="0.4", ls="--", lw=1.1)
ax.axhline(56.69, color="0.4", ls="--", lw=1.1)
ax.text(1158.6, 56.56, r"LJ/GA saturation limits $u_{\min}$/$u_{\max}$", color="0.3", fontsize=8)
ax.set_ylim(56.55, 56.85)
ax.legend(frameon=False, loc="upper right")
save(fig, "u_vs_tq")

# ---- 4. f* vs t_q ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_base_ax(ax, r"Query time $t_q$ (s)", r"Objective $f^*$")
_line(ax, ["LJ", "GA"], "f", ylog=True, ylim=(0.004, 120))
ax.legend(frameon=False, loc="best")
save(fig, "fopt_vs_tq")

# ---- 5. runtime vs t_q ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_base_ax(ax, r"Query time $t_q$ (s)", "Optimizer runtime (s)")
_line(ax, ["LJ", "GA"], "rt", ylim=(0, 560))
ax.axhline(120, color="0.45", ls=":", lw=1.1, label="120 s response-time reference")
ax.legend(frameon=False, loc="upper left")
save(fig, "runtime_vs_tq")

# ---- 6. variance boxplots at t_q = 1180 ----
fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.6))
labels = [r"$f^*$", "First action (A)", "Min predicted DSM (K)", "Optimizer runtime (s)"]
ylims = [(0.25, 0.55), (56.72, 56.76), (64.39, 64.43), (0, 380)]
for ax, key, lab, yl in zip(axes.ravel(), VAR, labels, ylims):
    bp = ax.boxplot([VAR[key][k] for k in ["LJ", "GA"]], labels=["LJ", "GA"],
                    widths=0.45, patch_artist=True, showmeans=True,
                    meanprops=dict(marker="D", markerfacecolor="0.2", markeredgecolor="0.2"),
                    medianprops=dict(color="black", linewidth=1.2))
    for patch, k in zip(bp["boxes"], ["LJ", "GA"]):
        patch.set_facecolor(COLORS[k]); patch.set_alpha(0.55)
    ax.set_title(lab, fontsize=10)
    ax.set_ylim(*yl)
    ax.grid(True, alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(r"20-seed distributions at fixed query state $t_q=1180$ s", fontsize=11)
fig.tight_layout()
save(fig, "variance_1180")

# ---- composed 3x2 overview, full-page portrait (thesis) ----
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9.5,
})
fig, axes = plt.subplots(3, 2, figsize=(7.4, 10.0))
panels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]

ax = axes[0, 0]; _base_ax(ax, r"$t_q$ (s)", "PCT (K)"); _acs_band(ax)
_line(ax, ["LJ", "GA", "PID"], "pct", ylim=(1135, 1185)); _target(ax)
ax.legend(frameon=False, loc="lower left", ncol=2, fontsize=8)

ax = axes[0, 1]; _base_ax(ax, r"$t_q$ (s)", "DSM (K)"); _acs_band(ax)
for s in ["LJ", "GA", "PID"]:
    d = data[s]
    ax.plot(d["tq"], 1244 - np.array(d["pct"]), color=COLORS[s], marker=MARK[s],
            markersize=4, linewidth=1.5, markevery=ME, label=s)
ax.set_ylim(45, 110); ax.legend(frameon=False, loc="lower right", fontsize=8)
ax.annotate("ACS activation transient", xy=(1160.5, 97.5), xytext=(1163.5, 103),
            arrowprops=dict(arrowstyle="->", color="0.3", lw=1.0), fontsize=7.5, color="0.25")

ax = axes[1, 0]; _base_ax(ax, r"$t_q$ (s)", "Applied ACS current (A)"); _acs_band(ax)
for s in ["LJ", "GA"]:
    d = data[s]
    ax.plot(d["tq"], d["u"], color=COLORS[s], marker=MARK[s], markersize=4,
            linewidth=1.5, markevery=ME, label=s)
_plot_pid_zoomed(ax)
ax.axhline(56.75, color="0.4", ls="--", lw=1.1); ax.axhline(56.69, color="0.4", ls="--", lw=1.1)
ax.text(1158.6, 56.56, r"LJ/GA bounds", color="0.3", fontsize=7)
ax.set_ylim(56.55, 56.85); ax.legend(frameon=False, loc="upper right", fontsize=8)

ax = axes[1, 1]; _base_ax(ax, r"$t_q$ (s)", r"Objective $f^*$")
_line(ax, ["LJ", "GA"], "f", ylog=True, ylim=(0.004, 120)); ax.legend(frameon=False, fontsize=8)

ax = axes[2, 0]; _base_ax(ax, r"$t_q$ (s)", "Optimizer runtime (s)")
_line(ax, ["LJ", "GA"], "rt", ylim=(0, 560))
ax.axhline(120, color="0.45", ls=":", lw=1.1, label="120 s reference")
ax.legend(frameon=False, loc="upper right", fontsize=8)

ax = axes[2, 1]
bp = ax.boxplot([VAR["f_opt"][k] for k in ["LJ", "GA"]], labels=["LJ", "GA"], widths=0.45,
                patch_artist=True, showmeans=True,
                meanprops=dict(marker="D", markerfacecolor="0.2", markeredgecolor="0.2"),
                medianprops=dict(color="black", linewidth=1.2))
for patch, k in zip(bp["boxes"], ["LJ", "GA"]):
    patch.set_facecolor(COLORS[k]); patch.set_alpha(0.55)
ax.set_title(r"$f^*$ at $t_q=1180$ s (20 seeds)", fontsize=10.5)
ax.set_ylabel(r"Objective $f^*$ (K$^2$)")
ax.grid(True, alpha=0.3); ax.spines[["top", "right"]].set_visible(False)

for ax, lab in zip(axes.ravel(), panels):
    ax.text(0.02, 0.98, lab, transform=ax.transAxes, va="top", ha="left",
            fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "phase2c_comparison.pdf")
fig.savefig(FIG / "phase2c_comparison.png")
plt.close(fig)
print("OK phase2c_comparison")
print("Figuras en:", FIG)
