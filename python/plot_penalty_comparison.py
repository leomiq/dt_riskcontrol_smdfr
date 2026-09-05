"""Penalised vs baseline LJ realization: publication figures.

Same visual language as plot_phase2c_trajectories.py (serif, per-indicator PDF +
PNG, composed overview).

  figures/pen_pct_vs_tq.{pdf,png}      PCT(t_q)            baseline vs penalised
  figures/pen_dsm_vs_tq.{pdf,png}      DSM(t_q) in K and normalised
  figures/pen_u_vs_tq.{pdf,png}        applied ACS current
  figures/pen_mdot_vs_tq.{pdf,png}     coolant mass flow + extra injected
  figures/pen_jterms_vs_tq.{pdf,png}   J_error / J_mdot / J_psm over time
  figures/pen_dominance.{pdf,png}      where J_psm outweighs J_error
  figures/penalty_comparison.{pdf,png} composed 3x2 overview

Data sources
  baseline : phase2c_fullloop_pilot_LJ_seed1.mat           (lambda = 0)
  penalised: penalised_thr0p00238_lam1e06_LJ_seed1.mat     (lambda = 1e6)
  J terms  : penalised_run3.log -- the diary records every cost evaluation as
             "J_error: ..., J_mdot: ..., J_psm: ..." and every LJ pass as
             "Pass k complete. f_best = ...", with "[LJ] t=... u=..." closing
             each query. The optimum of each query is taken as the evaluation
             whose total J is minimal within that query's block.
"""
import re
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
OUT = HERE / "results_penalty_saturation"
RES = OUT
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

TARGET = 1179.2711776
U_Y = 1244.0
SPAN = U_Y - TARGET
DSM_THR = 0.00237737
LAMBDA = 1e6
XTICKS = list(range(1160, 1191, 10))
XLIM = (1158, 1192)

C_BASE = "#1f77b4"     # baseline, lambda = 0
C_PEN = "#d62728"      # penalised, lambda = 1e6
C_PSM = "#ff7f0e"
C_ERR = "#1f77b4"
C_MDOT = "#2ca02c"

ANCHOR_U = np.array([56.69, 56.696, 56.702, 56.708, 56.714, 56.72,
                     56.726, 56.732, 56.738, 56.744, 56.75])
ANCHOR_M = np.array([0.5798, 1.8086, 3.0373, 4.2660, 5.4946, 6.7231,
                     7.9516, 9.1800, 10.4083, 11.6366, 12.8648])


def massflow(u):
    """CurrentToMassFlow over the admissible range, from measured anchors."""
    return np.interp(u, ANCHOR_U, ANCHOR_M)


def load_run(path):
    d = loadmat(path, squeeze_me=True, struct_as_record=False)
    return {
        "t": np.asarray(d["time_snapshots"], float).ravel(),
        "pct": np.asarray(d["PCT"], float).ravel(),
        "u": np.asarray(d["u_V2P_trajectory"], float).ravel(),
        "rt": np.asarray(d["optimizer_runtimes"], float).ravel(),
    }


def parse_j_terms(logpath):
    """Per-query optimal J decomposition, recovered from the diary.

    Each query block ends with an "[LJ] t=..." line, so evaluations are
    accumulated until that marker and the minimum-total evaluation of the block
    is reported as that query's optimum.
    """
    # fitness_cost prints "J_error: <v>, J_mdot: <v>, J_psm: <v>." -- the
    # trailing period must not be swallowed into the last number, so the
    # mantissa/exponent form is matched explicitly rather than as a char class.
    num = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    re_j = re.compile(rf"J_error:\s*({num}),\s*J_mdot:\s*({num}),\s*J_psm:\s*({num})")
    re_q = re.compile(r"^\[LJ\]\s+t=([\d.]+)\s+s\s+u=([\d.]+)")
    re_p = re.compile(rf"Pass\s+(\d+)\s+complete\.\s+f_best0?\s*=\s*({num})")

    out = {"t": [], "J_error": [], "J_mdot": [], "J_psm": [], "J": [],
           "n_eval": [], "n_pass": [], "f_first": [], "f_last": []}
    buf, passes = [], []
    with open(logpath, "r", errors="ignore") as fh:
        for line in fh:
            m = re_j.search(line)
            if m:
                je, jm, jp = (float(x) for x in m.groups())
                buf.append((je, jm, jp))
                continue
            m = re_p.search(line)
            if m:
                passes.append(float(m.group(2)))
                continue
            m = re_q.match(line)
            if m and buf:
                arr = np.array(buf)
                tot = arr.sum(axis=1)
                k = int(np.argmin(tot))
                out["t"].append(float(m.group(1)))
                out["J_error"].append(arr[k, 0])
                out["J_mdot"].append(arr[k, 1])
                out["J_psm"].append(arr[k, 2])
                out["J"].append(tot[k])
                out["n_eval"].append(len(buf))
                out["n_pass"].append(len(passes))
                out["f_first"].append(passes[0] if passes else np.nan)
                out["f_last"].append(passes[-1] if passes else np.nan)
                buf, passes = [], []
    return {k: np.asarray(v, float) for k, v in out.items()}


def _base_ax(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(*XLIM)
    ax.set_xticks(XTICKS)
    ax.grid(True, alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)


def _acs_band(ax):
    """The ACS activation transient, as in the existing figures."""
    ax.axvspan(1160, 1161.5, color="0.85", alpha=0.5, lw=0)


def _dominance_bands(ax, t, dominant, label=True):
    """Shade the windows where J_psm outweighs J_error."""
    first = True
    for i, dom in enumerate(dominant):
        if not dom:
            continue
        lo = t[i] - 0.25
        hi = t[i] + 0.25
        ax.axvspan(lo, hi, color=C_PSM, alpha=0.16, lw=0,
                   label="J$_{psm}$ > J$_{error}$" if (first and label) else None)
        first = False


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png")
    plt.close(fig)
    print("OK", name)


# ---------------------------------------------------------------- load
base = load_run(RES / "phase2c_fullloop_pilot_LJ_seed1.mat")
pen = load_run(RES / "penalised_thr0p00238_lam1e06_LJ_seed1.mat")
J = parse_j_terms(RES / "penalised_run3.log")

n = min(len(base["t"]), len(pen["t"]))
t = base["t"][:n]
mb, mp = massflow(base["u"][:n]), massflow(pen["u"][:n])
dsm_b = (base["pct"][:n] - TARGET) / SPAN
dsm_p = (pen["pct"][:n] - TARGET) / SPAN

# align the parsed J terms onto the query grid
Jm = {k: np.full(n, np.nan) for k in ["J_error", "J_mdot", "J_psm", "J", "n_eval"]}
for i, tq in enumerate(t):
    hit = np.isclose(J["t"], tq, atol=1e-6)
    if hit.any():
        k = int(np.argmax(hit))
        for key in Jm:
            Jm[key][i] = J[key][k]
dominant = Jm["J_psm"] > Jm["J_error"]

print(f"queries: {n}, J terms parsed for {int(np.isfinite(Jm['J']).sum())}")
print(f"J_psm dominates J_error in {int(np.nansum(dominant))} of {n} queries")

# ---- 1. PCT (full range + zoom, because the two curves nearly coincide) ----
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.9))
ax = axs[0]
_base_ax(ax, r"Query time $t_q$ (s)", "PCT (K)")
_acs_band(ax)
_dominance_bands(ax, t, dominant)
ax.plot(t, base["pct"][:n], color=C_BASE, marker="o", markersize=4, markevery=6,
        label=r"baseline ($\lambda=0$)")
ax.plot(t, pen["pct"][:n], color=C_PEN, marker="s", markersize=4, markevery=6,
        label=r"penalised ($\lambda=10^{6}$)")
ax.axhline(TARGET, color="0.45", ls="--", lw=1.1, label=r"setpoint $T_p^*=1179.27$ K")
ax.set_ylim(1135, 1185)
ax.legend(frameon=False, loc="lower right", ncol=1, fontsize=8)
ax.set_title("full range", fontsize=9.5, color="0.35")

ax = axs[1]
_base_ax(ax, r"Query time $t_q$ (s)", "PCT (K)")
_dominance_bands(ax, t, dominant)
ax.plot(t, base["pct"][:n], color=C_BASE, marker="o", markersize=4, markevery=4,
        label=r"baseline ($\lambda=0$)")
ax.plot(t, pen["pct"][:n], color=C_PEN, marker="s", markersize=4, markevery=4,
        label=r"penalised ($\lambda=10^{6}$)")
ax.axhline(TARGET, color="0.45", ls="--", lw=1.1)
ax.set_xlim(1161, 1192)
lo = min(base["pct"][2:n].min(), pen["pct"][2:n].min())
hi = max(base["pct"][2:n].max(), pen["pct"][2:n].max())
pad = 0.08 * (hi - lo)
ax.set_ylim(lo - pad, hi + pad)
ax.set_title(r"zoom, transient excluded (curves nearly coincide)",
             fontsize=9.5, color="0.35")
ax.legend(frameon=False, loc="upper right", fontsize=8)
fig.tight_layout()
save(fig, "pen_pct_vs_tq")

# ---- 2. DSM ----
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.9))
ax = axs[0]
_base_ax(ax, r"Query time $t_q$ (s)", "DSM (K)")
_acs_band(ax)
ax.plot(t, U_Y - base["pct"][:n], color=C_BASE, marker="o", markersize=4, markevery=6,
        label=r"baseline ($\lambda=0$)")
ax.plot(t, U_Y - pen["pct"][:n], color=C_PEN, marker="s", markersize=4, markevery=6,
        label=r"penalised ($\lambda=10^{6}$)")
ax.set_ylim(45, 110)
ax.annotate("ACS activation transient", xy=(1160.5, 97.5), xytext=(1164, 100),
            arrowprops=dict(arrowstyle="->", color="0.3", lw=1.0), fontsize=8, color="0.25")
ax.legend(frameon=False, loc="lower right", fontsize=8)
ax = axs[1]
_base_ax(ax, r"Query time $t_q$ (s)", r"Normalised DSM $=(T_p-T_p^*)/(U_y-T_p^*)$")
_dominance_bands(ax, t, dominant)
ax.plot(t, dsm_b, color=C_BASE, marker="o", markersize=4, markevery=4, label="baseline")
ax.plot(t, dsm_p, color=C_PEN, marker="s", markersize=4, markevery=4, label="penalised")
ax.axhline(DSM_THR, color=C_PSM, ls="-.", lw=1.2,
           label=rf"threshold $={DSM_THR:.5f}$")
ax.axhline(0, color="0.6", lw=0.8)
# the t=1160.5 transient reaches -0.5 and would flatten everything else
ax.set_xlim(1161, 1192)
lo = min(dsm_b[2:n].min(), dsm_p[2:n].min())
hi = max(dsm_b[2:n].max(), dsm_p[2:n].max())
pad = 0.12 * (hi - lo)
ax.set_ylim(lo - pad, hi + pad)
ax.set_title("zoom, transient excluded", fontsize=9.5, color="0.35")
ax.legend(frameon=False, loc="upper right", fontsize=8)
fig.tight_layout()
save(fig, "pen_dsm_vs_tq")

# ---- 3. applied control ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_base_ax(ax, r"Query time $t_q$ (s)", "Applied ACS current (A)")
_acs_band(ax)
_dominance_bands(ax, t, dominant)
ax.step(t, base["u"][:n], where="mid", color=C_BASE, label=r"baseline ($\lambda=0$)")
ax.step(t, pen["u"][:n], where="mid", color=C_PEN, label=r"penalised ($\lambda=10^{6}$)")
ax.axhline(56.75, color="0.4", ls="--", lw=1.1)
ax.axhline(56.69, color="0.4", ls="--", lw=1.1)
ax.text(1158.6, 56.6915, r"$u_{\min}$", color="0.3", fontsize=8)
ax.text(1158.6, 56.7455, r"$u_{\max}$", color="0.3", fontsize=8)
ax.set_ylim(56.682, 56.758)
ax.legend(frameon=False, loc="lower right", fontsize=8)
save(fig, "pen_u_vs_tq")

# ---- 4. coolant mass flow and extra injected ----
fig, axs = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True,
                        gridspec_kw=dict(height_ratios=[1.35, 1]))
ax = axs[0]
ax.set_ylabel(r"ACS coolant mass flow (kg s$^{-1}$)")
ax.set_xlim(*XLIM); ax.set_xticks(XTICKS)
ax.grid(True, alpha=0.3); ax.spines[["top", "right"]].set_visible(False)
_acs_band(ax)
_dominance_bands(ax, t, dominant)
ax.step(t, mb, where="mid", color=C_BASE, label=r"baseline ($\lambda=0$)")
ax.step(t, mp, where="mid", color=C_PEN, label=r"penalised ($\lambda=10^{6}$)")
ax.legend(frameon=False, loc="lower right", fontsize=8)
ax = axs[1]
_base_ax(ax, r"Query time $t_q$ (s)", r"Extra coolant injected (kg s$^{-1}$)")
_acs_band(ax)
d = mp - mb
ax.bar(t[d >= 0], d[d >= 0], width=0.42, color=C_MDOT, alpha=0.85, label="more than baseline")
ax.bar(t[d < 0], d[d < 0], width=0.42, color="0.55", alpha=0.85, label="less than baseline")
ax.axhline(0, color="0.3", lw=0.9)
ax.set_ylim(-13, 13)
ax.text(0.99, 0.04,
        rf"net $\Sigma\Delta\dot m = {d.sum():+.2f}$   mean $= {d.mean():+.3f}$",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color="0.25")
ax.legend(frameon=False, loc="upper right", fontsize=8, ncol=2)
fig.tight_layout()
save(fig, "pen_mdot_vs_tq")

# ---- 5. objective decomposition over time ----
fig, ax = plt.subplots(figsize=(7.2, 4.2))
_base_ax(ax, r"Query time $t_q$ (s)", r"Objective contribution")
_acs_band(ax)
ax.plot(t, Jm["J_error"], color=C_ERR, marker="o", markersize=3.6, markevery=5,
        label=r"$J_{error}$ (tracking)")
ax.plot(t, Jm["J_psm"], color=C_PSM, marker="s", markersize=3.6, markevery=5,
        label=r"$J_{psm}$ (DSM penalty)")
mdot_plot = np.where(Jm["J_mdot"] > 0, Jm["J_mdot"], np.nan)
ax.plot(t, mdot_plot, color=C_MDOT, marker="^", markersize=3.6, markevery=5,
        label=r"$J_{\dot m}$ (rate)")
ax.plot(t, Jm["J"], color="0.25", ls=":", lw=1.3, label=r"$J$ (total)")
ax.set_yscale("log")
ax.legend(frameon=False, loc="lower left", ncol=2, fontsize=8)
save(fig, "pen_jterms_vs_tq")

# ---- 6. dominance: share of the cost ----
fig, axs = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True)
ax = axs[0]
ax.set_ylabel(r"Share of total $J$")
ax.set_xlim(*XLIM); ax.set_xticks(XTICKS)
ax.grid(True, alpha=0.3); ax.spines[["top", "right"]].set_visible(False)
_acs_band(ax)
share_psm = Jm["J_psm"] / Jm["J"]
share_err = Jm["J_error"] / Jm["J"]
ax.fill_between(t, 0, share_psm, color=C_PSM, alpha=0.55, label=r"$J_{psm}$")
ax.fill_between(t, share_psm, share_psm + share_err, color=C_ERR, alpha=0.45,
                label=r"$J_{error}$")
ax.axhline(0.5, color="0.3", ls="--", lw=1.0)
ax.set_ylim(0, 1)
ax.legend(frameon=False, loc="center right", fontsize=8)
ax = axs[1]
_base_ax(ax, r"Query time $t_q$ (s)", r"$J_{psm}\,/\,J_{error}$")
_acs_band(ax)
ratio = Jm["J_psm"] / Jm["J_error"]
ax.plot(t, ratio, color=C_PSM, marker="s", markersize=3.6, markevery=4)
ax.axhline(1, color="0.3", ls="--", lw=1.1, label="parity")
ax.set_yscale("log")
ax.legend(frameon=False, loc="upper right", fontsize=8)
fig.tight_layout()
save(fig, "pen_dominance")

# ---- composed overview ----
plt.rcParams.update({"font.size": 11, "axes.labelsize": 11.5, "legend.fontsize": 8.5})
fig, axes = plt.subplots(3, 2, figsize=(7.6, 10.2))

# (a) and (b) are zoomed past the t=1160.5 transient: at full scale that single
# -0.5 excursion compresses everything else into a flat line.
ax = axes[0, 0]; _base_ax(ax, r"$t_q$ (s)", "PCT (K)")
_dominance_bands(ax, t, dominant, label=False)
ax.plot(t, base["pct"][:n], color=C_BASE, marker="o", markersize=3.5, markevery=4, label=r"$\lambda=0$")
ax.plot(t, pen["pct"][:n], color=C_PEN, marker="s", markersize=3.5, markevery=4, label=r"$\lambda=10^{6}$")
ax.axhline(TARGET, color="0.45", ls="--", lw=1.0, label="setpoint")
ax.set_xlim(1161, 1192)
_lo = min(base["pct"][2:n].min(), pen["pct"][2:n].min())
_hi = max(base["pct"][2:n].max(), pen["pct"][2:n].max())
ax.set_ylim(_lo - 0.08*(_hi-_lo), _hi + 0.08*(_hi-_lo))
ax.legend(frameon=False, loc="upper right", fontsize=8)

ax = axes[0, 1]; _base_ax(ax, r"$t_q$ (s)", "Normalised DSM")
_dominance_bands(ax, t, dominant, label=False)
ax.plot(t, dsm_b, color=C_BASE, marker="o", markersize=3.5, markevery=4, label=r"$\lambda=0$")
ax.plot(t, dsm_p, color=C_PEN, marker="s", markersize=3.5, markevery=4, label=r"$\lambda=10^{6}$")
ax.axhline(DSM_THR, color=C_PSM, ls="-.", lw=1.1, label="threshold")
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlim(1161, 1192)
_lo = min(dsm_b[2:n].min(), dsm_p[2:n].min())
_hi = max(dsm_b[2:n].max(), dsm_p[2:n].max())
ax.set_ylim(_lo - 0.12*(_hi-_lo), _hi + 0.12*(_hi-_lo))
ax.legend(frameon=False, loc="upper right", fontsize=8)

ax = axes[1, 0]; _base_ax(ax, r"$t_q$ (s)", "Applied current (A)"); _acs_band(ax)
_dominance_bands(ax, t, dominant, label=False)
ax.step(t, base["u"][:n], where="mid", color=C_BASE, label=r"$\lambda=0$")
ax.step(t, pen["u"][:n], where="mid", color=C_PEN, label=r"$\lambda=10^{6}$")
ax.axhline(56.75, color="0.4", ls="--", lw=1.0); ax.axhline(56.69, color="0.4", ls="--", lw=1.0)
ax.set_ylim(56.682, 56.758); ax.legend(frameon=False, loc="lower right", fontsize=8)

ax = axes[1, 1]; _base_ax(ax, r"$t_q$ (s)", r"$\Delta\dot m$ (kg s$^{-1}$)"); _acs_band(ax)
ax.bar(t[d >= 0], d[d >= 0], width=0.42, color=C_MDOT, alpha=0.85)
ax.bar(t[d < 0], d[d < 0], width=0.42, color="0.55", alpha=0.85)
ax.axhline(0, color="0.3", lw=0.9); ax.set_ylim(-13, 13)
ax.text(0.98, 0.05, rf"net ${d.sum():+.2f}$", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=8, color="0.25")

ax = axes[2, 0]; _base_ax(ax, r"$t_q$ (s)", "Objective contribution"); _acs_band(ax)
ax.plot(t, Jm["J_error"], color=C_ERR, marker="o", markersize=3.2, markevery=5, label=r"$J_{error}$")
ax.plot(t, Jm["J_psm"], color=C_PSM, marker="s", markersize=3.2, markevery=5, label=r"$J_{psm}$")
ax.plot(t, Jm["J"], color="0.25", ls=":", lw=1.2, label=r"$J$")
ax.set_yscale("log"); ax.legend(frameon=False, loc="lower left", fontsize=8, ncol=3)

ax = axes[2, 1]; _base_ax(ax, r"$t_q$ (s)", r"$J_{psm}/J_{error}$"); _acs_band(ax)
ax.plot(t, ratio, color=C_PSM, marker="s", markersize=3.2, markevery=4)
ax.axhline(1, color="0.3", ls="--", lw=1.0)
ax.set_yscale("log")

for ax, lab in zip(axes.ravel(), ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]):
    ax.text(0.02, 0.98, lab, transform=ax.transAxes, va="top", ha="left",
            fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "penalty_comparison.pdf")
fig.savefig(FIG / "penalty_comparison.png")
plt.close(fig)
print("OK penalty_comparison")

# ---------------------------------------------------------------- summary
print("\n" + "=" * 70)
print("PENALTY DOMINANCE WINDOWS (J_psm > J_error)")
print("=" * 70)
idx = np.where(dominant)[0]
if idx.size:
    runs, start = [], idx[0]
    for a, b in zip(idx, idx[1:]):
        if b != a + 1:
            runs.append((start, a)); start = b
    runs.append((start, idx[-1]))
    for a, b in runs:
        print(f"  t = {t[a]:7.1f} .. {t[b]:7.1f} s   "
              f"({b-a+1:2d} queries)  max ratio {np.nanmax(ratio[a:b+1]):9.1f}")
else:
    print("  none")
print(f"\ntotal: {int(np.nansum(dominant))} of {n} queries "
      f"({100*np.nansum(dominant)/n:.1f}%)")
print(f"J_psm/J_error: median {np.nanmedian(ratio):.2f}, "
      f"max {np.nanmax(ratio):.1f}, min {np.nanmin(ratio):.4f}")
print(f"\nfigures in {FIG}")
