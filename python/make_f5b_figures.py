"""Figures of F5b (plan §8.12): capacity map (step 1) and where the PCT peak comes from.

Reads only committed results (no MATLAB, no Julia):
  results/f5b1a_capacity/runs.csv, traj.csv.gz    F5b.1a capacity × ramp sweep
Usage: .venv\\Scripts\\python python/make_f5b_figures.py   -> docs/figures/fig_f5b_*.pdf/.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import csv
import gzip
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
RES = os.path.join(ROOT, "results")
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))
# reference palette (dataviz skill, light mode), same slots as make_closed_loop_figures.py:
# N2 = slot 1 (blue), PI = slot 2 (orange); legends + direct labels + markers as secondary encoding
BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df"
plt.rcParams.update({
    "font.family": "Inter Variable", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 8, "lines.linewidth": 1.5,
    "figure.dpi": 150, "savefig.bbox": "tight", "pdf.fonttype": 42,
})
W = 6.3
Y_REF = 1179.2711775986722
RAMP = 5.0


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)
    print("wrote", name)


def load(path):
    """CSV -> dict of numpy arrays (float where possible). pandas is avoided: its DLLs are
    blocked by Smart App Control on leo_acer."""
    rows = list(csv.DictReader(gzip.open(path, "rt") if path.endswith(".gz") else open(path)))
    out = {}
    for k in rows[0]:
        v = [x[k] for x in rows]
        try:
            out[k] = np.array([float(x) for x in v])
        except ValueError:
            out[k] = np.array(v)
    return out


def sel(d, m):
    return {k: v[m] for k, v in d.items()}


class G:   # per-controller curve along m_max (ramp RAMP)
    def __init__(self, d, ctrl):
        self.index = np.array(sorted(set(d["m_max"][(d["ctrl"] == ctrl) & np.isfinite(d["m_max"])])))
        grp = [sel(d, (d["ctrl"] == ctrl) & (d["m_max"] == m)) for m in self.index]
        self.n_lim = np.array([np.sum(x["limit"] == "true") for x in grp])
        self.pmean = np.array([x["pct_max"].mean() for x in grp]); self.pworst = np.array([x["pct_max"].max() for x in grp])
        self.hi = np.array([x["sat_hi"].mean() for x in grp]); self.free = np.array([x["free"].mean() for x in grp])


runs = load(os.path.join(RES, "f5b1a_capacity", "runs.csv"))
r = sel(runs, runs["ramp"] == RAMP)
g = {c: G(r, c) for c in ("N2", "PI")}
pi_inf = sel(runs, (runs["ctrl"] == "PI") & np.isinf(runs["m_max"]))

# ---- Fig 1: regime map along the capacity axis (ramp 5 kg/s per step) --------------------
fig, ax = plt.subplots(1, 3, figsize=(W, 2.6), layout="constrained")
for ctrl, col, mk, lab in (("N2", BLUE, "o", "N2+sesgo"), ("PI", ORANGE, "s", "PI del DET (misma capacidad)")):
    s = g[ctrl]
    ax[0].plot(s.index, s.n_lim, color=col, marker=mk, ms=4, label=lab)
    ax[1].plot(s.index, s.pworst, color=col, marker=mk, ms=4, label=lab)
    ax[1].plot(s.index, s.pmean, color=col, ls=(0, (3, 2)), lw=1)
ax[0].set(title="Realizaciones que llegan a 1220 K", xlabel=r"$\dot m_{\max}$ [kg/s]", ylabel="de 14")
ax[1].axhline(1220, color=MUTED, lw=0.9, ls=(0, (4, 2)))
ax[1].text(72, 1220.6, "fin del DET 1220 K", ha="right", va="bottom", fontsize=7, color=INK2)
ax[1].set(title="PCT máxima: peor (—) y media (- -)", xlabel=r"$\dot m_{\max}$ [kg/s]", ylabel="K")
s = g["N2"]
ax[2].plot(s.index, 100 * s.hi, color=BLUE, marker="o", ms=4)
ax[2].plot(s.index, 100 * s.free, color=BLUE, marker="o", ms=4, mfc="white", ls=(0, (3, 2)))
ax[2].text(42, 82, "en la cota", fontsize=7, color=INK2, va="bottom")
ax[2].text(60, 38, "interior", fontsize=7, color=INK2, ha="center", va="top")
ax[2].set(title="Primera acción de N2", xlabel=r"$\dot m_{\max}$ [kg/s]", ylabel="% de consultas", ylim=(0, 105))
for a in ax:
    a.axvspan(14, 22.5, color=GRID, alpha=0.5, lw=0)
    a.set_xlim(14, 72)
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc="outside lower center", ncol=2)
save(fig, "fig_f5b_capacity_map")

# ---- Fig 2: the PCT peak comes at the end of the ACS failure (N2, 40 kg/s, ramp 5) -------
tr = load(os.path.join(RES, "f5b1a_capacity", "traj.csv.gz"))
q = sel(r, (r["ctrl"] == "N2") & (r["m_max"] == 40))
fl, pm = q["fail_len"], q["pct_max"]
fig, ax = plt.subplots(1, 2, figsize=(W, 2.3), constrained_layout=True)
ax[0].scatter(fl, pm, s=22, color=BLUE, edgecolor="white", lw=0.6, zorder=3)
k, b = np.polyfit(fl, pm, 1)
xx = np.array([fl.min(), fl.max()])
ax[0].plot(xx, k * xx + b, color=MUTED, lw=0.9, ls=(0, (4, 2)))
ax[0].text(40, 1207.5, f"pendiente {k:.2f} K/s, r = {np.corrcoef(fl, pm)[0, 1]:.2f}",
           ha="left", va="center", fontsize=7, color=INK2)
ax[0].set(title="PCT máxima vs duración de la falla", xlabel="duración de la falla del ACS [s]", ylabel="K")
i = int(np.where(q["seed"] == 101)[0][0])
t = sel(tr, (tr["seed"] == 101) & (tr["ctrl"] == "N2") & (tr["m_max"] == 40) & (tr["ramp"] == RAMP))
fs = q["lead"][i] + q["fail_start"][i]; flen = q["fail_len"][i]
ax[1].axvspan(fs, fs + flen, color=GRID, alpha=0.7, lw=0)
ax[1].plot(t["t"], t["pct"], color=BLUE)
ax[1].axhline(Y_REF, color=INK2, lw=0.8, ls=":")
ax[1].text(t["t"].max(), Y_REF + 0.5, r"$y_{\rm ref}$", ha="right", va="bottom", fontsize=7, color=INK2)
ax[1].text(fs + flen / 2, t["pct"].max() + 0.5, "falla del ACS", ha="center", va="bottom", fontsize=7, color=INK2)
ax[1].set(title="Semilla 101: el pico llega al final de la falla", xlabel="t desde la activación [s]", ylabel="PCT [K]")
save(fig, "fig_f5b_peak_at_failure")
