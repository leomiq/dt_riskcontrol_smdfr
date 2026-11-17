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

# ---- Fig 3: failure-length survival P(D > t), F5a.3 extreme family vs DET-accelerated ------
import scipy.io as sio
db = sio.loadmat(os.path.join(ROOT, "models", "MCDET_conferencepaper_database.mat"))
lam = db["Lambda_t"].ravel() * 1e5; mu = db["Mu_t"].ravel() * 1e4; w = lam / lam.sum()
tt = np.linspace(0, 300, 601)
s_det = np.array([(w * np.exp(-mu * t)).sum() for t in tt])
s_old = np.clip((90 - tt) / 60, 0, 1)
fig, ax = plt.subplots(1, 1, figsize=(W * 0.62, 2.4), constrained_layout=True)
ax.axvspan(130, 150, color=GRID, alpha=0.7, lw=0)
ax.text(140, 0.86, "D que lleva de\n" r"$y_{\rm ref}$ a 1220 K", ha="center", va="center", fontsize=7, color=INK2)
ax.plot(tt, s_old, color=ORANGE, ls=(0, (4, 2)), label="F5a.3: una falla, D ~ U(30, 90) s")
ax.plot(tt, s_det, color=BLUE, label="DET acelerada: mezcla de exponenciales")
for a, q in ((0.05, 136.3), (0.01, 219.0)):
    ax.plot([q], [a], "o", color=BLUE, ms=4)
    ax.text(q + 4, a + 0.02, f"$q_{{{1 - a:.2f}}}$ = {q:.0f} s", fontsize=7, color=INK2, va="bottom")
ax.set(xlabel="duración de la falla t [s]", ylabel="P(D > t)", xlim=(0, 300), ylim=(0, 1.02),
       title="Duración de una falla del ACS")
ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.45))
save(fig, "fig_f5b_failure_duration")

# ---- Fig 4: cost-risk frontier of the risk bound (F5b.1c, DET-accelerated failures) ---------
rk = load(os.path.join(RES, "f5b1c_risk", "runs.csv"))
order = [("N2", "N2+sesgo", None), ("N2_risk_a0.1", "", 0.10), ("N2_risk_a0.05", r"$\alpha$=0.05", 0.05),
         ("N2_risk_a0.02", r"$\alpha$=0.02", 0.02), ("N2_risk_a0.01", r"$\alpha$=0.01", 0.01)]
base = sel(rk, rk["case"] == "N2"); bseed = {s: i for i, s in enumerate(base["seed"])}
fig, ax = plt.subplots(1, 2, figsize=(W, 2.4), constrained_layout=True)
xs, y_mean, y_worst, ms = [], [], [], []
for case, lab, a in order:
    q = sel(rk, rk["case"] == case)
    idx = [bseed[s] for s in q["seed"]]
    xs.append(np.mean(q["mass_plant"] - base["mass_plant"][idx]))
    y_mean.append(np.mean(q["pct_max"])); y_worst.append(np.max(q["pct_max"]))
    ms.append(np.mean(q["mean_dev"]))
ax[0].plot(xs, y_mean, color=BLUE, marker="o", ms=4, label="media")
ax[0].plot(xs, y_worst, color=BLUE, marker="o", ms=4, mfc="white", ls=(0, (3, 2)), label="peor")
for x, y, (case, lab, _) in zip(xs, y_worst, order):
    lab = r"N2 y $\alpha$=0.10" if case == "N2" else lab
    ax[0].text(x + 25, y - 1.0, lab, fontsize=6.5, color=INK2, ha="left", va="top")
ax[0].axhline(1220, color=MUTED, lw=0.9, ls=(0, (4, 2)))
ax[0].set(title="PCT máxima vs masa extra", xlabel="masa extra inyectada vs N2 [kg]", ylabel="K")
ax[0].legend(loc="center right")
ax[1].plot(ms, y_mean, color=BLUE, marker="o", ms=4)
for x, y, (case, lab, _) in zip(ms, y_mean, order):
    lab = r"N2 y $\alpha$=0.10" if case == "N2" else lab
    ax[1].text(x - 0.3, y, lab, fontsize=6.5, color=INK2, ha="right", va="center")
ax[1].set(title="Pre-enfriamiento vs PCT máx. media", xlim=(-14, 3), xlabel=r"desvío medio de la PCT respecto de $y_{\rm ref}$ [K]",
          ylabel="K")
save(fig, "fig_f5b_risk_frontier")

# ---- Fig 5: step 3.i.0, archived realisation, widened pump: N2, N2+DSM, thesis PID, DET PI -----
tj = load(os.path.join(RES, "f5b1b0_dsm_lambda", "traj_archived.csv"))
fig, ax = plt.subplots(1, 2, figsize=(W, 2.6), layout="constrained")
spec = [("N2", BLUE, "-", "N2+sesgo"), ("N2_DSM_U1220_l1000000", BLUE, (0, (1, 1.5)), r"N2+DSM, $\lambda$=10$^6$"),
        ("PID", ORANGE, "-", "PID de la tesis"), ("PI", MUTED, (0, (4, 2)), "PI del DET")]
for case, col, ls, lab in spec:
    q = sel(tj, (tj["act"] == "wide") & (tj["case"] == case) & (tj["n_req"] == 241))
    t = q["t"] - 1160.0
    ax[0].plot(t, q["pct"], color=col, ls=ls, label=lab, lw=1.3)
    ax[1].plot(t, np.clip(q["m_plant"], 0.1, None), color=col, ls=ls, label=lab, lw=1.3)
ax[0].axhline(Y_REF, color=INK2, lw=0.8, ls=":")
ax[0].text(120, Y_REF + 0.3, r"$y_{\rm ref}$", ha="right", va="bottom", fontsize=7, color=INK2)
ax[0].set(title="PCT (realización archivada, bomba 43.6 kg/s)", xlabel="t desde 1160 s [s]", ylabel="K",
          ylim=(1165, 1192))
ax[1].set_yscale("log")
ax[1].set(title="Caudal entregado", xlabel="t desde 1160 s [s]", ylabel="kg/s (log)")
h_, l_ = ax[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="outside lower center", ncol=4)
save(fig, "fig_f5b_dsm_lambda_archived")
