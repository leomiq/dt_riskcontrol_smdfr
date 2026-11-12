"""Figures of the closed-loop chapter of docs/estado_mpc_monolitico.tex (11-11-2026).

Reads only committed results (no MATLAB, no Julia):
  results/f1b_hrm_replica/validation_1180_run1.csv   F1b: HRM replica vs the archived run
  results/f5a_closed_loop/trajectories.csv            F5a.1: closed loop on the archived realisation
  results/f5a2_monte_carlo/runs_extreme*.csv          F5a.3: extreme family

Usage: .venv\\Scripts\\python python/make_closed_loop_figures.py   -> docs/figures/fig_cl_*.pdf
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
RES = os.path.join(ROOT, "results")
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))
# reference palette (dataviz skill, light mode); 4 slots validated (contrast WARN -> legends,
# direct labels and line styles as secondary encoding)
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
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


def save(fig, name):
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print("wrote", name)


def limits(ax, yref=True):
    for v, lab in ((1220, "1220 K: estado final del árbol DET (Simulink se detiene)"),):
        ax.axhline(v, color=MUTED, lw=0.9, ls=(0, (4, 2)))
        ax.text(ax.get_xlim()[1], v - 0.6, lab + " ", va="top", ha="right", fontsize=7, color=INK2)
    if yref:
        ax.axhline(Y_REF, color=INK2, lw=0.8, ls=":")


# --- F1b: HRM replica vs archive
d = pd.read_csv(os.path.join(RES, "f1b_hrm_replica", "validation_1180_run1.csv"))
fig, ax = plt.subplots(1, 2, figsize=(W, 2.2))
ax[0].plot(d.t, d.pct_archive, color=BLUE, label="Simulink (archivo)")
ax[0].plot(d.t, d.pct_replica, color=ORANGE, ls=(0, (3, 2)), label="réplica Julia")
ax[0].set(xlabel="t [s]", ylabel="PCT [K]", title="PCT del HRM, 0–1160 s")
ax[0].legend(loc="upper left")
err = np.abs(d.pct_replica - d.pct_archive).clip(lower=1e-12)
ax[1].semilogy(d.t, err, color=BLUE, lw=1)
ax[1].set(xlabel="t [s]", ylabel="|error| [K]", title="Error absoluto de la réplica")
act = d.t[d.act > 0]
if len(act):
    ax[1].axvline(act.iloc[0], color=MUTED, lw=0.8, ls=":")
    ax[1].text(act.iloc[0], ax[1].get_ylim()[1], " ACS activo", va="top", fontsize=7, color=INK2)
fig.tight_layout()
save(fig, "fig_cl_replica.pdf")

# --- F5a.1: closed loop on the archived realisation
t = pd.read_csv(os.path.join(RES, "f5a_closed_loop", "trajectories.csv"))
sel = [("arch", "zero", False, "sin ACS", MUTED, "-"),
       ("arch", "LJ", False, "LJ archivado", YELLOW, (0, (1, 1.5))),
       ("arch", "N2", True, "N2+sesgo, caja archivada", BLUE, "-"),
       ("phys", "N2", True, "N2+sesgo, actuador físico", ORANGE, (0, (4, 2)))]
fig, ax = plt.subplots(1, 2, figsize=(W, 2.4))
for a, c, b, lab, col, ls in sel:
    s = t[(t.act == a) & (t.ctrl == c) & (t.bias == b)]
    ax[0].plot(s.t_q, s.pct, color=col, ls=ls, label=lab)
    if c != "zero":
        ax[1].plot(s.t_q, s.m, color=col, ls=ls, label=lab)
ax[0].axhline(Y_REF, color=INK2, lw=0.8, ls=":")
ax[0].text(t.t_q.max(), Y_REF, "$y_{\\rm ref}$", va="bottom", ha="right", fontsize=7, color=INK2)
ax[0].set(xlabel="t [s]", ylabel="PCT del HRM [K]", title="PCT en lazo cerrado")
ax[1].set(xlabel="t [s]", ylabel="caudal ACS [kg/s]", title="Caudal decidido")
ax[0].legend(loc="upper right", fontsize=7)
fig.tight_layout()
save(fig, "fig_cl_loop.pdf")

# --- F5a.3: extreme family, one realisation
r = pd.read_csv(os.path.join(RES, "f5a2_monte_carlo", "runs_extreme.csv"))
x = pd.read_csv(os.path.join(RES, "f5a2_monte_carlo", "runs_extreme_traj.csv"))
seed = 101
info = r[r.seed == seed].iloc[0]
leak_act = info.drift - info.A * info.drift * (info.t_act - 500.0) ** (1 / info.B)  # |ṁ_leak|
fa, fb = info.fail_start - info.t_act, info.fail_start + info.fail_len - info.t_act
ext = [("none", "zero", False, "sin ACS", MUTED, "-"),
       ("arch", "N2", True, "N2+sesgo, caja archivada (= máx.)", BLUE, "-"),
       ("det", "PI", False, "PI del DET", AQUA, (0, (1, 1.5))),
       ("phys", "N2", True, "N2+sesgo, actuador físico", ORANGE, (0, (4, 2)))]
fig, ax = plt.subplots(2, 1, figsize=(W, 4.2), sharex=True, gridspec_kw={"height_ratios": [3, 2]})
for a, c, b, lab, col, ls in ext:
    s = x[(x.seed == seed) & (x.act == a) & (x.ctrl == c) & (x.bias == b)]
    ax[0].plot(s.t, s.pct, color=col, ls=ls, label=lab)
    ax[1].plot(s.t, s.m_plant.clip(lower=0.1), color=col, ls=ls, label=lab)
    if len(s) and s.pct.iloc[-1] > 1220:
        ax[0].plot(s.t.iloc[-1], s.pct.iloc[-1], "x", color=col, ms=7, mew=1.5)
for a_ in ax:
    a_.axvspan(fa, fb, color=GRID, alpha=0.7, lw=0)
ax[0].text((fa + fb) / 2, ax[0].get_ylim()[0], "falla del ACS", ha="center", va="bottom", fontsize=7, color=INK2)
limits(ax[0])
ax[0].set(ylabel="PCT del HRM [K]",
          title=f"Realización extrema {seed}: fuga de {leak_act:.1f} kg/s en la activación, "
                f"η = {info.eta:.2f}, falla de {info.fail_len:.0f} s")
ax[1].legend(loc="upper center", ncol=2, fontsize=7, bbox_to_anchor=(0.5, -0.32))
ax[1].set_yscale("log")
ax[1].set(xlabel="tiempo desde la activación del ACS [s]", ylabel="caudal recibido [kg/s]")
fig.tight_layout()
save(fig, "fig_cl_extreme.pdf")

# --- F5a.3: maximum PCT per realisation
r["lab"] = r.act + " " + r.ctrl + np.where(r.bias, "+b", "")
order = r[r.lab == "none zero"].sort_values("t_limit", ascending=False).seed.tolist()
pos = {s: i for i, s in enumerate(order)}
mk = [("none zero", "sin ACS", MUTED, "o"), ("arch N2+b", "N2+sesgo, caja archivada (= máx.)", BLUE, "s"),
      ("det PI", "PI del DET", AQUA, "^"), ("phys N2+b", "N2+sesgo, actuador físico", ORANGE, "D")]
fig, ax = plt.subplots(figsize=(W, 2.8))
for k, (lab, name, col, m) in enumerate(mk):
    s = r[r.lab == lab]
    xs = np.array([pos[q] for q in s.seed]) + (k - 1.5) * 0.12
    ax.scatter(xs, s.pct_max, s=22, color=col, marker=m, label=name, edgecolor="white", linewidth=0.6, zorder=3)
    bd = s[s.breakdown]
    if len(bd):
        ax.scatter(np.array([pos[q] for q in bd.seed]) + (k - 1.5) * 0.12, bd.pct_max, s=70,
                   facecolor="none", edgecolor=INK, linewidth=0.8, zorder=4, label="colapso del modelo")
ax.set_xticks(range(len(order)), [str(s) for s in order], fontsize=7)
ax.set_xlim(-0.6, len(order) - 0.4)
limits(ax, yref=False)
ax.set(xlabel="realización (de menos a más severa: tiempo hasta 1220 K sin ACS)", ylabel="PCT máxima [K]",
       title="Familia extrema: PCT máxima en 240 s (sobre 1220 K la corrida termina)")
ax.legend(loc="upper left", ncol=2, fontsize=7, bbox_to_anchor=(0, -0.22))
fig.tight_layout()
save(fig, "fig_cl_extreme_max.pdf")
