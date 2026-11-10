"""Figures of the status report docs/estado_mpc_monolitico.tex (6-11-2026).

Reads only committed results (no MATLAB, no Julia):
  results/f1_rom_replica/*.csv / rom_reference_*.mat   F1 replica, F2 NLP vs reference action
  results/e2_optimizers_1180/e2_summary.csv            E2 optimizers
  results/f4_risk/*                                    F4 ensemble, calibration, actuator sweep
  models/LJ_test_sampled_vel_ACS_2_loop_1180.mat       recorded leak (forecast window)

Usage: .venv\\Scripts\\python python/make_status_figures.py   -> docs/figures/*.pdf
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import scipy.io
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
os.makedirs(OUT, exist_ok=True)

font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))
# reference palette (dataviz skill, light mode): first three slots validate all-pairs
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df"
plt.rcParams.update({
    "font.family": "Inter Variable", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 8, "lines.linewidth": 1.5,
    "figure.dpi": 150, "savefig.bbox": "tight", "pdf.fonttype": 42,
})
W = 6.3  # text width in inches (A4, 25 mm margins)


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def num(a):
    return a.astype(float)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print("wrote", name)


# --- F1: replica vs Simulink (t_q = 1180 s window)
d = read_csv(os.path.join(ROOT, "results", "f1_rom_replica", "replica_vs_simulink_1180.csv"))
t, ys, yr = num(d["t"]), num(d["y_simulink"]), num(d["y_replica"])
fig, ax = plt.subplots(1, 2, figsize=(W, 2.2))
ax[0].plot(t, ys, color=BLUE, label="Simulink (ROM)")
ax[0].plot(t, yr, color=ORANGE, ls=(0, (4, 3)), label="Réplica algebraica (Julia)")
ax[0].axvline(1180, color=MUTED, lw=0.8, ls=":")
ax[0].text(1180.3, ys.max(), "$t_q$", color=INK2, va="top")
ax[0].set(xlabel="t [s]", ylabel="PCT estimada $\\hat y$ [K]", title="Trayectoria de la PCT")
ax[0].legend(loc="center right")
err = np.maximum(np.abs(yr - ys), 1e-9)
ax[1].semilogy(t, err, color=BLUE, marker="o", ms=3)
ax[1].axvline(1180, color=MUTED, lw=0.8, ls=":")
ax[1].set(xlabel="t [s]", ylabel="|réplica − Simulink| [K]", title="Error absoluto (escala log.)")
ax[1].text(1179.6, 2e-9, "calentamiento", color=INK2, fontsize=7.5, ha="right")
ax[1].text(1180.4, 2e-9, "horizonte", color=INK2, fontsize=7.5, ha="left")
save(fig, "fig_replica.pdf")

# --- F2: NLP vs reference action at both queries
fig, ax = plt.subplots(2, 2, figsize=(W, 3.9), sharex="col")
for j, tq in enumerate(("1161", "1180")):
    d = read_csv(os.path.join(ROOT, "results", "f1_rom_replica", f"f2_nlp_vs_lj_{tq}.csv"))
    r = scipy.io.loadmat(os.path.join(ROOT, "results", "f1_rom_replica", f"rom_reference_{tq}.mat"),
                         variable_names=["target"])
    yref = float(r["target"].item())
    t = num(d["t"])
    ref_lbl = "LJ archivado" if tq == "1180" else "Acción archivada del lazo"
    ax[0, j].step(t, num(d["mdot_LJ"]), where="post", color=ORANGE, label=ref_lbl)
    ax[0, j].step(t, num(d["mdot_NLP"]), where="post", color=BLUE, label="NLP monolítico")
    ax[0, j].set(title=f"$t_q$ = {tq} s", ylabel="$\\dot m_{\\rm ACS}$ [kg/s]" if j == 0 else None)
    ax[0, j].set_ylim(-0.5, 14.5)
    ax[1, j].plot(t, num(d["y_LJ"]), color=ORANGE, label=ref_lbl)
    ax[1, j].plot(t, num(d["y_NLP"]), color=BLUE, label="NLP monolítico")
    ax[1, j].axhline(yref, color=MUTED, lw=1, ls=(0, (4, 3)))
    ax[1, j].text(t[-1], yref, " $y_{\\rm ref}$", color=INK2, va="bottom", ha="right")
    ax[1, j].set(xlabel="t [s]", ylabel="$\\hat y$ [K]" if j == 0 else None)
    ax[1, j].ticklabel_format(axis="y", useOffset=False)
from matplotlib.lines import Line2D
fig.legend([Line2D([], [], color=ORANGE), Line2D([], [], color=BLUE)],
           ["Acción de referencia (1161 s: aplicada por el lazo; 1180 s: LJ archivado)", "NLP monolítico"],
           loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.04))
save(fig, "fig_f2_nlp.pdf")

# --- E2: optimizers at t_q = 1180 s (archived LJ off scale, stated in caption)
d = read_csv(os.path.join(ROOT, "results", "e2_optimizers_1180", "e2_summary.csv"))
labels, J, col = [], [], []
for m, s, j in zip(d["method"], d["seed"], num(d["J"])):
    if m == "LJ_archived":
        continue
    labels.append("NLP / QP" if m == "NLP" else f"{m}, semilla {int(float(s))}")
    J.append(j)
    col.append({"NLP": BLUE, "LJ": ORANGE, "GA": AQUA}[m])
fig, ax = plt.subplots(figsize=(W, 2.0))
y = np.arange(len(J))[::-1]
ax.barh(y, J, color=col, height=0.6)
for yi, ji in zip(y, J):
    ax.text(ji + 0.004, yi, f"{ji:.3f}", va="center", color=INK, fontsize=7.5)
ax.set_yticks(y, labels)
ax.set_xlabel("Costo P0 re-simulado en Simulink (menor es mejor)")
ax.grid(axis="y", visible=False)
ax.set_xlim(0, 0.43)
save(fig, "fig_e2.pdf")

# --- F4: forecast ensemble vs realized leak, t_q = 1161 s
e = scipy.io.loadmat(os.path.join(ROOT, "results", "f4_risk", "leak_ensemble_1161.mat"))
a = scipy.io.loadmat(os.path.join(ROOT, "models", "LJ_test_sampled_vel_ACS_2_loop_1180.mat"),
                     variable_names=["X_1_trajectory"])
X1 = a["X_1_trajectory"].ravel()
i = 1 + 2320 + 2           # 0-based index of t = 1161 s: [0; t = 0..1160 s (2321); 1160.5; 1161; ...]
H = 25
tf = 1161 + 0.5 * np.arange(1, H + 1)
th = 1161 - 0.5 * np.arange(13)[::-1]
L = e["Y"][:, :H, 0]
lo, med, hi = np.percentile(L, [5, 50, 95], axis=0)
real = e["real"].ravel()[:H]
fig, ax = plt.subplots(figsize=(W, 2.3))
ax.fill_between(tf, lo, hi, color=BLUE, alpha=0.18, lw=0, label="Banda 5–95 % del ensamble (1000)")
ax.plot(tf, med, color=BLUE, label="Mediana del pronóstico")
ax.plot(th, X1[i - 12:i + 1], color=INK2, lw=1.2, label="Fuga medida (fin de la ventana)")
assert abs(X1[i + 1] - e["real"].ravel()[0]) < 1e-12, "index of t_q in the record"
ax.plot(tf, real, color=ORANGE, marker="o", ms=3, label="Fuga real")
ax.axvline(1161, color=MUTED, lw=0.8, ls=":")
ax.set(xlabel="t [s]", ylabel="$\\dot m_{\\rm leak}$ [kg/s]")
ax.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.28))
save(fig, "fig_forecast.pdf")

# --- F4: calibration over 16 queries
d = read_csv(os.path.join(ROOT, "results", "f4_risk", "forecast_calibration.csv"))
tq, cov, z = num(d["t_q"]), num(d["coverage90"]), num(d["z"])
fig, ax = plt.subplots(1, 2, figsize=(W, 2.1))
ax[0].bar(tq, cov, width=0.35, color=BLUE)
ax[0].axhline(0.9, color=MUTED, ls=(0, (4, 3)), lw=1)
ax[0].text(tq[0] - 0.2, 0.92, "nominal 0.90", color=INK2, fontsize=7.5)
ax[0].set(ylim=(0, 1.05), xlabel="$t_q$ [s]", ylabel="Cobertura", title="Cobertura de la banda del 90 %")
ax[1].axhspan(-2, 2, color=GRID, lw=0)
ax[1].plot(tq, z, color=BLUE, marker="o", ms=4)
ax[1].text(tq[0], 1.0, "|z| ≤ 2 (calibrado)", color=INK2, fontsize=7.5)
ax[1].set(ylim=(-3, 10), xlabel="$t_q$ [s]", ylabel="z",
          title="Error de la fuga integrada (en desv. est.)")
save(fig, "fig_calibration.pdf")

# --- actuator sweep: optimal cost vs flow bound, soft vs hard ramp
d = read_csv(os.path.join(ROOT, "results", "f4_risk", "actuator_sweep.csv"))
fig, ax = plt.subplots(1, 2, figsize=(W, 2.3))
for j, tqs in enumerate(("1161", "1180")):
    for ramp, c, lab in (("soft", BLUE, "Rampa blanda (P0)"), ("hard", ORANGE, "Rampa dura, 5 kg/s por paso")):
        k = (d["t_q"] == tqs) & (d["ramp"] == ramp)
        ax[j].plot(num(d["m_max"][k]), num(d["J_p50"][k]), color=c, marker="o", ms=4, label=lab)
    ax[j].set_xscale("log"); ax[j].set_yscale("log")
    ax[j].yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax[j].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[j].set_yticks([3, 5, 10, 20, 50] if tqs == "1161" else [0.06, 0.1, 0.15, 0.25])
    ax[j].axvline(12.8648, color=MUTED, lw=0.8, ls=":")
    ax[j].text(13.5, ax[j].get_ylim()[0] * 1.15, "cota archivada\n(56.75 A)", color=INK2, fontsize=7)
    ax[j].set(title=f"$t_q$ = {tqs} s", xlabel="$\\dot m_{\\rm ACS,max}$ [kg/s]",
              ylabel="Costo óptimo P0" if j == 0 else None)
ax[0].legend(loc="upper right")
save(fig, "fig_actuator.pdf")

# --- F3b: NLP vs QP (N1, N2) control trajectories and their differences
fig, ax = plt.subplots(2, 2, figsize=(W, 3.9), sharex="col",
                       gridspec_kw={"height_ratios": [1.6, 1]})
for j, tq in enumerate(("1161", "1180")):
    d = read_csv(os.path.join(ROOT, "results", "f3_certificates", f"trajectories_{tq}.csv"))
    t, nlp = num(d["t"]), num(d["NLP"])
    ax[0, j].step(t, nlp, where="post", color=BLUE, lw=2.2, label="N3 (NLP)")
    ax[0, j].step(t, num(d["N1"]), where="post", color=ORANGE, ls=(0, (4, 2)), label="N1 (QP, LTI)")
    ax[0, j].step(t, num(d["N2_ode1be"]), where="post", color=AQUA, ls=(0, (1, 1.5)), lw=1.8,
                  label="N2 (QP, LTV, convergido)")
    ax[0, j].set(title=f"$t_q$ = {tq} s", ylabel=r"$\dot m_{\rm ACS}$ [kg/s]" if j == 0 else None,
                 ylim=(-0.5, 14.5))
    for name, col in (("N1", ORANGE), ("N2_ode1be", AQUA)):
        ax[1, j].semilogy(t, np.maximum(np.abs(num(d[name]) - nlp), 1e-12), color=col, marker="o", ms=2.5)
    ax[1, j].set(xlabel="t [s]", ylabel="|QP − NLP| [kg/s]" if j == 0 else None, ylim=(1e-12, 1))
h_, l_ = ax[0, 1].get_legend_handles_labels()
fig.legend(h_, l_, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.05))
save(fig, "fig_nlp_vs_qp.pdf")

# --- F3b: SQP convergence of N2 at t_q = 1180 s
d = read_csv(os.path.join(ROOT, "results", "f3_certificates", "n2_iterations_1180.csv"))
k = d["nominal_scheme"] == "ode1be"
it = num(d["iter"][k])
fig, ax = plt.subplots(figsize=(W * 0.62, 2.1))
ax.semilogy(it, num(d["step"][k]), color=BLUE, marker="o", ms=4, label=r"$\max_k|u^{(i)}_k-u^{(i-1)}_k|$ (paso entre iteraciones)")
ax.semilogy(it, num(d["max_du_vs_NLP"][k]), color=ORANGE, marker="s", ms=4, label=r"$\max_k|u^{(i)}_k-u^{\rm NLP}_k|$ (distancia al NLP)")
ax.set(xlabel="iteración $i$", ylabel="kg/s", xticks=it)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.45))
save(fig, "fig_sqp.pdf")

# --- F3b: structure of the sparse monolithic QP vs the condensed QP (N1, 1180 s)
import h5py  # MAT.jl writes MATLAB v7.3 (HDF5)
with h5py.File(os.path.join(ROOT, "results", "f3_certificates", "sparsity_1180.mat"), "r") as s:
    r = np.array(s["rows"]).ravel().astype(int) - 1
    cidx = np.array(s["cols"]).ravel().astype(int) - 1
    nr, nc = int(np.array(s["nrows"]).item()), int(np.array(s["ncols"]).item())
    G = np.array(s["G"]).T          # HDF5 stores column-major
fig, ax = plt.subplots(1, 2, figsize=(W, 2.9), gridspec_kw={"width_ratios": [1.5, 1]})
ax[0].scatter(cidx, r, s=0.15, color=BLUE, marker="s", lw=0)
ax[0].set(xlim=(-5, nc + 5), ylim=(nr + 5, -5), title=f"Monolítico disperso: {nr} × {nc}",
          xlabel="variables (estados ζ, luego 26 controles)", ylabel="restricciones de dinámica")
ax[0].grid(False)
im = ax[1].imshow(np.log10(np.maximum(np.abs(G), 1e-12)), cmap="Blues", vmin=-6, vmax=np.log10(np.abs(G).max()))
ax[1].set(title="Condensado: $G$ (26 × 26)", xlabel="control $u_j$", ylabel=r"salida $\hat y_k$")
ax[1].grid(False)
cb = fig.colorbar(im, ax=ax[1], fraction=0.046, pad=0.04)
cb.set_label("log$_{10}|G_{kj}|$ [K/(kg/s)]")
save(fig, "fig_structure.pdf")
