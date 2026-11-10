"""Didactic figures of docs/estado_mpc_monolitico.tex (introductory chapters).

Uses only stored data:
  models/LJ_test_sampled_vel_ACS_2_loop_1180.mat  archived accident transient (HRM)
  results/didactic/*.csv                          from julia/didactic_analysis.jl
plus two purely illustrative toy examples (stiff ODE, local vs global optimum), which are
labelled as such in the captions.

Usage: .venv\\Scripts\\python python/make_didactic_figures.py   -> docs/figures/did_*.pdf
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
DID = os.path.join(ROOT, "results", "didactic")
os.makedirs(OUT, exist_ok=True)
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))
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
W = 6.3


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def save(fig, name):
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)
    print("wrote", name)


# --- 1. the accident in time (archived HRM transient)
a = scipy.io.loadmat(os.path.join(ROOT, "models", "LJ_test_sampled_vel_ACS_2_loop_1180.mat"),
                     variable_names=["X_1_trajectory", "tw_max_V2P_trajectory", "target"])
leak = a["X_1_trajectory"].ravel()[1:2322]            # 0..1160 s every 0.5 s
t = np.arange(len(leak)) * 0.5
tw = a["tw_max_V2P_trajectory"]
m = tw[:, 0] > 0
tp, pct = tw[m, 0], tw[m, 1]
# one isolated recording artefact (t = 1160.5 s, ~35 K below its neighbours, at the HRM restart
# of the first MPC query) is removed: keep samples within 5 K of the median of their neighbours
med = np.array([np.median(pct[max(0, i - 3):i + 4]) for i in range(len(pct))])
keep = (np.abs(pct - med) < 5.0) & (tp >= 5.0)        # also skip the first 5 s (simulation start-up)
print("PCT samples removed as recording artefacts:", tp[~keep])
tp, pct = tp[keep], pct[keep]
dM = np.cumsum(leak) * 0.5                            # inventory lost (no ACS before 1160 s)
yref = float(a["target"].item())
fig, ax = plt.subplots(3, 1, figsize=(W, 5.0), sharex=True, gridspec_kw={"height_ratios": [1.3, 1, 1]})
ax[0].plot(tp, pct, color=BLUE)
for yv, lab, ls in ((1244, "límite de seguridad $U_y$ = 1244 K", (0, (4, 3))),
                    (1180, "umbral de demanda del ACS = 1180 K", (0, (1, 1.5))),
                    (yref, "", None)):
    if ls:
        ax[0].axhline(yv, color=MUTED, lw=1, ls=ls)
        ax[0].text(20, yv + 1.2, lab, color=INK2, fontsize=7.5)
ax[0].set(ylabel="PCT [K]", ylim=(1135, 1250), title="Temperatura pico de la vaina (HRM)")
ax[1].plot(t, leak, color=ORANGE)
ax[1].set(ylabel=r"$\dot m_{\rm leak}$ [kg/s]", title="Caudal de la fuga")
ax[2].plot(t, -dM, color=AQUA)
ax[2].set(ylabel="kg", xlabel="t [s]", title="Inventario de plomo perdido (fuga acumulada)")
for axi in ax:
    axi.axvline(500, color=MUTED, lw=0.8, ls=":")
    axi.axvspan(1160, 1180, color=BLUE, alpha=0.08, lw=0)
ax[0].text(510, 1230, "ruptura (t = 500 s)", color=INK2, fontsize=7.5)
ax[0].text(1000, 1215, "ventana del MPC\n(1160–1190 s)", color=INK2, fontsize=7.5, ha="right")
save(fig, "did_transient.pdf")

# --- 2. ACS pump characteristic
def pump_mdot(I):
    """ECM of the ACS pump, as matlab/CurrentToMassFlow.m and julia/SMDFRROM.jl:
    ΔP = C I − B Q − A Q², solved for the flow Q ≥ 0; ṁ = ρ Q."""
    rho, Hd, Wd, L, d_, Bm, dP, f_, q0lm, q0w, Kk = 10449.0, 0.515, 0.515, 1.0, 0.001, 1.0, 15000.0, 0.02, 0.118, 1.2, 1.0
    Sx = Wd * Hd
    Rlm = q0lm * Wd / (Hd * L); Rw = q0w * Wd / (2 * d_ * L); Rb = q0lm * Kk**2 * Hd / d_
    Rver = Rw * Rb / (Rw + Rb); Rtot = Rlm + Rver
    C = Bm * Rver * Hd**2 / (Rtot * d_); B = Bm**2 * Rtot * Hd**2 / d_
    A = f_ * rho * L * (Wd + Hd) / (4 * Sx**3)
    Q = (-B + np.sqrt(B**2 - 4 * A * (dP - C * np.asarray(I)))) / (2 * A)
    return rho * np.maximum(Q, 0.0)


d = read_csv(os.path.join(DID, "pump_curve.csv"))            # Julia's map, used as a check
assert np.allclose(pump_mdot(d["current_A"]), d["mdot_kgs"], rtol=1e-9, atol=1e-9)
Ifull = np.linspace(56.6, 87.5, 600)
Idet = np.linspace(56.66, 56.95, 600)
fig, ax = plt.subplots(1, 2, figsize=(W, 2.5))
ax[0].plot(Ifull, pump_mdot(Ifull), color=BLUE)
ax[0].axvspan(56.69, 56.75, color=ORANGE, alpha=0.35, lw=0)
ax[0].set(xlabel="corriente $u$ [A]", ylabel=r"caudal $\dot m$ [kg/s]", title="Curva completa del modelo ECM")
ax[0].annotate("ventana archivada\n56.69–56.75 A\n(0.6–12.9 kg/s)", xy=(56.8, 60), xytext=(63, 2600),
               color=INK2, fontsize=7.5, arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
ax[0].annotate("caudal nominal del\ncircuito: 5550 kg/s", xy=(86.9, 5550), xytext=(68, 4700), color=INK2,
               fontsize=7.5, arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
ax[1].plot(Idet, pump_mdot(Idet), color=BLUE)
ax[1].axvspan(56.69, 56.75, color=ORANGE, alpha=0.35, lw=0)
ax[1].axvline(56.687, color=MUTED, lw=0.8, ls=":")
ax[1].text(56.76, 4, "arranque: 56.687 A\npendiente ≈ 205 kg/s por A", color=INK2, fontsize=7.5)
ax[1].set(xlabel="corriente $u$ [A]", title="Detalle: el pie de la curva")
save(fig, "did_pump.pdf")

# --- 3. stiffness: eigenvalues of the ROM Jacobian
fig, ax = plt.subplots(figsize=(W, 1.9))
for tq, col, yy in (("1180", BLUE, 1.0), ("1161", ORANGE, 0.0)):
    e = read_csv(os.path.join(DID, f"eigenvalues_{tq}.csv"))
    re = -e["re"]
    nz = re > 1e-9                                   # λ = 0 is the inventory integrator
    ax.scatter(re[nz], np.full(nz.sum(), yy), color=col, s=22, zorder=3, label=f"$t_q$ = {tq} s")
    ax.text(4e-4, yy, r"y $\lambda=0$" + "\n(inventario)", color=INK2, fontsize=7, va="center")
ax.set_xscale("log")
ax.set_xlim(3e-4, 5e3)
ax.set(yticks=[0, 1], yticklabels=["1161 s", "1180 s"], ylim=(-0.6, 1.6),
       xlabel=r"$-\mathrm{Re}\,\lambda$ [1/s]  (velocidad de decaimiento de cada modo)")
ax.axvline(2 / 0.5, color=MUTED, lw=1, ls=(0, (4, 3)))
ax.text(2 / 0.5 * 1.15, 1.45, r"$2/\Delta t$ = 4 s$^{-1}$: a la derecha, Euler explícito con $\Delta t$ = 0.5 s es inestable",
        color=INK2, fontsize=7.5)
ax.grid(axis="y", visible=False)
save(fig, "did_eigs.pdf")

# --- 4. stiff toy problem: explicit vs implicit Euler (illustrative)
lam, h, T = 50.0, 0.05, 1.0
f = lambda tt, y: -lam * (y - np.cos(tt))
ts = np.arange(0, T + 1e-9, h)
ye, yi = [0.0], [0.0]                                    # start off the slow solution y ≈ cos t
for n in range(len(ts) - 1):
    ye.append(ye[-1] + h * f(ts[n], ye[-1]))
    yi.append((yi[-1] + h * lam * np.cos(ts[n + 1])) / (1 + h * lam))
tf = np.linspace(0, T, 400)
fig, ax = plt.subplots(figsize=(W, 2.2))
yex = (lam**2 * np.cos(tf) + lam * np.sin(tf)) / (lam**2 + 1) - lam**2 / (lam**2 + 1) * np.exp(-lam * tf)
ax.plot(tf, yex, color=MUTED, lw=1, label="solución exacta")
ax.plot(ts, yi, color=BLUE, marker="o", ms=3, label="Euler implícito, $h\\lambda$ = 2.5")
ax.plot(ts, ye, color=ORANGE, marker="s", ms=3, label="Euler explícito, $h\\lambda$ = 2.5")
ax.set(ylim=(-3, 3), xlabel="t", ylabel="y",
       title=r"Ejemplo ilustrativo: $\dot y=-50\,(y-\cos t)$, $y(0)=0$, paso $h$ = 0.05")
ax.legend(loc="lower left", ncol=3, fontsize=7.5)
save(fig, "did_stiff.pdf")

# --- 5. frequency response of the linearized plant (integrator)
fig, ax = plt.subplots(1, 2, figsize=(W, 2.3))
for tq, col in (("1180", BLUE), ("1161", ORANGE)):
    b = read_csv(os.path.join(DID, f"bode_{tq}.csv"))
    ax[0].loglog(b["omega"], b["mag"], color=col, label=f"$t_q$ = {tq} s")
    ax[1].semilogx(b["omega"], b["phase_deg"], color=col)
w = np.logspace(-4, -1, 10)
ax[0].loglog(w, 1.6e-2 / w, color=MUTED, ls=(0, (4, 3)), lw=1)     # |G_M / (iω)|, G_M = ∂ŷ/∂ΔM
ax[0].text(2e-4, 1.6e-2 / 2e-4 * 0.35, r"pendiente $-1$: integrador", color=INK2, fontsize=7.5)
ax[0].set(xlabel=r"$\omega$ [rad/s]", ylabel=r"$|G(i\omega)|$ [K/(kg/s)]", title="Ganancia")
ax[0].text(3, 6e-4, r"meseta: término directo $D=p_4'$", color=INK2, fontsize=7.5)
ax[1].set(xlabel=r"$\omega$ [rad/s]", ylabel="fase [°]", title="Fase")
ax[0].legend(loc="lower left")
save(fig, "did_bode.pdf")

# --- 6. actuator authority over the horizon (nonlinear replica)
fig, ax = plt.subplots(1, 2, figsize=(W, 2.4), sharey=False)
for j, tq in enumerate(("1161", "1180")):
    d = read_csv(os.path.join(DID, f"authority_{tq}.csv"))
    keys = [k for k in d if k != "t"]
    cols = plt.cm.Blues(np.linspace(0.35, 0.95, len(keys)))
    for k_, c_ in zip(keys, cols):
        ax[j].plot(d["t"], d[k_], color=c_, label=f"{float(k_[1:]):.1f} kg/s")
    ax[j].axhline(yref, color=MUTED, lw=1, ls=(0, (4, 3)))
    ax[j].text(d["t"][0], yref + 0.05, r"$y_{\rm ref}$", color=INK2, fontsize=7.5)
    ax[j].set(title=f"$t_q$ = {tq} s", xlabel="t [s]", ylabel="PCT estimada [K]" if j == 0 else None)
    ax[j].ticklabel_format(axis="y", useOffset=False)
h_, l_ = ax[1].get_legend_handles_labels()
fig.legend(h_, l_, title="caudal constante del ACS durante el horizonte", fontsize=7.5, title_fontsize=7.5,
           loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.20))
save(fig, "did_authority.pdf")

# --- 7. local vs global optimum and the optimality gap (illustrative)
x = np.linspace(-2, 2.2, 600)
J = 0.6 * x**4 - 0.4 * x**3 - 1.5 * x**2 + 0.3 * x + 2
xl = x[(x < 0)][np.argmin(J[(x < 0)])]                   # local (worse) minimum, left basin
xg = x[np.argmin(J)]                                     # global minimum, right basin


def lower_convex_hull(xs, ys):
    """Lower convex envelope of the points (monotone chain), evaluated on xs."""
    hull = []
    for p in zip(xs, ys):
        while len(hull) >= 2 and ((hull[-1][0] - hull[-2][0]) * (p[1] - hull[-2][1])
                                  - (hull[-1][1] - hull[-2][1]) * (p[0] - hull[-2][0])) <= 0:
            hull.pop()
        hull.append(p)
    hx, hy = zip(*hull)
    return np.interp(xs, hx, hy)


L = lower_convex_hull(x, J) - 0.25                       # convex and below J on the whole interval
assert np.all(L <= J + 1e-12)
fig, ax = plt.subplots(figsize=(W, 2.4))
ax.plot(x, J, color=BLUE, label="costo $J(u)$ (no convexo)")
ax.plot(x, L, color=ORANGE, ls=(0, (4, 3)), label="relajación convexa $L(u)\\leq J(u)$")
ax.scatter([xl], [J[x == xl]], color=INK, zorder=3, s=18)
ax.scatter([xg], [J.min()], color=AQUA, zorder=3, s=22)
ax.axhline(J[x == xl][0], color=MUTED, lw=0.8, ls=":")
ax.axhline(L.min(), color=MUTED, lw=0.8, ls=":")
ax.text(-1.70, J[x == xl][0] + 0.10, "UB", color=INK2, fontsize=7.5)
ax.text(-1.70, L.min() - 0.42, "LB", color=INK2, fontsize=7.5)
ax.annotate("óptimo local", xy=(xl, J[x == xl][0]), xytext=(-0.55, 3.2), color=INK2, fontsize=7.5,
            arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
ax.text(xg + 0.08, J.min() + 0.25, "óptimo global", color=INK2, fontsize=7.5)
ax.annotate("", xy=(-1.78, L.min()), xytext=(-1.78, J[x == xl][0]),
            arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.8))
ax.text(-1.74, (L.min() + J[x == xl][0]) / 2, "brecha", color=INK2, fontsize=7.5, va="center")
ax.set(xlabel="decisión $u$", ylabel="costo", title="Ejemplo ilustrativo en una dimensión")
ax.legend(loc="upper center")
save(fig, "did_gap.pdf")
