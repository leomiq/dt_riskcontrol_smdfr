"""Figures of F9 (results/f9_rozas/): the one-scenario result of Homem-de-Mello et al. (2026).

  fig_f9_theorem   ECDF of the decision gap |u − u_SAA|∞ per query: the Theorem-1 scenario ξ*
                   (exact), the best PHYSICAL one-scenario (one failure duration, bilevel) and the
                   mean-duration point forecast; (a) ρ = 100, (b) ρ = 1000
  fig_f9_learning  (a) CPU vs sample size N (SAA QP vs scenario search, paper Table 3);
                   (b) contextual forecast Φ(θ, x): share of the VSS captured vs training size
  fig_f9_closed    closed loop: risk vs pre-cooling for N2, SAA and the one-scenario MPCs
Style as make_f6_figures.py (reference palette of the dataviz skill, light mode).
Usage: python python/make_f9_figures.py
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
RES = os.path.join(ROOT, "results", "f9_rozas")
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"   # categorical slots 1–4
RAMP = {"100": "#6da7ec", "1000": "#104281"}                         # ordinal blue: ρ
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df"
plt.rcParams.update({
    "font.family": "Inter Variable", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 7, "lines.linewidth": 1.5,
    "figure.dpi": 150, "savefig.bbox": "tight", "pdf.fonttype": 42,
})
W = 6.3


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)
    print("wrote", name)


def rows(name):
    p = os.path.join(RES, name + ".csv")
    return list(csv.DictReader(open(p, encoding="utf-8"))) if os.path.isfile(p) else []


def col(R, k):
    v = np.array([float(r.get(k, "nan") or "nan") for r in R])
    return v[np.isfinite(v)]


def ecdf(ax, v, color, label, lw=1.6, ls="-"):
    v = np.sort(np.maximum(v, 1e-9))
    ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=color, lw=lw, ls=ls, label=label)


A = rows("A_queries")
# ---- Fig 1: Theorem 1 vs physical one-scenario -----------------------------------------------
if A:
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.6), layout="constrained", sharey=True)
    for ax, ρ in zip(axs, ("100", "1000")):
        ecdf(ax, col(A, f"r{ρ}_T1_gap"), BLUE, "escenario ξ* del teorema 1\n(un margen por paso)")
        ecdf(ax, col(A, f"r{ρ}_BL_gap"), ORANGE, "mejor falla única $\hat{D}^*$ (bilevel):\nplan completo")
        ecdf(ax, col(A, f"r{ρ}_BL_gap1"), ORANGE, "mejor falla única $\hat{D}^*$:\nacción aplicada u₁", ls="--", lw=1.2)
        ecdf(ax, col(A, f"r{ρ}_mean_gap"), MUTED, "duración media E[D]\n(pronóstico puntual)")
        ax.set_xscale("log"); ax.set_xlim(1e-9, 1e2)
        ax.axvline(1e-2, color=INK2, lw=0.6, ls=":")
        ax.set_title(f"({'a' if ρ == '100' else 'b'}) pérdida esperada, ρ = {ρ}", loc="left")
        ax.set_xlabel("|u − u_SAA|∞ [kg/s] (log)")
    axs[0].set_ylabel("fracción de consultas")
    axs[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    save(fig, "fig_f9_theorem")

# ---- Fig 2: CPU vs N and learning curve -----------------------------------------------------
S = [r for r in A if r.get("cpu_saa_N10", "nan") not in ("", "nan", "NaN")]
Lr = rows("B_learning")
if S:
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.5), layout="constrained")
    Ns = np.array([10, 40, 100, 400, 1000, 4000])
    for k, c, lab in (("saa", BLUE, "SAA: un QP con N escenarios"), ("bl", ORANGE, "búsqueda del escenario único\n(301 QP de un escenario + evaluar N)")):
        med = [np.median(col(S, f"cpu_{k}_N{N}")) for N in Ns]
        lo = [np.percentile(col(S, f"cpu_{k}_N{N}"), 10) for N in Ns]; hi = [np.percentile(col(S, f"cpu_{k}_N{N}"), 90) for N in Ns]
        axs[0].fill_between(Ns, lo, hi, color=c, alpha=0.15, lw=0)
        axs[0].plot(Ns, med, "o-", color=c, ms=4, label=lab)
    axs[0].set_xscale("log"); axs[0].set_yscale("log")
    axs[0].set_xlabel("muestras N de la duración"); axs[0].set_ylabel("tiempo por consulta [s]")
    axs[0].set_title("(a) costo computacional (tabla 3 del artículo)", loc="left"); axs[0].set_ylim(0.01, 300); axs[0].legend(loc="upper left")
    if Lr:
        for ρ in sorted({r["rho"] for r in Lr}, key=float):
            key = "100" if float(ρ) == 100 else "1000"
            Rr = [r for r in Lr if r["rho"] == ρ]
            Nt = sorted({int(r["N_train"]) for r in Rr})
            m = [col([r for r in Rr if int(r["N_train"]) == N], "VSS_captured") for N in Nt]
            axs[1].plot(Nt, [v.mean() for v in m], "o-", color=RAMP[key], ms=4, label=f"ρ = {key}")
            for N, v in zip(Nt, m):
                axs[1].plot([N] * len(v), v, ".", color=RAMP[key], ms=3, alpha=0.5)
        axs[1].set_xscale("log"); axs[1].axhline(1, color=INK2, lw=0.6, ls=":"); axs[1].set_ylim(-1.0, 1.08)
        axs[1].set_xlabel("consultas de entrenamiento"); axs[1].set_ylabel("fracción del VSS capturada (prueba)")
        axs[1].set_title("(b) pronóstico Φ(θ, x) aprendido", loc="left"); axs[1].legend(loc="lower right")
    save(fig, "fig_f9_learning")

# ---- Fig 3: closed loop ----------------------------------------------------------------------
C = rows("C_rb")
if C:
    fig, ax = plt.subplots(figsize=(W * 0.62, 2.8), layout="constrained")
    style = {"SAA": (BLUE, "SAA (40 escenarios)"), "scen_Dconst": (ORANGE, "un escenario: $\hat{D}$ fijo"),
             "scen_ADlinE": (AQUA, "un escenario: Φ(θ,x) lineal,\ncontexto completo (oscila)"),
             "scen_ADslowE": (YELLOW, "un escenario: Φ(θ,x) lineal,\nsólo variables lentas")}
    def pt(R):
        p = col(R, "p_cond"); m = -col(R, "mean_dev")
        return m.mean(), p.mean(), 1.96 * m.std(ddof=1) / np.sqrt(len(m)), 1.96 * p.std(ddof=1) / np.sqrt(len(p))
    x, y, ex, ey = pt([r for r in C if r["case"] == "N2"])
    ax.errorbar(x, y, xerr=ex, yerr=ey, fmt="o", color=INK2, ms=6, capsize=0, lw=0.8)
    ax.annotate("N2 sin riesgo", (x, y), textcoords="offset points", xytext=(6, 4), fontsize=7, color=INK2)
    for pre, (c, lab) in style.items():
        for ρ, mk in (("100", "o"), ("1000", "s")):
            R = [r for r in C if r["case"] == f"{pre}_r{ρ}"]
            if not R:
                continue
            x, y, ex, ey = pt(R)
            ax.errorbar(x, y, xerr=ex, yerr=ey, fmt=mk, color=c, mfc=c, mec="white", mew=1.0, ms=7, capsize=0, lw=0.8,
                        label=f"{lab}" if ρ == "100" else None)
            if ρ == "1000":                       # 4 series: direct labels as well (relief rule)
                if pre != "scen_ADlinE":          # Φ completo and Φ lento coincide: one label
                    short = {"SAA": "SAA", "scen_Dconst": "D̂ fijo", "scen_ADslowE": "Φ lento ≈ Φ completo"}[pre]
                    off = {"SAA": (7, -6), "scen_Dconst": (7, 3), "scen_ADslowE": (-30, 9)}[pre]
                    ax.annotate(short, (x, y), textcoords="offset points", xytext=off, fontsize=6.5, color=INK2)
    ax.plot([], [], "o", color=MUTED, label="ρ = 100"); ax.plot([], [], "s", color=MUTED, label="ρ = 1000")
    ax.set_xlabel("pre-enfriamiento medio bajo y_ref [K]"); ax.set_ylabel("P(llegar a 1220 K | falla)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))
    save(fig, "fig_f9_closed")
