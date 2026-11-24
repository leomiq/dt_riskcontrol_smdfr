"""Didactic figure of the risk metrics (status report, section «Las métricas de riesgo desde cero»).

  fig_risk_tutorial  (a) one-step toy: heat-up Z = 0.3 K/s x D, D ~ Exp(mean 46 s), so Z ~ Exp(13.8 K);
                         P(Z > m) and the margin m that each risk metric asks for
                     (b) the SMDFR case: survival S(D) of the accelerated DET mixture, quantiles
                         q_{1-a}(D) used by the bound, N2 operating point and the 0.002 floor
                     (c) why Rao-Blackwell: 95 % half-width of the estimate of P(cross | failure)
                         vs number of closed-loop realisations, counting vs Rao-Blackwell
All numbers are closed-form (toy) or come from models/MCDET_conferencepaper_database.mat
(rates x1e5 / x1e4, as julia/F6Tasks.jl). Style as make_f9_figures.py (dataviz palette).
Usage: python python/make_risk_tutorial_figures.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.io import loadmat

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
VIOLET = "#7c5cd6"
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

# ---- toy -------------------------------------------------------------------------------------
THETA = 0.3 * 46.0                  # mean heat-up [K]
P = lambda m: np.exp(-m / THETA)    # P(Z > m)
M0 = 1220.0 - 1179.0                # margin at the set-point [K]


def wilks_margin(level):            # S(D_max) = level  ->  margin
    return THETA * np.log(1.0 / level)


toy = [  # label, margin, colour
    ("N2 determinista", M0, MUTED),
    ("cuantil α=0.05", THETA * np.log(20), BLUE),
    ("pérdida esperada ρ=100", M0 + 2.2, AQUA),
    ("pérdida esperada ρ=1000", M0 + 11.3, AQUA),
    ("cuantil α=0.02", THETA * np.log(50), BLUE),
    ("CVaR$_{0.95}$", THETA * (np.log(20) + 1), ORANGE),
]

# ---- SMDFR mixture -----------------------------------------------------------------------------
db = loadmat(os.path.join(ROOT, "models", "MCDET_conferencepaper_database.mat"))
LAM = np.ravel(db["Lambda_t"]).astype(float) * 1e5
MU = np.ravel(db["Mu_t"]).astype(float) * 1e4
WGT = LAM / LAM.sum()
S = lambda t: np.sum(WGT[:, None] * np.exp(-MU[:, None] * np.atleast_1d(t)[None, :]), axis=0)


def quant(p):
    lo, hi = 0.0, 1e7
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if 1 - S(m)[0] < p else (lo, m)
    return 0.5 * (lo + hi)


def main():
    fig, ax = plt.subplots(1, 3, figsize=(W, 2.9), gridspec_kw=dict(wspace=0.55))

    # (a) toy
    a = ax[0]
    m = np.linspace(0, 95, 400)
    a.semilogy(m, P(m), color=INK2, lw=1.3)
    lo, med, hi = wilks_margin(1 - 0.1 ** (1 / 59)), wilks_margin(1 - 0.5 ** (1 / 59)), wilks_margin(1 - 0.9 ** (1 / 59))
    a.axvspan(lo, hi, color=VIOLET, alpha=0.10, lw=0)
    a.plot([med], [P(med)], "o", color=VIOLET, ms=4)
    a.annotate("Wilks 59", (med, P(med)), xytext=(70, 0.3), fontsize=6.5, color=VIOLET, ha="center")
    seen = set()
    for lab, mm, c in toy:
        key = lab.split(" ")[0] + (" " + lab.split(" ")[1] if lab.startswith("pérdida") else "")
        a.plot([mm], [P(mm)], "o", color=c, ms=3.6, zorder=3)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=c, ms=3.6, label=l) for l, c in
         [("N2 determinista (m=41 K)", MUTED), ("cuantil α=0.05 / 0.02", BLUE),
          ("pérdida esperada ρ=100 / 1000", AQUA), ("CVaR$_{0.95}$", ORANGE)]]
    a.legend(handles=h, loc="lower left", fontsize=5.8, handletextpad=0.2, borderaxespad=0.1)
    a.set_xlabel("margen $m$ = 1220 − PCT [K]")
    a.set_ylabel("P(cruce | falla) = P(Z > m)")
    a.set_title("(a) Juguete", loc="left")
    a.set_ylim(1e-3, 1.2); a.set_xlim(0, 95)

    # (b) SMDFR mixture
    b = ax[1]
    t = np.linspace(0, 450, 900)
    b.semilogy(t, S(t), color=INK2, lw=1.3)
    b.axhline(WGT[-1], color=MUTED, ls=":", lw=1)
    b.text(445, WGT[-1] * 1.25, "piso 0.002", fontsize=6.2, color=MUTED, ha="right")
    for al, c in [(0.05, BLUE), (0.02, BLUE), (0.01, BLUE)]:
        q = quant(1 - al)
        b.plot([q, q], [1e-3, al], color=c, lw=0.8, ls="--")
        b.plot([q], [al], "o", color=c, ms=3.4)
        b.text(q + 5, al * 1.2, f"{q:.0f} s", fontsize=6.2, color=c, ha="left") if al == 0.01 else b.text(q - 5, al * 0.55, f"{q:.0f} s", fontsize=6.2, color=c, ha="right")
    tn = 159.6
    b.plot([tn], [S(tn)[0]], "o", color=ORANGE, ms=4.2, zorder=4)
    b.annotate("N2: τ≈160 s\n→ 0.031",(tn, S(tn)[0]), xytext=(14, 10), textcoords="offset points",
               fontsize=6.2, color=ORANGE)
    b.set_xlabel("duración $D$ [s]")
    b.set_ylabel("S(D) = P(duración > D)")
    b.set_title("(b) Caso SMDFR", loc="left")
    b.set_ylim(1e-3, 1.2); b.set_xlim(0, 450)

    # (c) counting vs Rao-Blackwell
    c = ax[2]
    n = np.logspace(1, 4.3, 200)
    p = 0.031
    hw_count = 1.96 * np.sqrt(p * (1 - p) / n)
    s_rb = 0.004 / 1.96 * np.sqrt(50)        # per-realisation sd implied by +-0.004 with 50 runs
    hw_rb = 1.96 * s_rb / np.sqrt(n)
    c.loglog(n, hw_count, color=MUTED, label="conteo de cruces")
    c.loglog(n, hw_rb, color=BLUE, label="Rao–Blackwell")
    c.axhline(0.1 * p, color=ORANGE, ls="--", lw=1)
    c.text(2000, 0.1 * p * 0.55, "10 % de P", fontsize=6.2, color=ORANGE, ha="center")
    n_c = (1.96 / 0.1) ** 2 * (1 - p) / p
    n_r = (1.96 * s_rb / (0.1 * p)) ** 2
    for nn, col in [(n_c, MUTED), (n_r, BLUE)]:
        c.plot([nn], [0.1 * p], "o", color=col, ms=4)
        c.annotate(f"n≈{nn:,.0f}".replace(",", " "), (nn, 0.1 * p), xytext=(0, 8), textcoords="offset points",
                   fontsize=6.4, color=col, ha="center")
    c.plot([50], [0.004], "s", color=BLUE, ms=3.6)
    c.annotate("N2: ±0.004 (50 real.)", (50, 0.004), xytext=(-6, -34), textcoords="offset points", fontsize=6.2,
               color=BLUE)
    c.set_xlabel("realizaciones")
    c.set_ylabel("semiancho del IC 95 %")
    c.set_title("(c) Estimadores", loc="left")
    c.legend(loc="upper right")
    fig.savefig(os.path.join(OUT, "fig_risk_tutorial.pdf"))
    fig.savefig(os.path.join(OUT, "fig_risk_tutorial.png"))
    print("wrote fig_risk_tutorial;", f"n_count={n_c:.0f} n_rb={n_r:.0f}",
          "q:", {a: round(quant(1 - a), 1) for a in (0.05, 0.02, 0.01)})


if __name__ == "__main__":
    main()
