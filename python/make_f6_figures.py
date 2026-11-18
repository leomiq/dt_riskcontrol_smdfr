"""Figures of the F6 paper campaign (plan §10.2; results/f6_paper/).

  fig_f6_regimes   risk vs pre-cooling for the three failure accelerations μ× (40 kg/s):
                   the risk layer only works in the intermediate regime
  fig_f6_frontier  (a) every risk formulation on one frontier (50 realisations, μ×1e4);
                   (b) partial failures: N2 risk = p_0 × P(total) for the four severity vectors
Style as make_f5b_figures.py (reference palette of the dataviz skill, light mode):
N2 grey, risk bounds on an ordinal blue ramp (α), PI = slot 2 orange, PID = slot 3 aqua.
Usage: python python/make_f6_figures.py
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
RES = os.path.join(ROOT, "results", "f6_paper")
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
A_RAMP = {0.05: "#6da7ec", 0.02: "#2a78d6", 0.01: "#104281"}      # ordinal blue (steps 300/450/650)
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
    return list(csv.DictReader(open(os.path.join(RES, name + ".csv"))))


def stat(R, k, sign=1.0):
    v = sign * np.array([float(r[k]) for r in R])
    v = v[np.isfinite(v)]
    return v.mean(), 1.96 * v.std(ddof=1) / np.sqrt(len(v))


# ---- Fig 1: three regimes of the failure acceleration ------------------------------------------
S = [r for r in rows("sens") if r["act"] == "cap40"]
fig, ax = plt.subplots(1, 3, figsize=(W, 2.55), layout="constrained", sharey=True)
titles = {"3e3": "Fallas largas: media ≈ 150 s\n(μ×3·10³)", "1e4": "Intermedias: media ≈ 46 s\n(μ×10⁴, supuesto)",
          "3e4": "Cortas: media ≈ 15 s\n(μ×3·10⁴)"}
for a, mu in zip(ax, ("3e3", "1e4", "3e4")):
    k = f"p_cond_mu{mu}"

    def pt(case, col, mk, lab=None, hollow=False, ms=5):
        R = [r for r in S if r["case"] == case]
        x, xe = stat(R, "mean_dev", -1); y, ye = stat(R, k)
        a.errorbar(x, y, xerr=xe, yerr=ye, color=col, marker=mk, ms=ms, ls="none", lw=1, capsize=0,
                   mfc="white" if hollow else col, mew=1.3, label=lab)
        return x, y
    pt("N2", MUTED, "D", "N2+sesgo")
    for al in (0.05, 0.02, 0.01):
        x, y = pt(f"N2_var_a{al}_mu{mu}", A_RAMP[al], "o", rf"cota $\alpha$={al} (μ correcto)")
        a.axhline(al, color=A_RAMP[al], lw=0.8, ls=(0, (2, 2)))
        if mu != "1e4":
            pt(f"N2_var_a{al}_mu1e4", A_RAMP[al], "o", rf"cota $\alpha$={al} (diseñada con μ×10⁴)", hollow=True)
    pt("PI", ORANGE, "s", "PI del DET")
    pt("PID", AQUA, "^", "PID de la tesis")
    a.set_yscale("log"); a.set_ylim(1e-3, 0.8); a.set_xlim(-14, 26)
    a.set(title=titles[mu], xlabel=r"K bajo $y_{\rm ref}$ (pre-enfriamiento)")
ax[0].set_ylabel("P(llegar a 1220 K | falla)")
ax[0].text(15, 0.075, "cotas inalcanzables:\nACS saturado", ha="center", va="center", fontsize=6.5, color=INK2)
ax[2].text(6, 0.0024, "piso: componente lenta", ha="center", va="bottom", fontsize=6.5, color=INK2)
ax[1].text(25, 0.05, r"metas $\alpha$ (trazos)", ha="right", va="bottom", fontsize=6.5, color=INK2)
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc="outside lower center", ncol=3, fontsize=6.5)
save(fig, "fig_f6_regimes")

# ---- Fig 2: (a) one frontier, 50 realisations; (b) partial failures = p0 × total ---------------
F = rows("form") + [r for r in S if r["case"] in ("N2", "PI", "PID") or r["case"].endswith("_mu1e4")]
fam = [("N2", ["N2"], MUTED, "D", "N2+sesgo"),
       ("cuantil", ["N2_var_a0.1", "N2_var_a0.05_mu1e4", "N2_var_a0.02_mu1e4", "N2_var_a0.01_mu1e4"], BLUE, "o",
        r"cuantil $\alpha$ = 0.1 / 0.05 / 0.02 / 0.01"),
       ("cvar", ["N2_cvar_b0.9", "N2_cvar_b0.95_mu1e4", "N2_cvar_b0.98"], "#104281", "s", r"CVaR$_\beta$, β = 0.90 / 0.95 / 0.98"),
       ("wilks", ["N2_wilks59"], ORANGE, "*", "PSM de Wilks (59)"),
       ("exp", ["N2_exp_r10", "N2_exp_r100", "N2_exp_r1000"], AQUA, "^", "pérdida esperada, ρ = 10 / 100 / 1000")]
fig, ax = plt.subplots(1, 2, figsize=(W, 2.6), layout="constrained", gridspec_kw={"width_ratios": [1.45, 1]})
for _, cases, col, mk, lab in fam:
    pts = []
    for c in cases:
        R = [r for r in F if r["case"] == c]
        x, xe = stat(R, "mean_dev", -1); y, ye = stat(R, "p_cond_mu1e4")
        pts.append((x, y))
        ax[0].errorbar(x, 1e2 * y, xerr=xe, yerr=1e2 * ye, color=col, marker=mk, ms=7 if mk == "*" else 5, ls="none",
                       lw=1, mew=1.2, capsize=0)
    pts.sort()
    ax[0].plot(*zip(*[(x, 1e2 * y) for x, y in pts]), color=col, lw=1, label=lab, marker=mk, ms=7 if mk == "*" else 4.5, ls="-")
ax[0].set(title="Una sola frontera (50 realizaciones, 40 kg/s)", xlabel=r"K bajo $y_{\rm ref}$ (pre-enfriamiento)",
          ylabel="P(llegar a 1220 K | falla) [%]")
ax[0].legend(loc="upper right", fontsize=6.3)
E = [r for r in rows("eta") if r["case"] == "N2"]
p0 = {"base": 0.1, "uniform": 0.25, "twotrain": 0.3, "pess": 0.4}
lab = {"base": "base", "uniform": "uniforme", "twotrain": "dos trenes", "pess": "pesimista"}
ptot, _ = stat(E, "p_total_only")
xx = np.linspace(0, 0.45, 10)
ax[1].plot(xx, 1e2 * ptot * xx, color=MUTED, lw=1, ls=(0, (3, 2)), label=r"$p_0\times$P(cruce | falla total)")
for k, v in p0.items():
    y, ye = stat(E, f"p_cond_{k}")
    ax[1].errorbar(v, 1e2 * y, yerr=1e2 * ye, color=BLUE, marker="o", ms=5, ls="none", lw=1, capsize=0)
    ax[1].annotate(lab[k], (v, 1e2 * y), xytext=(-4, 6), textcoords="offset points", ha="right", fontsize=6.5, color=INK2)
ax[1].set(title="Fallas parciales: sólo cuenta $p_0$", xlabel=r"fracción de fallas totales $p_0$",
          ylabel="P(cruce | falla) de N2 [%]", xlim=(0, 0.45), ylim=(0, 1.5))
ax[1].legend(loc="upper left", fontsize=6.5)
save(fig, "fig_f6_frontier")
