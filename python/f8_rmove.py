"""F8: move suppression R Σ (ṁ_k − ṁ_{k−1})² in the closed loop (results/f8_rmove/, Olab).

Campaigns (julia/F6Tasks.jl, `rrisk` and `rmove`; R grid 0, 1e-3 … 100 [K²/(kg/s)²]):
  rrisk  failure-free loop + Rao–Blackwell P(cross | failure), N2 and N2 + VaR α = 0.02 (μ×1e4),
         50 realisations of the det family, 40 kg/s
  rmove  closed loop with DET-accelerated failures (det family, N2 for every R, N3 for R ∈
         {0, 0.1, 10}) and the nominal family (N2), 50 realisations each, 40 kg/s

Writes results/f8_rmove/summary.md and the figures
  fig_f8_tradeoff  IAE vs total variation of the commanded flow (det and nominal), one point per R
  fig_f8_risk      P(cross | failure) and total variation of the flow vs R, without and with the
                   risk layer (R ≥ 0.03 removes the period-2 chattering of the risk layer)
Style as make_f6_figures.py. Usage: python python/f8_rmove.py
"""
import os, csv, gzip
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
RES = os.path.join(ROOT, "results", "f8_rmove")
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df"
GREY = "#6f6e69"
plt.rcParams.update({
    "font.family": "Inter Variable", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 7, "lines.linewidth": 1.5,
    "figure.dpi": 150, "savefig.bbox": "tight", "pdf.fonttype": 42,
})
W = 6.3
R0X = 3e-4            # where R = 0 is drawn on the log axes


def rows(name):
    p = os.path.join(RES, name + ".csv")
    return list(csv.DictReader(open(p))) if os.path.isfile(p) else []


def f(r, k):
    try:
        return float(r[k])
    except (KeyError, ValueError):
        return float("nan")


def base(case):
    return case.rsplit("_R", 1)[0]


def ci(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 2:
        return (v.mean() if len(v) else np.nan), np.nan, len(v)
    return v.mean(), 1.96 * v.std(ddof=1) / np.sqrt(len(v)), len(v)


def group(Rows, fam=None):
    """{(base case, R): {seed: row}}"""
    G = defaultdict(dict)
    for r in Rows:
        if fam and r["family"] != fam:
            continue
        G[(base(r["case"]), f(r, "R"))][r["seed"]] = r
    return G


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)
    print("wrote", name)


def xr(R):
    return R if R > 0 else R0X


rr, rm = rows("rrisk"), rows("rmove")
RS = sorted({f(r, "R") for r in rr + rm})
L = ["# F8 — move suppression R Σ(ṁ_k − ṁ_{k−1})² (Olab)\n",
     f"R grid [K²/(kg/s)²]: {', '.join(f'{x:g}' for x in RS)}. 40 kg/s, ramp 5 kg/s per step (hard). "
     "Means over realisations ± 95 % CI half-width; paired differences against R = 0 on the same seeds.\n"]

# ---- 1. exactness of the soft margin: largest predicted PCT over every query --------------------
ym = np.array([f(r, "yhat_max") for r in rr + rm]); ym = ym[np.isfinite(ym)]
pm = np.array([f(r, "pct_max") for r in rr + rm]); pm = pm[np.isfinite(pm)]
L.append("## 1. Is the soft margin λ_pct[ŷ − 1244]₊² ever active?\n")
L.append(f"- runs: {len(ym)}; largest predicted PCT over all queries (ROM, applied plan, no bias "
         f"correction): **{ym.max():.2f} K**; largest plant PCT: {pm.max():.2f} K; U_y = 1244 K → "
         f"penalty active in {int((ym > 1244).sum())} runs.\n")
r0 = [f(r, "yhat_max") for r in rr + rm if f(r, "R") == 0]
r0 = np.array([x for x in r0 if np.isfinite(x)])
if len(r0):
    L.append(f"- R = 0 runs (the F5–F7 formulation): {len(r0)}, max predicted PCT {r0.max():.2f} K.\n")

# ---- 2. closed loop with failures and nominal family: tracking vs smoothness -------------------
MET = [("iae", "IAE [K·s]"), ("total_variation", "TV ṁ [kg/s]"), ("mass_cmd", "mass [kg]"),
       ("pct_max", "PCT max [K]"), ("min_overshoot", "min(PCT − y_ref) [K]")]
TR = {}
for fam, title in (("detacc", "det family (DET-accelerated failures)"), ("nominal", "nominal family")):
    G = group(rm, fam)
    if not G:
        continue
    L.append(f"\n## 2{'a' if fam == 'detacc' else 'b'}. Closed loop, {title}: N2 + bias\n")
    L.append("| R | n | " + " | ".join(m[1] for m in MET) + " | limit hits | ΔIAE vs R=0 | ΔTV vs R=0 |")
    L.append("|---|---|" + "---|" * (len(MET) + 3))
    ref = G.get(("N2b", 0.0), {})
    pts = []
    for R in RS:
        g = G.get(("N2b", R))
        if not g:
            continue
        cells = []
        for k, _ in MET:
            m, h, n = ci([f(r, k) for r in g.values()])
            cells.append(f"{m:.4g} ± {h:.2g}")
        lim = sum(r["limit"] in ("true", "1") for r in g.values())
        sd = [s for s in g if s in ref]
        dI = ci([f(g[s], "iae") - f(ref[s], "iae") for s in sd])
        dT = ci([f(g[s], "total_variation") - f(ref[s], "total_variation") for s in sd])
        L.append(f"| {R:g} | {len(g)} | " + " | ".join(cells) +
                 f" | {lim} | {dI[0]:+.3g} ± {dI[1]:.2g} | {dT[0]:+.3g} ± {dT[1]:.2g} |")
        pts.append((R, ci([f(r, "total_variation") for r in g.values()])[0],
                    ci([f(r, "iae") for r in g.values()])[0],
                    ci([f(r, "mass_cmd") for r in g.values()])[0]))
    TR[fam] = pts

# N3 vs N2 at the same R (det family)
G = group(rm, "detacc")
traj = defaultdict(dict)
tp = os.path.join(RES, "rmove_traj.csv")
tz = tp + ".gz"
if os.path.isfile(tp) or os.path.isfile(tz):
    fh = gzip.open(tz, "rt") if (os.path.isfile(tz) and not os.path.isfile(tp)) else open(tp)
    for r in csv.DictReader(fh):
        traj[r["key"]][r["t"]] = (float(r["pct"]), float(r["m"]))
L.append("\n## 2c. N3 (exact NLP) vs N2 (sequential QP) at the same R, det family\n")
L.append("| R | pairs | max abs ΔIAE [K·s] | max rel ΔIAE | max abs ΔPCT along the loop [K] | max abs Δṁ [kg/s] |")
L.append("|---|---|---|---|---|---|")
for R in RS:
    a, b = G.get(("N2b", R)), G.get(("N3b", R))
    if not a or not b:
        continue
    sd = [s for s in a if s in b]
    dI = [abs(f(a[s], "iae") - f(b[s], "iae")) for s in sd]
    rI = [abs(f(a[s], "iae") - f(b[s], "iae")) / max(f(a[s], "iae"), 1e-9) for s in sd]
    dP, dM = [], []
    for s in sd:
        ta, tb = traj.get(a[s]["key"], {}), traj.get(b[s]["key"], {})
        for t in ta.keys() & tb.keys():
            dP.append(abs(ta[t][0] - tb[t][0])); dM.append(abs(ta[t][1] - tb[t][1]))
    L.append(f"| {R:g} | {len(sd)} | {max(dI):.3g} | {max(rI):.2g} | "
             f"{(max(dP) if dP else float('nan')):.3g} | {(max(dM) if dM else float('nan')):.3g} |")

# ---- 3. risk: P(cross | failure) vs R, without and with the risk layer -------------------------
G = group(rr)
RK = {}
if G:
    L.append("\n## 3. Risk (failure-free loop + Rao–Blackwell, det family): P(cross | failure)\n")
    L.append("| case | R | n | P μ×3e3 | P μ×1e4 | P μ×3e4 | pre-cooling [K] | mass [kg] | IAE [K·s] | TV ṁ |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for c in ("N2", "N2_var_a0.02"):
        pts = []
        for R in RS:
            g = G.get((c, R))
            if not g:
                continue
            P = [ci([f(r, f"p_cond_mu{mu}") for r in g.values()]) for mu in ("3e3", "1e4", "3e4")]
            pc = ci([-f(r, "mean_dev") for r in g.values()])
            ms = ci([f(r, "mass_cmd") for r in g.values()])
            ia = ci([f(r, "iae") for r in g.values()])
            tv = ci([f(r, "total_variation") for r in g.values()])
            L.append(f"| {c} | {R:g} | {len(g)} | " + " | ".join(f"{p[0]:.4g} ± {p[1]:.2g}" for p in P) +
                     f" | {pc[0]:.3g} | {ms[0]:.4g} | {ia[0]:.4g} | {tv[0]:.4g} |")
            pts.append((R, P[1], ms[0], pc[0], tv[0]))
        RK[c] = pts

open(os.path.join(RES, "summary.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))

# ---- figures ------------------------------------------------------------------------------------
if TR:
    fig, ax = plt.subplots(1, 2, figsize=(W, 2.6), layout="constrained")
    for a, fam, title in ((ax[0], "detacc", "Con fallas del ACS (familia det)"),
                          (ax[1], "nominal", "Familia nominal")):
        pts = TR.get(fam, [])
        if not pts:
            continue
        x = [p[1] for p in pts]; y = [p[2] for p in pts]
        a.plot(x, y, "-", color=BLUE, lw=1.2, zorder=2)
        a.scatter(x, y, s=16, color=BLUE, zorder=3)
        for (R, tv, ia, _) in pts:
            if R in (0, 0.01, 0.1, 1, 10, 100):
                a.annotate("R = 0" if R == 0 else f"{R:g}", (tv, ia), textcoords="offset points",
                           xytext=(4, 3), fontsize=6.5, color=INK2)
        a.set_yscale("log")
        a.set_xlabel("Variación total del caudal comandado [kg/s]")
        a.set_title(title)
    ax[0].set_ylabel("IAE de la PCT [K·s]")
    save(fig, "fig_f8_tradeoff")

if RK:
    fig, ax = plt.subplots(1, 2, figsize=(W, 2.6), layout="constrained")
    sty = {"N2": (GREY, "N2 (sin capa de riesgo)"), "N2_var_a0.02": (BLUE, "N2 + cota de riesgo, α = 0.02")}
    for c, pts in RK.items():
        col, lab = sty[c]
        x = [xr(p[0]) for p in pts]
        m = np.array([p[1][0] for p in pts]); h = np.array([p[1][1] for p in pts])
        ax[0].errorbar(x, m, yerr=h, fmt="o-", ms=3, color=col, label=lab, capsize=2, lw=1.2)
        ax[1].plot(x, [p[4] for p in pts], "o-", ms=3, color=col, label=lab, lw=1.2)
    for a in ax:
        a.set_xscale("log")
        ticks = [R0X, 1e-3, 1e-2, 1e-1, 1, 10, 100]
        a.set_xticks(ticks, ["0", "10⁻³", "10⁻²", "10⁻¹", "1", "10", "10²"])
        a.set_xlabel("Peso de supresión de movimiento R [K²/(kg/s)²]")
    ax[0].set_ylabel("P(cruzar 1220 K | falla), μ×10⁴")
    ax[1].set_ylabel("Variación total del caudal [kg/s]")
    ax[1].set_yscale("log")
    ax[0].legend(loc="upper left")
    ax[0].set_title("Riesgo")
    ax[1].set_title("Oscilación del caudal (llega a la planta)")
    save(fig, "fig_f8_risk")
