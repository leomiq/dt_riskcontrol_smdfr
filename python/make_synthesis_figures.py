"""Synthesis figures of the consolidated status report (docs/estado_mpc_monolitico.tex, 18-11-2026).

Condense the integrated results into a few figures meant to be read first (and to seed the
manuscript). Reads only committed results (no MATLAB, no Julia):
  results/e2_optimizers_1180/e2_summary.csv        open loop: LJ / GA wall time
  results/f3_certificates/*_1180.csv               open loop: N1 / N2 / N3 solve times
  results/f3c_optimality/decisions_compared.csv    open loop: excess over the best plan, DE
  results/f6_paper/*.csv(.gz)                      closed loop, 50 realisations per family
  models/MCDET_conferencepaper_database.mat        DET failure / repair rates (μ×1e4 case)

  syn_overview        graphical abstract: solve well / control well / manage risk well
  syn_openloop        time per decision vs cost excess over the best plan (t_q = 1180 s)
  syn_cl_compare      representative closed loops: PCT and delivered flow, three families
  syn_cl_fan          50 closed loops with ACS failures aligned at the failure start (fan charts)
  syn_scorecard       thesis metrics by controller and family (mean ± 95 % CI, all runs)
  syn_mechanism       the margin construction: ΔT(D) does not depend on the controller, and
                      every risk measure on D maps through ΔT to one margin m = 1220 − U
  syn_audit           11 400 audited closed-loop queries: distance to the best plan and solve time

Style as the other report scripts (reference palette of the dataviz skill, light mode).
Colour follows the controller in every synthesis figure:
  MPC N1 / N2 / N3 = ordinal blue ramp (fidelity ladder), N2 + risk bound = navy,
  PI of the DET = orange (slot 2), thesis PID = aqua (slot 3), metaheuristics = grey.
Usage: .venv\\Scripts\\python python/make_synthesis_figures.py   -> docs/figures/syn_*.pdf/.png
"""
import csv
import gzip
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from matplotlib import font_manager
from matplotlib import patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
RES = os.path.join(ROOT, "results")
F6 = os.path.join(RES, "f6_paper")
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
B300, B650, B700 = "#6da7ec", "#104281", "#0d366b"
INK, INK2, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df", "#ecebe7"
COL = {"PI": ORANGE, "PID": AQUA, "N1b": B300, "N2b": BLUE, "N2": BLUE, "N3b": B650,
       "N2_var_a0.02": B700, "meta": MUTED}
LAB = {"PI": "PI del DET", "PID": "PID de la tesis", "N1b": "N1 (QP-LTI)", "N2b": "N2 (QP-LTV)",
       "N2": "N2", "N3b": "N3 (NLP)", "N2_var_a0.02": r"N2 + cota de riesgo ($\alpha$ = 0.02)"}
plt.rcParams.update({
    "font.family": "Inter Variable", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 7, "lines.linewidth": 1.5,
    "figure.dpi": 150, "savefig.bbox": "tight", "pdf.fonttype": 42,
})
W = 6.3                      # text width [in] (A4, 25 mm margins)
Y_REF = 1179.2711775986722   # PCT set-point [K]
LIM = 1220.0                 # end state of the DET tree [K]
K_H = 26                     # MPC horizon (audit ΔRMS)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)
    print("wrote", name)


def rows(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def f6(name):
    return rows(os.path.join(F6, name + ".csv"))


def ci(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return v.mean(), 1.96 * v.std(ddof=1) / np.sqrt(len(v))


def load_traj(name, keep=None):
    """key -> dict of arrays (t, pct, m, m_plant, solve); keep(key) filters while streaming."""
    acc = defaultdict(lambda: defaultdict(list))
    with gzip.open(os.path.join(F6, name + "_traj.csv.gz"), "rt") as f:
        for r in csv.DictReader(f):
            if keep and not keep(r["key"]):
                continue
            a = acc[r["key"]]
            for c in ("t", "pct", "m", "m_plant", "solve"):
                a[c].append(float(r[c]))
    return {k: {c: np.array(v) for c, v in a.items()} for k, a in acc.items()}


def parts(key):          # campaign|family|seed|act|case
    c, fam, s, act, case = key.split("|")
    return int(s), act, case


def limit_lines(a, x1, label=True, ref=True, below=True):
    a.axhline(LIM, color=INK2, lw=0.9, ls=(0, (4, 2)))
    if label:
        a.text(x1, LIM - 0.8 if below else LIM + 0.5, "1220 K ", ha="right", va="top" if below else "bottom",
               fontsize=6.3, color=INK2)
    if ref:
        a.axhline(Y_REF, color=MUTED, lw=0.8, ls=":")


def bold(c=INK, w=0.55):
    """The variable Inter font exposes only one weight to matplotlib: fake bold with a stroke."""
    return dict(path_effects=[pe.withStroke(linewidth=w, foreground=c)], color=c)


def tag(a, s):
    a.set_title(f"({s})  " + a.get_title())


# =================================================================================================
# data shared by several figures
# =================================================================================================
DET = f6("det")
DETROW = {(int(r["seed"]), r["case"]): r for r in DET}
DET_CASES = ("PI", "PID", "N2b", "N2_var_a0.02")
det_traj = load_traj("det", keep=lambda k: k.split("|")[-1] in DET_CASES + ("N1b", "N3b", "N2_wilks59"))

# failure windows (first failure; plant flow is 0 while the ACS is down)
def failure_window(seed, case):
    r = DETROW[(seed, case)]
    return float(r["first_fail"]), float(r["max_fail_len"]) if r["n_fail"] == "1" else np.nan


# aligned fans: s = t − t_fail on a common grid; after a crossing the loop stops (≥ 1220 K)
S_GRID = np.arange(-60, 180.01, 0.5)


def aligned(case, col="pct"):
    M, crossed, D = [], [], []
    for (s, a, c), tr in ((parts(k), v) for k, v in det_traj.items()):
        if c != case:
            continue
        tf = float(DETROW[(s, c)]["first_fail"])
        x = tr["t"] - tf
        y = np.interp(S_GRID, x, tr[col], left=np.nan, right=np.nan)
        end = x[-1]
        hit = float(DETROW[(s, c)]["pct_max"]) >= LIM - 0.1
        if hit and col == "pct":
            y[S_GRID > end] = 1226.0          # beyond the axis: the run has crossed
        M.append(y); crossed.append(hit); D.append(float(DETROW[(s, c)]["max_fail_len"]))
    return np.array(M), np.array(crossed), np.array(D)


# open-loop figures (t_q = 1180 s)
e2 = rows(os.path.join(RES, "e2_optimizers_1180", "e2_summary.csv"))
dec = [r for r in rows(os.path.join(RES, "f3c_optimality", "decisions_compared.csv")) if r["t_q"] == "1180"]
qp = rows(os.path.join(RES, "f3_certificates", "qp_solvers_1180.csv"))
n2it = [r for r in rows(os.path.join(RES, "f3_certificates", "n2_iterations_1180.csv")) if r["nominal_scheme"] == "ode1be"]
t_lj = np.mean([float(r["wall_s"]) for r in e2 if r["method"] == "LJ"])
t_ga = np.mean([float(r["wall_s"]) for r in e2 if r["method"] == "GA"])
t_n1 = float([r for r in qp if r["model"] == "N1" and r["form"] == "condensed" and r["solver"] == "clarabel"][0]["time_s"])
t_n2 = sum(float(r["time_s"]) for r in n2it)                 # 4 SQP iterations, condensed QP
t_n3 = float([r for r in qp if r["model"] == "N3_NLP"][0]["time_s"])
t_de = 91.6                                                  # f3c summary_1180.txt (90 060 evaluations)


def excess(pattern):
    v = [float(r["excess_vs_best_pct"]) for r in dec if pattern in r["method"]]
    return np.mean(v), np.min(v), np.max(v)


# =================================================================================================
# 1. syn_openloop: time per decision vs cost excess (open loop, t_q = 1180 s)
# =================================================================================================
FLOOR = 1e-9
pts = [  # label, time [s], excess (mean, min, max) [%], colour, marker
    ("Luus–Jaakola\n(3 semillas)", t_lj, excess("LJ corregido"), MUTED, "o"),
    ("Algoritmo genético\n(2 semillas)", t_ga, excess("GA, semilla"), MUTED, "s"),
    ("Evolución diferencial\n(90 000 eval.)", t_de, excess("Evolución diferencial"), MUTED, "D"),
    ("N1 (QP-LTI)", t_n1, excess("N1 QP-LTI"), B300, "o"),
    ("N2 (QP-LTV iterado)", t_n2, excess("N2 QP-LTV"), BLUE, "o"),
    ("N3 (NLP, Ipopt)", t_n3, excess("N3 NLP (Ipopt"), B650, "o"),
]
fig, ax = plt.subplots(1, 2, figsize=(W, 2.9), layout="constrained", gridspec_kw={"width_ratios": [1.5, 1]})
a = ax[0]
a.axvspan(1e-3, 0.5, color="#eaf2fc", lw=0, zorder=0)
a.text(0.45, 2e3, "cabe en el período\nde muestreo (0.5 s)", ha="right", va="top", fontsize=6.5, color=B650)
a.axhspan(FLOOR * 0.3, 2.39, color=BAND, lw=0, zorder=0)
a.text(2.5e3, 1.3, "≤ 2.39 %: distancia certificada al\nóptimo global (cota inferior, F3c)", ha="right",
       va="top", fontsize=6.2, color=INK2)
offs = {"Luus": (-8, 9, "right"), "Algo": (-8, -9, "right"), "Evol": (-8, 0, "right")}
for lab, t, (m, lo, hi), c, mk in pts:
    y = max(m, FLOOR)
    yerr = None if hi - lo < 1e-12 else [[y - max(lo, FLOOR)], [hi - y]]
    a.errorbar(t, y, yerr=yerr, color=c, marker=mk, ms=6, mec="white", mew=0.8, lw=1, capsize=0, zorder=3)
    dx, dy, ha = offs.get(lab[:4], (8, 0, "left"))
    a.annotate(lab.replace("\n", " "), (t, y), xytext=(dx, dy), textcoords="offset points", ha=ha, va="center",
               fontsize=6.3, color=INK)
a.set(xscale="log", yscale="log", xlim=(3e-3, 3e3), ylim=(3e-10, 5e3),
      xlabel="tiempo de cálculo por decisión [s]", ylabel="costo sobre el mejor plan [%]",
      title="Lazo abierto: calidad y tiempo de cada método")
a.set_yticks([1e-9, 1e-6, 1e-3, 1, 1e3])
a.set_yticklabels(["≤10⁻⁹", "10⁻⁶", "10⁻³", "1", "10³"])
tag(a, "a")
# (b) the same decision quality, ranked bar of time on a log axis
b = ax[1]
order = sorted(pts, key=lambda p: p[1])
y = np.arange(len(order))
b.barh(y, [p[1] for p in order], color=[p[3] for p in order], height=0.62, left=1e-3)
b.axvline(0.5, color=B650, lw=0.9, ls=(0, (3, 2)))
b.text(0.55, len(order) - 0.45, "0.5 s", color=B650, fontsize=6.5, va="center")
for yi, p in zip(y, order):
    s = f"{p[1] * 1e3:.0f} ms" if p[1] < 1 else f"{p[1]:.0f} s"
    b.text(max(p[1] * 1.25, 0.75), yi, s, va="center", fontsize=6.5, color=INK2)
b.text(1.2e4, 1.5, f"N2 ≈ {round(t_lj / t_n2, -3):,.0f}×".replace(",", " ") + "\nmás rápido que LJ",
       fontsize=7.5, ha="right", va="center", **bold(INK, 0.4))
b.set_yticks(y)
b.set_yticklabels([p[0].replace("\n", " ") for p in order], fontsize=6.4)
b.set(xscale="log", xlim=(1e-3, 2e4), xlabel="segundos por decisión", title="Tiempo por decisión")
b.grid(axis="y", visible=False)
tag(b, "b")
save(fig, "syn_openloop")

# =================================================================================================
# 2. syn_cl_compare: representative closed loops in three families
# =================================================================================================
nom_traj = load_traj("nom", keep=lambda k: k.split("|")[2] == "1" and k.split("|")[3] == "arch")
ext_rows = f6("ext")
worst = max((r for r in ext_rows if r["case"] == "N2b" and r["act"] == "phys"), key=lambda r: float(r["pct_max"]))
SEED_EXT = int(worst["seed"])
ext_traj = load_traj("ext", keep=lambda k: k.split("|")[2] == str(SEED_EXT))
SEED_DET = 208     # PI crosses 1220 K; one failure of 118 s starting at t = 37 s
panels = [
    ("Nominal (escenario de la tesis)\nbomba archivada, semilla 1", nom_traj, lambda c: f"nom|nominal|1|arch|{c}",
     ("PI", "PID", "N2b"), None),
    (f"Fallas del ACS como en el DET\nACS de 40 kg/s, semilla {SEED_DET}", det_traj,
     lambda c: f"det|detacc|{SEED_DET}|cap40|{c}", ("PI", "PID", "N2b", "N2_var_a0.02"),
     failure_window(SEED_DET, "N2b")),
    (f"Familia extrema, actuador físico\npeor caso de N2 (semilla {SEED_EXT})", ext_traj,
     lambda c: f"ext|extreme|{SEED_EXT}|phys|{c}", ("PI", "PID", "N2b"), None),
]
fig, ax = plt.subplots(2, 3, figsize=(W, 4.3), layout="constrained", sharex="col",
                       gridspec_kw={"height_ratios": [1.55, 1]})
for j, (title, T, key, cases, win) in enumerate(panels):
    a, b = ax[0, j], ax[1, j]
    if win is not None:
        for z in (a, b):
            z.axvspan(win[0], win[0] + win[1], color=BAND, lw=0, zorder=0)
        a.text(win[0] + win[1] / 2, 1174, "ACS fallado", ha="center", va="bottom", fontsize=6.3, color=INK2)
    tmax = 0
    for c in cases:
        tr = T[key(c)]
        lw = 1.7 if c.startswith("N2") else 1.2
        a.plot(tr["t"], tr["pct"], color=COL[c], lw=lw, label=LAB[c], zorder=3 if c.startswith("N2") else 2)
        b.plot(tr["t"], tr["m_plant"], color=COL[c], lw=lw if c != "PID" else 0.6,
               alpha=1 if c != "PID" else 0.75, zorder=3 if c.startswith("N2") else 2)
        tmax = max(tmax, tr["t"][-1])
        if tr["pct"].max() >= LIM - 0.1:
            a.plot(tr["t"][-1], tr["pct"][-1], marker="X", color=COL[c], ms=7, mec="white", mew=0.8, zorder=4)
    limit_lines(a, tmax, label=(j > 0))
    a.set_title(title, fontsize=7.4)
    b.set_xlabel("t desde la activación del ACS [s]", fontsize=7.4)
    a.text(tmax, Y_REF - 0.3, r"$y_{\rm ref}$", ha="right", va="top", fontsize=6.5, color=INK2)
    if j == 0:
        a.set_ylabel("PCT [K]")
        b.set_ylabel("caudal entregado [kg/s]")
nmin = min(T_["pct"].min() for k_, T_ in nom_traj.items() if k_.split("|")[-1] in ("PI", "PID", "N2b"))
ax[0, 0].set_ylim(nmin - 0.4, 1188)
ax[0, 1].set_ylim(1172, 1224)
pmin = min(T_["pct"].min() for k_, T_ in ext_traj.items() if k_.split("|")[-1] in ("PI", "PID", "N2b"))
ax[0, 2].set_ylim(np.floor(pmin / 10) * 10 - 2, 1224)
mx = {c: ext_traj[f"ext|extreme|{SEED_EXT}|phys|{c}"]["m_plant"].max() for c in ("PI", "PID", "N2b")}
ax[1, 2].set_yscale("symlog", linthresh=50)
ax[1, 2].set_ylim(-3, 1e4)
ax[1, 2].set_yticks([0, 25, 50, 500, 5000])
ax[1, 2].set_yticklabels(["0", "25", "50", "500", "5000"])
ax[1, 2].text(238, 3000, f"escala lineal hasta 50 kg/s,\nlog. encima; PID: hasta {mx['PID']:,.0f} kg/s".replace(",", " "),
              fontsize=5.9, color=INK2, ha="right", va="top")
ax[1, 0].text(58, 9, "PI: zona muerta,\nno inyecta (0 kg/s)", fontsize=6.0, color=INK2, ha="right", va="center")
ax[1, 1].annotate("pre-enfría antes\nde la falla", (12, 40), xytext=(60, 30), fontsize=6.0, color=INK2,
                  arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6), va="center")
X = Line2D([], [], marker="X", color=INK2, ls="none", ms=6, mec="white", label="cruce de 1220 K (la corrida se detiene)")
h = [Line2D([], [], color=COL[c], lw=1.7, label=LAB[c]) for c in ("N2b", "N2_var_a0.02", "PI", "PID")] + [X]
fig.legend(handles=h, loc="outside lower center", ncol=3, fontsize=6.8)
save(fig, "syn_cl_compare")

# =================================================================================================
# 3. syn_cl_fan: 50 closed loops aligned at the failure start (DET family, 40 kg/s)
# =================================================================================================
fig, ax = plt.subplots(2, 4, figsize=(W, 3.9), layout="constrained", sharex=True, sharey="row",
                       gridspec_kw={"height_ratios": [1.6, 1]})
for j, c in enumerate(("PI", "PID", "N2b", "N2_var_a0.02")):
    a, b = ax[0, j], ax[1, j]
    P, crossed, D = aligned(c)
    Q = np.nanpercentile(P, [5, 25, 50, 75, 95], axis=0)
    a.fill_between(S_GRID, Q[0], Q[4], color=COL[c], alpha=0.16, lw=0)
    a.fill_between(S_GRID, Q[1], Q[3], color=COL[c], alpha=0.32, lw=0)
    a.plot(S_GRID, Q[2], color=COL[c], lw=1.6)
    for y in P[crossed]:
        a.plot(S_GRID, np.where(y > LIM + 1, np.nan, y), color=COL[c], lw=0.6, alpha=0.9)
    a.axvline(0, color=INK2, lw=0.7)
    limit_lines(a, 178, label=(j == 2))
    n = int(crossed.sum())
    a.text(-55, 1215.5, f"{n}/50 cruzan", fontsize=7.2, va="top", **(bold() if n else dict(color=INK2)))
    a.set_title(LAB[c].replace(r" ($\alpha$ = 0.02)", "\n" + r"($\alpha$ = 0.02)"), fontsize=7.6)
    M, _, _ = aligned(c, "m_plant")
    Qm = np.nanpercentile(M, [5, 25, 50, 75, 95], axis=0)
    b.fill_between(S_GRID, Qm[0], Qm[4], color=COL[c], alpha=0.16, lw=0)
    b.fill_between(S_GRID, Qm[1], Qm[3], color=COL[c], alpha=0.32, lw=0)
    b.plot(S_GRID, Qm[2], color=COL[c], lw=1.4)
    b.axvline(0, color=INK2, lw=0.7)
    b.set_xlabel("s desde el inicio de la falla [s]", fontsize=7.2)
ax[0, 0].set_ylabel("PCT [K]")
ax[1, 0].set_ylabel("caudal entregado\n[kg/s]")
ax[0, 0].set_ylim(1168, 1224)
ax[0, 0].text(2, 1170, "falla", fontsize=6.3, color=INK2)
h = [Patch(color=INK2, alpha=0.32, label="50 % central"), Patch(color=INK2, alpha=0.16, label="90 % central"),
     Line2D([], [], color=INK2, lw=1.6, label="mediana"),
     Line2D([], [], color=INK2, lw=0.6, label="corridas que cruzan (se detienen en 1220 K)")]
fig.legend(handles=h, loc="outside lower center", ncol=4, fontsize=6.6)
save(fig, "syn_cl_fan")

# =================================================================================================
# 4. syn_scorecard: thesis metrics by controller and family
# =================================================================================================
fams = [("Nominal (tesis)\nbomba archivada", [r for r in f6("nom") if r["act"] == "arch"]),
        ("Fallas del DET\nACS de 40 kg/s", DET),
        ("Familia extrema\nactuador físico", [r for r in ext_rows if r["act"] == "phys"])]
ctrls = ("PI", "PID", "N1b", "N2b", "N3b")
mets = [("IAE [K·s]", "IAE", True), ("PCT pico [K]", "peak_PCT", False),
        ("variación total del\ncaudal entregado [kg/s]", "tv_plant", True)]
fig, ax = plt.subplots(len(mets), 3, figsize=(W, 4.7), layout="constrained", sharey=True)
for j, (fname, R) in enumerate(fams):
    for i, (mname, k, logx) in enumerate(mets):
        a = ax[i, j]
        for yi, c in enumerate(ctrls):
            v = np.array([float(r[k]) for r in R if r["case"] == c])
            if not len(v):
                continue
            m, e = ci(v)
            if k == "peak_PCT":
                jit = (np.random.default_rng(yi).random(len(v)) - 0.5) * 0.45
                a.scatter(v, yi + jit, s=5, color=COL[c], alpha=0.45, lw=0)
                a.plot([m], [yi], marker="|", ms=11, mew=2, color=INK)
                n = int((v >= LIM - 0.1).sum())
                if n:
                    a.text(LIM + 2, yi, f"{n}/50 cruzan", va="center", fontsize=6.0, linespacing=0.9, **bold())
            else:
                a.barh(yi, m, xerr=e, color=COL[c], height=0.62, error_kw=dict(lw=0.8, ecolor=INK2),
                       left=1e-1 if logx else 0)
                a.text(m * 1.5 if logx else m, yi, f"{m:,.0f}".replace(",", " ") if m >= 10 else f"{m:.1f}",
                       va="center", fontsize=6.2, color=INK2)
        if k == "peak_PCT":
            a.axvline(LIM, color=INK2, lw=0.9, ls=(0, (4, 2)))
            a.axvline(Y_REF, color=MUTED, lw=0.8, ls=":")
            a.set_xlim(1176, 1236)
        if logx:
            a.set_xscale("log")
            a.set_xlim(1, 3e5 if k == "tv_plant" else 3e4)
        a.grid(axis="y", visible=False)
        if j == 0:
            a.set_ylabel(mname, fontsize=7.4)
        if i == 0:
            a.set_title(fname, fontsize=7.8)
ax[0, 0].set_yticks(range(len(ctrls)))
ax[0, 0].invert_yaxis()
for i in range(len(mets)):
    ax[i, 0].set_yticklabels([LAB[c] for c in ctrls], fontsize=6.8)
ax[2, 0].text(1.3, 0, "0 (no inyecta)", va="center", fontsize=6.2, color=INK2)
fig.text(0.5, -0.01, "Barras: media de 50 realizaciones ± IC 95 %. Puntos: cada realización; | = media. "
         "Escala log. en IAE y variación total.", ha="center", fontsize=6.6, color=INK2)
save(fig, "syn_scorecard")

# =================================================================================================
# 5. syn_mechanism: (a) the heat-up during a failure does not depend on the controller (A1);
#    (b) failure-length survival S(D); (c) time-to-limit τ that each controller keeps in the plant.
#    Risk = S evaluated at τ (Rao–Blackwell): a bound on the PCT is a floor on τ, i.e. one margin.
# =================================================================================================
db = scipy.io.loadmat(os.path.join(ROOT, "models", "MCDET_conferencepaper_database.mat"),
                      variable_names=["Lambda_t", "Mu_t"])
LAM = db["Lambda_t"].ravel() * 1e5
MU = db["Mu_t"].ravel() * 1e4
WGT = LAM / LAM.sum()
surv = lambda d: (WGT[:, None] * np.exp(-MU[:, None] * np.atleast_1d(d)[None, :])).sum(0)
XMAX = 340
Dg = np.linspace(0, XMAX, 1361)
Sg = surv(Dg)

curves = defaultdict(list)          # ΔT(s) = PCT(t_f + s) − PCT(t_f) while the ACS is down
for k, tr in det_traj.items():
    s, act, c = parts(k)
    tf, Dlen = failure_window(s, c)
    if not np.isfinite(Dlen):
        continue
    x = tr["t"] - tf
    sel = (x >= 0) & (x <= Dlen)
    if sel.sum() < 4:
        continue
    curves[c].append((x[sel], tr["pct"][sel] - np.interp(0.0, x, tr["pct"])))
# same seed, different controllers: spread of ΔT at equal s (the A1 test)
by_seed = defaultdict(dict)
for k, tr in det_traj.items():
    s, act, c = parts(k)
    tf, Dlen = failure_window(s, c)
    if np.isfinite(Dlen) and Dlen >= 40:
        x = tr["t"] - tf
        if x[-1] >= 40:
            by_seed[s][c] = np.interp(40.0, x, tr["pct"]) - np.interp(0.0, x, tr["pct"])
spread = [max(v.values()) - min(v.values()) for v in by_seed.values() if len(v) >= 5]
level = [np.mean(list(v.values())) for v in by_seed.values() if len(v) >= 5]

sens = [r for r in f6("sens") if r["act"] == "cap40"]
tau_rows = defaultdict(list)
with gzip.open(os.path.join(F6, "sens_tau.csv.gz"), "rt") as f:
    for r in csv.DictReader(f):
        kk = r["key"].split("|")
        if kk[3] == "cap40" and float(r["t"]) >= 60:
            tau_rows[kk[4]].append(float(r["tau"]))
TAU_CASES = [("PI", "PI del DET", ORANGE), ("PID", "PID de la tesis", AQUA), ("N2", "N2 sin cota", BLUE),
             ("N2_var_a0.05_mu1e4", r"cota $\alpha$ = 0.05", B300),
             ("N2_var_a0.02_mu1e4", r"cota $\alpha$ = 0.02", B650),
             ("N2_var_a0.01_mu1e4", r"cota $\alpha$ = 0.01", B700)]

fig = plt.figure(figsize=(W, 5.6), layout="constrained")
gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 0.85, 1.25])
a = fig.add_subplot(gs[0])
for c in ("PI", "PID", "N1b", "N2b", "N3b", "N2_var_a0.02", "N2_wilks59"):
    for xs, dT in curves.get(c, []):
        a.plot(xs, dT, color=COL.get(c, B700), lw=0.5, alpha=0.4)
a.set(xlim=(0, XMAX), ylim=(-2, 50), ylabel="ΔT durante la falla [K]",
      title="El calentamiento durante una falla depende poco del controlador (A1, aproximado)")
d_pre = [v["N2_var_a0.02"] - v["N2b"] for v in by_seed.values() if "N2_var_a0.02" in v and "N2b" in v]
d_wil = [v["N2_wilks59"] - v["N2b"] for v in by_seed.values() if "N2_wilks59" in v and "N2b" in v]
d_n13 = [max(abs(v["N1b"] - v["N2b"]), abs(v["N3b"] - v["N2b"])) for v in by_seed.values()
         if all(c in v for c in ("N1b", "N2b", "N3b"))]
a.text(0.99, 0.95, f"{sum(len(v) for v in curves.values())} fallas observadas, 7 controladores.\n"
       f"Misma semilla, 40 s después del inicio: ΔT = {np.median(level):.1f} K (mediana).\n"
       f"N1, N2, N3: idénticos (|dif.| ≤ {max(d_n13):.0e} K). Con pre-enfriamiento se calienta menos:\n"
       f"cota α = 0.02: {np.median(d_pre):+.1f} K; Wilks: {np.median(d_wil):+.1f} K (medianas), del lado seguro.",
       transform=a.transAxes, ha="right", va="top", fontsize=6.3, color=INK2)
tag(a, "a")
bx = fig.add_subplot(gs[1], sharex=a)
bx.semilogy(Dg, Sg, color=INK, lw=1.6)
bx.set(ylim=(2e-3, 1.3), ylabel="S(d) = P(D > d)",
       title=r"Duración de la falla: mezcla de exponenciales del DET ($\mu\times10^4$)")
for al, c in ((0.05, B300), (0.02, B650), (0.01, B700)):
    q = np.interp(-al, -Sg, Dg)
    bx.plot([0, q, q], [al, al, 2e-3], color=c, lw=0.9, ls=(0, (2, 2)))
    bx.text(q + 3, al * 1.15, rf"$q_{{1-\alpha}}$ = {q:.0f} s ($\alpha$ = {al})", fontsize=6.2, color=c,
            va="bottom")
tag(bx, "b")
cx = fig.add_subplot(gs[2], sharex=a)
for i, (c, lab, col) in enumerate(TAU_CASES):
    v = np.array(tau_rows[c])
    p5, p25, p50, p75, p95 = np.percentile(v, [5, 25, 50, 75, 95])
    yi = len(TAU_CASES) - 1 - i
    cx.plot([p5, p95], [yi, yi], color=col, lw=1.2, solid_capstyle="butt")
    cx.plot([p25, p75], [yi, yi], color=col, lw=6, solid_capstyle="butt")
    cx.plot([p50], [yi], marker="|", ms=10, mew=2, color=INK)
    R = [r for r in sens if r["case"] == c]
    p, _ = ci([float(r["p_cond_mu1e4"]) for r in R])
    pre, _ = ci([-float(r["mean_dev"]) for r in R])
    cx.text(XMAX - 2, yi, f"P = {100 * p:.1f} %  ·  pre-enfría {max(pre, 0):.0f} K", ha="right", va="center",
            fontsize=6.3, color=INK, bbox=dict(fc="white", ec="none", pad=1.0))
cx.set_yticks(range(len(TAU_CASES)))
cx.set_yticklabels([t[1] for t in TAU_CASES][::-1], fontsize=6.8)
cx.grid(axis="y", visible=False)
cx.set(xlabel="tiempo [s]: duración de la falla D en (a)–(b); τ = tiempo hasta 1220 K si el ACS falla ahora en (c)",
       title="Lo que cada controlador garantiza: τ en la planta (t ≥ 60 s, 50 realizaciones)")
cx.text(0.01, 0.02, "riesgo = S(τ) promediado sobre el instante de la falla (Rao–Blackwell):\n"
        r"exigir P $\leq\alpha$ equivale a exigir τ $\geq q_{1-\alpha}$ $\Rightarrow$ un solo margen sobre la PCT",
        transform=cx.transAxes, fontsize=6.3, color=INK2, va="bottom")
tag(cx, "c")
for z in (a, bx):
    plt.setp(z.get_xticklabels(), visible=False)
save(fig, "syn_mechanism")

# =================================================================================================
# 6. syn_audit: 11 400 audited closed-loop queries — distance to the best plan and solve time
# =================================================================================================
aud = {}
for fam in ("nom", "det", "ext"):
    for r in f6(fam + "_audit"):
        c = r["key"].split("|")[-1]
        jb, g = float(r["Jbest"]), float(r["gap_loop"])
        d = np.sqrt(max(jb * (1 + g), 0) / K_H) - np.sqrt(max(jb, 0) / K_H)
        aud.setdefault(c, {"d": [], "t": []})
        aud[c]["d"].append(max(d, 1e-12))
        aud[c]["t"].append(float(r["t_" + c[:2]]))
risk_t = np.concatenate([tr["solve"] for k, tr in det_traj.items()
                         if k.split("|")[-1] in ("N2_var_a0.02", "N2_wilks59") and np.isfinite(tr["solve"]).all()])
fig, ax = plt.subplots(1, 2, figsize=(W, 3.0), layout="constrained")


def ecdf(a, v, c, lab, ls="-"):
    v = np.sort(np.asarray(v))
    a.plot(v, np.arange(1, len(v) + 1) / len(v), color=c, lw=1.6, ls=ls, label=lab, drawstyle="steps-post")


for c in ("N1b", "N2b", "N3b"):
    ecdf(ax[0], aud[c]["d"], COL[c], f"{LAB[c]}  (máx. {max(aud[c]['d']):.0e} K)".replace("e-0", "e-"))
ax[0].axvline(1e-3, color=INK2, lw=0.8, ls=(0, (3, 2)))
ax[0].text(1.15e-3, 0.06, "1 mK", fontsize=6.5, color=INK2)
ax[0].set(xscale="log", xlim=(1e-12, 1), ylim=(0, 1.02), xlabel="exceso de error RMS sobre el mejor plan [K]",
          ylabel="fracción de consultas", title=f"Calidad: {sum(len(v['d']) for v in aud.values()):,} consultas".replace(",", " "))
ax[0].set_xticks([1e-12, 1e-9, 1e-6, 1e-3, 1])
ax[0].set_xticklabels(["≤10⁻¹²", "10⁻⁹", "10⁻⁶", "10⁻³", "1"])
ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=6.2, ncol=2, columnspacing=1.0)
tag(ax[0], "a")
T99 = max(np.percentile(aud[c]["t"], 99) for c in aud)
TMX = max(max(aud[c]["t"]) for c in aud)
for c in ("N1b", "N2b", "N3b"):
    ecdf(ax[1], aud[c]["t"], COL[c], LAB[c])
ecdf(ax[1], risk_t, B700, "N2 + cota de riesgo (en el lazo)", ls=(0, (3, 1.5)))
ax[1].axvline(0.5, color=B650, lw=0.9, ls=(0, (3, 2)))
ax[1].text(0.45, 0.06, "período de\nmuestreo 0.5 s", fontsize=6.3, color=B650, ha="right")
ax[1].set(xscale="log", xlim=(5e-3, 1), ylim=(0, 1.02), xlabel="tiempo de resolución por consulta [s]",
          title=f"Tiempo: p99 ≤ {T99:.2f} s (máx. {TMX:.2f} s)")
ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), fontsize=6.2, ncol=2, columnspacing=1.0)
tag(ax[1], "b")
save(fig, "syn_audit")
for c in ("N1b", "N2b", "N3b"):
    print(c, "p99 t =", np.percentile(aud[c]["t"], 99), "max t =", max(aud[c]["t"]), "max d =", max(aud[c]["d"]))
print("risk p99 t =", np.percentile(risk_t, 99), "max", risk_t.max())

# =================================================================================================
# 7. syn_overview: graphical abstract (three columns, a headline number each)
# =================================================================================================
fig = plt.figure(figsize=(W, 3.7), layout="constrained")
gs = fig.add_gridspec(2, 3, height_ratios=[0.6, 1])
p_n2, _ = ci([float(r["p_cond_mu1e4"]) for r in sens if r["case"] == "N2"])
heads = [("1 · Resolver bien", f"≈ {round(t_lj / t_n2, -3):,.0f}×".replace(",", " "),
          "más rápido que Luus–Jaakola (N2: 42 ms),\ny ≥ 28 % más barato (cota certificada)"),
         ("2 · Controlar bien", "0 / 50", "cruces de 1220 K con N2 bajo fallas del\nACS (40 kg/s); PI del DET: 7 / 50"),
         ("3 · Gestionar el riesgo", "1 margen",
          "cuantil, CVaR, Wilks y pérdida esperada\ncaen sobre una sola frontera")]
for j, (h1, big, sub) in enumerate(heads):
    t = fig.add_subplot(gs[0, j])
    t.axis("off")
    t.text(0, 1.0, h1, fontsize=8.2, va="top", **bold(B650, 0.4))
    t.text(0, 0.74, big, fontsize=17, va="top", **bold(INK, 0.6))
    t.text(0, 0.16, sub, fontsize=6.6, color=INK2, va="bottom")
# (1) time per decision
a = fig.add_subplot(gs[1, 0])
bars = [("LJ", t_lj, MUTED), ("GA", t_ga, MUTED), ("N3", t_n3, B650), ("N2", t_n2, BLUE), ("N1", t_n1, B300)]
a.barh(range(5), [b[1] for b in bars], color=[b[2] for b in bars], height=0.6, left=1e-3)
a.set_yticks(range(5))
a.set_yticklabels([b[0] for b in bars])
a.invert_yaxis()
a.axvline(0.5, color=B650, lw=0.9, ls=(0, (3, 2)))
a.text(0.7, -0.62, "0.5 s", color=B650, fontsize=6.3, ha="left", va="center")
for i, b in enumerate(bars):
    a.text(max(b[1] * 1.3, 0.8), i, f"{b[1] * 1e3:.0f} ms" if b[1] < 1 else f"{b[1]:.0f} s", va="center",
           fontsize=6.2, color=INK2)
a.set(xscale="log", xlim=(1e-3, 3e4), xlabel="s por decisión (lazo abierto)")
a.grid(axis="y", visible=False)
# (2) peak PCT of every realisation under ACS failures (DET family, 40 kg/s)
a = fig.add_subplot(gs[1, 1])
OVC = ("PI", "PID", "N2b", "N2_var_a0.02")
for yi, c in enumerate(OVC):
    v = np.array([float(DETROW[(s, c)]["pct_max"]) for s in sorted({k[0] for k in DETROW})])
    jit = (np.random.default_rng(yi).random(len(v)) - 0.5) * 0.5
    a.scatter(v, yi + jit, s=7, color=COL[c], alpha=0.6, lw=0)
    n = int((v >= LIM - 0.1).sum())
    a.text(1226, yi, f"{n}/50", va="center", fontsize=6.6, **(bold() if n else dict(color=INK2)))
a.axvline(LIM, color=INK2, lw=0.9, ls=(0, (4, 2)))
a.axvline(Y_REF, color=MUTED, lw=0.8, ls=":")
a.set_yticks(range(len(OVC)))
a.set_yticklabels(["PI DET", "PID tesis", "N2", "N2 + cota"], fontsize=7)
a.invert_yaxis()
a.grid(axis="y", visible=False)
a.text(1226, -0.62, "cruces", fontsize=6.2, color=INK2, va="center")
a.set(xlim=(1172, 1236), xlabel="PCT pico por realización [K]")
# (3) one frontier
a = fig.add_subplot(gs[1, 2])
form = f6("form") + [r for r in sens if r["case"].endswith("_mu1e4") or r["case"] == "N2"]
fam_ov = [(["N2_var_a0.1", "N2_var_a0.05_mu1e4", "N2_var_a0.02_mu1e4", "N2_var_a0.01_mu1e4"], BLUE, "o", "cuantil"),
          (["N2_cvar_b0.9", "N2_cvar_b0.95_mu1e4", "N2_cvar_b0.98"], B650, "s", "CVaR"),
          (["N2_wilks59"], ORANGE, "*", "Wilks 59"),
          (["N2_exp_r10", "N2_exp_r100", "N2_exp_r1000"], AQUA, "^", "pérdida esperada")]
allp = []
for cases, c, mk, lab in fam_ov:
    P = []
    for cs in cases:
        R = [r for r in form if r["case"] == cs]
        x, _ = ci([-float(r["mean_dev"]) for r in R])
        y, _ = ci([float(r["p_cond_mu1e4"]) for r in R])
        P.append((x, 100 * y))
    allp += P
    a.plot(*zip(*P), marker=mk, color=c, ls="none", ms=7 if mk == "*" else 4.5, mec="white", mew=0.6, label=lab)
allp.sort()
a.plot(*zip(*allp), color=MUTED, lw=0.8, zorder=0)
x0, _ = ci([-float(r["mean_dev"]) for r in sens if r["case"] == "N2"])
a.plot([x0], [100 * p_n2], marker="D", color=INK2, ms=4.5, ls="none", label="N2 sin cota")
a.set(xlabel="pre-enfriamiento [K]", ylabel="P(cruce | falla) [%]", xlim=(-1, 16), ylim=(0.5, 3.6))
a.legend(loc="lower left", fontsize=5.9, handletextpad=0.2)
save(fig, "syn_overview")
