"""Closed-loop trajectory dynamics compared across formulations (status report, 20-11-2026).

The synthesis figures (make_synthesis_figures.py) compare N2 with the thesis PID and the DET PI.
This script compares the MPC *formulations* with each other, trajectory by trajectory, on the
F6 campaign (50 realisations per family, common random numbers):

  traj_formulations   N1 / N2 / N3: PCT and delivered flow (median, 10–90 % band) and the
                      largest pointwise PCT difference to N3 over the 50 realisations
  traj_extreme_n1     extreme family: where and why N1 leaves N3 (worst realisation + scatter)
  traj_risk_dynamics  DET family aligned at the failure start: the risk formulations pre-cool
                      by different amounts and then heat up by the same ΔT
  traj_hierarchy      what changes a closed-loop trajectory: per-realisation max |ΔPCT| of
                      every pairwise comparison, on one log axis

Also writes results/f6_paper/trajectory_distances.csv (the numbers quoted in the report).
Reads only results/f6_paper (committed). Same style as make_synthesis_figures.py.
Usage: .venv\\Scripts\\python python/make_trajectory_figures.py
"""
import csv
import gzip
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
F6 = os.path.join(ROOT, "results", "f6_paper")
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
B300, B650, B700 = "#6da7ec", "#104281", "#0d366b"
VIOLET, MAG, YEL = "#7b5fd6", "#c2449b", "#c49a00"
INK, INK2, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df", "#ecebe7"
COL = {"N1b": B300, "N2b": BLUE, "N3b": B650, "N2": "#9aa9bb", "PI": ORANGE, "PID": AQUA,
       "N2_var_a0.02": B700, "N2_cvar_b0.95": VIOLET, "N2_wilks59": MAG, "N2_exp_r100": YEL,
       "N2_DSM_U1244_l1e4": AQUA}
LAB = {"N1b": "N1 (QP-LTI)", "N2b": "N2 (QP-LTV)", "N3b": "N3 (NLP)", "N2": "N2 sin sesgo",
       "PI": "PI del DET", "PID": "PID de la tesis",
       "N2_var_a0.02": r"cuantil $\alpha$=0.02", "N2_cvar_b0.95": r"CVaR$_{0.95}$",
       "N2_wilks59": "Wilks 59", "N2_exp_r100": r"pérdida esperada $\rho$=100",
       "N2_DSM_U1244_l1e4": r"DSM $\lambda$=10$^4$"}
plt.rcParams.update({
    "font.family": "Inter Variable", "font.size": 8.5, "axes.titlesize": 8.6, "axes.labelsize": 8.3,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 7, "lines.linewidth": 1.5,
    "figure.dpi": 150, "savefig.bbox": "tight", "pdf.fonttype": 42,
})
W = 6.3
Y_REF = 1179.2711775986722
LIM = 1220.0
RES = 1e-4               # PCT is logged with 4 decimals: differences below this are not resolved


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)
    print("wrote", name)


def rows(name):
    with open(os.path.join(F6, name + ".csv"), newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def load_traj(name):
    """(act, case) -> seed -> dict of arrays."""
    acc = defaultdict(lambda: defaultdict(list))
    with gzip.open(os.path.join(F6, name + "_traj.csv.gz"), "rt") as f:
        for r in csv.DictReader(f):
            a = acc[r["key"]]
            for c in ("t", "pct", "m", "m_plant", "ubound"):
                a[c].append(float(r[c]))
    out = defaultdict(dict)
    for k, a in acc.items():
        _, _, s, act, case = k.split("|")
        out[(act, case)][int(s)] = {c: np.array(v) for c, v in a.items()}
    return out


def tag(a, s):
    a.set_title(f"({s})  " + a.get_title())


def stack(T, act, case, col, n):
    """seeds × n matrix (NaN after a loop stops at 1220 K)."""
    M = []
    for s in sorted(T[(act, case)]):
        v = T[(act, case)][s][col]
        y = np.full(n, np.nan)
        y[:min(n, len(v))] = v[:n]
        M.append(y)
    return np.array(M)


def absdiff(T, act, a, b, n):
    A, B = stack(T, act, a, "pct", n), stack(T, act, b, "pct", n)
    return np.abs(A - B)


def per_real(T, act, a, b, col="pct"):
    """Per-realisation max |Δ| between controllers a and b (common seeds, common length)."""
    out = {}
    for s in sorted(set(T[(act, a)]) & set(T[(act, b)])):
        A, B = T[(act, a)][s], T[(act, b)][s]
        n = min(len(A["t"]), len(B["t"]))
        d = np.abs(A[col][:n] - B[col][:n])
        out[s] = (d.max(), np.sqrt((d ** 2).mean()), A["t"][np.argmax(d > 0.01)] if (d > 0.01).any() else np.inf)
    return out


NOM, EXT, DET = load_traj("nom"), load_traj("ext"), load_traj("det")
EXTROW = {(int(r["seed"]), r["case"]): r for r in rows("ext") if r["act"] == "phys"}
DETROW = {(int(r["seed"]), r["case"]): r for r in rows("det")}

# =================================================================================================
# numbers: trajectory distances (CSV, quoted in docs/estado/11k_trayectorias.tex)
# =================================================================================================
COMPARE = [  # (family label, data, act, a, b)
    ("nominal, bomba archivada", NOM, "arch", "N1b", "N3b"),
    ("nominal, bomba archivada", NOM, "arch", "N2b", "N3b"),
    ("nominal, 40 kg/s", NOM, "cap40", "N1b", "N3b"),
    ("nominal, 40 kg/s", NOM, "cap40", "N2b", "N3b"),
    ("fallas del DET, 40 kg/s", DET, "cap40", "N1b", "N3b"),
    ("fallas del DET, 40 kg/s", DET, "cap40", "N2b", "N3b"),
    ("extrema, actuador físico", EXT, "phys", "N1b", "N3b"),
    ("extrema, actuador físico", EXT, "phys", "N2b", "N3b"),
    ("fallas del DET, 40 kg/s", DET, "cap40", "N2_DSM_U1244_l1e4", "N2b"),
    ("extrema, actuador físico", EXT, "phys", "N2", "N2b"),
    ("fallas del DET, 40 kg/s", DET, "cap40", "N2_exp_r100", "N2b"),
    ("fallas del DET, 40 kg/s", DET, "cap40", "N2_var_a0.02", "N2b"),
    ("fallas del DET, 40 kg/s", DET, "cap40", "N2_cvar_b0.95", "N2b"),
    ("fallas del DET, 40 kg/s", DET, "cap40", "N2_wilks59", "N2b"),
    ("nominal, 40 kg/s", NOM, "cap40", "PID", "N2b"),
    ("extrema, actuador físico", EXT, "phys", "PI", "N2b"),
]
dist_rows, DIST = [], {}
for fam, T, act, a, b in COMPARE:
    pr = per_real(T, act, a, b)
    pm = per_real(T, act, a, b, col="m_plant")
    mx = np.array([v[0] for v in pr.values()])
    rms = np.array([v[1] for v in pr.values()])
    first = np.array([v[2] for v in pr.values()])
    dm = np.array([v[0] for v in pm.values()])
    DIST[(fam, a, b)] = mx
    fin = first[np.isfinite(first)]
    dist_rows.append(dict(family=fam, act=act, a=a, b=b, n=len(mx),
                          max_dpct_med=np.median(mx), max_dpct_p90=np.quantile(mx, 0.9), max_dpct_max=mx.max(),
                          rms_dpct_med=np.median(rms), max_dm_med=np.median(dm), max_dm_max=dm.max(),
                          frac_above_0p01K=np.mean(np.isfinite(first)),
                          t_first_0p01K_med=np.median(fin) if len(fin) else np.nan))
with open(os.path.join(F6, "trajectory_distances.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(dist_rows[0]))
    w.writeheader()
    for r in dist_rows:
        w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
for r in dist_rows:
    print(f"{r['family']:26s} {r['a']:>18s} vs {r['b']:4s}  max|dPCT| med {r['max_dpct_med']:.3g}"
          f"  p90 {r['max_dpct_p90']:.3g}  max {r['max_dpct_max']:.3g}  |dm| med {r['max_dm_med']:.3g}"
          f"  >0.01K {r['frac_above_0p01K']:.2f}  t1 {r['t_first_0p01K_med']:.3g}")

# =================================================================================================
# 1. traj_formulations: N1 / N2 / N3 in the three families
# =================================================================================================
FAMS = [("Nominal, 40 kg/s", NOM, "cap40", 121), ("Fallas del DET, 40 kg/s", DET, "cap40", 481),
        ("Extrema, actuador físico", EXT, "phys", 481)]
FORMS = ("N3b", "N2b", "N1b")
fig, ax = plt.subplots(3, 3, figsize=(W, 6.0), layout="constrained", sharex="col",
                       gridspec_kw={"height_ratios": [1.3, 1, 1]})
for j, (title, T, act, n) in enumerate(FAMS):
    t = np.arange(n) * 0.5
    a, b, c = ax[0, j], ax[1, j], ax[2, j]
    for case in FORMS:
        P, M = stack(T, act, case, "pct", n), stack(T, act, case, "m_plant", n)
        q = lambda X, p: np.nanquantile(X, p, axis=0)
        ls = {"N3b": "-", "N2b": (0, (5, 2)), "N1b": (0, (1.5, 1.5))}[case]
        if case == "N3b":
            a.fill_between(t, q(P, 0.1), q(P, 0.9), color=COL[case], alpha=0.15, lw=0)
            b.fill_between(t, q(M, 0.1), q(M, 0.9), color=COL[case], alpha=0.15, lw=0)
        else:   # N1/N2 band edges as thin lines: visible if they leave N3's band
            for p in (0.1, 0.9):
                a.plot(t, q(P, p), color=COL[case], lw=0.5, ls=ls)
                b.plot(t, q(M, p), color=COL[case], lw=0.5, ls=ls)
        a.plot(t, q(P, 0.5), color=COL[case], lw=1.4 if case == "N3b" else 1.1, ls=ls, label=LAB[case])
        b.plot(t, q(M, 0.5), color=COL[case], lw=1.4 if case == "N3b" else 1.1, ls=ls)
    a.axhline(Y_REF, color=MUTED, lw=0.8, ls=":")
    if j > 0:
        a.axhline(LIM, color=INK2, lw=0.8, ls=(0, (4, 2)))
    for case, mk in (("N1b", B300), ("N2b", BLUE)):
        D = absdiff(T, act, case, "N3b", n)
        c.plot(t, np.maximum(np.nanmax(D, axis=0), RES), color=COL[case], lw=1.1,
               label=f"{LAB[case].split()[0]} − N3: máximo")
        c.plot(t, np.maximum(np.nanmedian(D, axis=0), RES), color=COL[case], lw=0.8, ls=(0, (2, 1.5)),
               label=f"{LAB[case].split()[0]} − N3: mediana")
    c.set_yscale("log"); c.set_ylim(5e-5, 50)
    c.axhspan(5e-5, RES, color=BAND, lw=0)
    a.set_title(title)
    c.set_xlabel("t desde el inicio del lazo [s]")
    if j == 0:
        a.set_ylabel("PCT [K]"); b.set_ylabel("caudal entregado [kg/s]")
        c.set_ylabel(r"$|\Delta$PCT$|$ frente a N3 [K]")
        c.text(1, 1.3e-4, "bajo la resolución del registro", fontsize=6, color=INK2, va="bottom")
    if j == 2:
        b.set_yscale("symlog", linthresh=50); b.set_ylim(0, 3000)
for k, s in enumerate("abc"):
    tag(ax[0, k], s)
ax[0, 0].legend(loc="upper right", handlelength=2.4)
ax[2, 0].legend(loc="upper right", fontsize=6.2, handlelength=2.2)
save(fig, "traj_formulations")

# =================================================================================================
# 2. traj_extreme_n1: where N1 leaves N3 in the extreme family
# =================================================================================================
pr = per_real(EXT, "phys", "N1b", "N3b")
sw = max(pr, key=lambda s: pr[s][0])
fig = plt.figure(figsize=(W, 3.6), layout="constrained")
gs = fig.add_gridspec(2, 2, width_ratios=[1.6, 1], height_ratios=[1.4, 1])
a, b, c = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[:, 1])
r = EXTROW[(sw, "N3b")]
tf, D = float(r["first_fail"]), float(r["max_fail_len"])
for z in (a, b):
    if np.isfinite(tf) and D > 0:
        z.axvspan(tf, tf + D, color=BAND, lw=0, zorder=0)
for case in ("N2", "N1b", "N2b", "N3b"):
    tr = EXT[("phys", case)][sw]
    ls = {"N3b": "-", "N2b": (0, (5, 2)), "N1b": "-", "N2": "-"}[case]
    a.plot(tr["t"], tr["pct"], color=COL[case], lw=1.3, ls=ls, label=LAB[case])
    b.plot(tr["t"], tr["m_plant"], color=COL[case], lw=1.1, ls=ls)
a.axhline(Y_REF, color=MUTED, lw=0.8, ls=":")
a.text(tf + D / 2, a.get_ylim()[0] + 0.5, f"ACS fallado\n({D:.0f} s, η = {float(r['eta']):.2f})",
       ha="center", va="bottom", fontsize=6.2, color=INK2)
a.set_title(f"Realización {sw}: la de mayor diferencia N1 − N3")
a.set_ylabel("PCT [K]"); b.set_ylabel("caudal [kg/s]"); b.set_xlabel("t [s]")
b.set_yscale("symlog", linthresh=50)
a.legend(loc="upper right", fontsize=6.4)
plt.setp(a.get_xticklabels(), visible=False)
# scatter: per realisation max |ΔPCT| vs the failure severity η·D (mass of cooling lost)
xs, ys, cs = [], [], []
for s, v in pr.items():
    rr = EXTROW[(s, "N3b")]
    xs.append(float(rr["max_fail_len"]) * (1 - float(rr["eta"])) if rr["n_fail"] != "0" else 0.0)
    ys.append(max(v[0], RES)); cs.append(float(rr["pct_max"]))
sc = c.scatter(xs, ys, c=cs, cmap="Blues", s=16, edgecolor=INK2, lw=0.3, vmin=1180, vmax=1215)
c.set_yscale("log"); c.set_ylim(5e-5, 60)
c.set_xlabel(r"falla efectiva $(1-\eta)\,D$ [s]")
c.set_ylabel(r"máx $|$PCT$_{N1}-$PCT$_{N3}|$ [K]")
c.set_title("50 realizaciones")
cb = fig.colorbar(sc, ax=c, shrink=0.8, pad=0.01)
cb.set_label("PCT pico de N3 [K]", fontsize=7)
tag(a, "a"); tag(c, "b")
save(fig, "traj_extreme_n1")
cor = np.corrcoef(np.array(xs), np.log10(np.array(ys)))[0, 1]
print(f"extreme: worst seed {sw}, max dPCT {pr[sw][0]:.2f} K; corr(log dPCT, (1-eta)D) = {cor:.2f}")

# =================================================================================================
# 3. traj_risk_dynamics: risk formulations aligned at the failure start (DET family)
# =================================================================================================
RISK = ("N2b", "N2_DSM_U1244_l1e4", "N2_exp_r100", "N2_cvar_b0.95", "N2_var_a0.02", "N2_wilks59")
S = np.arange(-60, 150.01, 0.5)


def aligned(case, col):
    M, at0, peak, dur = [], [], [], []
    for s, tr in DET[("cap40", case)].items():
        rr = DETROW[(s, case)]
        if rr["n_fail"] != "1":      # single-failure realisations (D is that failure)
            continue
        t0 = float(rr["first_fail"])
        x = tr["t"] - t0
        M.append(np.interp(S, x, tr[col], left=np.nan, right=np.nan))
        at0.append(np.interp(0.0, x, tr["pct"]))
        # peak of the failure episode (pct_max of the run also holds the initial transient)
        w = (x >= 0) & (x <= float(rr["max_fail_len"]) + 30.0)
        # failures cut by the end of the run (240 s) have no complete episode: D set to NaN
        complete = t0 + float(rr["max_fail_len"]) + 5.0 <= tr["t"][-1]
        peak.append(tr["pct"][w].max()); dur.append(float(rr["max_fail_len"]) if complete else np.nan)
    return np.array(M), np.array(at0), np.array(peak), np.array(dur)


fig, ax = plt.subplots(1, 3, figsize=(W, 2.95), layout="constrained", gridspec_kw={"width_ratios": [1.15, 1.0, 1.0]})
a, b, c = ax
for case in RISK:
    P, at0, pk, D = aligned(case, "pct")
    M, _, _, _ = aligned(case, "m_plant")
    lw = 1.5 if case in ("N2b", "N2_var_a0.02") else 1.0
    a.plot(S, np.nanmedian(P, axis=0), color=COL[case], lw=lw, label=LAB[case] if case != "N2b" else "N2 (sin capa de riesgo)")
    b.plot(S, np.nanmedian(M, axis=0), color=COL[case], lw=lw)
    c.scatter(D, pk - at0, s=7, color=COL[case], alpha=0.65, lw=0)
for z in (a, b):
    z.axvline(0, color=INK2, lw=0.7, ls=(0, (3, 2)))
a.axhline(Y_REF, color=MUTED, lw=0.8, ls=":")
a.set_title("PCT (mediana)"); a.set_xlabel("t − inicio de la falla [s]"); a.set_ylabel("PCT [K]")
b.set_title("Caudal (mediana)"); b.set_xlabel("t − inicio de la falla [s]"); b.set_ylabel("kg/s")
c.set_title("Calentamiento en la falla")
c.set_xlabel("duración de la falla D [s]")
c.set_ylabel("pico − PCT al inicio de la falla [K]")
h_, l_ = a.get_legend_handles_labels()
fig.legend(h_, l_, loc="outside lower center", ncol=3, fontsize=6.6, handlelength=1.8)
for z, s in zip(ax, "abc"):
    tag(z, s)
save(fig, "traj_risk_dynamics")
for case in RISK:
    _, at0, pk, _ = aligned(case, "pct")
    print(f"{case:20s} pre-cool med {np.median(Y_REF - at0):6.2f} K  heat-up med {np.median(pk - at0):6.2f} K"
          f"  peak med {np.median(pk):7.2f}")

# =================================================================================================
# 4. traj_hierarchy: what changes a closed-loop trajectory
# =================================================================================================
ORDER = [
    ("N2 − N3", "nominal, 40 kg/s", "N2b", "N3b"),
    ("N2 − N3", "fallas del DET, 40 kg/s", "N2b", "N3b"),
    ("N2 − N3", "extrema, actuador físico", "N2b", "N3b"),
    ("N1 − N3", "nominal, 40 kg/s", "N1b", "N3b"),
    ("N1 − N3", "fallas del DET, 40 kg/s", "N1b", "N3b"),
    ("N1 − N3", "extrema, actuador físico", "N1b", "N3b"),
    ("margen DSM (λ=10⁴) − N2", "fallas del DET, 40 kg/s", "N2_DSM_U1244_l1e4", "N2b"),
    ("sin sesgo − con sesgo (N2)", "extrema, actuador físico", "N2", "N2b"),
    ("pérdida esperada ρ=100 − N2", "fallas del DET, 40 kg/s", "N2_exp_r100", "N2b"),
    ("cuantil α=0.02 − N2", "fallas del DET, 40 kg/s", "N2_var_a0.02", "N2b"),
    ("Wilks 59 − N2", "fallas del DET, 40 kg/s", "N2_wilks59", "N2b"),
    ("PID de la tesis − N2", "nominal, 40 kg/s", "PID", "N2b"),
    ("PI del DET − N2", "extrema, actuador físico", "PI", "N2b"),
]
GROUPC = [BLUE] * 3 + [B300] * 3 + [AQUA] + ["#9aa9bb"] + [B700] * 3 + [ORANGE] * 2
fig, a = plt.subplots(figsize=(W, 3.9), layout="constrained")
rng = np.random.default_rng(0)
for i, ((lab, fam, x, y), col) in enumerate(zip(ORDER, GROUPC)):
    v = np.maximum(DIST[(fam, x, y)], RES)
    yi = len(ORDER) - 1 - i
    a.scatter(v, yi + rng.uniform(-0.22, 0.22, len(v)), s=6, color=col, alpha=0.55, lw=0)
    a.plot([np.median(v)] * 2, [yi - 0.32, yi + 0.32], color=INK, lw=1.4)
    a.text(1.6e2, yi, f"{np.median(v):.2g}", va="center", fontsize=6.4, color=INK2)
a.set_yticks(range(len(ORDER)))
a.set_yticklabels([f"{lab}  ·  {fam}" for lab, fam, _, _ in ORDER][::-1], fontsize=6.6)
a.set_xscale("log"); a.set_xlim(5e-5, 4e2)
a.axvspan(5e-5, RES, color=BAND, lw=0)
a.axvline(0.30, color=INK2, lw=0.8, ls=(0, (4, 2)))
a.text(0.32, len(ORDER) - 0.4, "banda de\nasentamiento\n0.30 K", fontsize=6, color=INK2, va="top")
a.text(1.6e2, len(ORDER) - 0.3, "mediana", fontsize=6, color=INK2, va="bottom")
a.set_xlabel(r"máximo en el tiempo de $|$PCT$_A-$PCT$_B|$ por realización [K]  (50 realizaciones; marca = mediana)")
a.set_title("Qué cambia una trayectoria de lazo cerrado")
for yb in (2.5, 5.5, 6.5, 7.5, 10.5):
    a.axhline(len(ORDER) - 1 - yb, color=GRID, lw=0.8)
save(fig, "traj_hierarchy")
