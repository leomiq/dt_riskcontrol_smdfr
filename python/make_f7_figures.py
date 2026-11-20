"""F7 figures: certificate of global optimality of closed-loop MPC queries (status report).

Reads results/f7_cert/<campaign>_<seed>.csv and *_hist.csv (julia/f7_certify.jl) and writes
  docs/figures/fig_f7_certificate   (a) certified gap (UB − LB)/UB per query, root vs final;
                                    (b) certified excess over the global optimum of the loop's
                                        N3, N2 and N1 decisions
  docs/figures/fig_f7_convergence   gap vs B&B time for every query
and results/f7_cert/summary.md (the numbers quoted in docs/estado/11l_certificacion.tex).
Usage: .venv\\Scripts\\python python/make_f7_figures.py
"""
import csv
import glob
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "figures")
RES = os.path.join(ROOT, "results", "f7_cert")
font_manager.fontManager.addfont(os.path.join(ROOT, "assets", "fonts", "InterVariable.ttf"))
BLUE, B300, B650, ORANGE, AQUA = "#2a78d6", "#6da7ec", "#104281", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df", "#ecebe7"
FAMC = {"nom": AQUA, "det": BLUE, "ext": ORANGE}
FAML = {"nom": "nominal (40 kg/s)", "det": "fallas del DET (40 kg/s)", "ext": "extrema (actuador físico)"}
plt.rcParams.update({
    "font.family": "Inter Variable", "font.size": 8.5, "axes.titlesize": 8.8, "axes.labelsize": 8.3,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 7, "figure.dpi": 150, "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
W = 6.3
FLOOR = 1e-7


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)
    print("wrote", name)


rows, hist = [], defaultdict(list)
for f in sorted(glob.glob(os.path.join(RES, "*.csv"))):
    b = os.path.basename(f)
    if b.startswith("knitro") or b == "summary.csv":
        continue
    with open(f, newline="") as fh:
        R = list(csv.DictReader(fh))
    if b.endswith("_hist.csv"):
        for r in R:
            hist[(r["key"], int(r["k"]))].append((float(r["t"]), float(r["LB"]), float(r["UB"])))
    elif not b.endswith("_local.csv"):
        rows += R
for r in rows:
    for c in ("UB", "LB", "root_LB", "gap_rel", "J_loopN3", "J_N1", "J_N2", "J_ipopt_ms", "t_bnb", "rms_gap_K",
              "cert_loopN3_rel", "cert_N1_rel", "cert_N2_rel", "raw_violation"):
        r[c] = float(r[c])
    r["k"] = int(r["k"])
    r["root_gap"] = (r["UB"] - r["root_LB"]) / r["UB"]
rows.sort(key=lambda r: ({"nom": 0, "det": 1, "ext": 2}[r["campaign"]], r["UB"]))
n = len(rows)
print(n, "queries")

# ---------------------------------------------------------------- summary
lines = ["# F7 — certificate of global optimality, closed-loop queries", "",
         f"{n} queries (julia/f7_certify.jl). gap = (UB − LB)/UB; excess of a decision = (J − LB)/UB.", ""]
lines.append("| familia | consultas | certificadas (≤0.01 %) | brecha mediana | brecha máx. | brecha raíz mediana | "
             "exceso N3 máx. | exceso N2 máx. | exceso N1 máx. | RMS equiv. máx. [K] | multistart mejora UB |")
lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
for fam in ("nom", "det", "ext"):
    R = [r for r in rows if r["campaign"] == fam]
    if not R:
        continue
    g = np.array([r["gap_rel"] for r in R])
    better = sum(r["J_ipopt_ms"] < r["J_loopN3"] - 1e-9 * max(1, r["J_loopN3"]) for r in R)
    lines.append(f"| {fam} | {len(R)} | {sum(r['bnb_status'].startswith('certified') or r['gap_rel'] <= 1e-4 for r in R)} | "
                 f"{np.median(g):.2e} | {g.max():.2e} | {np.median([r['root_gap'] for r in R]):.2e} | "
                 f"{max(r['cert_loopN3_rel'] for r in R):.2e} | {max(r['cert_N2_rel'] for r in R):.2e} | "
                 f"{max(r['cert_N1_rel'] for r in R):.2e} | {max(r['rms_gap_K'] for r in R):.2e} | {better} |")
with open(os.path.join(RES, "summary.md"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))

# ---------------------------------------------------------------- fig 1
fig, ax = plt.subplots(1, 2, figsize=(W, 3.25), layout="constrained", sharey=False)
a, b = ax
x = np.arange(n)
for i, r in enumerate(rows):
    col = FAMC[r["campaign"]]
    a.plot([i, i], [max(r["gap_rel"], FLOOR), max(r["root_gap"], FLOOR)], color=col, lw=0.6, alpha=0.5)
    a.scatter(i, max(r["root_gap"], FLOOR), s=9, facecolor="white", edgecolor=col, lw=0.7, zorder=3)
    a.scatter(i, max(r["gap_rel"], FLOOR), s=11, color=col, zorder=4)
a.axhline(1e-4, color=INK2, lw=0.8, ls=(0, (4, 2)))
a.text(n - 0.5, 1.15e-4, "tolerancia 0.01 %", ha="right", va="bottom", fontsize=6.3, color=INK2)
a.set_yscale("log"); a.set_xticks([])
a.set_xlabel("consultas (ordenadas por familia y costo)")
a.set_ylabel("brecha certificada (UB − LB)/UB")
a.set_title("(a)  Raíz (círculo) y final del B&B (punto)")
for fam in ("nom", "det", "ext"):
    a.scatter([], [], color=FAMC[fam], s=11, label=FAML[fam])
fig.legend(*a.get_legend_handles_labels(), loc="outside lower center", ncol=3, fontsize=6.6)
for dx, key, lab, col in ((-0.25, "cert_N1_rel", "N1", B300), (0.0, "cert_N2_rel", "N2", BLUE),
                          (0.25, "cert_loopN3_rel", "N3 (lazo)", B650)):
    v = np.maximum([r[key] for r in rows], FLOOR)
    b.scatter(x + dx, v, s=7, color=col, label=lab, alpha=0.85, lw=0)
b.set_yscale("log"); b.set_xticks([])
b.axhline(1e-4, color=INK2, lw=0.8, ls=(0, (4, 2)))
b.set_xlabel("consultas (mismo orden)")
b.set_ylabel(r"exceso garantizado $(J(u)-LB)/UB$")
b.set_title("(b)  Cuánto puede distar cada decisión del óptimo global")
b.set_ylim(top=6e-2); b.legend(loc="upper center", ncol=3, fontsize=6.5)
save(fig, "fig_f7_certificate")

# ---------------------------------------------------------------- fig 2
fig, a = plt.subplots(figsize=(W, 2.6), layout="constrained")
for r in rows:
    H = hist.get((r["key"], r["k"]))
    if not H:
        continue
    H = np.array(H)
    g = np.maximum((H[:, 2] - H[:, 1]) / H[:, 2], FLOOR)
    a.plot(np.maximum(H[:, 0], 1.0), g, color=FAMC[r["campaign"]], lw=0.7, alpha=0.7)
a.axhline(1e-4, color=INK2, lw=0.8, ls=(0, (4, 2)))
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel("tiempo de B&B [s]"); a.set_ylabel("(UB − LB)/UB")
a.set_title("Cierre de la brecha en el branch and bound espacial")
for fam in ("nom", "det", "ext"):
    a.plot([], [], color=FAMC[fam], label=FAML[fam])
a.legend(loc="lower left", fontsize=6.5)
save(fig, "fig_f7_convergence")
