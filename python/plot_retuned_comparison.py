"""Retuned LJ with a moderated penalty: three-configuration comparison.

    A  baseline          u_max = 56.75,  lambda = 0,     LJ original
    C  over-penalised    u_max = 56.90,  lambda = 1e6,   LJ original
    D  retuned           u_max = 56.90,  lambda = 1e4,   LJ retuned

C is kept in the figures as the intermediate step: it is what exposed both the
optimizer defect and the over-weighted penalty, so showing A/C/D together makes
the correction legible rather than presenting D as if it had come first.

  figures/ret_pct_vs_tq          PCT, full range + zoom
  figures/ret_dsm_vs_tq          normalised DSM with the threshold
  figures/ret_u_vs_tq            applied ACS current with all bounds
  figures/ret_mdot_vs_tq         coolant mass flow (log)
  figures/ret_jterms_vs_tq       J_error / J_mdot / J_psm for C and D
  figures/ret_kpi_bars           closed-loop KPIs
  figures/ret_lj_quality         LJ optimality gap, original vs retuned
  figures/retuned_comparison     composed 4x2 overview
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

plt.rcParams.update({
    "font.family": "serif", "font.size": 10, "axes.labelsize": 11,
    "axes.titlesize": 10.5, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9, "lines.linewidth": 1.5, "figure.dpi": 150,
})

HERE = Path(r"external revision/Response to reviewers/evidence/phase2c_fixed_state_experiment")
OUT = HERE / "results_retuned_lambda1e4"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

TARGET = 1179.2711776
U_Y = 1244.0
SPAN = U_Y - TARGET
DSM_THR = 0.00237737
XTICKS = list(range(1160, 1191, 10))
XLIM = (1158, 1192)
ME = 5

CFG = {
    "A": dict(file="phase2c_fullloop_pilot_LJ_seed1.mat", log=None,
              lbl=r"A  baseline: $u_{\max}$=56.75, $\lambda$=0",
              short=r"A $\lambda$=0", col="#1f77b4", mrk="o"),
    "C": dict(file="widepen_umax56p90_LJ_seed1.mat", log="widepen_only.log",
              lbl=r"C  over-penalised: $u_{\max}$=56.90, $\lambda$=$10^{6}$, LJ original",
              short=r"C $\lambda$=$10^6$", col="#d62728", mrk="s"),
    "D": dict(file="retuned_lam1e4_umax56p90_LJ_seed1.mat", log="retuned_lam1e4.log",
              lbl=r"D  retuned: $u_{\max}$=56.90, $\lambda$=$10^{4}$, LJ retuned",
              short=r"D $\lambda$=$10^4$", col="#2ca02c", mrk="^"),
}
C_PSM, C_MDOT = "#ff7f0e", "#8c564b"

# ECM pump model (CurrentToMassFlow.m)
RHO, HD, WD, L, DTH = 10449.0, 0.515, 0.515, 1.0, 0.001
S_ = WD * HD
BM, DELTAP, FR = 1.0, 15000.0, 0.02
Q0_LM, Q0_W, KF = 0.118, 1.2, 1.0
R_LM = Q0_LM * WD / (HD * L); R_W = Q0_W * WD / (2 * DTH * L)
R_B = Q0_LM * (KF**2 * HD) / DTH
R_VER = (R_W * R_B) / (R_W + R_B); R_TOT = R_LM + R_VER
C_ECM = (BM * R_VER * HD**2) / (R_TOT * DTH)
B_ECM = (BM**2 * R_TOT * HD**2) / DTH
A_ECM = (FR * RHO * L * (WD + HD)) / (4 * S_**3)


def massflow(I):
    I = np.atleast_1d(np.asarray(I, float))
    disc = B_ECM**2 - 4 * A_ECM * (DELTAP - C_ECM * I)
    q = np.where(disc < 0, 0.0, (-B_ECM + np.sqrt(np.maximum(disc, 0))) / (2 * A_ECM))
    return RHO * np.maximum(q, 0.0)


def load_run(path):
    d = loadmat(path, squeeze_me=True, struct_as_record=False)
    m = d.get("metrics")
    kpi = {}
    if m is not None and hasattr(m, "_fieldnames"):
        for f in m._fieldnames:
            v = np.asarray(getattr(m, f))
            if v.size == 1:
                kpi[f] = float(v.ravel()[0])
    return dict(t=np.asarray(d["time_snapshots"], float).ravel(),
                pct=np.asarray(d["PCT"], float).ravel(),
                u=np.asarray(d["u_V2P_trajectory"], float).ravel(),
                rt=np.asarray(d["optimizer_runtimes"], float).ravel(), kpi=kpi)


def parse_j(path, n, t):
    num = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    re_j = re.compile(rf"J_error:\s*({num}),\s*J_mdot:\s*({num}),\s*J_psm:\s*({num})")
    re_q = re.compile(r"^\[LJ\]\s+t=([\d.]+)")
    tt, Je, Jm, Jp = [], [], [], []
    buf = []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            m = re_j.search(line)
            if m:
                buf.append([float(x) for x in m.groups()]); continue
            m = re_q.match(line)
            if m and buf:
                a = np.array(buf); k = int(np.argmin(a.sum(axis=1)))
                tt.append(float(m.group(1))); Je.append(a[k, 0])
                Jm.append(a[k, 1]); Jp.append(a[k, 2]); buf = []
    tt = np.asarray(tt); out = {k: np.full(n, np.nan) for k in ("J_error", "J_mdot", "J_psm")}
    for i, tq in enumerate(t):
        hit = np.isclose(tt, tq, atol=1e-6)
        if hit.any():
            k = int(np.argmax(hit))
            out["J_error"][i] = Je[k]; out["J_mdot"][i] = Jm[k]; out["J_psm"][i] = Jp[k]
    return out


def _ax(ax, xl, yl):
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_xlim(*XLIM); ax.set_xticks(XTICKS)
    ax.grid(True, alpha=0.3); ax.spines[["top", "right"]].set_visible(False)


def _bands(ax, t, mask, label=None):
    first = True
    for i, m in enumerate(mask):
        if not m: continue
        ax.axvspan(t[i] - 0.25, t[i] + 0.25, color=C_PSM, alpha=0.15, lw=0,
                   label=label if first else None)
        first = False


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf"); fig.savefig(FIG / f"{name}.png")
    plt.close(fig); print("OK", name)


# ------------------------------------------------------------------ data
R = {k: load_run(OUT / v["file"]) for k, v in CFG.items()}
n = min(len(R[k]["t"]) for k in R)
t = R["A"]["t"][:n]
MD = {k: massflow(R[k]["u"][:n]) for k in R}
DSMn = {k: (R[k]["pct"][:n] - TARGET) / SPAN for k in R}
J = {k: parse_j(OUT / CFG[k]["log"], n, t) for k in ("C", "D")}
actD = J["D"]["J_psm"] > 0
actC = J["C"]["J_psm"] > 0
print(f"queries {n}; penalty active: C {int(actC.sum())}, D {int(actD.sum())}")

# ---- 1. PCT ----
fig, axs = plt.subplots(1, 2, figsize=(9.8, 3.9))
for ax, zoom in zip(axs, (False, True)):
    _ax(ax, r"Query time $t_q$ (s)", "PCT (K)")
    _bands(ax, t, actD, r"D: $J_{psm}>0$" if not zoom else None)
    for k in ("A", "C", "D"):
        ax.plot(t, R[k]["pct"][:n], color=CFG[k]["col"], marker=CFG[k]["mrk"],
                markersize=3.8, markevery=ME, label=CFG[k]["lbl"] if not zoom else None)
    ax.axhline(TARGET, color="0.45", ls="--", lw=1.1,
               label=r"setpoint $T_p^*$" if not zoom else None)
    if zoom:
        ax.set_xlim(1161, 1192)
        lo = min(R[k]["pct"][2:n].min() for k in R); hi = max(R[k]["pct"][2:n].max() for k in R)
        ax.set_ylim(lo - 0.08*(hi-lo), hi + 0.08*(hi-lo))
        ax.set_title("zoom, transient excluded", fontsize=9.5, color="0.35")
    else:
        ax.set_ylim(1140, 1185); ax.legend(frameon=False, loc="lower right", fontsize=7.5)
        ax.set_title("full range", fontsize=9.5, color="0.35")
fig.tight_layout(); save(fig, "ret_pct_vs_tq")

# ---- 2. normalised DSM ----
fig, ax = plt.subplots(figsize=(7.4, 4.0))
_ax(ax, r"Query time $t_q$ (s)", r"Normalised DSM $=(T_p-T_p^*)/(U_y-T_p^*)$")
_bands(ax, t, actD, r"D: $J_{psm}>0$")
for k in ("A", "C", "D"):
    ax.plot(t, DSMn[k], color=CFG[k]["col"], marker=CFG[k]["mrk"], markersize=3.8,
            markevery=4, label=CFG[k]["short"])
ax.axhline(DSM_THR, color=C_PSM, ls="-.", lw=1.3, label=rf"threshold ${DSM_THR:.5f}$")
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlim(1161, 1192)
lo = min(DSMn[k][2:].min() for k in DSMn); hi = max(DSMn[k][2:].max() for k in DSMn)
ax.set_ylim(lo - 0.12*(hi-lo), hi + 0.12*(hi-lo))
ax.legend(frameon=False, loc="lower left", fontsize=8)
save(fig, "ret_dsm_vs_tq")

# ---- 3. applied current ----
fig, ax = plt.subplots(figsize=(7.4, 4.0))
_ax(ax, r"Query time $t_q$ (s)", "Applied ACS current (A)")
_bands(ax, t, actD, r"D: $J_{psm}>0$")
for k in ("A", "C", "D"):
    ax.step(t, R[k]["u"][:n], where="mid", color=CFG[k]["col"], label=CFG[k]["short"],
            alpha=0.9 if k != "C" else 0.55)
for y, lab in [(56.69, r"$u_{\min}$ (0.58 kg s$^{-1}$)"),
               (56.75, r"archived $u_{\max}$"), (56.90, r"widened $u_{\max}$")]:
    ax.axhline(y, color="0.45", ls="--", lw=1.0)
    ax.text(1158.3, y + 0.004, lab, fontsize=7.5, color="0.35")
ax.set_ylim(56.675, 56.928)
ax.legend(frameon=False, loc="center right", fontsize=8)
save(fig, "ret_u_vs_tq")

# ---- 4. mass flow ----
fig, ax = plt.subplots(figsize=(7.4, 4.0))
_ax(ax, r"Query time $t_q$ (s)", r"ACS coolant mass flow (kg s$^{-1}$)")
_bands(ax, t, actD, r"D: $J_{psm}>0$")
for k in ("A", "C", "D"):
    ax.step(t, MD[k], where="mid", color=CFG[k]["col"], label=CFG[k]["short"],
            alpha=0.9 if k != "C" else 0.55)
ax.set_yscale("log")
txt = "  ".join(rf"{CFG[k]['short']}: {MD[k].sum():.0f}" for k in ("A", "C", "D"))
ax.text(0.99, 0.04, r"$\Sigma\dot m$  " + txt, transform=ax.transAxes,
        ha="right", va="bottom", fontsize=8, color="0.25")
ax.legend(frameon=False, loc="upper right", fontsize=8)
save(fig, "ret_mdot_vs_tq")

# ---- 5. objective terms, C vs D ----
fig, axs = plt.subplots(1, 2, figsize=(9.8, 4.0), sharey=True)
for ax, k in zip(axs, ("C", "D")):
    _ax(ax, r"Query time $t_q$ (s)", "Objective contribution" if k == "C" else "")
    _bands(ax, t, J[k]["J_psm"] > 0, r"$J_{psm}>0$")
    ax.plot(t, J[k]["J_error"], color="#1f77b4", marker="o", markersize=3.4,
            markevery=ME, label=r"$J_{error}$")
    ax.plot(t, np.where(J[k]["J_psm"] > 0, J[k]["J_psm"], np.nan), color=C_PSM,
            marker="s", markersize=3.4, markevery=3, label=r"$J_{psm}$")
    ax.plot(t, np.where(J[k]["J_mdot"] > 0, J[k]["J_mdot"], np.nan), color=C_MDOT,
            marker="^", markersize=3.4, markevery=ME, label=r"$J_{\dot m}$")
    ax.set_yscale("log")
    ax.set_title(CFG[k]["lbl"], fontsize=9, color="0.3")
    ax.legend(frameon=False, loc="lower right", fontsize=8, ncol=3)
fig.tight_layout(); save(fig, "ret_jterms_vs_tq")

# ---- 6. KPI bars ----
KEYS = [("IAE", r"IAE (K$\cdot$s)", "{:.2f}"),
        ("steady_state_error", "Steady-state error (K)", "{:.4f}"),
        ("total_variation", "Control TV (A)", "{:.3f}"),
        ("min_safety_margin", "Min safety margin (K)", "{:.4f}")]
fig, axes = plt.subplots(1, 4, figsize=(11.4, 3.5))
for ax, (key, lab, fmt) in zip(axes, KEYS):
    vals = [R[k]["kpi"].get(key, np.nan) for k in ("A", "C", "D")]
    ax.bar([0, 1, 2], vals, color=[CFG[k]["col"] for k in ("A", "C", "D")],
           alpha=0.78, width=0.62, edgecolor="0.25", linewidth=0.7)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["A", "C", "D"])
    ax.set_ylabel(lab); ax.grid(True, axis="y", alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    if key == "min_safety_margin":
        ax.set_ylim(62.0, 63.2)
    for x, v in zip([0, 1, 2], vals):
        if np.isfinite(v):
            ax.annotate(fmt.format(v), xy=(x, v), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=8.5, color="0.2")
fig.tight_layout(); save(fig, "ret_kpi_bars")

# ---- 7. LJ optimality gap ----
try:
    vs = loadmat(next(OUT.glob("verify_lj_suboptimality_*.mat")), squeeze_me=True,
                 struct_as_record=False)["V"]
    vr = loadmat(next(OUT.glob("verify_lj_retuned_*.mat")), squeeze_me=True,
                 struct_as_record=False)["V"]
    tq = np.asarray(vs.tq, float).ravel()
    gap_o = 100*np.asarray(vs.gap_rel, float).ravel()
    gap_r = 100*np.asarray(vr.gap_rel, float).ravel()
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    x = np.arange(len(tq)); w = 0.36
    ax.bar(x - w/2, np.abs(gap_o), w, color="#d62728", alpha=0.8,
           edgecolor="0.25", linewidth=0.6, label="LJ original")
    ax.bar(x + w/2, np.abs(gap_r), w, color="#2ca02c", alpha=0.8,
           edgecolor="0.25", linewidth=0.6, label="LJ retuned")
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels([f"{v:.0f}" for v in tq])
    ax.set_xlabel(r"Query time $t_q$ (s)")
    ax.set_ylabel(r"$|$relative cost gap vs grid$|$  (%)")
    ax.axhline(100, color="0.45", ls=":", lw=1.0)
    ax.grid(True, axis="y", alpha=0.3); ax.spines[["top", "right"]].set_visible(False)
    for xi, v in zip(x - w/2, np.abs(gap_o)):
        ax.annotate(f"{v:.0f}", xy=(xi, v), xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=7.5, color="0.25")
    for xi, v in zip(x + w/2, np.abs(gap_r)):
        ax.annotate(f"{v:.0f}", xy=(xi, v), xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=7.5, color="0.25")
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("Optimizer quality on the queries where the drift appeared",
                 fontsize=9.5, color="0.35")
    save(fig, "ret_lj_quality")
except (StopIteration, KeyError, AttributeError) as e:
    print("  (LJ quality figure skipped:", e, ")")

# ---- composed 4x2 ----
plt.rcParams.update({"font.size": 10, "axes.labelsize": 10.5, "legend.fontsize": 7.5})
fig, axes = plt.subplots(4, 2, figsize=(7.8, 12.4))

ax = axes[0, 0]; _ax(ax, r"$t_q$ (s)", "PCT (K)"); _bands(ax, t, actD)
for k in ("A", "C", "D"):
    ax.plot(t, R[k]["pct"][:n], color=CFG[k]["col"], marker=CFG[k]["mrk"],
            markersize=3.2, markevery=ME, label=CFG[k]["short"])
ax.axhline(TARGET, color="0.45", ls="--", lw=1.0); ax.set_ylim(1140, 1185)
ax.legend(frameon=False, loc="lower right")

ax = axes[0, 1]; _ax(ax, r"$t_q$ (s)", "PCT (K), zoom"); _bands(ax, t, actD)
for k in ("A", "C", "D"):
    ax.plot(t, R[k]["pct"][:n], color=CFG[k]["col"], marker=CFG[k]["mrk"],
            markersize=3.2, markevery=4)
ax.axhline(TARGET, color="0.45", ls="--", lw=1.0); ax.set_xlim(1161, 1192)
_lo = min(R[k]["pct"][2:n].min() for k in R); _hi = max(R[k]["pct"][2:n].max() for k in R)
ax.set_ylim(_lo - 0.08*(_hi-_lo), _hi + 0.08*(_hi-_lo))

ax = axes[1, 0]; _ax(ax, r"$t_q$ (s)", "Normalised DSM"); _bands(ax, t, actD)
for k in ("A", "C", "D"):
    ax.plot(t, DSMn[k], color=CFG[k]["col"], marker=CFG[k]["mrk"], markersize=3.2, markevery=4)
ax.axhline(DSM_THR, color=C_PSM, ls="-.", lw=1.1); ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlim(1161, 1192)
_lo = min(DSMn[k][2:].min() for k in DSMn); _hi = max(DSMn[k][2:].max() for k in DSMn)
ax.set_ylim(_lo - 0.12*(_hi-_lo), _hi + 0.12*(_hi-_lo))

ax = axes[1, 1]; _ax(ax, r"$t_q$ (s)", "Applied current (A)"); _bands(ax, t, actD)
for k in ("A", "C", "D"):
    ax.step(t, R[k]["u"][:n], where="mid", color=CFG[k]["col"], alpha=0.9 if k != "C" else 0.5)
for y in (56.69, 56.75, 56.90): ax.axhline(y, color="0.45", ls="--", lw=0.9)
ax.set_ylim(56.675, 56.928)

ax = axes[2, 0]; _ax(ax, r"$t_q$ (s)", r"$\dot m$ (kg s$^{-1}$)"); _bands(ax, t, actD)
for k in ("A", "C", "D"):
    ax.step(t, MD[k], where="mid", color=CFG[k]["col"], alpha=0.9 if k != "C" else 0.5)
ax.set_yscale("log")

ax = axes[2, 1]; _ax(ax, r"$t_q$ (s)", "Objective terms, D"); _bands(ax, t, actD)
ax.plot(t, J["D"]["J_error"], color="#1f77b4", marker="o", markersize=3, markevery=ME, label=r"$J_{error}$")
ax.plot(t, np.where(J["D"]["J_psm"] > 0, J["D"]["J_psm"], np.nan), color=C_PSM,
        marker="s", markersize=3, markevery=3, label=r"$J_{psm}$")
ax.plot(t, np.where(J["D"]["J_mdot"] > 0, J["D"]["J_mdot"], np.nan), color=C_MDOT,
        marker="^", markersize=3, markevery=ME, label=r"$J_{\dot m}$")
ax.set_yscale("log"); ax.legend(frameon=False, loc="lower right", ncol=3)

ax = axes[3, 0]; _ax(ax, r"$t_q$ (s)", "Optimizer runtime (s)"); _bands(ax, t, actD)
for k in ("A", "C", "D"):
    ax.plot(t, R[k]["rt"][:n], color=CFG[k]["col"], marker=CFG[k]["mrk"],
            markersize=3.2, markevery=ME)
ax.axhline(120, color="0.45", ls=":", lw=1.0); ax.set_ylim(0, 900)

ax = axes[3, 1]
vals_iae = [R[k]["kpi"].get("IAE", np.nan) for k in ("A", "C", "D")]
vals_tv = [R[k]["kpi"].get("total_variation", np.nan) for k in ("A", "C", "D")]
x = np.arange(3); w = 0.36
ax.bar(x - w/2, vals_iae, w, color="0.55", alpha=0.85, edgecolor="0.25",
       linewidth=0.6, label=r"IAE (K$\cdot$s)")
ax.bar(x + w/2, vals_tv, w, color=C_MDOT, alpha=0.85, edgecolor="0.25",
       linewidth=0.6, label="Control TV (A)")
ax.set_xticks(x); ax.set_xticklabels(["A", "C", "D"])
ax.set_ylabel("KPI value"); ax.grid(True, axis="y", alpha=0.3)
ax.spines[["top", "right"]].set_visible(False)
for xi, v in zip(x - w/2, vals_iae):
    ax.annotate(f"{v:.1f}", xy=(xi, v), xytext=(0, 3), textcoords="offset points",
                ha="center", fontsize=7.5)
for xi, v in zip(x + w/2, vals_tv):
    ax.annotate(f"{v:.2f}", xy=(xi, v), xytext=(0, 3), textcoords="offset points",
                ha="center", fontsize=7.5)
ax.legend(frameon=False, fontsize=7.5)

for ax, lab in zip(axes.ravel(), ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]):
    ax.text(0.02, 0.98, lab, transform=ax.transAxes, va="top", ha="left",
            fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "retuned_comparison.pdf"); fig.savefig(FIG / "retuned_comparison.png")
plt.close(fig); print("OK retuned_comparison")

# ---- summary ----
def slope(x, y):
    return np.polyfit(x, y, 1)[0] * 1000.0


print("\n" + "=" * 82)
print("THREE-CONFIGURATION COMPARISON")
print("=" * 82)
print(f"  {'metric':<30}{'A lam=0':>16}{'C lam=1e6':>16}{'D lam=1e4':>16}")
print("  " + "-" * 78)
for key, lab, fmt in [("peak_PCT", "Peak PCT (K)", "{:.4f}"),
                      ("min_safety_margin", "Min safety margin (K)", "{:.4f}"),
                      ("IAE", "IAE (K.s)", "{:.4f}"),
                      ("steady_state_error", "Steady-state error (K)", "{:+.4f}"),
                      ("total_variation", "Control TV (A)", "{:.4f}"),
                      ("max_control_rate", "Peak rate (A/s)", "{:.4f}"),
                      ("settling_time", "Settling time (s)", "{:.1f}")]:
    row = f"  {lab:<30}"
    for k in ("A", "C", "D"):
        v = R[k]["kpi"].get(key, np.nan)
        row += f"{fmt.format(v):>16}" if np.isfinite(v) else f"{'inf':>16}"
    print(row)
row = f"  {'final PCT (K)':<30}"
for k in ("A", "C", "D"): row += f"{R[k]['pct'][n-1]:>16.4f}"
print(row)
row = f"  {'final error (K)':<30}"
for k in ("A", "C", "D"): row += f"{R[k]['pct'][n-1]-TARGET:>+16.4f}"
print(row)
row = f"  {'final slope (mK/s)':<30}"
for k in ("A", "C", "D"): row += f"{slope(t[-10:], R[k]['pct'][n-10:n]):>+16.3f}"
print(row)
row = f"  {'total coolant (kg/s)':<30}"
for k in ("A", "C", "D"): row += f"{MD[k].sum():>16.2f}"
print(row)
row = f"  {'interior control (%)':<30}"
for k in ("A", "C", "D"):
    u = R[k]["u"][:n]
    row += f"{100*np.mean((u > 56.6901) & (u < 56.8999)):>16.1f}"
print(row)
row = f"  {'penalty-active queries':<30}{'--':>16}"
row += f"{int(actC.sum()):>16d}{int(actD.sum()):>16d}"
print(row)

np.savez(OUT / "retuned_comparison.npz",
         t=t, **{f"pct_{k}": R[k]["pct"][:n] for k in R},
         **{f"u_{k}": R[k]["u"][:n] for k in R},
         **{f"mdot_{k}": MD[k] for k in MD},
         **{f"dsm_{k}": DSMn[k] for k in DSMn},
         **{f"rt_{k}": R[k]["rt"][:n] for k in R},
         J_error_C=J["C"]["J_error"], J_psm_C=J["C"]["J_psm"], J_mdot_C=J["C"]["J_mdot"],
         J_error_D=J["D"]["J_error"], J_psm_D=J["D"]["J_psm"], J_mdot_D=J["D"]["J_mdot"],
         dsm_thr=DSM_THR, target=TARGET, U_y=U_Y)
print(f"\n  saved {OUT / 'retuned_comparison.npz'}")
print(f"  figures in {FIG}")
