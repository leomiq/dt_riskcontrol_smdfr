"""Full per-indicator figure set for the DSM-penalty comparison.

Mirrors the archived Phase 2C set (plot_phase2c_trajectories.py) one-for-one,
but the two series are now CONFIGURATIONS rather than solvers:

    BASE  archived bounds  u_max = 56.75 A,  lambda = 0
    PEN   widened pump     u_max = 56.90 A,  lambda = 1e6,  dsm_thr = 0.00237737

  figures/set_pct_vs_tq.{pdf,png}       PCT(t_q)
  figures/set_dsm_vs_tq.{pdf,png}       DSM(t_q) in K
  figures/set_dsmn_vs_tq.{pdf,png}      normalised DSM with the threshold
  figures/set_u_vs_tq.{pdf,png}         applied current with all three bounds
  figures/set_mdot_vs_tq.{pdf,png}      coolant mass flow (log)
  figures/set_jterms_vs_tq.{pdf,png}    J_error / J_mdot / J_psm / J
  figures/set_fopt_vs_tq.{pdf,png}      objective f*(t_q)
  figures/set_runtime_vs_tq.{pdf,png}   optimizer runtime
  figures/set_kpi_bars.{pdf,png}        closed-loop KPI bars
  figures/penalty_full_comparison.{pdf,png}  composed 4x2 overview

Caveat carried in the captions: the comparison varies TWO factors at once (pump
saturation limit and penalty weight). The fixed-state check
(verify_widened_bounds.m) is what licenses reading the trajectory change as the
penalty acting: at u_max = 56.90 the penalised optimum becomes interior
(56.7425 A) instead of pinned at the bound, which the archived bounds prevented.
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
OUT = HERE / "results_penalty_saturation"    # baseline vs penalty + new saturation
RES = OUT                                    # realization MATs and diaries live here
WORK = OUT
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

TARGET = 1179.2711776
U_Y = 1244.0
SPAN = U_Y - TARGET
DSM_THR = 0.00237737
XTICKS = list(range(1160, 1191, 10))
XLIM = (1158, 1192)
ME = 5

LBL = {"BASE": r"baseline: $u_{\max}$=56.75, $\lambda$=0",
       "PEN": r"penalised: $u_{\max}$=56.90, $\lambda$=$10^{6}$"}
COL = {"BASE": "#1f77b4", "PEN": "#d62728"}
MRK = {"BASE": "o", "PEN": "s"}
C_PSM, C_MDOT = "#ff7f0e", "#2ca02c"

# ECM pump model, evaluated directly (the widened run leaves the old table)
RHO, HD, WD, L, DTH = 10449.0, 0.515, 0.515, 1.0, 0.001
S_ = WD * HD
BM, DELTAP, FR = 1.0, 15000.0, 0.02
Q0_LM, Q0_W, KF = 0.118, 1.2, 1.0
R_LM = Q0_LM * WD / (HD * L)
R_W = Q0_W * WD / (2 * DTH * L)
R_B = Q0_LM * (KF**2 * HD) / DTH
R_VER = (R_W * R_B) / (R_W + R_B)
R_TOT = R_LM + R_VER
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
    m = d.get("metrics", None)
    kpi = {}
    if m is not None and hasattr(m, "_fieldnames"):
        for f in m._fieldnames:
            v = getattr(m, f)
            if np.isscalar(v) or (hasattr(v, "size") and v.size == 1):
                kpi[f] = float(np.asarray(v).ravel()[0])
    return {
        "t": np.asarray(d["time_snapshots"], float).ravel(),
        "pct": np.asarray(d["PCT"], float).ravel(),
        "u": np.asarray(d["u_V2P_trajectory"], float).ravel(),
        "rt": np.asarray(d["optimizer_runtimes"], float).ravel(),
        "kpi": kpi,
    }


def parse_j(logpath):
    """Per-query optimal J decomposition and the per-pass f_best sequence."""
    num = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    re_j = re.compile(rf"J_error:\s*({num}),\s*J_mdot:\s*({num}),\s*J_psm:\s*({num})")
    re_q = re.compile(r"^\[LJ\]\s+t=([\d.]+)")
    re_p = re.compile(rf"Pass\s+\d+\s+complete\.\s+f_best0?\s*=\s*({num})")
    out = {k: [] for k in ["t", "J_error", "J_mdot", "J_psm", "J", "f_opt"]}
    buf, passes = [], []
    with open(logpath, "r", errors="ignore") as fh:
        for line in fh:
            m = re_j.search(line)
            if m:
                buf.append([float(x) for x in m.groups()])
                continue
            m = re_p.search(line)
            if m:
                passes.append(float(m.group(1)))
                continue
            m = re_q.match(line)
            if m and buf:
                a = np.array(buf)
                k = int(np.argmin(a.sum(axis=1)))
                out["t"].append(float(m.group(1)))
                out["J_error"].append(a[k, 0])
                out["J_mdot"].append(a[k, 1])
                out["J_psm"].append(a[k, 2])
                out["J"].append(a[k].sum())
                out["f_opt"].append(passes[-1] if passes else np.nan)
                buf, passes = [], []
    return {k: np.asarray(v, float) for k, v in out.items()}


def _ax(ax, xl, yl):
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    ax.set_xlim(*XLIM); ax.set_xticks(XTICKS)
    ax.grid(True, alpha=0.3); ax.spines[["top", "right"]].set_visible(False)


def _bands(ax, t, mask, label=None):
    first = True
    for i, m in enumerate(mask):
        if not m:
            continue
        ax.axvspan(t[i] - 0.25, t[i] + 0.25, color=C_PSM, alpha=0.15, lw=0,
                   label=label if first else None)
        first = False


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf"); fig.savefig(FIG / f"{name}.png")
    plt.close(fig); print("OK", name)


# ------------------------------------------------------------------ data
base = load_run(RES / "phase2c_fullloop_pilot_LJ_seed1.mat")
pen = load_run(WORK / "widepen_umax56p90_LJ_seed1.mat")
Jb = parse_j(RES / "penalised_run3.log")     # archived-bounds penalised diary
Jp = parse_j(WORK / "widepen_only.log")      # widened diary

n = min(len(base["t"]), len(pen["t"]))
t = base["t"][:n]
D = {"BASE": base, "PEN": pen}
MD = {k: massflow(D[k]["u"][:n]) for k in D}
DSMn = {k: (D[k]["pct"][:n] - TARGET) / SPAN for k in D}

Jm = {k: np.full(n, np.nan) for k in ["J_error", "J_mdot", "J_psm", "J", "f_opt"]}
for i, tq in enumerate(t):
    hit = np.isclose(Jp["t"], tq, atol=1e-6)
    if hit.any():
        k = int(np.argmax(hit))
        for key in Jm:
            Jm[key][i] = Jp[key][k]
active = Jm["J_psm"] > 0

fb = np.full(n, np.nan)
for i, tq in enumerate(t):
    hit = np.isclose(Jb["t"], tq, atol=1e-6)
    if hit.any():
        fb[i] = Jb["f_opt"][int(np.argmax(hit))]

print(f"queries {n}; J_psm active in {int(np.nansum(active))}")

# ---- 1. PCT ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_ax(ax, r"Query time $t_q$ (s)", "PCT (K)")
_bands(ax, t, active, r"$J_{psm}>0$")
for k in D:
    ax.plot(t, D[k]["pct"][:n], color=COL[k], marker=MRK[k], markersize=4,
            markevery=ME, label=LBL[k])
ax.axhline(TARGET, color="0.45", ls="--", lw=1.1, label=r"setpoint $T_p^*=1179.27$ K")
ax.set_ylim(1140, 1185)
ax.legend(frameon=False, loc="lower right", fontsize=8)
save(fig, "set_pct_vs_tq")

# ---- 2. DSM in K ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_ax(ax, r"Query time $t_q$ (s)", "DSM (K)")
_bands(ax, t, active, r"$J_{psm}>0$")
for k in D:
    ax.plot(t, U_Y - D[k]["pct"][:n], color=COL[k], marker=MRK[k], markersize=4,
            markevery=ME, label=LBL[k])
ax.set_ylim(60, 100)
ax.annotate("ACS activation transient", xy=(1160.5, 97.5), xytext=(1165, 93),
            arrowprops=dict(arrowstyle="->", color="0.3", lw=1.0), fontsize=8.5, color="0.25")
ax.legend(frameon=False, loc="lower right", fontsize=8)
save(fig, "set_dsm_vs_tq")

# ---- 3. normalised DSM (zoom past the transient) ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_ax(ax, r"Query time $t_q$ (s)", r"Normalised DSM $=(T_p-T_p^*)/(U_y-T_p^*)$")
_bands(ax, t, active, r"$J_{psm}>0$")
for k in D:
    ax.plot(t, DSMn[k], color=COL[k], marker=MRK[k], markersize=4, markevery=4, label=LBL[k])
ax.axhline(DSM_THR, color=C_PSM, ls="-.", lw=1.3,
           label=rf"threshold ${DSM_THR:.5f}$")
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlim(1161, 1192)
lo = min(DSMn["BASE"][2:].min(), DSMn["PEN"][2:].min())
hi = max(DSMn["BASE"][2:].max(), DSMn["PEN"][2:].max())
ax.set_ylim(lo - 0.12*(hi-lo), hi + 0.12*(hi-lo))
ax.legend(frameon=False, loc="lower left", fontsize=8)
save(fig, "set_dsmn_vs_tq")

# ---- 4. applied current ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_ax(ax, r"Query time $t_q$ (s)", "Applied ACS current (A)")
_bands(ax, t, active, r"$J_{psm}>0$")
for k in D:
    ax.step(t, D[k]["u"][:n], where="mid", color=COL[k], label=LBL[k])
for y, lab in [(56.69, r"$u_{\min}$ (pump onset 56.6872 A)"),
               (56.75, r"archived $u_{\max}$"),
               (56.90, r"widened $u_{\max}$")]:
    ax.axhline(y, color="0.45", ls="--", lw=1.0)
    ax.text(1158.3, y + 0.004, lab, fontsize=7.5, color="0.35")
ax.set_ylim(56.675, 56.928)
ax.legend(frameon=False, loc="center right", fontsize=8)
save(fig, "set_u_vs_tq")

# ---- 5. coolant mass flow ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_ax(ax, r"Query time $t_q$ (s)", r"ACS coolant mass flow (kg s$^{-1}$)")
_bands(ax, t, active, r"$J_{psm}>0$")
for k in D:
    ax.step(t, MD[k], where="mid", color=COL[k], label=LBL[k])
ax.set_yscale("log")
ax.text(0.99, 0.05,
        rf"$\Sigma\dot m$: {MD['BASE'].sum():.1f} $\to$ {MD['PEN'].sum():.1f} kg s$^{{-1}}$"
        rf"  (${MD['PEN'].sum()/MD['BASE'].sum():.2f}\times$)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color="0.25")
ax.legend(frameon=False, loc="lower right", fontsize=8)
save(fig, "set_mdot_vs_tq")

# ---- 6. objective decomposition (penalised run) ----
fig, ax = plt.subplots(figsize=(7.2, 4.2))
_ax(ax, r"Query time $t_q$ (s)", "Objective contribution")
_bands(ax, t, active, r"$J_{psm}>0$")
ax.plot(t, Jm["J_error"], color=COL["BASE"], marker="o", markersize=3.6, markevery=ME,
        label=r"$J_{error}$ (tracking)")
ax.plot(t, np.where(Jm["J_psm"] > 0, Jm["J_psm"], np.nan), color=C_PSM, marker="s",
        markersize=3.6, markevery=3, label=r"$J_{psm}$ (DSM penalty)")
ax.plot(t, np.where(Jm["J_mdot"] > 0, Jm["J_mdot"], np.nan), color=C_MDOT, marker="^",
        markersize=3.6, markevery=ME, label=r"$J_{\dot m}$ (flow rate)")
ax.plot(t, Jm["J"], color="0.25", ls=":", lw=1.3, label=r"$J$ (total)")
ax.set_yscale("log")
ax.legend(frameon=False, loc="lower right", ncol=2, fontsize=8)
save(fig, "set_jterms_vs_tq")

# ---- 7. objective f* ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_ax(ax, r"Query time $t_q$ (s)", r"Objective $f^*$")
_bands(ax, t, active, r"$J_{psm}>0$")
ax.plot(t, fb, color=COL["BASE"], marker="o", markersize=4, markevery=ME,
        label=r"archived bounds, $\lambda=10^{6}$ (for reference)")
ax.plot(t, Jm["f_opt"], color=COL["PEN"], marker="s", markersize=4, markevery=ME,
        label=r"widened pump, $\lambda=10^{6}$")
ax.set_yscale("log")
ax.legend(frameon=False, loc="upper right", fontsize=8)
ax.set_title(r"$f^*$ includes the penalty, so it is not comparable to the $\lambda=0$ baseline",
             fontsize=8.5, color="0.35")
save(fig, "set_fopt_vs_tq")

# ---- 8. optimizer runtime ----
fig, ax = plt.subplots(figsize=(7.2, 4.0))
_ax(ax, r"Query time $t_q$ (s)", "Optimizer runtime (s)")
_bands(ax, t, active, r"$J_{psm}>0$")
for k in D:
    ax.plot(t, D[k]["rt"][:n], color=COL[k], marker=MRK[k], markersize=4,
            markevery=ME, label=LBL[k])
ax.axhline(120, color="0.45", ls=":", lw=1.1, label="120 s response-time reference")
ax.set_ylim(0, 900)
ax.legend(frameon=False, loc="upper left", fontsize=8)
save(fig, "set_runtime_vs_tq")

# ---- 9. closed-loop KPI bars ----
KPIS = [("peak_PCT", "Peak PCT (K)", "{:.4f}", None),
        ("min_safety_margin", "Min safety margin (K)", "{:.4f}", (62.0, 63.2)),
        ("IAE", r"IAE (K$\cdot$s)", "{:.2f}", None),
        ("total_variation", "Control TV (A)", "{:.3f}", None)]
fig, axes = plt.subplots(1, 4, figsize=(11.0, 3.4))
for ax, (key, lab, fmt, ylim) in zip(axes, KPIS):
    vals = [D[k]["kpi"].get(key, np.nan) for k in ["BASE", "PEN"]]
    bars = ax.bar([0, 1], vals, color=[COL["BASE"], COL["PEN"]], alpha=0.75,
                  width=0.6, edgecolor="0.25", linewidth=0.7)
    ax.set_xticks([0, 1]); ax.set_xticklabels([r"$\lambda$=0" + "\n56.75", r"$\lambda$=$10^6$" + "\n56.90"],
                                              fontsize=8)
    ax.set_ylabel(lab); ax.grid(True, axis="y", alpha=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    if ylim: ax.set_ylim(*ylim)
    for x, v in zip([0, 1], vals):
        if np.isfinite(v):
            ax.annotate(fmt.format(v), xy=(x, v), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=8.5, color="0.2")
fig.tight_layout()
save(fig, "set_kpi_bars")

# ---- composed 4x2 ----
plt.rcParams.update({"font.size": 10, "axes.labelsize": 10.5, "legend.fontsize": 7.5})
fig, axes = plt.subplots(4, 2, figsize=(7.8, 12.4))

ax = axes[0, 0]; _ax(ax, r"$t_q$ (s)", "PCT (K)"); _bands(ax, t, active)
for k in D: ax.plot(t, D[k]["pct"][:n], color=COL[k], marker=MRK[k], markersize=3.2, markevery=ME, label=LBL[k])
ax.axhline(TARGET, color="0.45", ls="--", lw=1.0); ax.set_ylim(1140, 1185)
ax.legend(frameon=False, loc="lower right")

ax = axes[0, 1]; _ax(ax, r"$t_q$ (s)", "PCT (K), zoom"); _bands(ax, t, active)
for k in D: ax.plot(t, D[k]["pct"][:n], color=COL[k], marker=MRK[k], markersize=3.2, markevery=4)
ax.axhline(TARGET, color="0.45", ls="--", lw=1.0); ax.set_xlim(1161, 1192)
_lo = min(D['BASE']["pct"][2:n].min(), D['PEN']["pct"][2:n].min())
_hi = max(D['BASE']["pct"][2:n].max(), D['PEN']["pct"][2:n].max())
ax.set_ylim(_lo - 0.08*(_hi-_lo), _hi + 0.08*(_hi-_lo))

ax = axes[1, 0]; _ax(ax, r"$t_q$ (s)", "Normalised DSM"); _bands(ax, t, active)
for k in D: ax.plot(t, DSMn[k], color=COL[k], marker=MRK[k], markersize=3.2, markevery=4)
ax.axhline(DSM_THR, color=C_PSM, ls="-.", lw=1.1); ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlim(1161, 1192); ax.set_ylim(lo - 0.12*(hi-lo), hi + 0.12*(hi-lo))

ax = axes[1, 1]; _ax(ax, r"$t_q$ (s)", "Applied current (A)"); _bands(ax, t, active)
for k in D: ax.step(t, D[k]["u"][:n], where="mid", color=COL[k])
for y in (56.69, 56.75, 56.90): ax.axhline(y, color="0.45", ls="--", lw=0.9)
ax.set_ylim(56.675, 56.928)

ax = axes[2, 0]; _ax(ax, r"$t_q$ (s)", r"$\dot m$ (kg s$^{-1}$)"); _bands(ax, t, active)
for k in D: ax.step(t, MD[k], where="mid", color=COL[k])
ax.set_yscale("log")

ax = axes[2, 1]; _ax(ax, r"$t_q$ (s)", "Objective contribution"); _bands(ax, t, active)
ax.plot(t, Jm["J_error"], color=COL["BASE"], marker="o", markersize=3, markevery=ME, label=r"$J_{error}$")
ax.plot(t, np.where(Jm["J_psm"] > 0, Jm["J_psm"], np.nan), color=C_PSM, marker="s", markersize=3, markevery=3, label=r"$J_{psm}$")
ax.plot(t, np.where(Jm["J_mdot"] > 0, Jm["J_mdot"], np.nan), color=C_MDOT, marker="^", markersize=3, markevery=ME, label=r"$J_{\dot m}$")
ax.set_yscale("log"); ax.legend(frameon=False, loc="lower right", ncol=3)

ax = axes[3, 0]; _ax(ax, r"$t_q$ (s)", "Optimizer runtime (s)"); _bands(ax, t, active)
for k in D: ax.plot(t, D[k]["rt"][:n], color=COL[k], marker=MRK[k], markersize=3.2, markevery=ME)
ax.axhline(120, color="0.45", ls=":", lw=1.0); ax.set_ylim(0, 900)

ax = axes[3, 1]
vals_tv = [D[k]["kpi"].get("total_variation", np.nan) for k in ["BASE", "PEN"]]
vals_iae = [D[k]["kpi"].get("IAE", np.nan) for k in ["BASE", "PEN"]]
x = np.arange(2); w = 0.36
ax.bar(x - w/2, vals_iae, w, color="0.55", alpha=0.85, edgecolor="0.25", linewidth=0.6, label=r"IAE (K$\cdot$s)")
ax.bar(x + w/2, vals_tv, w, color=C_MDOT, alpha=0.85, edgecolor="0.25", linewidth=0.6, label="Control TV (A)")
ax.set_xticks(x); ax.set_xticklabels([r"$\lambda$=0, 56.75", r"$\lambda$=$10^6$, 56.90"], fontsize=8)
ax.set_ylabel("KPI value"); ax.grid(True, axis="y", alpha=0.3)
ax.spines[["top", "right"]].set_visible(False)
for xi, v in zip(x - w/2, vals_iae):
    ax.annotate(f"{v:.2f}", xy=(xi, v), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=7.5)
for xi, v in zip(x + w/2, vals_tv):
    ax.annotate(f"{v:.2f}", xy=(xi, v), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=7.5)
ax.legend(frameon=False, fontsize=7.5)

for ax, lab in zip(axes.ravel(), ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]):
    ax.text(0.02, 0.98, lab, transform=ax.transAxes, va="top", ha="left",
            fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "penalty_full_comparison.pdf")
fig.savefig(FIG / "penalty_full_comparison.png")
plt.close(fig); print("OK penalty_full_comparison")

# ---- KPI table ----
print("\n" + "=" * 74)
print("CLOSED-LOOP KPIs")
print("=" * 74)
print(f"  {'metric':<32}{'lam=0/56.75':>14}{'lam=1e6/56.90':>16}{'delta':>13}")
print("  " + "-" * 72)
for key, lab in [("peak_PCT", "Peak PCT (K)"), ("min_safety_margin", "Min safety margin (K)"),
                 ("percent_overshoot", "Overshoot (%)"), ("settling_time", "Settling time (s)"),
                 ("IAE", "IAE (K.s)"), ("steady_state_error", "Steady-state error (K)"),
                 ("total_variation", "Control TV (A)"), ("max_control_rate", "Peak rate (A/s)")]:
    b = D["BASE"]["kpi"].get(key, np.nan); p = D["PEN"]["kpi"].get(key, np.nan)
    if np.isfinite(b) and np.isfinite(p):
        print(f"  {lab:<32}{b:>14.4f}{p:>16.4f}{p-b:>+13.4f}")
    else:
        bs = f"{b:.4f}" if np.isfinite(b) else "inf"
        ps = f"{p:.4f}" if np.isfinite(p) else "inf"
        print(f"  {lab:<32}{bs:>14}{ps:>16}{'--':>13}")
np.savez(OUT / "penalty_full_set.npz",
         t=t, pct_base=D["BASE"]["pct"][:n], pct_pen=D["PEN"]["pct"][:n],
         u_base=D["BASE"]["u"][:n], u_pen=D["PEN"]["u"][:n],
         mdot_base=MD["BASE"], mdot_pen=MD["PEN"],
         dsm_base=DSMn["BASE"], dsm_pen=DSMn["PEN"],
         rt_base=D["BASE"]["rt"][:n], rt_pen=D["PEN"]["rt"][:n],
         J_error=Jm["J_error"], J_psm=Jm["J_psm"], J_mdot=Jm["J_mdot"], J=Jm["J"],
         f_opt=Jm["f_opt"], dsm_thr=DSM_THR, target=TARGET, U_y=U_Y)

print(f"\n  total coolant   {MD['BASE'].sum():>12.2f}{MD['PEN'].sum():>16.2f}"
      f"{MD['PEN'].sum()-MD['BASE'].sum():>+13.2f}  kg/s")
print(f"  J_psm active in {int(np.nansum(active))} of {n} queries; "
      f"J_mdot active in {int(np.nansum(Jm['J_mdot']>0))}")
print(f"\n  figures in {FIG}")
