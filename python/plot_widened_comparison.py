"""Widened-pump penalised realization vs the archived-bounds baseline.

  baseline : phase2c_fullloop_pilot_LJ_seed1.mat        lambda = 0,   u_max = 56.75
  widened  : widepen_umax56p90_LJ_seed1.mat            lambda = 1e6, u_max = 56.90

CAVEAT stated up front: this comparison varies TWO things at once (the pump
saturation limit and the penalty). It is the shortcut Leo authorised instead of
also running a widened baseline. What licenses reading a trajectory change as
"the penalty acting" is the fixed-state verification (verify_widened_bounds.m):
at u_max = 56.90 the penalised optimum becomes INTERIOR (56.7425 A) rather than
pinned at the bound, which is what the archived bounds prevented.

  figures/wide_pct_vs_tq.{pdf,png}
  figures/wide_u_mdot_vs_tq.{pdf,png}
  figures/wide_jterms_vs_tq.{pdf,png}
  figures/widened_comparison.{pdf,png}   composed 3x2
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
OUT = HERE / "results_penalty_saturation"
RES = OUT
WORK = OUT
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

TARGET = 1179.2711776
U_Y = 1244.0
SPAN = U_Y - TARGET
DSM_THR = 0.00237737
XTICKS = list(range(1160, 1191, 10))
XLIM = (1158, 1192)

C_BASE = "#1f77b4"   # archived bounds, lambda = 0
C_WIDE = "#d62728"   # widened bounds, lambda = 1e6
C_PSM = "#ff7f0e"
C_MDOT = "#2ca02c"

# ECM pump model (CurrentToMassFlow.m) -- evaluated directly, not interpolated,
# because the widened run reaches currents outside the old anchor table.
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
    return {
        "t": np.asarray(d["time_snapshots"], float).ravel(),
        "pct": np.asarray(d["PCT"], float).ravel(),
        "u": np.asarray(d["u_V2P_trajectory"], float).ravel(),
        "rt": np.asarray(d["optimizer_runtimes"], float).ravel(),
    }


def parse_j(logpath):
    num = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    re_j = re.compile(rf"J_error:\s*({num}),\s*J_mdot:\s*({num}),\s*J_psm:\s*({num})")
    re_q = re.compile(r"^\[LJ\]\s+t=([\d.]+)")
    out = {k: [] for k in ["t", "J_error", "J_mdot", "J_psm", "J"]}
    buf = []
    with open(logpath, "r", errors="ignore") as fh:
        for line in fh:
            m = re_j.search(line)
            if m:
                buf.append([float(x) for x in m.groups()])
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
                buf = []
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
        ax.axvspan(t[i] - 0.25, t[i] + 0.25, color=C_PSM, alpha=0.16, lw=0,
                   label=label if first else None)
        first = False


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf"); fig.savefig(FIG / f"{name}.png")
    plt.close(fig); print("OK", name)


base = load_run(RES / "phase2c_fullloop_pilot_LJ_seed1.mat")
wide = load_run(WORK / "widepen_umax56p90_LJ_seed1.mat")
J = parse_j(WORK / "widepen_only.log")

n = min(len(base["t"]), len(wide["t"]))
t = base["t"][:n]
mb, mw = massflow(base["u"][:n]), massflow(wide["u"][:n])
dsm_b = (base["pct"][:n] - TARGET) / SPAN
dsm_w = (wide["pct"][:n] - TARGET) / SPAN

Jm = {k: np.full(n, np.nan) for k in ["J_error", "J_mdot", "J_psm", "J"]}
for i, tq in enumerate(t):
    hit = np.isclose(J["t"], tq, atol=1e-6)
    if hit.any():
        k = int(np.argmax(hit))
        for key in Jm:
            Jm[key][i] = J[key][k]
active = Jm["J_psm"] > 0

print(f"queries {n}; J_psm active in {int(np.nansum(active))} of {n}")

# ---- 1. PCT ----
fig, axs = plt.subplots(1, 2, figsize=(9.6, 3.9))
for ax, zoom in zip(axs, [False, True]):
    _ax(ax, r"Query time $t_q$ (s)", "PCT (K)")
    _bands(ax, t, active, r"$J_{psm}>0$" if not zoom else None)
    ax.plot(t, base["pct"][:n], color=C_BASE, marker="o", markersize=4, markevery=5,
            label=r"archived bounds, $\lambda=0$")
    ax.plot(t, wide["pct"][:n], color=C_WIDE, marker="s", markersize=4, markevery=5,
            label=r"widened pump, $\lambda=10^{6}$")
    ax.axhline(TARGET, color="0.45", ls="--", lw=1.1, label=r"setpoint $T_p^*$")
    if zoom:
        ax.set_xlim(1161, 1192)
        lo = min(base["pct"][2:n].min(), wide["pct"][2:n].min())
        hi = max(base["pct"][2:n].max(), wide["pct"][2:n].max())
        ax.set_ylim(lo - 0.08*(hi-lo), hi + 0.08*(hi-lo))
        ax.set_title("zoom, transient excluded", fontsize=9.5, color="0.35")
    else:
        ax.set_ylim(1140, 1185)
        ax.set_title("full range", fontsize=9.5, color="0.35")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
fig.tight_layout()
save(fig, "wide_pct_vs_tq")

# ---- 2. control and coolant ----
fig, axs = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True)
ax = axs[0]
ax.set_ylabel("Applied ACS current (A)")
ax.set_xlim(*XLIM); ax.set_xticks(XTICKS)
ax.grid(True, alpha=0.3); ax.spines[["top", "right"]].set_visible(False)
_bands(ax, t, active, r"$J_{psm}>0$")
ax.step(t, base["u"][:n], where="mid", color=C_BASE, label=r"archived, $\lambda=0$")
ax.step(t, wide["u"][:n], where="mid", color=C_WIDE, label=r"widened, $\lambda=10^{6}$")
for y, lab in [(56.69, r"$u_{\min}$"), (56.75, r"archived $u_{\max}$"),
               (56.90, r"widened $u_{\max}$")]:
    ax.axhline(y, color="0.45", ls="--", lw=1.0)
    ax.text(1158.4, y + 0.004, lab, fontsize=7.5, color="0.35")
ax.set_ylim(56.675, 56.925)
ax.legend(frameon=False, loc="center right", fontsize=8)
ax = axs[1]
_ax(ax, r"Query time $t_q$ (s)", r"ACS coolant mass flow (kg s$^{-1}$)")
_bands(ax, t, active)
ax.step(t, mb, where="mid", color=C_BASE, label=r"archived, $\lambda=0$")
ax.step(t, mw, where="mid", color=C_WIDE, label=r"widened, $\lambda=10^{6}$")
ax.set_yscale("log")
ax.text(0.99, 0.04, rf"$\Sigma\dot m$: {mb.sum():.1f} $\to$ {mw.sum():.1f} kg s$^{{-1}}$"
                    rf"  ({mw.sum()/mb.sum():.2f}$\times$)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color="0.25")
ax.legend(frameon=False, loc="upper right", fontsize=8)
fig.tight_layout()
save(fig, "wide_u_mdot_vs_tq")

# ---- 3. J terms ----
fig, ax = plt.subplots(figsize=(7.2, 4.2))
_ax(ax, r"Query time $t_q$ (s)", "Objective contribution")
_bands(ax, t, active, r"$J_{psm}>0$")
ax.plot(t, Jm["J_error"], color=C_BASE, marker="o", markersize=3.6, markevery=5,
        label=r"$J_{error}$ (tracking)")
ax.plot(t, np.where(Jm["J_psm"] > 0, Jm["J_psm"], np.nan), color=C_PSM,
        marker="s", markersize=3.6, markevery=5, label=r"$J_{psm}$ (DSM penalty)")
ax.plot(t, np.where(Jm["J_mdot"] > 0, Jm["J_mdot"], np.nan), color=C_MDOT,
        marker="^", markersize=3.6, markevery=5, label=r"$J_{\dot m}$ (flow rate)")
ax.plot(t, Jm["J"], color="0.25", ls=":", lw=1.3, label=r"$J$ (total)")
ax.set_yscale("log")
ax.legend(frameon=False, loc="lower right", ncol=2, fontsize=8)
save(fig, "wide_jterms_vs_tq")

# ---- composed ----
plt.rcParams.update({"font.size": 10.5, "legend.fontsize": 8})
fig, axes = plt.subplots(3, 2, figsize=(7.8, 10.2))

ax = axes[0, 0]; _ax(ax, r"$t_q$ (s)", "PCT (K)")
_bands(ax, t, active)
ax.plot(t, base["pct"][:n], color=C_BASE, marker="o", markersize=3.4, markevery=5, label=r"$\lambda=0$, 56.75")
ax.plot(t, wide["pct"][:n], color=C_WIDE, marker="s", markersize=3.4, markevery=5, label=r"$\lambda=10^6$, 56.90")
ax.axhline(TARGET, color="0.45", ls="--", lw=1.0)
ax.set_ylim(1140, 1185); ax.legend(frameon=False, loc="lower right")

ax = axes[0, 1]; _ax(ax, r"$t_q$ (s)", "PCT (K), zoom")
_bands(ax, t, active)
ax.plot(t, base["pct"][:n], color=C_BASE, marker="o", markersize=3.4, markevery=4)
ax.plot(t, wide["pct"][:n], color=C_WIDE, marker="s", markersize=3.4, markevery=4)
ax.axhline(TARGET, color="0.45", ls="--", lw=1.0)
ax.set_xlim(1161, 1192)
_lo = min(base["pct"][2:n].min(), wide["pct"][2:n].min())
_hi = max(base["pct"][2:n].max(), wide["pct"][2:n].max())
ax.set_ylim(_lo - 0.08*(_hi-_lo), _hi + 0.08*(_hi-_lo))

ax = axes[1, 0]; _ax(ax, r"$t_q$ (s)", "Applied current (A)")
_bands(ax, t, active)
ax.step(t, base["u"][:n], where="mid", color=C_BASE)
ax.step(t, wide["u"][:n], where="mid", color=C_WIDE)
for y in (56.69, 56.75, 56.90):
    ax.axhline(y, color="0.45", ls="--", lw=0.9)
ax.set_ylim(56.675, 56.925)

ax = axes[1, 1]; _ax(ax, r"$t_q$ (s)", r"$\dot m$ (kg s$^{-1}$)")
_bands(ax, t, active)
ax.step(t, mb, where="mid", color=C_BASE)
ax.step(t, mw, where="mid", color=C_WIDE)
ax.set_yscale("log")

ax = axes[2, 0]; _ax(ax, r"$t_q$ (s)", "Objective contribution")
_bands(ax, t, active)
ax.plot(t, Jm["J_error"], color=C_BASE, marker="o", markersize=3.2, markevery=5, label=r"$J_{error}$")
ax.plot(t, np.where(Jm["J_psm"] > 0, Jm["J_psm"], np.nan), color=C_PSM, marker="s",
        markersize=3.2, markevery=5, label=r"$J_{psm}$")
ax.plot(t, np.where(Jm["J_mdot"] > 0, Jm["J_mdot"], np.nan), color=C_MDOT, marker="^",
        markersize=3.2, markevery=5, label=r"$J_{\dot m}$")
ax.set_yscale("log"); ax.legend(frameon=False, loc="lower right", ncol=3)

ax = axes[2, 1]; _ax(ax, r"$t_q$ (s)", "Normalised DSM")
_bands(ax, t, active)
ax.plot(t, dsm_b, color=C_BASE, marker="o", markersize=3.2, markevery=4)
ax.plot(t, dsm_w, color=C_WIDE, marker="s", markersize=3.2, markevery=4)
ax.axhline(DSM_THR, color=C_PSM, ls="-.", lw=1.1)
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlim(1161, 1192)
_lo = min(dsm_b[2:n].min(), dsm_w[2:n].min()); _hi = max(dsm_b[2:n].max(), dsm_w[2:n].max())
ax.set_ylim(_lo - 0.12*(_hi-_lo), _hi + 0.12*(_hi-_lo))

for ax, lab in zip(axes.ravel(), ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]):
    ax.text(0.02, 0.98, lab, transform=ax.transAxes, va="top", ha="left",
            fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig(FIG / "widened_comparison.pdf"); fig.savefig(FIG / "widened_comparison.png")
plt.close(fig); print("OK widened_comparison")

# ---- summary ----
print("\n" + "=" * 74)
print("ARCHIVED BOUNDS (lambda=0)  vs  WIDENED PUMP (lambda=1e6)")
print("=" * 74)
print(f"  {'metric':<34}{'archived':>13}{'widened':>13}{'delta':>13}")
print("  " + "-" * 71)
for name, b, w in [
    ("peak PCT [K]", base["pct"][:n].max(), wide["pct"][:n].max()),
    ("min PCT [K]", base["pct"][:n].min(), wide["pct"][:n].min()),
    ("final PCT [K]", base["pct"][n-1], wide["pct"][n-1]),
    ("max DSM", dsm_b.max(), dsm_w.max()),
    ("min safety margin [K]", (U_Y-base["pct"][:n]).min(), (U_Y-wide["pct"][:n]).min()),
    ("mean |PCT-target| [K]", np.abs(base["pct"][:n]-TARGET).mean(),
     np.abs(wide["pct"][:n]-TARGET).mean()),
    ("mean u [A]", base["u"][:n].mean(), wide["u"][:n].mean()),
    ("max u [A]", base["u"][:n].max(), wide["u"][:n].max()),
    ("mean mdot [kg/s]", mb.mean(), mw.mean()),
    ("max mdot [kg/s]", mb.max(), mw.max()),
    ("total coolant [kg/s]", mb.sum(), mw.sum()),
    ("control TV [A]", np.abs(np.diff(base["u"][:n])).sum(),
     np.abs(np.diff(wide["u"][:n])).sum()),
]:
    print(f"  {name:<34}{b:>13.4f}{w:>13.4f}{w-b:>+13.4f}")

print(f"\n  J_psm active in {int(np.nansum(active))} of {n} queries")
idx = np.where(active)[0]
if idx.size:
    print(f"  active window: t = {t[idx[0]]:.1f} .. {t[idx[-1]]:.1f} s")
    print(f"  J_psm: max {np.nanmax(Jm['J_psm']):.4e}, "
          f"median over active {np.nanmedian(Jm['J_psm'][active]):.4e}")
mdot_act = int(np.nansum(Jm["J_mdot"] > 0))
print(f"  J_mdot active in {mdot_act} of {n} queries "
      f"(it was identically 0 with the archived bounds)")
np.savez(OUT / "figure_widened_comparison.npz",
         t=t, pct_base=base["pct"][:n], pct_wide=wide["pct"][:n],
         u_base=base["u"][:n], u_wide=wide["u"][:n], mdot_base=mb, mdot_wide=mw,
         dsm_base=dsm_b, dsm_wide=dsm_w,
         J_error=Jm["J_error"], J_psm=Jm["J_psm"], J_mdot=Jm["J_mdot"], J=Jm["J"],
         dsm_thr=DSM_THR, target=TARGET, U_y=U_Y)
print(f"\n  saved {OUT / 'figure_widened_comparison.npz'}")
print(f"  figures in {FIG}")
