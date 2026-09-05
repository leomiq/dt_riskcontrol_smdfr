"""Why does the widened penalised PCT drift below the setpoint instead of converging?

Three candidate causes, and what would distinguish them in the data we already have:

  H1  lambda too large / threshold too tight
      The penalty keeps charging even when the DSM is already negative, so the
      optimizer keeps buying margin it does not need.
      SIGNATURE: J_psm > 0 while DSM < 0, i.e. active in the drift region.

  H2  horizon too short
      The controller is still recovering from the early over-cooling and would
      converge given more time.
      SIGNATURE: the PCT trend flattens toward t=1190 (decreasing |slope|), and
      the applied control moves toward the value that holds the setpoint.

  H3  LJ solves the late queries badly (sub-optimal solutions)
      The drift is an artefact of the optimizer, not of the cost function.
      SIGNATURES: (a) the objective f* stops improving / LJ hits its early-stop
      rule more often late; (b) the applied u is far from the u that minimises
      J_error alone; (c) the per-query cost decomposition shows J dominated by a
      term the optimizer could have reduced.

The three are not mutually exclusive; the point is to see which dominates before
spending 4 h on a re-run.
"""
from pathlib import Path
import re

import numpy as np

HERE = Path(r"external revision/Response to reviewers/evidence/phase2c_fixed_state_experiment")
OUT = HERE / "results_penalty_saturation"

TARGET = 1179.2711776
U_Y = 1244.0
SPAN = U_Y - TARGET
DSM_THR = 0.00237737
LAMBDA = 1e6

d = np.load(OUT / "penalty_full_set.npz")
t = d["t"]
pct_b, pct_p = d["pct_base"], d["pct_pen"]
u_b, u_p = d["u_base"], d["u_pen"]
md_b, md_p = d["mdot_base"], d["mdot_pen"]
dsm_b, dsm_p = d["dsm_base"], d["dsm_pen"]
Je, Jp_, Jm_, Jt = d["J_error"], d["J_psm"], d["J_mdot"], d["J"]
rt_b, rt_p = d["rt_base"], d["rt_pen"]
n = len(t)

print("=" * 78)
print("PCT TRAJECTORY, penalised widened run")
print("=" * 78)
print(f"{'t_q':>8}{'PCT':>11}{'DSM':>10}{'u':>10}{'mdot':>9}"
      f"{'J_error':>11}{'J_psm':>11}{'J_mdot':>9}{'runtime':>9}")
print("-" * 78)
for i in range(n):
    flag = ""
    if Jp_[i] > 0:
        flag = " psm"
    if dsm_p[i] < 0 and Jp_[i] > 0:
        flag = " PSM+DSM<0 !!"
    print(f"{t[i]:8.1f}{pct_p[i]:11.4f}{dsm_p[i]:+10.5f}{u_p[i]:10.5f}{md_p[i]:9.3f}"
          f"{Je[i]:11.3e}{Jp_[i]:11.3e}{Jm_[i]:9.3f}{rt_p[i]:9.1f}{flag}")

# ---------------------------------------------------------------- H1
print("\n" + "=" * 78)
print("H1  lambda too large / threshold too tight")
print("=" * 78)
act = Jp_ > 0
neg = dsm_p < 0
both = act & neg
print(f"queries with J_psm > 0                : {int(act.sum())} of {n}")
print(f"queries with DSM < 0 (below setpoint) : {int(neg.sum())}")
print(f"queries with BOTH (penalty on a safe state): {int(both.sum())}")
if both.any():
    print(f"  -> t = {t[both][0]:.1f} .. {t[both][-1]:.1f}")
    print("  H1 SUPPORTED: the penalty charges while the margin is already safe.")
else:
    print("  H1 NOT SUPPORTED: once DSM goes negative the penalty is exactly zero,")
    print("  so lambda cannot be what drives the late drift.")

last = act.nonzero()[0][-1] if act.any() else -1
if last >= 0:
    print(f"\nlast query with an active penalty: t = {t[last]:.1f} "
          f"(index {last} of {n-1})")
    print(f"drift region t > {t[last]:.1f}: PCT {pct_p[last]:.4f} -> {pct_p[-1]:.4f} K "
          f"({pct_p[-1]-pct_p[last]:+.4f} K over {n-1-last} queries)")
    print("  the drift happens with J_psm == 0 throughout, i.e. AFTER the penalty")
    print("  has switched off -- it is not the penalty pushing the PCT down.")

# ---------------------------------------------------------------- H2
print("\n" + "=" * 78)
print("H2  horizon too short (still recovering)")
print("=" * 78)


def slope(x, y):
    if len(x) < 3:
        return np.nan
    return np.polyfit(x, y, 1)[0] * 1000.0     # mK per s


win = 8
print(f"{'window':>18}{'PCT slope (mK/s)':>20}{'mean u':>10}{'mean mdot':>11}")
print("-" * 60)
for a in range(0, n - win + 1, win):
    b = min(a + win, n)
    print(f"  t {t[a]:6.1f}-{t[b-1]:6.1f}{slope(t[a:b], pct_p[a:b]):20.3f}"
          f"{u_p[a:b].mean():10.5f}{md_p[a:b].mean():11.3f}")
s_first = slope(t[last+1:last+1+10], pct_p[last+1:last+1+10]) if last + 11 <= n else np.nan
s_last = slope(t[-10:], pct_p[-10:])
print(f"\nslope just after the penalty switches off : {s_first:+.3f} mK/s")
print(f"slope over the final 10 queries           : {s_last:+.3f} mK/s")
if np.isfinite(s_first) and np.isfinite(s_last):
    if abs(s_last) < 0.6 * abs(s_first):
        print("  H2 SUPPORTED: the drift is decelerating -- a longer horizon would")
        print("  plausibly converge.")
    elif abs(s_last) > 1.2 * abs(s_first):
        print("  H2 CONTRADICTED: the drift is ACCELERATING, so more time would make")
        print("  it worse, not better.")
    else:
        print("  H2 INCONCLUSIVE: the drift rate is roughly constant -- extending the")
        print("  horizon would extend the drift, not resolve it.")

# where does the baseline sit for comparison?
print(f"\nbaseline over the same final window: slope {slope(t[-10:], pct_b[-10:]):+.3f} mK/s, "
      f"final PCT {pct_b[-1]:.4f} K (error {pct_b[-1]-TARGET:+.4f} K)")
print(f"penalised final PCT {pct_p[-1]:.4f} K (error {pct_p[-1]-TARGET:+.4f} K)")

# ---------------------------------------------------------------- H3
print("\n" + "=" * 78)
print("H3  LJ sub-optimal on the late queries")
print("=" * 78)

# (a) is J_error alone larger late than a naive hold would give?
#     With the penalty off, J = J_error + J_mdot, so a well-solved query should
#     drive J_error down. Track it against the DSM it produced.
print(f"{'window':>18}{'mean J_error':>14}{'mean J_mdot':>13}{'J_mdot share':>14}")
print("-" * 60)
for a in range(0, n - win + 1, win):
    b = min(a + win, n)
    je, jm = np.nanmean(Je[a:b]), np.nanmean(Jm_[a:b])
    print(f"  t {t[a]:6.1f}-{t[b-1]:6.1f}{je:14.3f}{jm:13.3f}{jm/max(je+jm,1e-12):14.3f}")

print(f"\nJ_error: first 8 queries {np.nanmean(Je[:8]):.3f}, "
      f"last 8 {np.nanmean(Je[-8:]):.3f}  "
      f"({np.nanmean(Je[-8:])/max(np.nanmean(Je[:8]),1e-12):.2f}x)")

# (b) control effort late: is the optimizer using the range or thrashing?
print("\ncontrol behaviour in the drift region:")
dr = slice(last + 1, n)
print(f"  u range   : [{u_p[dr].min():.5f}, {u_p[dr].max():.5f}]  "
      f"(full band available: [56.69, 56.90])")
print(f"  at u_min   : {int((u_p[dr] < 56.6901).sum())} of {n-last-1}")
print(f"  at u_max   : {int((u_p[dr] > 56.8999).sum())} of {n-last-1}")
print(f"  |du| mean  : {np.abs(np.diff(u_p[dr])).mean():.5f} A, "
      f"max {np.abs(np.diff(u_p[dr])).max():.5f} A")
print(f"  mdot range : [{md_p[dr].min():.3f}, {md_p[dr].max():.3f}] kg/s")

# (c) does the optimizer keep improving? runtime is a proxy for early-stopping
print(f"\noptimizer runtime: first 8 {rt_p[:8].mean():.1f} s, last 8 {rt_p[-8:].mean():.1f} s")
print(f"  baseline for comparison: first 8 {rt_b[:8].mean():.1f} s, last 8 {rt_b[-8:].mean():.1f} s")

# ---------------------------------------------------------------- the key test
print("\n" + "=" * 78)
print("KEY TEST: is the applied control consistent with minimising J_error?")
print("=" * 78)
print("With the penalty off (DSM < 0), the cost is J_error + J_mdot. J_error is")
print("minimised by holding the PCT at the setpoint. The PCT is BELOW setpoint, so")
print("the correct action is LESS coolant, i.e. u -> u_min. Check what LJ chose:")
print()
below = (dsm_p < 0) & (Jp_ == 0)
if below.any():
    idx = below.nonzero()[0]
    at_min = (u_p[idx] < 56.6901).sum()
    above_mid = (u_p[idx] > 56.75).sum()
    print(f"  queries below setpoint with no penalty: {len(idx)}")
    print(f"    LJ chose u_min (correct)            : {int(at_min)} "
          f"({100*at_min/len(idx):.0f}%)")
    print(f"    LJ chose u > archived u_max 56.75   : {int(above_mid)} "
          f"({100*above_mid/len(idx):.0f}%)  <- MORE cooling while already too cold")
    if above_mid > 0.2 * len(idx):
        print("\n  >> H3 SUPPORTED: in a sizeable share of the drift region LJ commands")
        print("     MORE coolant than the archived maximum while the PCT is already")
        print("     below the setpoint. That is inconsistent with minimising J_error,")
        print("     so those are sub-optimal solutions, not a cost-function effect.")
    else:
        print("\n  >> H3 NOT the main driver: LJ mostly commands minimum coolant, which")
        print("     is what J_error asks for. The drift then reflects plant dynamics")
        print("     (thermal inertia), not optimizer failure.")

print("\n" + "=" * 78)
print("SUMMARY OF EVIDENCE")
print("=" * 78)
print(f"  penalty active window      : t = {t[act][0]:.1f} .. {t[act][-1]:.1f} "
      f"({int(act.sum())} queries)")
print(f"  penalty active while safe  : {int(both.sum())} queries")
print(f"  drift after penalty off    : {pct_p[-1]-pct_p[last]:+.4f} K over "
      f"{n-1-last} queries ({slope(t[last+1:], pct_p[last+1:]):+.3f} mK/s)")
print(f"  final tracking error       : {pct_p[-1]-TARGET:+.4f} K "
      f"(baseline {pct_b[-1]-TARGET:+.4f} K)")
print(f"  J_mdot share, last 8       : "
      f"{np.nanmean(Jm_[-8:])/max(np.nanmean(Je[-8:]+Jm_[-8:]),1e-12):.3f}")
