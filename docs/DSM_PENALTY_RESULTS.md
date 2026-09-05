# Deterministic safety margin (DSM) penalty — experimental results

Purpose: demonstrate that the MPC optimisation problem, as formulated, actually reshapes the
PCT trajectory during the accident — that the safety-margin penalty activates in critical
situations, drives more ACS coolant injection, and smooths the control action.

`DSM = (Tp − target) / (U_y − target)`, with `U_y = 1244 K` (absolute safety limit) and
`target = 1179.2712 K`, so the DSM span is 64.7288 K. `DSM = 1` means the PCT sits at the limit.

Penalty form used for the sweeps (`fitness_cost_dsm.m`):

```
J_psm = lambda_pct * sum( max(0, DSM − dsm_thr)^2 )
```

written on the DSM itself rather than on `Tp − U_y` in Kelvin, so `lambda` is dimensionless and
means the same thing across thresholds and operating points. The production `fitness_cost.m`
keeps its Kelvin form and was not modified.

---

## Preconditions verified before any sweep

### 1. The penalty applies to FUTURE instants only

This was a known past bug. `audit_psm_future_only.m` confirms it is fixed:

```
decision instant t_root         = 1180.0 s
warm-up (PAST) window           = [1170.0 .. 1179.5]  (20 samples)
reconstructed range             = [1180.00 .. 1192.50]
samples strictly before t_root  = 0
samples inside the warm-up window = 0
>> CLEAN
```

The `keep = t_signal > time_delayed(end)` filter guarantees this for the root and every child,
so the optimizer is never charged for a trajectory it cannot influence. **Sweeps are valid.**

### 2. Lambda alone is powerless

With the Kelvin form and `peak Tp = 1179.6 K` against `U_y = 1244 K`, the violation term is
identically zero, so `lambda × 0 = 0`. Raising `lambda_pct` from 10 to 1e6 changes nothing.
**The threshold must come down first** — which is why the sweep is joint in `(dsm_thr, lambda)`.

### 3. The actuator has real authority

`test_control_authority.m` at t = 1180:

| `u` (A) | coolant mass flow | peak Tp (K) | peak DSM |
|---|---|---|---|
| 56.690 | 0.5798 | 1180.5539 | 0.019816 |
| 56.720 | 6.7231 | 1179.6054 | 0.005163 |
| 56.750 | 12.8648 | 1179.5789 | 0.004755 |

The 0.06 A range maps to **12.28 units of coolant mass flow** (a 22× factor — the current→flow
relation is strongly nonlinear) and moves the PCT by **0.975 K**. Less coolant → higher PCT →
more margin consumed: the expected physics.

---

## Result A — t = 1180 s (convergence phase): the penalty DOES reorient the decision

Control response over the grid (11 DET evaluations; the sweep itself is analytic because the
trajectory does not depend on the cost parameters):

| `u` (A) | mass flow | peak Tp (K) | peak DSM | `J_error` |
|---|---|---|---|---|
| 56.690 | 0.5798 | 1180.5539 | 0.019816 | 19.53 |
| 56.702 | 3.0373 | 1180.0582 | 0.012158 | 8.94 |
| 56.714 | 5.4946 | 1179.6126 | 0.005274 | 2.73 |
| **56.726** | **7.9516** | **1179.6001** | **0.005082** | **0.8195** ← unpenalised optimum |
| 56.738 | 10.4083 | 1179.5895 | 0.004918 | 3.16 |
| 56.750 | 12.8648 | 1179.5789 | 0.004755 | 9.69 |

Reachable peak-DSM band `[0.004755, 0.019816]`; PCT authority 0.9749 K.

**Reorientation achieved** with `dsm_thr = 0.002377`, `lambda = 1e6`:

| | unpenalised | penalised |
|---|---|---|
| `u*` | 56.72600 A | **56.74400 A** (+0.018 A) |
| coolant mass flow | 7.9516 | **11.6366** (+46%) |
| peak DSM | 0.005082 | **0.004837** |
| `J_psm` share of `J` | 0 | **65.7%** |

**Graded activation** is visible along that threshold row as lambda rises 1e1 → 1e6:

| lambda | 1e1 | 1e2 | 1e3 | 1e4 | 1e5 | 1e6 |
|---|---|---|---|---|---|---|
| `J_psm` share | 0.03% | 0.34% | 3.3% | 25.7% | 56.7% | 65.7% |

The margin engages progressively rather than switching on — which is the intended "activate in
critical situations and smooth the control action" behaviour.

### Caveat on Result A

Only **2 of 49** combinations move the optimum, and both need `dsm_thr = 0.002377`, which sits
*below* the reachable band and is therefore violated at every `u`. The penalty then acts as an
always-on gradient rather than a selectively activated constraint, and at `lambda = 1e6` it owns
66% of the cost (the controller is minimising DSM, not tracking the setpoint).

Every threshold *inside* the reachable band (0.00476–0.0198) moves the optimum at no lambda,
because `J_error` has a sharp interior minimum at 56.726 (0.82 versus 9.69 at `u_max`) while the
DSM barely varies (0.00476–0.00508) in that neighbourhood — the penalty has no leverage there.

---

## Result B — t = 1161 s (acute phase): the penalty is redundant, and that is informative

| `u` (A) | mass flow | peak Tp (K) | peak DSM | `J_error` |
|---|---|---|---|---|
| 56.690 | 0.5798 | 1182.2828 | 0.046526 | 170.29 |
| 56.714 | 5.4946 | 1181.3371 | 0.031917 | 107.79 |
| **56.750** | **12.8648** | **1181.3047** | **0.031416** | **48.28** ← unpenalised optimum |

**0 of 49** combinations move the optimum. The reason is structural, not a failure: at t = 1161
the PCT sits **2.03 K above the setpoint**, so `J_error` is monotonically decreasing in `u` —
more cooling is better for tracking *and* for safety simultaneously. The unpenalised optimum is
already saturated at `u_max`, so there is nowhere for the penalty to move it.

Worth noting: at `u_min` the DSM reaches 0.0465, i.e. **without coolant injection the margin
would be 47% consumed**. The safety mechanism is real; the tracking term simply already covers it.

---

## The two regimes compared

| | t = 1180 (convergence) | t = 1161 (acute) |
|---|---|---|
| PCT vs setpoint | +0.33 K | **+2.03 K** |
| `J_error` vs `u` | non-monotonic, interior minimum | **monotonically decreasing** |
| Unpenalised optimum | 56.726 A (interior) | **56.750 A (saturated)** |
| Reachable DSM band | [0.0048, 0.0198] | [0.0314, 0.0465] |
| Combinations moving the optimum | 2 of 49 | **0 of 49** |

**Central finding:** the DSM penalty can only reorient the decision when the tracking optimum is
**interior** — i.e. when the PCT is near the setpoint and over-cooling carries a cost. When the
PCT is well above the setpoint, tracking and safety agree and the penalty is redundant.

This is a defensible statement about the formulation: **the safety constraint never contradicts
the tracking objective.** In the acute phase both demand maximum coolant; in the convergence
phase, where tracking alone would tolerate a slightly high PCT, the DSM penalty is what prevents
it.

Honest limitation: in this scenario the DSM penalty is *not* what saves the margin — `J_error`
is, because the setpoint (1179.27 K) already sits 64 K below the limit (1244 K). The DSM term
adds value only where tracking and safety **compete**.

---

## Result C — severe scenario regenerated from the HRM

Varying ROM parameters cannot make the scenario severe: `vel_ACS` is an **HRM** parameter and the
ROM receives the ACS action through its `u_ACS` inport. Verified — writing `vel_ACS` straight into
the ROM's own model workspace and sweeping `[0.1, 0.05, 0.02, 0.01, 0.002]` (a 50× degradation)
changed the predicted trajectory by **nothing at all**, and no ROM block references `vel_ACS`.
Separately, `DET_2.m` lines 164/449 set `vel_ACS` with `'Workspace', name_Wks` where `name_Wks`
is a **`.m` init script**, not a Simulink model — those calls are inert.

So severity was generated at the source: **HRM → BiLSTM MC-dropout → ROM → DET**
(`generate_severe_scenario.m`). Environment: Python 3.11, torch with CUDA on an RTX 3050.
Per scenario: HRM 12 s, inference 0.2 s, ROM 5 s.

### Note on a retracted claim

An earlier draft called these cases invalid because `delta_tw_max ≈ 72 K` "versus ~2 K nominal",
implying BiLSTM extrapolation. **That was wrong.** The archived nominal workspace also carries
~72.7 K (`det_params.Y_denorm` col 2: 72.6221 … 72.6999), and the regenerated input signals match
the archived ones (`mdot_c_leak` last value −4.5593 vs −4.5594, both 2361 samples). The network is
in its normal range and the severe cases are valid.

### `B` controls severity monotonically; `A` barely matters

With `A = 0.2834`, `vel_ACS = 0.09`, query t = 1180:

| `B` | future DSM band | authority |
|---|---|---|
| 10 (nominal) | [0.056, 0.072] | 1.030 K |
| 12 | [0.135, 0.153] | 1.215 K |
| **15** | **[0.210, 0.232]** | **1.381 K** |

Raising `A` from 0.2834 to 0.36 moves `delta_tw_max` by only 1.7 K and *reduces* the future DSM
(0.232 → 0.026) — a counterintuitive inverse relation.

**Metric note:** the global PCT peak always lands at `t = t_q − delayed` (= 1170), the ROM
initialisation transient in the warm-up window, identical for every `u` by construction and hence
uncontrollable. The controllable quantity is the peak over the **future horizon**, which
`fitness_cost_dsm.m` isolates.

### Sweep at `B = 15` (the most severe valid case): 0 of 49

| `u` (A) | mass flow | `J_error` | future peak DSM |
|---|---|---|---|
| 56.690 | 0.5798 | 5358.6 | 0.231789 |
| 56.720 | 6.7231 | 4824.9 | 0.211196 |
| **56.750** | **12.8648** | **4333.8** ← min | **0.210442** ← min |

`J_error` and the DSM are **both monotonically decreasing in `u` and minimise at the same bound**.
The unpenalised optimum is already at maximum coolant, so no `(dsm_thr, lambda)` pair can displace
it; at `u_max` the DSM (0.2104) sits below every threshold inside the band, so `J_psm` share is 0
throughout.

---

## Central structural finding

| Case | PCT vs setpoint | unpenalised optimum | combos moving `u*` |
|---|---|---|---|
| t=1180 nominal | +0.33 K | 56.726 A **interior** | **2 of 49** |
| t=1161 acute | +2.03 K | 56.750 A saturated | 0 of 49 |
| severe B=10 | — | 56.750 A saturated | 0 of 49 |
| severe B=15 | — | 56.750 A saturated | 0 of 49 |

**The DSM penalty can only reorient the decision when the tracking optimum is interior.** Whenever
the PCT sits above the setpoint, tracking and safety are the same objective — more coolant improves
both — and the optimum saturates at `u_max`.

This is a property of the formulation, not a failure: **the safety constraint never contradicts the
control objective.** In the acute phase both demand maximum cooling; in the convergence phase, where
tracking alone would tolerate a slightly high PCT, the DSM penalty is what prevents it, and there it
demonstrably acts (Result A).

For the DSM term to bind harder, either over-cooling must carry a cost `J_error` does not capture,
or the setpoint must sit closer to the limit — with setpoint 1179.27 K against a 1244 K limit there
are 64.7 K of headroom that tracking covers unaided. That is a formulation change, not an experiment.

## Superseded: Result C placeholder

To make tracking and safety compete, the archive provides `LJ_test_sampled_vel_ACS_09.mat`:

| | nominal series | degraded series |
|---|---|---|
| `sampled_vel_ACS` | 2.0 | **0.09559** (≈1/21 of capacity) |
| `time_snapshot_root` | 1161 / 1180 | **1300** (late accident) |
| `tsim` | 1173.5 / 1192.5 | 1325 |
| `sampled_A` | 0.2834 | 0.2834 (same leak) |
| archived `outRoot` peak tw_max | 1181.18 K | **1181.64 K** |

Same leak, same drift, but the auxiliary cooling system is at a twentieth of its capacity — the
condition under which the PCT can climb *despite* maximum coolant demand. Sweep running.

---

---

## Result D — closed-loop realization, archived bounds (61 queries, 4.01 h)

The decisive experiment: LJ, seed 1, same accident stream, `dsm_thr = 0.00237737`,
`lambda = 1e6`, **the cost function the only difference** from the archived pilot.

Files: `results_penalised_run/penalised_thr0p00238_lam1e06_LJ_seed1.mat` versus
`phase2c_fullloop_pilot_LJ_seed1.mat`. Figures: `figures/pen_*.{pdf,png}` and
`figures/penalty_comparison.{pdf,png}`.

### The penalty was unambiguously active

`J_psm` dominated `J_error` over **36 of 61 queries (59%)**, in a single contiguous window:

```
t = 1160.0 .. 1177.5 s      ratio J_psm/J_error: median 56.7, max 197.3
```

The crossing at t = 1177.5 is abrupt (vertical on a log axis) and coincides exactly with the
normalised DSM falling below the threshold: after that instant the PCT sits *below* the setpoint,
so `J_psm` switches off on its own. **The constraint engages while the margin is compromised and
withdraws when it is not** — the intended design behaviour, visible without ambiguity.

### But the trajectory barely moved

| metric | baseline | penalised | delta |
|---|---|---|---|
| peak PCT [K] | 1181.2132 | 1181.2132 | **+0.0000** |
| max DSM | 0.0300 | 0.0300 | **+0.0000** |
| min safety margin [K] | 62.7868 | 62.7868 | **+0.0000** |
| mean \|PCT − target\| [K] | 1.1698 | 1.1750 | +0.0052 |
| total coolant (Σ) | 559.4855 | 557.3571 | −2.1284 |
| control total variation [A] | 0.8913 | 0.8625 | **−0.0288 (−3.2%)** |
| queries at `u_max` | 23 | 21 | −2 |

Two structural reasons, not a weak penalty:

1. **The PCT peak occurs at t = 1160.0, the first query** — before any control action has had
   an effect. No controller can move it.
2. **The control was SATURATED at `u_max` through most of the dominance window.** The safety term
   had already obtained everything it could ask for; saturation, not the penalty, was the binding
   constraint.

Per-query differences alternate sign with amplitudes up to ±12 kg s⁻¹ — an order of magnitude
above the net effect — so with one realization per configuration the penalty's contribution
cannot be separated from optimizer noise (LJ is stochastic and `parfor` reorders evaluations).
Distinguishing them would need a baseline rerun with a different seed.

---

## Result E — the saturation limit was the binding constraint

The MHD pump model (`CurrentToMassFlow.m`, ECM of Nowak et al. 2022, DOI
10.1016/j.anucene.2022.109142) explains why the archived bounds leave so little room:

```
zero-flow current  DeltaP/C          = 56.6872 A
archived bounds    [56.69, 56.75]    = only 0.0028 .. 0.0628 A above onset
slope              dmdot/dI          ~ 204.7 kg/s per A  (essentially linear)
```

The window sits **right at the pump's onset**, capping the coolant at 12.86 kg s⁻¹ — not because
the pump is weak, but because of where the bounds are placed:

| `u_max` | max mass flow | factor |
|---|---|---|
| 56.75 (archived) | 12.86 kg s⁻¹ | 1× |
| 56.90 | 43.55 | **3.4×** |
| 57.00 | 63.99 | 5.0× |

The archived **PID series already drove this same plant model to 57.26 A** (and down to 30 A), so
widening the MPC bound toward ~57.0 A uses a region the model has already been through — a
correction of an artificially narrow MPC constraint, not invented physics.

### Fixed-state verification at `u_max = 56.90` (`verify_widened_bounds.m`)

| `u` [A] | mass flow | `J_error` | `J_psm` | violations |
|---|---|---|---|---|
| 56.7250 | 7.747 | 0.82 | **30.72** | 10 ← unpenalised optimum |
| **56.7425** | **11.330** | **5.12** | **11.96** | **4** ← penalised optimum |
| 56.7600 | 14.912 | 18.28 | 6.95 | 3 |
| 56.9000 | 43.550 | 421.87 | 0.11 | 1 |

```
shift          = +0.01750 A    extra coolant = +3.583 kg/s (1.46x)
peak DSM       : +0.005095 -> +0.004857
J_psm share    : 70.0%
penalised optimum at a bound: NO -> INTERIOR
PCT authority  : 1.1077 K  (was 0.9749 K with u_max 56.75)
```

`J_error` rises monotonically and `J_psm` falls monotonically, so the penalised optimum is a
genuine balance between them, and threshold violations drop from 10 to 4. **With the pump given
room, the safety term produces an interior optimum instead of pinning the actuator at its bound.**

---

## Status

- [x] Past-instant contamination: verified absent
- [x] Control authority: confirmed (0.06 A → 12.28 kg s⁻¹ → 0.975 K PCT)
- [x] Penalty reformulated on the DSM (`fitness_cost_dsm.m`), production cost untouched
- [x] Result A (t=1180 fixed state): reorientation demonstrated, with a stated caveat
- [x] Result B (t=1161): redundancy explained structurally
- [x] Result C (severe HRM scenario): valid, `B` controls severity; still saturates at `u_max`
- [x] Result D (61-query closed loop, archived bounds): penalty dominates 36/61 queries but the
      trajectory is unchanged — saturation is the binding constraint
- [x] Result E (fixed state, `u_max = 56.90`): penalised optimum becomes interior
- [ ] **Running**: 61-query closed loop with `u_max = 56.90` and `lambda = 1e6`
      (`widepen_umax56p90_LJ_seed1.mat`), to be compared against the archived-bounds baseline.
      Note this comparison varies two things at once (bound + penalty); Result E is what licenses
      attributing a trajectory change to the penalty having room to act.
- [ ] Optional: baseline rerun with a different seed, to size optimizer noise against the
      penalty's effect

## Figures produced

| file | content |
|---|---|
| `pen_pct_vs_tq` | PCT, full range + zoom past the transient |
| `pen_dsm_vs_tq` | DSM in K and normalised, with the threshold |
| `pen_u_vs_tq` | applied ACS current with `u_min`/`u_max` |
| `pen_mdot_vs_tq` | coolant mass flow + extra injected per query |
| `pen_jterms_vs_tq` | `J_error` / `J_mdot` / `J_psm` / `J` over time |
| `pen_dominance` | share of the cost and `J_psm`/`J_error` ratio |
| `penalty_comparison` | composed 3×2 overview |
| `solver_kpi_bars`, `solver_actuation`, `solver_comparison` | LJ/GA/PID, archived KPIs only |

Orange bands in the penalty figures mark the windows where `J_psm > J_error`.

## Solver comparison (separate axis)

Kept apart on purpose: this asks *which optimizer solves the same problem better*, while the
penalty figures ask *what the constraint does to the trajectory*. Mixing them would confound the
two effects.

| metric | LJ | GA | PID |
|---|---|---|---|
| min safety margin (K) | 62.786805 | 62.786805 | 62.554817 |
| peak PCT (K) | 1181.213195 | 1181.213195 | 1181.445183 |
| IAE (K·s) | 36.51 | 54.47 | **17.76** |
| settling time (s) | 1186.5 | inf | 1177.5 |
| control total variation (A) | 0.858 | **0.433** | **28.404** |
| peak control rate (A s⁻¹) | 0.111 | **0.054** | **53.360** |

LJ and GA agree to the 15th digit on both safety metrics — with a single-node DET they
reconstruct the same root trajectory, so the optimizer choice cannot change the thermal outcome.
PID tracks best in IAE but pays **33× the total variation and 481× the peak rate**, and it cannot
carry a `J_psm` term at all: a Kp/Ki/Kd law does not minimise a cost.

---

## Closing position: what the DSM penalty is claimed to do

This is the conclusion carried into the defense deck (slides 51–53 of
`phd defense/main_editable_v3.pptx`). The defensible comparison is **A baseline (λ=0) against
D retuned (λ=10⁴)**. Configuration C (λ=10⁶) is retained in this record because it is what
exposed both defects, but it is *not* presented as a result: it is a diagnostic step, and
showing three series saturated the panels.

| | A baseline | D retuned |
|---|---|---|
| peak PCT | 1181.21 K | 1181.21 K |
| minimum DSM margin | 62.79 K | 62.79 K |
| IAE | 35.16 K·s | **28.65 K·s** |
| steady-state error | **+0.07 K** | +0.70 K |
| control total variation | **0.86 A** | 2.54 A |
| median optimizer time / query | 351.8 s | **76.0 s** |

**Three independent claims, and no more:**

1. *The envelope does not move.* Peak PCT and minimum margin are identical. The soft
   constraint changes how the loop reaches the setpoint, not what it protects. This is a
   preservation claim, not an improvement claim.
2. *Tracking improves.* IAE falls 18.5 %, at the cost of a busier actuator (TV 0.86 → 2.54 A).
   The trade is stated, not hidden.
3. *Real-time feasibility improves.* The corrected LJ solves 4.6× faster, inside the
   per-query budget.

### An integrated-margin KPI was considered and rejected

The obvious way to argue "D is safer on average" is to accumulate the margin,

```
IM = ∫ (U_y − PCT) dt        A: 1941.2 K·s     D: 1964.0 K·s   (+22.8, +1.17 %)
```

It was implemented, plotted, and then removed. Two reasons, both disqualifying:

- **It is not independent of the IAE.** `U_y` is constant, so `IM = U_y·T − ∫PCT dt` is an
  affine transform of the mean PCT. With D sitting below the setpoint for most of the window,
  IM and IAE are the same observation in different units — not corroboration.
- **It rewards the known defect.** D accumulates more margin *because* it settles below the
  setpoint, which is the residual −34 mK/s over-cooling that the actuator limitation section
  above explains. Presenting that as a safety benefit invites the obvious rebuttal.

If a genuinely independent companion to the IAE is wanted later, the candidate is the time the
margin spends compromised,

```
T_exp = ∫ 1[DSM > DSM_thr] dt
```

which measures what the soft constraint actually targets and cannot be improved by
over-cooling. It has not been computed.

### Cost function as implemented

Slide 41 of the deck previously showed only the tracking term. The implemented objective
(`fitness_cost.m`, lines 102–161) is:

```
J       = J_error + J_psm + J_mdot
J_error = sum_j ( yhat(t+j|t) − w )^2
J_psm   = lambda_DSM  * sum_j max( 0, DSM(t+j) − DSM_thr )^2      DSM = (Tp − w)/(U_y − w)
J_mdot  = lambda_mdot * sum_j max( 0, |dm(t+j)| − dm_lim )^2
```

with `lambda_mdot = 1e-4`, `dm_lim = 5 kg/s`, and for run D `lambda_DSM = 1e4`,
`DSM_thr = 0.00237737`. The deck now states all three terms with this nomenclature.
