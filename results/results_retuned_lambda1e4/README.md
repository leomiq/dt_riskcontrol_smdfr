# Retuned LJ with a moderated penalty — final configuration

**Question:** can the MPC protect the safety margin *and* converge toward the
setpoint, instead of over-cooling?

Three configurations, all LJ, seed 1, same accident stream:

| | `u_max` | `lambda` | optimizer |
|---|---|---|---|
| **A** baseline | 56.75 | 0 | LJ original |
| **C** over-penalised | 56.90 | 1e6 | LJ original |
| **D** retuned | 56.90 | **1e4** | **LJ retuned** |

`dsm_thr = 0.00237737` in C and D. C is kept in the figures because it is what
exposed both defects; showing A/C/D together makes the correction legible.

## Results

| metric | A (`λ`=0) | C (`λ`=1e6) | **D (`λ`=1e4)** |
|---|---|---|---|
| **IAE (K·s)** | 35.16 | 34.59 | **28.65** ← best |
| final PCT (K) | 1179.156 | 1177.545 | 1178.426 |
| final error (K) | −0.115 | −1.726 | **−0.845** |
| **final slope (mK/s)** | −3.9 | **−152.0** | **−34.1** |
| steady-state error (K) | +0.074 | +1.270 | +0.678 |
| control TV (A) | 0.864 | 3.157 | 2.543 |
| peak rate (A s⁻¹) | 0.120 | 0.420 | 0.330 |
| interior control (%) | 88.5 | 50.8 | 57.4 |
| total coolant (kg s⁻¹) | 559.5 | 792.0 | 650.5 |
| penalty-active queries | — | 17 | 12 |
| peak PCT / min margin | identical across all three | | |
| wall-clock | — | 4.01 h | **1.55 h** |

**D improves every failure mode of C**: the drift is 4.5× smaller, the final
error halves, actuation is smoother, and it computes in a third of the time —
while achieving the best IAE of all three.

## Two independent defects were found and fixed

### 1. LJ was returning sub-optimal actions

Its pass reset was `r = |u_best - old_u_best|`, which discarded the within-pass
contraction and let the radius grow again, so it never resolved the narrow
interior minimum the penalised cost has (~0.05 A wide, cost ×179 at its edges).
The early-stop test was absolute (`improvement < 1`) while late-horizon costs
reach ~3000.

Measured against a 22-point grid reference on the four queries where the drift
appeared:

| `t_q` | LJ original | LJ retuned |
|---|---|---|
| 1176 | +161% | −66% |
| 1180 | **+18933%** | −29% |
| 1184 | +584% | −73% |
| 1188 | +17832% | −91% |

(Negative gaps mean LJ beat the discrete grid — it searches in continuum, so the
true minimum is deeper and narrower than 22 points could resolve.)

Fixes in `luus_jaakola_det.m`: geometric ceiling on the radius
(`r = min(|Δu_best|, r_init·γ^(pass−1))`) and a relative early stop (0.1 % of
`f_best`). In `run_phase2c_full_loop.m`: `γ` 0.95 → 0.90, `R` 15 → 24, to suit a
range that grew 3.5× (0.06 → 0.21 A).

### 2. `lambda = 1e6` was too large

With the retuned LJ, the true optimum at the late queries was still `u_max` —
maximum coolant even with the PCT already below setpoint. So the over-cooling was
**not** only an optimizer artefact. Fixed-state check with the retuned optimizer:

| `t_q` | `λ`=1e4 | `λ`=1e5 |
|---|---|---|
| 1180 | `u*` = 56.762, `J_psm` 0.3 % | `u*` = 56.751, `J_psm` 7.7 % |
| 1188 | `u*` = 56.690 (`u_min`, correct) | `u*` = 56.863, **`J_error` = 266** |

`λ`=1e5 pays a large tracking error to remove a negligible penalty — worse
calibrated than 1e6. `λ`=1e4 keeps the penalty significant (`J_psm` and
`J_error` within the same order of magnitude) without forcing saturation.

## Residual drift: an actuator limitation, not a tuning one

D still drifts at −34 mK/s and `settling_time` remains infinite. The cause is
physical, not algorithmic:

```
u_min = 56.69 A  →  0.580 kg/s
```

**`u_min` is not zero coolant.** Even at the minimum command the ACS keeps
injecting 0.58 kg s⁻¹, and **the ACS can only cool** — there is no actuator that
heats. Once the PCT is driven below the setpoint during the acute phase, the loop
cannot recover in that direction. A never over-cools that far, so 0.58 kg s⁻¹
suffices to hold it.

A minor residual remains: LJ occasionally selects ~56.81 (22 kg s⁻¹) while the
PCT is already 0.6 K below setpoint, which `J_error` alone does not justify.
Small compared with what was corrected, but it suggests `γ = 0.85` could tighten
things further.

## Figures (`figures/`, PDF vector + PNG)

Orange bands mark queries where D's penalty is active.

| file | content |
|---|---|
| `ret_pct_vs_tq` | PCT, full range + zoom past the transient |
| `ret_dsm_vs_tq` | normalised DSM with the threshold |
| `ret_u_vs_tq` | applied current with `u_min` / both `u_max` |
| `ret_mdot_vs_tq` | coolant mass flow (log) |
| `ret_jterms_vs_tq` | objective decomposition, C vs D side by side |
| `ret_kpi_bars` | IAE, steady-state error, TV, safety margin |
| `ret_lj_quality` | optimizer gap vs grid, original vs retuned (log) |
| `retuned_comparison` | composed 4×2 overview |

## Data

- `phase2c_fullloop_pilot_LJ_seed1.mat` (A), `widepen_umax56p90_LJ_seed1.mat` (C),
  `retuned_lam1e4_umax56p90_LJ_seed1.mat` (D) + D's per-query checkpoint
- `widepen_only.log`, `retuned_lam1e4.log` — diaries with every cost evaluation
- `verify_lj_suboptimality_*.mat`, `verify_lj_retuned_*.mat` — the grid-reference
  evidence for the optimizer fix
- `retuned_comparison.npz` — plotting arrays

## Regenerate

```bash
python plot_retuned_comparison.py
```

Related: `../results_penalty_saturation/` (how A and C were obtained and why the
pump bound was widened), `../results_solver_comparison/` (LJ/GA/PID, a separate
axis).
