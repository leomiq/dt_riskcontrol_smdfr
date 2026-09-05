# DSM penalty and pump saturation — LJ baseline vs penalised

**Question:** does the deterministic safety-margin penalty actually reshape the
PCT trajectory, and does it need room in the actuator to do so?

All runs use the **same solver (LJ), same seed (1), same accident stream**. Kept
separate from `../results_solver_comparison/` so *the effect of the constraint* is
never confounded with *the effect of the search algorithm*.

Penalty form (on the DSM itself, so `lambda` is dimensionless):

```
J_psm = lambda * sum( max(0, DSM - dsm_thr)^2 )
DSM   = (Tp - target) / (U_y - target),   U_y = 1244 K,  target = 1179.2712 K
```

## The three configurations

| # | run | `u_max` | `lambda` | file |
|---|---|---|---|---|
| A | baseline | 56.75 | 0 | `phase2c_fullloop_pilot_LJ_seed1.mat` |
| B | penalised, archived bounds | 56.75 | 1e6 | `penalised_thr0p00238_lam1e06_LJ_seed1.mat` |
| C | penalised, widened pump | **56.90** | 1e6 | `widepen_umax56p90_LJ_seed1.mat` |

`dsm_thr = 0.00237737` in B and C.

## Finding 1 (B vs A): the penalty dominates the cost but cannot move the trajectory

`J_psm` outweighed `J_error` over **36 of 61 queries**, in one contiguous window
`t = 1160.0 … 1177.5 s`, with the ratio reaching **197×**. The crossing at 1177.5
is abrupt and coincides exactly with the normalised DSM falling below the
threshold — the constraint engages while the margin is compromised and withdraws
on its own.

Yet peak PCT, max DSM and min safety margin came out **identical to four
decimals**. Two structural reasons:

1. The PCT peak occurs at `t = 1160.0`, the **first** query, before any control
   action has had an effect.
2. The control was **saturated at `u_max`** through most of the dominance window,
   so the safety term had already obtained everything it could ask for.
   Saturation, not the penalty, was the binding constraint.

## Finding 2: the archived bounds sit at the pump's onset

From the ECM pump model (`CurrentToMassFlow.m`, Nowak et al. 2022,
DOI 10.1016/j.anucene.2022.109142):

```
zero-flow current  DeltaP/C       = 56.6872 A
archived bounds    [56.69, 56.75] = only 0.0028 .. 0.0628 A above onset
slope              dmdot/dI       ~ 204.7 kg/s per A
```

| `u_max` | max mass flow | factor |
|---|---|---|
| 56.75 (archived) | 12.86 kg s⁻¹ | 1× |
| **56.90** (used in C) | **43.55** | **3.4×** |
| 57.00 | 63.99 | 5.0× |

The archived **PID series already drove this same plant model to 57.26 A**, so
widening toward ~57 A uses a region the model has been through — a correction of
an artificially narrow MPC constraint, not invented physics.

Fixed-state check at `u_max = 56.90` (`verify_widened_bounds.m`): the penalised
optimum becomes **interior** (56.7425 A, +3.58 kg s⁻¹, `J_psm` 70% of the cost)
instead of pinned at the bound. **This is what licenses reading C's trajectory
change as the penalty acting**, since C varies two factors at once.

## Finding 3 (C vs A): with room in the pump, the trajectory does change

| metric | A (`λ`=0, 56.75) | C (`λ`=1e6, 56.90) | delta |
|---|---|---|---|
| peak PCT (K) | 1181.2132 | 1181.2132 | +0.0000 |
| min safety margin (K) | 62.7868 | 62.7868 | +0.0000 |
| final PCT (K) | 1179.156 | **1177.545** | **−1.611** |
| IAE (K·s) | 35.16 | 34.59 | −0.58 |
| steady-state error (K) | 0.074 | **1.270** | +1.196 |
| settling time (s) | 1176.0 | **inf** | — |
| control total variation (A) | 0.864 | **3.157** | +2.29 |
| peak control rate (A s⁻¹) | 0.120 | 0.420 | +0.30 |
| total coolant (kg s⁻¹, Σ) | 559.5 | **792.0** | **+232.5 (+42%)** |
| max mass flow (kg s⁻¹) | 12.86 | **43.55** | 3.4× |

The normalised DSM falls below the threshold at **t ≈ 1163** in C versus
**t ≈ 1177** in A — the margin is secured ~14 s earlier. `J_psm` is then active in
only **17 of 61** queries (`t = 1160.0 … 1169.5`) and switches off.

**Interpretation:** the margin is prioritised over setpoint-tracking accuracy, and
the price is over-cooling. The final PCT sits 1.6 K *below* the setpoint, the
steady-state error is 17× larger, the loop never re-enters the tolerance band, and
actuation cost rises 3.6×. Conservative for safety, expensive in operation.

Two things did **not** change: peak PCT and the minimum safety margin — the peak is
set before any control action can influence it.

Side effect: `J_mdot` (flow-rate soft constraint) is active in **61 of 61** queries
in C, having been identically zero with the archived bounds — the wider range makes
step-to-step flow changes exceed `delta_limit = 5 kg s⁻¹`.

## Figures (`figures/`, PDF vector + PNG preview)

Orange bands mark the queries where the penalty is active.

**Full per-indicator set (A vs C), mirroring the archived layout:**
`set_pct_vs_tq`, `set_dsm_vs_tq`, `set_dsmn_vs_tq`, `set_u_vs_tq`,
`set_mdot_vs_tq`, `set_jterms_vs_tq`, `set_fopt_vs_tq`, `set_runtime_vs_tq`,
`set_kpi_bars`, and the composed **`penalty_full_comparison`** (4×2).

**A vs B (archived bounds):** `pen_pct_vs_tq`, `pen_dsm_vs_tq`, `pen_u_vs_tq`,
`pen_mdot_vs_tq`, `pen_jterms_vs_tq`, `pen_dominance`, `penalty_comparison`.

**A vs C (focused):** `wide_pct_vs_tq`, `wide_u_mdot_vs_tq`, `wide_jterms_vs_tq`,
`widened_comparison`.

## Data

- Realization MATs (A, B, C) and their per-query checkpoints
- `penalised_run3.log`, `widepen_only.log` — diaries carrying every cost
  evaluation, from which the `J_i` decompositions are parsed
- `figure1_penalty_effect.npz`, `figure_widened_comparison.npz`,
  `penalty_full_set.npz` — plotting arrays

## Regenerate

```bash
python plot_penalty_full_set.py      # A vs C, full indicator set (main)
python plot_penalty_comparison.py    # A vs B
python plot_widened_comparison.py    # A vs C, focused
```

## Open

- A widened **baseline** (`λ`=0, `u_max`=56.90) would separate the bound's effect
  from the penalty's cleanly; not run (~4 h).
- Per-query differences alternate sign with amplitudes up to ±12 kg s⁻¹, so with
  one realization per configuration the penalty's contribution cannot be
  separated from LJ's stochastic noise. A baseline rerun with a different seed
  would size that noise.
