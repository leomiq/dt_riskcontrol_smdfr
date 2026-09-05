# Solver comparison — LJ / GA / PID (archived realization)

**Question:** which optimizer solves the same MPC problem better?

Single archived closed-loop realization per solver, 61 query instants
`t_q = 1160 … 1190 s`, identical scenario. No DSM penalty anywhere: this axis is
purely about the search algorithm.

Kept separate from `../results_penalty_saturation/` on purpose. Mixing them would
confound *the effect of the search algorithm* with *the effect of the safety
constraint*.

## Results

| metric | LJ | GA | PID |
|---|---|---|---|
| min safety margin (K) | 62.786805 | 62.786805 | 62.554817 |
| peak PCT (K) | 1181.213195 | 1181.213195 | 1181.445183 |
| overshoot (%) | −2.7784 | −2.7739 | −1.3341 |
| settling time (s) | 1186.5 | inf | 1177.5 |
| IAE (K·s) | 36.51 | 54.47 | **17.76** |
| steady-state error (K) | −0.2799 | −0.5986 | +0.0024 |
| control total variation (A) | 0.858 | **0.433** | **28.404** |
| peak control rate (A s⁻¹) | 0.111 | **0.054** | **53.360** |

Two facts the figures make legible:

1. **LJ and GA agree to the 15th digit** on `min_safety_margin`
   (62.786805032494385 K) and `peak_PCT` (1181.2131949675056 K). With a
   single-node DET both reconstruct the same root trajectory, so the optimizer
   choice cannot change the thermal outcome in this scenario. They differ only in
   control smoothness (TV 0.858 vs 0.433) and in convergence of the objective.
2. **PID is the non-optimised reference.** It tracks best in IAE but pays **33×
   the total variation and 481× the peak actuation rate**, and it cannot carry a
   `J_psm` term at all — a Kp/Ki/Kd law does not minimise a cost, so there is
   nowhere to insert a safety penalty.

## Figures (`figures/`, PDF vector + PNG preview)

| file | content |
|---|---|
| `pct_vs_tq` | PCT(t_q), all three solvers |
| `dsm_vs_tq` | DSM(t_q) in K |
| `u_vs_tq` | applied ACS current (PID off-scale values annotated) |
| `fopt_vs_tq` | objective `f*(t_q)`, LJ vs GA |
| `runtime_vs_tq` | optimizer runtime with the 120 s reference |
| `variance_1180` | 20-seed distributions at the fixed state `t_q = 1180 s` |
| `phase2c_comparison` | composed 3×2 overview (thesis) |
| `solver_kpi_bars` | safety / peak / tracking KPI bars |
| `solver_actuation` | total variation and peak rate (log scale) |
| `solver_comparison` | composed 2×2 KPI overview |

## Data

- `phase2c_archived_trajectory_kpis.csv` — closed-loop KPIs per solver
- `phase2c_archived_trajectory_per_instant.csv` — per-query `u`, PCT, `f_opt`, runtime

## Regenerate

```bash
python plot_phase2c_trajectories.py    # per-indicator set + composed overview
python plot_solver_comparison.py       # KPI bars
```

Both write into this folder's `figures/`. Sources are the archived series in
`Models for training (v01-07)\Version5\GBM_RMM_CINECA\{LJ,GA} test 11_10` plus the
root PID series; `extract_archived_trajectory.py` rebuilds the two CSVs from them.
