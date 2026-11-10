# Grey-Box Digital Twin for Risk Monitoring and Predictive Control of an SMDFR

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22407139.svg)](https://doi.org/10.5281/zenodo.22407139)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Reference implementation accompanying the PhD thesis *A Grey-Box Digital
Twin-based Approach for Condition-informed Risk Monitoring and Control of
New-Generation Nuclear Power Plants* (Politecnico di Milano).

The code covers the full loop studied in the thesis: a reduced-order model of a
Small Modular Dual Fluid Reactor (SMDFR), a grey-box virtual sensor for the peak
cladding temperature (PCT), an Importance-Sampling Monte Carlo Dynamic Event
Tree (IS-MCDET) for dynamic risk assessment, and a model predictive controller
whose objective carries a deterministic-safety-margin (DSM) soft constraint.

## What the loop does

The case study is a small-break loss-of-coolant accident (LOCA). The lead coolant
leaks, the PCT rises, and above 1180 K the auxiliary cooling system (ACS) is
demanded. The ACS is driven by a magnetohydrodynamic (MHD) pump whose input
current `u` [A] sets the injected coolant flow. The safety limit is
`U_y = 1244 K`; the controller tracks the ACS setpoint `y_ref = 1179.27 K`.

Every `Δt = 0.5 s`, between t = 1160 s and 1190 s (61 decisions), the controller:

1. reads the plant state, here simulated by the high-resolution model (HRM);
2. forecasts the leak flow and the ROM's PCT error with a BiLSTM run under MC
   dropout (59 trajectories; the median is used);
3. propagates the reduced-order model (ROM) over a 25-step (12.5 s) horizon
   inside a dynamic event tree that branches when the ACS state changes
   (demanded, available, failed);
4. searches the 26 future currents with Luus–Jaakola (a genetic algorithm and a
   PID controller serve as baselines), applies only the first one, and repeats.

```
 HRM plant ──► BiLSTM (MC dropout) ──► leak forecast d̂, ROM correction Δ̂
                                           │
 candidate u ──► MHD pump map ──► ROM + DET ──► PCT estimate ŷ = h(x) + Δ̂ ──► cost J(u)
                                                                                  │
                                                  Luus–Jaakola / GA ◄────────────┘
```

## The control problem as implemented (P0)

$$
\min_{u_0,\dots,u_N}\;
\underbrace{\sum_{j}(\hat y_j-y_{\rm ref})^2}_{\text{tracking}}
+\lambda_m\underbrace{\sum_k\big[\,|\dot m_{k+1}-\dot m_k|-\Delta\dot m\,\big]_+^2}_{\text{actuator smoothness}}
+\lambda_s\underbrace{\sum_j\big[\,\mathrm{DSM}_j-\tau\,\big]_+^2}_{\text{margin erosion}}
\quad\text{s.t.}\quad u_{\min}\le u_k\le u_{\max}
$$

with `DSM = (ŷ − y_ref)/(U_y − y_ref)`, `Δṁ = 5 kg/s`, `λ_m = 1e-4`,
`τ = 0.00238` and `λ_s ∈ {0, 1e4}` (see `matlab/fitness_cost.m`). The
probabilistic safety margin (PSM) from IS-MCDET is computed after the decision,
in `matlab/DET_scheduler_postprocessing_perc_online.m`. It does not enter the cost.

## Results at a glance

Closed loop, 61 queries, Luus–Jaakola with seed 1, same accident stream
(`results/results_retuned_lambda1e4/`):

| | A: baseline (`λ_s = 0`, `u_max = 56.75 A`) | D: retuned (`λ_s = 1e4`, `u_max = 56.90 A`) |
|---|---:|---:|
| Peak PCT / minimum margin | 1181.21 K / 62.79 K | 1181.21 K / 62.79 K |
| IAE | 35.16 K·s | 28.65 K·s |
| Steady-state error | +0.07 K | +0.70 K |
| Control total variation | 0.86 A | 2.54 A |
| Median optimizer time per query | 351.8 s | 76.0 s |

The safety envelope does not move, because the PCT peak occurs at the first
query, before any control action takes effect. Tracking improves at the cost of
a busier actuator, and the corrected optimizer is 4.6× faster. The margin
penalty can only move the decision when the tracking optimum is interior: in
the acute phase both objectives ask for maximum coolant and the action
saturates. Details are in `docs/DSM_PENALTY_RESULTS.md`.

## Research extension (2026): risk-informed stochastic MPC

The continuation asks whether the *probability* that the ACS fails should change
the control action. P0 cannot answer that, for four reasons:

1. Its cost concatenates the nodes of one event tree without path probabilities.
2. The PSM is computed after the decision.
3. In a 12.5 s window the tree rarely branches, so ACS failures are invisible to the optimizer.
4. The DSM describes one (median) trajectory.

The target formulation (P1) replaces the cost with an expectation over weighted
scenarios plus a risk measure ρ (exceedance probability, quantile PSM or CVaR),
optionally with a chance constraint, over a risk horizon longer than the
tracking horizon:

$$
\min_{\pi}\ \mathbb E_p\Big[\sum_k \ell(\hat y_k,\dot m_k)\Big]+\lambda_R\,\rho(Z^\pi)
\quad\text{s.t.}\quad \mathbb P_p\big(\max_{k\le N_R}\hat y_k>U_y\big)\le\varepsilon,
\quad u_k=\pi_k(\text{observed information}).
$$

A technical discussion on 22 October 2026 changed the order of work. Even the
deterministic MPC is not solved with margin: the median is 76 s of a 120 s
response window. A probabilistic margin evaluated for every candidate would
multiply that cost by about 59, the Wilks sample size for 95/95. So the problem
must first be made tractable, and only then made stochastic.

Exploration avenues, in execution order:

| | Avenue | Idea | Priority |
|---|---|---|---|
| V1 | Monolithic formulation | Write P0 as a standard optimization problem: dynamics as equality constraints, slack variables for the squared hinges, the margin as a hard linear constraint `ŷ_k ≤ y_ref + τ(U_y − y_ref)`. Solve it with a conventional (QP, conic or NLP) solver instead of a metaheuristic. | high |
| V2 | Linearized reactor model | Identify `δx⁺ = A δx + B δṁ + E δd̂`, `ŷ = ȳ + C δx + Δ̂` from ROM responses around one or more operating points. The result is a linear-quadratic MPC (a QP) that supports robust, distributionally robust and contextual MPC, plus chance constraints with closed-form tightening. | high |
| V3 | Exact reference | Solve the same query states with Luus–Jaakola, the genetic algorithm and the solver from V1–V2, re-simulate the actions in the ROM, and measure the metaheuristics' optimality gap and runtime. | high |
| V4 | Risk in every candidate | Chance constraint, or `Ĵ + λ_R ρ̂` over common weighted scenarios. | high |
| V5 | Extended risk horizon | Evaluate ρ over `N_R ≫ N` so that ACS failures become visible. | high |
| V6 | Decision-focused scenario | Learn offline a map from context to a single scenario such that the deterministic MPC fed with it reproduces the stochastic decision. The hard part moves offline. Conceptual for now. | medium |
| V7 | Nonlinear-model compute | Move blocking and warm starts, if the ROM has to stay in the loop. | medium |
| V8 | Multi-stage MPC | Different actions after observing ACS success or failure. | medium |
| V9 | Several protection systems | With one active system a failure has no backup; with two, activation order matters. | medium |
| V10 | Bounds for the nonlinear model | Rigorous ROM enclosures over control boxes. | low |

The actuator part of V1–V2 is already linear. Over `u ∈ [56.69, 56.90] A`,
the MHD pump map is `ṁ ≈ 204.6 (u − 56.687)` with an error below 0.007 kg/s,
so the coolant flow can be the decision variable. The remaining nonlinearity is
the thermal dynamics of the ROM, which includes the heat-transfer lag.
The current bound itself has no physical meaning: the pump ECM is sized as a
main-circuit pump, and the ACS make-up flow uses the foot of its curve. The bound
is therefore set in kg/s, and an actuator ramp is added (plan §8.3, §9.6).

Stages:

- **E0**: freeze a reproducible baseline.
- **E1**: build the monolithic problem and the linear model (V1–V2). Close when
  the linear prediction error is well below 1 K, the control authority over the
  horizon.
- **E2**: build the exact reference (V3).
- **E3**: put risk into the decision (V4–V5).
- **E4**: validate in closed loop under 120 s per query, with V6 as the online
  option.
- **E5**: extend to multi-stage decisions and several protection systems (V8–V9).

Certified bounds for the nonlinear model (V10) stay secondary.

The groundwork for V1–V3 is laid out in a separate plan. The ROM used as the MPC
prediction model is a 0-D model with 16 states: 6-group point kinetics with circulating
precursors, fuel, wall and coolant temperatures, and coolant inventory. The plan writes it
once, algebraically, and derives four levels from it: steady state (N0), linearization at
one point (N1), linearization along the trajectory (N2) and the full nonlinear model as
constraints (N3). It then checks when the optimal-scenario result of Homem-de-Mello et al.
(2026, *Math. Program.*, doi:10.1007/s10107-026-02359-9) applies to the forecast
uncertainty.

Status (9 Nov 2026, plan §8.1–8.5; summary with figures in
[docs/estado_mpc_monolitico.pdf](docs/estado_mpc_monolitico.pdf), Spanish):

- **F1–F3 closed.** The algebraic ROM replica in `julia/` reproduces Simulink to
  about 1e-6 K. The monolithic NLP solves P0 in about 0.05 s at t_q = 1161 and 1180 s.
  The linearized QP is written in sparse and condensed form (77 variables, ~10 ms);
  Clarabel certifies the QP optimum with a dual gap of ~1e-12.
- **F3c, optimality of the nonlinear P0.** The NLP solution is a strict local minimum
  (KKT and second-order conditions). It is the only optimum found in ~590 multistarts
  and in a 90k-evaluation differential evolution. A certified lower bound (Taylor model
  with an interval remainder) puts the global optimum within 2.3 % (1180 s) and 11.4 %
  (1161 s) of the NLP solution; the sampled, semi-rigorous gap is 0.3 % and 1.2 %.
  With every decision evaluated on the same nonlinear cost, the NLP and the iterated
  QP (N2) tie as best, the one-point QP (N1) is 0.0026 % behind, and the corrected
  Luus–Jaakola and GA are at least 28 % and 43 % above the global optimum
  (`julia/f3c_*.jl`, `results/f3c_optimality/`).
- **F4, first result.** The BiLSTM MC-dropout ensemble was regenerated without MATLAB
  (`python/f4_leak_ensemble.py`). The forecast is biased and overconfident: 90 % band
  coverage is 0.44. Its effect on the PCT is still only about 0.025 K, so neither
  quantile decisions nor a two-stage SAA change the first action, with or without
  actuator saturation (`julia/f4_risk_decisions.jl`, `julia/f4_actuator_sweep.jl`).
- **Next steps.** Set the injection ramp and flow limits from the ACS specification.
  Run ACS-failure scenarios and the ROM–HRM mismatch as the risk sources that can move
  the decision. Then run the closed loop on the 61 queries. Tightening the certified
  gap (plan §8.5, step 4b.5) is scheduled for the next round.

- [Monolithic formulation plan](docs/monolithic_formulation_plan.md) and
  [its presentation](presentations/monolithic_formulation.pdf) (Spanish)
- [Improvement plan and exploration avenues](docs/implementation_route.md) (Spanish)
- [Presentation, USM, October 2026](presentations/usm_risk_control.pdf) (Spanish)
- [Technical report on the cost lower bound](docs/technical_report.pdf) and
  [its results](docs/m1_cost_bound_results.md) (Spanish)

`matlab/certified_cost_lower_bound_box.m` and `matlab/m1_branch_and_bound.m`
implement a *conditional* lower bound on the P0 cost and a box branch-and-bound.
They close the gap on a two-mode test problem. For the real ROM/DET the only
lower bound available is `J ≥ 0`, because rigorous thermal enclosures are not
built yet.

## Layout

```
matlab/          DET engine, MPC loop, optimizers (Luus-Jaakola, GA, PID),
                 physics/parameters, verification and sweep scripts
matlab/mex_sources/   C sources for the seven Level-2 S-functions
models/          Simulink models, trained BiLSTM weights, normalisation and
                 fixed-state inputs needed to run the pilot
python/          MC-dropout inference bridge, training pipeline, plotting and
                 KPI-compilation scripts; F4 ensemble and calibration (f4_*.py)
julia/           Algebraic ROM replica, monolithic NLP (F2), linearized QP (F3),
                 risk and actuator studies (F4)
results/         Arrays and tables that regenerate the thesis figures without
                 re-running the simulations
docs/            DSM-penalty record, improvement plan, reproducibility notes,
                 cost-bound report
presentations/   Beamer source and PDF of the research presentation
```

## Requirements

Step-by-step setup on a new Windows PC: `docs/setup_entorno.md`; check with
`tools/check_env.ps1`.

- MATLAB R2024a/R2024b with Simulink (Parallel Computing Toolbox recommended), and a
  C compiler for `mex` only to rebuild the S-functions (verified with MinGW64)
- Julia 1.10.4 for the monolithic NLP/QP work (`julia/`, pinned by `Manifest.toml`)
- Python 3.11 (3.9–3.12 for MATLAB's `pyenv`) with `python/requirements.txt`;
  `python/requirements-bilstm.txt` (torch) additionally for the grey-box inference path
- Windows x64 for the pre-compiled `.mexw64` binaries. On any other platform,
  rebuild them from `matlab/mex_sources/` first (see below).

The ROM initialization reads the five archived polynomial coefficients from
`matlab/fit_poly_coefficients.m` instead of deserializing the `cfit` object in
`models/fit_poly.mat`. The pilot therefore no longer needs the Curve Fitting
Toolbox. To check the values against the original MAT archive, run
`python python/verify_fit_poly_coefficients.py`.

**Known reproducibility gap.** Under MATLAB R2024b the pilot returns
`f = 0.362026`, against the archived `0.326151` from R2024a (+11.0 %). The cause
is not isolated yet, so new methods should not be compared against the archived
value until it is. See `docs/pilot_reproducibility.md`.

## Building the S-functions

The Simulink model loads seven Level-2 S-functions as compiled MEX. Pre-built
Windows binaries are included; to rebuild:

```matlab
cd matlab/mex_sources
mex -largeArrayDims msfcn_indirectps_hrm_c.c
mex -largeArrayDims msfcn_indirectps_rom_v2.c
mex -largeArrayDims msfcn_indirectps_rom_v3.c
mex -largeArrayDims msfcn_limintm_hrm_c.c
mex -largeArrayDims msfcn_limintm_rom_v2.c
mex -largeArrayDims msfcn_schedule_hrm_c.c
mex -largeArrayDims msfcn_schedule_rom_v2.c
```

Move the resulting `.mexw64` (or platform equivalent) next to the model.
Rebuilt binaries were checked against the archived ones: the PCT trajectories
agree to 9.3e-06 K (relative 7.9e-09) and the first cost evaluation is
bit-identical. See `matlab/mex_sources/README.md`.

## Running

A single fixed-state optimisation at the query instant `t_q = 1180 s`
(~2 minutes, the quickest way to confirm the installation works):

```matlab
addpath('matlab'); addpath('models');
pilot_single_run
```

A full closed-loop realization over the 61-instant grid:

```matlab
run_phase2c_full_loop
```

Be aware of the cost: roughly 350 s per query for the baseline configuration
and about 6 hours for a complete realization. The retuned optimizer brings the
median down to ~76 s per query.

To regenerate the thesis figures from the stored results, without simulating:

```bash
python python/plot_retuned_comparison.py
python python/compile_thesis_kpis.py
```

## Configuration

Paths that point at bulk data not redistributed here are read from environment
variables, each falling back to a repository-relative default:

| variable | used for |
|---|---|
| `P2C_RUN_DIR` | directory holding run artefacts to compare |
| `P2C_OUT_DIR` | where comparison output is written |
| `ARCHIVE_DIR` | archived per-solver sweep directories |
| `TRAINING_DB_DIR` | failure-injection training databases (`.xls` batches) |
| `FAILURE_INJECTION_DIR` | failure-injection database for the `Tp`/`mdot_c` fit |

## What is not included

- The failure-injection training databases and the archived per-solver sweeps.
  These are bulk raw output; the trained weights in `models/` and the summary
  arrays in `results/` are shipped instead.
- Large intermediate `.mat` files from exploratory runs.

## Scientific basis

The SMDFR plant model builds on the Dual Fluid Reactor concept as analysed by
X. Wang (*Analysis and Evaluation of the Dual Fluid Reactor Concept*, PhD
thesis, Technische Universität München) and on the control scheme of C. Liu,
R. Luo and R. Macián-Juan (*A New Uncertainty-Based Control Scheme of the SMDFR
and Its Optimization*, Energies 14(20), 2021, doi:10.3390/en14206708). The
grey-box modelling, IS-MCDET integration, DSM-penalised MPC formulation and all
code in this repository are original contributions of the thesis.

## Citing

This software is archived on Zenodo and citable by DOI:

> Miqueles Leiva, L. A. (2026). *Grey-Box Digital Twin for Risk Monitoring and
> Predictive Control of a Small Modular Dual Fluid Reactor* (v1.0.0) [Software].
> Zenodo. https://doi.org/10.5281/zenodo.22407139

```bibtex
@software{miqueles2026dtriskcontrol,
  author    = {Miqueles Leiva, Leonardo Andr\'{e}s},
  title     = {Grey-Box Digital Twin for Risk Monitoring and Predictive
               Control of a Small Modular Dual Fluid Reactor},
  version   = {1.0.0},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22407139},
  url       = {https://doi.org/10.5281/zenodo.22407139}
}
```

The DOI above is the *concept* DOI: it always resolves to the most recent
archived version. Zenodo also mints a version-specific DOI for each release if
you need to pin one exactly. Machine-readable metadata is in `CITATION.cff`.

## Licence

MIT — see `LICENSE`. Released open-access in line with the intellectual-property
and open-science requirements applying to EU Horizon-funded research.

## Acknowledgements

Developed at the Laboratory of Analysis of Systems for the Assessment of
Reliability, Risk and Resilience (LASAR), Politecnico di Milano, under the
supervision of Prof. Enrico Zio and Prof. Francesco Di Maio.

Thanks to Giovanni Floreale, Stefano Marchetti and Ibrahim Ahmed for helpful
discussions during this work.
