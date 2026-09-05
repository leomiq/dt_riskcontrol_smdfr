# Grey-Box Digital Twin for Risk Monitoring and Predictive Control of an SMDFR

Reference implementation accompanying the PhD thesis *A Grey-Box Digital
Twin-based Approach for Condition-informed Risk Monitoring and Control of
New-Generation Nuclear Power Plants* (Politecnico di Milano).

The code covers the full loop studied in the thesis: a reduced-order model of a
Small Modular Dual Fluid Reactor (SMDFR), a grey-box virtual sensor for the peak
cladding temperature (PCT), an Importance-Sampling Monte Carlo Dynamic Event
Tree (IS-MCDET) for dynamic risk assessment, and a model predictive controller
whose objective carries a deterministic-safety-margin (DSM) soft constraint.

## Layout

```
matlab/          DET engine, MPC loop, optimizers (Luus-Jaakola, GA, PID),
                 physics/parameters, verification and sweep scripts
matlab/mex_sources/   C sources for the seven Level-2 S-functions
models/          Simulink models, trained BiLSTM weights, normalisation and
                 fixed-state inputs needed to run the pilot
python/          MC-dropout inference bridge, training pipeline, plotting and
                 KPI-compilation scripts
results/         Arrays and tables that regenerate the thesis figures without
                 re-running the simulations
docs/            Experimental record for the DSM-penalty study
```

## Requirements

- MATLAB R2024a with Simulink, and a C compiler for `mex` (verified with MinGW64)
- Python 3.10+ with `numpy`, `scipy`, `matplotlib`; `torch` additionally for the
  grey-box inference path
- Windows x64 for the pre-compiled `.mexw64` binaries. On any other platform,
  rebuild them from `matlab/mex_sources/` first (see below).

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

Please cite the thesis and this repository. A `CITATION.cff` is provided; a DOI
will be minted on archival release.

## Licence

MIT — see `LICENSE`. Released open-access in line with the intellectual-property
and open-science requirements applying to EU Horizon-funded research.

## Acknowledgements

Developed at the Laboratory of Analysis of Systems for the Assessment of
Reliability, Risk and Resilience (LASAR), Politecnico di Milano, under the
supervision of Prof. Enrico Zio and Prof. Francesco Di Maio.

Thanks to Giovanni Floreale, Stefano Marchetti and Ibrahim Ahmed for helpful
discussions during this work.
