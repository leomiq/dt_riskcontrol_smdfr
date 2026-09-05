# Level-2 C S-function sources

The C sources for the seven Level-2 S-functions the Simulink model
(`rom_MCDET_v2_2022b_control_CMEX.slx`) loads as compiled `.mexw64` binaries,
plus the MathWorks `sfuntmpl_basic.c` template they derive from.

Recovered on 2026-09-05 from
`PhD - Polimi/Grey-Box Modelling/Models for training (v01-07)/Version6/GBM_RMM_CINECA/`.
The Version5 and Version6 copies are byte-identical, so there is no ambiguity
about which revision the shipped binaries correspond to.

They are kept here so the model can be rebuilt on another machine — the
`.mexw64` binaries are Windows x64 only, and without source the experiment
cannot be reproduced on any other platform.

## Rebuilding

```matlab
mex -largeArrayDims msfcn_indirectps_hrm_c.c
mex -largeArrayDims msfcn_indirectps_rom_v2.c
mex -largeArrayDims msfcn_indirectps_rom_v3.c
mex -largeArrayDims msfcn_limintm_hrm_c.c
mex -largeArrayDims msfcn_limintm_rom_v2.c
mex -largeArrayDims msfcn_schedule_hrm_c.c
mex -largeArrayDims msfcn_schedule_rom_v2.c
```

Verified with MATLAB R2024a and the MinGW64 C compiler.

## Equivalence check against the shipped binaries (2026-09-05)

Two identical copies of the experiment folder were prepared, one keeping the
archived `.mexw64` and one with freshly compiled ones, and `pilot_single_run.m`
(the fixed-state LJ run at `t_q = 1180 s`, seed 1) was executed in each:

| quantity | agreement |
|---|---|
| PCT trajectory, max abs. difference | 9.3e-06 K  (relative 7.9e-09) |
| `f_opt` | 1.4e-05 relative |
| first cost evaluation | bit-identical (`J_error = 1.2387350000`) |

The first cost evaluation matching bit-for-bit is the decisive result: it is
the deterministic evaluation before the stochastic search begins, so the
recompiled physics is identical. The residual differences are at
floating-point level and of the same order as the run-to-run variation of the
`parfor` search itself, which is not bit-reproducible between two runs of the
*same* binary.

**Method note.** An earlier comparison appeared to show a 9.5 % gap in
`f_opt`. That was an artefact: `diary` appends, so a second run's log was read
as if it were one run of 902 evaluations instead of two of 451. With fresh
logs the two trees agree as tabulated above. If this check is repeated, delete
`pilot.log` between runs.

## One difference worth noting

`msfcn_indirectps_rom_v2.mexw64` as archived is 261 KB against 40 KB when
rebuilt. It carries debug sections (`.debug_info`, `.debug_line`), i.e. it was
built with `mex -g`; the other six were not. This is a build-flag difference
only — the object code is otherwise equivalent, as the run above confirms.
