# N=80 matched baseline/candidate decomposition pair

**Date:** 10 October 2026  
**Purpose:** Test whether the N=40 signed component-attribution pattern persists at the next spatial resolution. Diagnostic only; no production physics or defaults are changed.

## Launch configuration

- Marker: `[run-n80-decomposition-pair]`
- Resolution: `N=80`, `r_max=40`
- Amplitude: `0.01`; width: `7`; D amplitude: `1e-10`
- Radiation ON; CFL `0.03`; final coordinate time `3.0`
- Sample targets: `0, 0.75, 1.5, 2.25, 3.0`
- Expected fixed-step count from `dr=0.5` and `dt=0.015`: 200, if the trajectory reaches the target normally
- Tests: run `tests/test_discrete_consistent_initial.py` before either trajectory

The workflow runs baseline and candidate sequentially from the same commit and uploads separate JSON reports. Baseline uses `V55ProductionKernel.initialize(...)`. Candidate uses `discrete_consistent_state(...)` with the existing candidate-associated CMC lapse recomputation. The test distinguishes the two initialization procedures as packages; it does not isolate individual substeps within candidate initialization.

## Evidence motivating the run

At N=40, baseline central H moved from `+1.007883465926e-3` to `+1.370796034671e-3`; candidate moved from approximately roundoff to `-6.539932812662e-5`. The cumulative metric contribution changed from `+2.183573122244e-4` (baseline) to `-2.095316855657e-4` (candidate), while the cumulative matter term was nearly the same in both.

The N=80 pair asks whether this signed metric-sector difference persists with resolution. It does not presuppose that any equation is wrong.

## Admission and interpretation gates

1. Tests must pass before the paired run.
2. Both reports must state mode, resolution, settings, commit, completion status, accepted steps, sample history, and no evolution failure.
3. Confirm per-step all-cell decomposition closure and central cumulative-sum closure from both reports.
4. Compare central and all-cell maximum H at identical samples and signed component sums.
5. Do not infer physical causality from the largest contribution alone. No production admission, equation edits, fitted coefficients, damping, floors, clamps, or invented physical laws follow from this diagnostic alone.
