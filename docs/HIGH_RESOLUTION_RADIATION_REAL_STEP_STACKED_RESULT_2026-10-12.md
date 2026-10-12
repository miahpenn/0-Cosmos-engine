# Stacked radiation admissibility repair — high-resolution real-step result

**Date:** 2026-10-12  
**Classification:** diagnostic evidence only; not production admission  
**Status:** `diagnostic_capture_complete`

## Question

On the combined candidate containing (1) the radiation RHS/admissibility stage-order repair and (2) the zero radiation rest-current representation repair, does the production kernel complete one ordinary step at (N=240) and (N=320) while every captured radiation admissibility check passes?

## Reproducibility

- Diagnostic workflow: [run 38320321651](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38320321651)
- Artifact: [high-resolution-radiation-real-step](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38320321651/artifacts/11778526198)
- Artifact ID: `11778526198`
- Artifact ZIP SHA-256 (workflow upload digest): `64a52661aa9fa5bfd2a78475559d9306f4e4f538417e44beb8db42fba01869f1`
- Embedded `report.json` SHA-256, independently rechecked after extraction: `db501090b7e2100ab5bf2a744910b8d8a4646050519e532d572c09871cb87b4d`
- Diagnostic source branch: `diag/high-resolution-radiation-real-step-stacked-2026-10-12`
- Diagnostic source commit: `8ac3d19a80caf0dac5dccafdd64ea4301a3cbea2`
- Candidate base tested: `cdcec8282bff0443f2fe4674b15929822ff38736`
- Geometry/PIRK vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- Run result: successful job, no captured step error.

## Results

| Resolution | Δr | Δt | Step | Final max |S|/E | Minimum recorded cone margin |
|---|---:|---:|---|---:|---:|
| (N=240), (r_{max}=80) | 1/3 | 0.01 | completed | 0.009919261711863676 | (2.4402507110609815\times10^{-10}) |
| (N=320), (r_{max}=80) | 0.25 | 0.0075 | completed | 0.00991880271985033 | (1.3726509266587036\times10^{-10}) |

For each resolution, all initial and final conservative values were finite, all metric cells were valid, no nonfinite ratio cells were observed, and all cells passed the unchanged (E\geq0), |S|\leq E admissibility test. The minimum recorded cone margin occurs at cell 0, so the state passes but the centre remains the limiting location.

There were eight recorded validator calls per resolution: two each for `initial`, `explicit-predictor`, `primary-predictor`, and `final`. Each recorded call reports `passed`. The duplicated stage labels correspond to the pre-consumer and post-consumer checks added to ensure primitive recovery/CMC/boundary consumers cannot run before the applicable actual-state guard.

### Workflow checks

Run 38320321651 passed:
- pinned geometry/PIRK vendor check;
- diagnostic/production-module compile check;
- focused radiation RHS/admissibility regressions: **7 passed**;
- `python -m pytest -q tests`: **33 passed** (the `tests/` subtree, not the whole repository);
- exactly one actual production step at each of (N=240) and (N=320);
- JSON report and checksum preservation.

Separately, the combined PR #10 + PR #11 regression gate [run 38320154231](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38320154231) passed **4 rest-current tests, 7 radiation RHS/admissibility tests, and 90 tests under the repository-wide `python -m pytest -q` command**, plus vendor and compile checks. Keep these as two distinct pieces of evidence; do not describe the one-step diagnostic's `tests/` command as repository-wide.

## What this establishes

The two combined software repairs remove the artificial unit-step RHS-extraction rejection and ensure admissibility is checked before downstream stage consumers as well as after them. The zero-rest-current representation invariant also passes its direct regression. With both repairs stacked together, the actual production kernel completes one step on each tested high resolution without a recorded radiation-cone rejection.

## What it does not establish

- It is a one-step diagnostic, not a multi-step trajectory or campaign convergence result.
- The very small positive centre cone margins warrant continued stage-by-stage monitoring; passing this step is not evidence of robust long-time margins.
- It does not resolve the distinct late-time (N=80) outer-cell radiation predictor cone crossing previously captured near accepted (t\approx46.98), attempted (t\approx47.01). That failure was associated with the explicit metric update consuming the matter cone margin; its underlying geometric cause remains unidentified.
- It does not close centre-series regularity, semi-discrete momentum-constraint growth, Misner–Sharp conservation/current balance, or multi-resolution long-horizon convergence.
- It does not establish or support a bounce, turnaround, cycle, ejection event, or any specific physical interpretation.

## Disposition and next gate

Keep PR #10 and stacked PR #11 draft and unmerged pending independent review. Keep this branch diagnostic-only. Do not loosen the admissibility condition, clamp or project the state, add damping, change a source coefficient, or infer a geometric law from this pass.

The next justified work is a bounded diagnostic of the previously captured late-time (N=80) crossing on the corrected candidate: preserve the exact accepted state and attempted predictor, replay the actual stage order, and resolve the cone-margin budget between matter update and metric update. Compare the matched stage snapshots and source/operator contributions without introducing a repair hypothesis. Only after that attribution should a new long-horizon run be considered. The separate momentum/centre-series and Misner–Sharp admission gates remain required.

**Production status:** no production `main` change or physics promotion was made by this work.
