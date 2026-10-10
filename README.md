# 0-Cosmos Engine

GEAR-derived cyclic-cosmos numerical engine and diagnostic archive.


## Current engineering checkpoint — 2026-10-10

**Active diagnostic branch:** `physics/spatial-beta-covariant-source-closure`  
**Verified source baseline:** `9eed1cf95c52f001c2b87afbbe6da37faa56d50e` (isolated; not merged to `main`)  
**Latest regression CI:** [Run #576](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938171) — **127 tests passed in 59.58 s**.  
**Radiation diagnostic:** [Trace run #7](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938182) reproduced and captured the late outer predictor failure; see [the trace report](docs/RADIATION_PREDICTOR_TRACE_2026-10-10.md).

The current branch verifies the covariant scalar/DM source pair, a deterministic homogeneous-COSMOS limit, PIRK stage synchronization, and the production CMC residual at multiple initial configurations and on the actual predictor/final slices of a controlled one-step update. Its CMC proper-volume gauge now rejects non-finite or non-positive spatial metric factors before metric/matter operators are called; no metric clamp is used.

These regression checks are software/numerical consistency gates, **not long-run physical validation**. A separate diagnostic-only N=80 radiation run now reproduces the late outer predictor realizability failure at t=46.98 accepted / t=47.01 attempted; its uploaded artifact and ordered cone-margin budget are documented in the trace report. The result localizes the final crossing to the predictor-metric evaluation but does not establish whether the upstream driver is the geometry RHS, boundary treatment, or coupled stage admissibility. No bounce, turnaround, or completed cycle is established.

**Governance note:** GitHub currently reports `main` and the active experimental branch as unprotected. The connected integration could not change repository protection settings. A repository administrator should review branch protection and required checks.

For exact commit/run chronology, preserved failed attempts, equations and next gates, see [the engineering audit](docs/ENGINEERING_AUDIT_2026-10-10.md). The older campaign notes below are retained as historical records and should not be read as proof that those jobs are running now.

## Experimental status

The active experimental work is on **0star-central-clock** and is not merged into `main`.

The current validation sequence is:

1. repository regression/smoke tests;
2. spherical-center regularity validation;
3. resolution convergence;
4. CFL/domain independence;
5. synchronized D-mode phase collapse;
6. only then, extended turnaround/bounce searches.

The project follows a strict rule: repair numerical defects in the numerical layer before interpreting strong-field behavior physically. No fitted feedback coefficients, damping/clamps, lapse floors, manufactured bounce laws, or silent physical terms are introduced.

## Current numerical repair

The production kernel now reapplies the pinned spherical-center algebraic regularity projection during the PIRK/CMC stages. This is a numerical consistency projection already defined by the reference BSSN kernel; it does not alter the physical equations.

See `docs/center-regularity-repair.md` for the finding and validation gates.

## Experimental branch

`0star-central-clock`

Draft PR #1: archive-derived central proper-time CMC gauge diagnostic.

`main` remains untouched by the experimental repair.

## Campaigns

Manual GitHub Actions workflows under `.github/workflows/` are the canonical expensive-run entry points. Artifacts should be preserved with each campaign and compared before any subsequent physics interpretation.

## Current validation status

The production branch is one shared spherical spacetime/state graph containing the local S/D fields, corrected COSMOS scalar, conservative dark matter, baryons, and radiation. There is **not yet** a reciprocal homogeneous COSMOS evolution/source channel driven by the local production solution.

The completed full space-time D-lock campaign established the current diagnostic picture: strong D produces a moving D-centered spacetime/clock structure with ordered geometric responses, while the completed runs showed no turnaround/re-expansion/handoff event through t=22.5. This is a diagnostic result, not evidence for a bounce or cycle.

The current production validation sequence is therefore:

1. repository regression/smoke tests;
2. spherical-center regularity validation;
3. resolution/CFL/domain checks;
4. D-field source-order and propagating-ridge diagnostics;
5. only then, any reciprocal COSMOS coupling derived from existing physical ledger quantities.

No fitted feedback coefficient, free volume normalization, manufactured bounce law, lapse floor, or physical stop condition is permitted.


## Historical launch notes (retained record)

The entries below preserve prior trigger messages and launch context; they are not a current run-status indicator.

### Earlier launch notes
The full bidirectional production diagnostic is being launched now: D-amplitude sensitivity plus an independent CMC-operator probe. This launch changes no physical equations.


The D-mode causality/source-order diagnostic pair is also being launched to compare the response across initial D amplitudes and source-order observables.


A coordinate-CFL central-clock control is being launched alongside the coupled and D-mode campaigns; it is a numerical control, not a physics modification.


Coupled launch retry: the workflow is now present on the branch; this commit exists solely to trigger the full bidirectional D/CMC diagnostic.


The coupled D-amplitude/CFL discriminator matrix is now active: D=0, 5e-11, 1e-10 crossed with CFL=0.06 and 0.03 at N=160, Rmax=80, target t=50.


The coupled D phase-lock discriminator is now active: fixed D amplitude with phase rotations at -pi/2, 0, pi/4, pi/2 on the full bidirectional machine.


Phase-lock retry uses immutable ScalarFields replacement; no physical equations changed.


The full bidirectional machine is also being tested at Rmax=160, N=320 (fixed dr=0.5) for D=0 versus D=1e-10 through t=80.


[run-coupled-r160-mid] Launch high-domain intermediate D=5e-11 control.


[run-coupled-r160-roff] Launch high-domain radiation-off D=1e-10 control.
