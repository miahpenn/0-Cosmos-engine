# 0-Cosmos Engine

GEAR-derived cyclic-cosmos numerical engine and diagnostic archive.

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

## Campaign progression

The True-CMC stage-aware strong-field gate completed successfully at N=40,60,80 through t=24.

The matched-spacing enlarged-domain turnaround gate then completed successfully at N=80,120,160 through t=40, with one turnaround event and one handoff in each resolution and no re-expansion before t=40.

The next gate is the existing long post-turnaround campaign `[run-true-cmc-post]`.


## Launch note
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
