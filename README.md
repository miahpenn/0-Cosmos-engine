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

## Campaign trigger

The True-CMC stage-aware strong-field gate is being launched from the experimental branch using the repository's existing `[run-true-cmc-gates]` push trigger.
