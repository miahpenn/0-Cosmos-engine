# Radiation rest-current channel audit — diagnostic-only handoff

**Date:** 12 October 2026  
**Scope:** One-step, fixed-initial-state comparison; no production source or physics change.  
**Base:** PR #10 radiation RHS guard-repair candidate `5a5a9657e5d0e45c899addd5851ab2bd8db1cfe0`.  
**Question:** Does the generic all-species HLL path populate the otherwise-zero radiation `rest` slot, and can a diagnostic dispatch to existing radiation-specific helpers change E/S or geometry?

## Source-derived identity

In `engine/matter_system.py`:

- `initialize_radiation()` initializes `ConservedSpecies.rest` to exact zeros.
- `_radiation_conserved()` and `_radiation_flux()` return zero as their first component, while their energy/momentum components match the generic Valencia helpers.
- `_hll_flux()` unconditionally calls `_valencia_flux()` and `_valencia_conserved()`, even for `Species.RADIATION`.
- `evolve_species()` updates the rest slot from the first face-flux component and also fills the outer face with generic `_valencia_flux()`.
- `recover_radiation()` consumes the stored Eulerian energy and momentum and does not use the rest slot; the radiation stress-energy projection is likewise based on recovered primitives.

This suggests a passive rest-current representation inconsistency. It is not yet attributed to the late outer radiation-cone failure. Do not change the production path from static inspection alone.

## Diagnostic design

At N=80 and N=160, run the same default radiation-enabled one-step state twice:
1. existing generic dispatch;
2. diagnostic-process-only dispatch through the existing radiation-specific `_radiation_flux` and `_radiation_conserved` helpers for radiation.

The wrapper restores all monkeypatches after the one-step call. All other source equations, stage ordering, boundaries, cone validators, and vendor remain unchanged. The JSON records full radiation profiles, rest maxima, E/S cone checks, geometry profiles, constraint summaries, candidate-minus-baseline deltas and SHA-256.

Expected isolating signature: the generic arm may produce a nonzero radiation rest slot; the radiation-specific arm should retain exact zero. Energy/momentum and geometry should be identical if this is indeed a passive channel. The runner fails closed if the corrected-rest arm is nonzero or if E/S or geometry changes.

## Admission rules

A green workflow validates the diagnostic and this one-step identity only. It does not justify a production patch unless the intended invariant is confirmed against the fluid/state contract and a regression test is accepted. No clamp, physical source, radiation cone change, physics-law addition, or long evolution is part of this audit. Keep the late-time outer-cell predictor-cone failure as a separate issue.
