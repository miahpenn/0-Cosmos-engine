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


## Result admitted — 12 October 2026

The generic-versus-radiation-specific one-step comparison completed on run [38317797540](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317797540): success. Artifact [11777312270](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317797540/artifacts/11777312270), ZIP SHA-256 `e34c99e6b575222f0f8ce5c134f9ac2bdb9a958ece0dc235f593a8a63e38e090`, JSON SHA-256 `566b51387ca3272fa76716b19d5e1128450ef02c6ae784ee851f2759697a785b`. Vendor pin, compile, 1 helper test, 4 guard tests and 27 full-suite tests passed; all four arms passed.

Result:
- N=80 generic rest max `1.0085731565667868e-11`; radiation-specific rest max `0.0`.
- N=160 generic rest max `5.012616800265632e-12`; radiation-specific rest max `0.0`.
- Full radiation E/S profiles, all derived cone profiles, all geometry fields and H/M/C_Lambda witnesses are bit-identical between arms on both grids.

The isolated source-repair candidate is at [`fix/radiation-rest-current-channel-2026-10-12`](https://github.com/miahpenn/0-Cosmos-engine/tree/fix/radiation-rest-current-channel-2026-10-12). It dispatches radiation HLL and outer-face flux through the existing radiation-specific helpers and adds an actual N=80/N=160 rest-invariant regression. Its dedicated repair gate [38317961532](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317961532) was queued/running at the time of this handoff update; await its result before calling the code candidate admitted.

This identifies a narrow representation defect and supports the minimum dispatch repair, but its E/S and geometry effects are exactly zero in the captured A/B. It does not explain the late-time cone crossing and is not a physical-law change or full-machine admission.


## Follow-on fix candidate

The diagnostic is complete. Its exact A/B showed generic dispatch creates rest drift, while the existing radiation-specific helpers remove the drift with E/S, geometry and H/M/C_Lambda bit-identical. A separate implementation candidate has passed its dedicated gate (3 focused tests, 4 radiation guard tests, 29 full-suite tests): draft stacked [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11), based on [PR #10](https://github.com/miahpenn/0-Cosmos-engine/pull/10). This branch remains diagnostic-only; the repair candidate branch is separate. No merge to main occurred.
