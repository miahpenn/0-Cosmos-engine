# V5.5 Driver Audit — Derived Clock and E-fold Integration

**Date:** 2026-10-13 UTC  
**Refresh base:** current `main` at `a0d610d0f6c0de5587b8331bc3703225bbac425e`  
**Original candidate base:** `c01f2467347e8aebd57d220374fd7bc4e93045bc`  
**Candidate branch:** `diag/driver-v55-clock-ledger-2026-10-13`  
**Scope:** numerical bookkeeping only; no physics/source or vendor changes.

## Driver boundary

`engine/coupled_engine.py` is a historical CMC/reference driver and is explicitly marked non-production by the README. The canonical campaign runner in `engine/run_production.py` defaults to `V55ProductionKernel`; `engine/campaign.py::prepare_campaign` applies the capability contract before a campaign begins. The current evidence therefore does not justify wiring the historical driver into the production path.

## Findings on the production V5.5 step

1. Geometry, scalar and conservative-matter variables use stage-paired / trapezoidal evolution, but the derived accumulators used the endpoint-only rectangle rules `tau += dt * tau_rate_end` and `e_folds += dt * H_eff_end`. This made their integration order lower than the step's endpoint pairing.
2. Diagnostics were generated before `candidate.e_folds` was advanced. The history record therefore persisted the incoming e-fold value on the new endpoint unless explicitly overwritten.
3. On the first accepted step, event detection substituted endpoint H for incoming H when history was empty, which could hide a genuine first-step sign crossing.

## Isolated candidate design

- Keep `coupled_engine.py` historical; do not add a second production owner or duplicate COSMOS RHS.
- Compute the left-end rates from the saved accepted slice before the step-local CMC lapse refresh: volume-weighted `H_eff_start` and the recorded central-lapse rate `tau_rate_start`. This makes the accumulated proper time reproducible from accepted-endpoint rates. Record the refreshed start-stage lapse rate and its difference from the accepted rate separately, exposing the existing CMC re-solve discrepancy instead of silently substituting one for the other.
- Trapezoid-integrate H and proper-time rates against the final accepted state's diagnostics.
- Detect crossings using incoming and outgoing H even when the history list is empty.
- Persist the updated `e_folds` alongside the accepted endpoint's `t` and `tau`.
- Preserve the original three-argument `append_efolds` behavior for existing callers; its explicit start-rate option enables endpoint trapezoidal integration.

## Verification plan and limits

**Latest refresh note:** this branch is refreshed from current `main=a0d610d0f6c0de5587b8331bc3703225bbac425e`. Earlier workflow results validate their own pinned heads only; fresh exact-head gates are required for this tree.

**Refresh note:** the workflow results and numerical witness quoted above were produced on the original candidate base `c01f2467347e8aebd57d220374fd7bc4e93045bc`. They are retained as historical evidence only; they do not validate this refreshed tree. Fresh exact-head software gates are required before promotion.

Focused tests exercise a linear-rate quadrature, backward-compatible helper use, a synthetic one-step crossing, both accumulators, and ledger endpoint consistency. Exact-head repository CI and its bounded N=40 one-step smoke are the acceptance checks for this software candidate.

This update changes reported derived coordinates and event labels, not evolved field equations. It does not repair or claim closure of centre regularity, constraints, the radiation predictor boundary, fixed-areal-radius Misner–Sharp balance, scalar–DM source-pair conservation, or the full production campaign. It is diagnostic/software work and must remain isolated until the exact-head evidence is reviewed.
