# Strict campaign admission + radiation retry integration candidate

**Date:** 12 October 2026  
**Disposition:** draft engineering integration; not merged; not scientific admission.

## Objective

Reconcile the tested strict campaign runner and the typed radiation-cone retry mechanism without changing physical equations, coefficients, guards, boundary conditions, or the pinned vendor kernel.

## Source provenance

- Strict runner/diagnostic base: `fix/strict-campaign-admission-2026-10-12`, head `d5a5b97eadcd5c565722e61262a22195baf3c1d5`.
- Campaign preflight candidate: `fix/campaign-config-preflight-2026-10-12`, head `c9041799ece39cf64ae64ddf1cb45ae94a4b665c`.
- Typed retry helper and radiation state-contract candidate: `fix/radiation-adaptive-step-on-cone-2026-10-12`, head `651dc124484e7c056b784033b69bedba6da0b5a5`.
- PR #3's strict constraint/centre diagnostic additions in `engine/production_kernel.py` are preserved; the diverged adaptive branch's kernel was not copied wholesale.

## Integration invariants

1. `prepare_campaign()` validates campaign parameters, the resolution set, timestep progress, checkpoint spacing, and production capability contract before evolution.
2. The retry helper catches only the typed `RadiationAdmissibilityError`; unrelated geometry, gauge, and runtime errors are not retried. It halves only the trial timestep and fails closed when coordinate-time progress is no longer representable.
3. Radiation admissibility is checked at initialization and the actual initial, explicit predictor, primary predictor, and completed radiation states before source recovery can consume inadmissible values.
4. On a successful retry, whole-state finiteness, positive geometry, diagnostic finiteness, endpoint canonicalization, and checkpoint logic still apply to the accepted state.
5. Per-resolution `radiation_step_retries.json` and `ledgers.json` preserve retry attempts. `campaign_summary.json` records exact resolution/horizon admission and remains fail-closed on numerical failure or an incomplete horizon.
6. No clipping, projection, retry of non-radiation errors, cone relaxation, new source, fitted coefficient, or physics promotion is introduced.

## Required regression gates

- Retry one injected typed cone rejection, then take real production steps to the full short horizon; require an accepted retry ledger and `admission.passed=true`.
- Inject persistent typed cone rejection, let the actual retry helper reach its time-progress limit, and require a persisted failure ledger plus `admission.passed=false`.
- Run the complete repository suite and bounded N=40 production smoke in the integration-specific Actions workflow.
- Preserve the stage-level radiation-admissibility and zero-rest-current regression suites from the adaptive candidate.
- After independent integration review and green CI, a bounded N=40/60/80 to t=12 run may be considered. Do not launch a long-horizon campaign until the separate centre-series/Hamiltonian/momentum/connection, radiation N=240/320, Misner–Sharp balance, and convergence gates are reviewed.

## Interpretation boundary

Even if the integration workflow passes, it establishes software-level retry and strict campaign-admission compatibility only. It does not demonstrate scientific convergence, conservation, physical turnaround/re-expansion, a bounce, a completed cycle, or production readiness.
