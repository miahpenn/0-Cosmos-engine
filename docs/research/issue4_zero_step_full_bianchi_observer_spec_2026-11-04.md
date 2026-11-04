# Issue #4 — Zero-step full-Bianchi observer specification

Date: 2026-11-04  
Status: authorized bounded tooling task; diagnostic/research only  
Coordination: Issue #4 (scientific lane), Issue #5 (coordination), PR #427  
Numerical execution: **not authorized by this specification alone**; implement and review the observer first.

## 1. Purpose

Provide one reviewed, zero-production-step observer that applies the already-recorded full matter-coupled Bianchi predictor to two additional state classes:

1. manufactured regular profiles;
2. the frozen production initial state at (t=0).

The observer must use the same predictor, source/vendor pins, conventions, signed accounting, and report schema as the accepted-state comparison already recorded in PR #427. Its purpose is consistent measurement, not new continuum derivation or a closure claim.

The accepted-state (t=0.60) comparison remains the reference class. The three state classes must remain separately identified; their generation histories are not interchangeable.

## 2. Locked provenance contract

The observer must fail closed unless all required pins match the accepted-state comparison:

- production source anchor: `4aec88d0e3ff028e0a65a8ac36b30948625df792`;
- vendor source: `d6052d605673ce9d82cdc99b0df79855fa2ea215`;
- tangent helper blob: `3f00a4e0a2063b86d2e007c5a825f9edfff61a4c`;
- total-stress helper blob: `c10f9e8b9e654d5195fd55c493c3fe06ba5c7082`.

The total-stress helper's historical source-reference commit `7f47669a1c0edfafef9fc5f27b2db80840d3f035` is metadata for its protected helper lineage only. It is **not** an alternate accepted production source pin and must never be used to satisfy the source-anchor gate.

At startup and in the output manifest, record the resolved repository commit, exact source/vendor/helper blob identities, Python and NumPy versions, runner/OS, CPU model/features, BLAS implementation/core dispatch, NumPy SIMD dispatch where available, and all input/output digests. Missing hardware metadata must be reported as unavailable, not inferred.

No fallback, auto-discovery substitution, or "close enough" pin matching is permitted. A mismatch or missing required identity terminates before state evaluation and produces a machine-readable refusal.

## 3. Zero-step and immutability contract

- No call to `ProductionKernel.step`, no full-step or half-step consistency check, no PIRK advance, and no hidden evolution through a helper.
- The observer may evaluate the already-pinned instantaneous RHS components and directional perturbation states only where those operations are part of the accepted-state observer contract. It must not advance the supplied production state.
- Record the canonical input-state digest immediately before evaluation and immediately after evaluation; require exact equality.
- Record `production_steps_run = 0` and an explicit `state_advanced = false`.
- Add a test that makes any production-step call fail immediately, plus a test that detects input-state mutation.
- If state construction for a class would require a production step, refuse that class. Do not approximate or reconstruct a substitute state.

## 4. Predictor contract

Reuse the accepted-state full-Bianchi predictor implementation and operand definitions. Do not fork or independently re-derive a second predictor.

The report must preserve the signed contributions and ownership of:

1. shift transport (zero on the currently specified states; still report the measured shift and the evaluated zero contribution);
2. normal geometric trace pair, including (alpha K M_{r,mathrm{tot}}), (-	frac13(D_ralpha)H_{mathrm{tot}}), and (+	frac16alpha D_rH_{mathrm{tot}});
3. evolved-connection / tracefree spatial-residual response (D_j(alpha W^j{}_r));
4. independently reconstructed total-stress Bianchi source (-8pialphamathcal D_r^{mathrm{tot}}).

Preserve the accepted-state convention split:

- metric-only Ricci for the geometric constraints;
- evolved-(Lambda) Ricci offset retained separately;
- (chi=X^2) as the conformal factor in the Ricci offset;
- (chi_r=-D_1X/(2X)) only in the spherical tensor-divergence geometry;
- (C^r=Lambda-Lambda_{mathrm{metric}}) odd at the origin and differentiated with the odd-parity operator (D_1^{(-)});
- even metric fields use their declared even-parity operators;
- (F_C=alpha W^r{}_r) and the pinned spherical tensor-divergence operator are reused unchanged.

No sign, parity, index-position, density-weight, or operator substitution is allowed to be selected by residual size.

## 5. Total-stress accounting

Independently reconstruct the total stress projections ((ho,p_r,p_t,j)) from the same six sectors and the same pinned conventions used by the accepted-state comparison.

Retain separate channels for:

- pointwise covariant total-stress divergence;
- finite-volume face-flux total-stress divergence;
- expanded six-sector divergence reconstruction.

Keep the per-sector density-rate accounting explicit. The raw sum of sector divergences must not be described as passing if the total-projection density-rate conversion is required. Report both the raw result and the explicitly accounted result, including the unexplained signed remainder.

Do not substitute an individual sector-exchange covector for the total-stress Bianchi source. Do not count internal exchange twice.

## 6. State-class adapters

The observer may add the smallest research-only input adapters needed to expose each state class in the common observer's canonical state interface. These adapters may load or construct frozen diagnostic inputs; they may not alter production equations, the initializer, vendor code, or the state-generation equations.

### A. Manufactured regular profiles

- Use the already-declared manufactured regular profile definition and parameters; record them fully.
- Record resolution, grid, parity, regularity conditions, and complete state digest.
- Do not promote a manufactured profile to a production-generated state.
- If its fields do not supply every required predictor operand under the pinned definitions, refuse with the exact missing operand list.

### B. Frozen production initial state, (t=0)

- Load or invoke only the pinned production initializer to obtain the initial state; do not call a production evolution step.
- Record initializer source/blob identity, all initialization parameters, resolution/grid, and complete state digest.
- Verify that the state time is exactly the declared initial time and that no state advancement occurred.
- If the initializer cannot be invoked without stepping or cannot reproduce the pinned state identity, refuse rather than synthesize a replacement.

The two adapters must not introduce class-specific equations, coefficients, source terms, or corrections.

## 7. Signed output and comparability

Use the accepted-state report's term names, signs, units, normalization, masks, and summary definitions wherever they apply. The output must include:

- each signed predictor operand and each signed predictor term;
- measured semi-discrete momentum tangent, where the accepted-state contract supports it;
- pointwise and finite-volume divergence channels separately;
- raw and density-rate-accounted six-sector reconstructions;
- signed residual profiles, not only absolute maxima;
- cells 0 and 1 reported individually (cell 2 retained as a useful adjacent-cell diagnostic);
- the same fixed physical-radius masks and annuli as the accepted-state report;
- resolution and state-class labels on every result;
- input/output digests and complete provenance manifest.

Do not compare different resolutions as though they were the same grid. Do not omit a term because it is small or unavailable. Missing terms are explicit refusals or named unavailable operands, not zeros.

## 8. Inherited assumptions and limits

The observer inherits, without independently proving:

- the locked projected-Bianchi continuum identity and its sign/index conventions;
- the pinned BSSN mapping of the normal geometric and tracefree connection blocks;
- the declared metric-only versus evolved-connection Ricci split;
- the six-sector stress definitions and canonical projection rules;
- the pointwise and finite-volume discrete divergence definitions.

The report must label these as inherited conventions/assumptions and point to their derivation/source records. The observer tests their application to the state classes; it does not independently prove the continuum derivation, physical validity of each sector model, or conservation closure.

The continuum predictor, semi-discrete tangent, and finite-step/PIRK secant remain distinct. This task evaluates the first against the same instantaneous tangent contract as the accepted-state report; it does not claim PIRK finite-step closure.

## 9. Required tests before any numerical run

1. Exact accepted source/vendor/helper pins accepted.
2. Each single-pin mismatch independently rejected before state evaluation.
3. Historical `7f47669…` source reference cannot satisfy the accepted source-anchor gate.
4. Production-step guard proves zero calls.
5. Input-state digest unchanged.
6. (chi) versus (chi_r) and (C^r) parity conventions match the accepted-state implementation.
7. Predictor term-by-term outputs match the accepted-state observer on its existing archived accepted-state fixture within the already-declared numerical tolerances.
8. Both state adapters either pass schema/operand validation or refuse with precise missing-input reasons.
9. No production, initializer, vendor, dependency, or workflow file is modified.

Do not start the manufactured or (t=0) numerical campaign until these tests pass and the observer implementation has received review.

## 10. Explicit non-goals and authorization boundary

This task does not authorize:

- production evolution or any production step;
- changes to production equations, initializer, CMC centre row, vendor implementation, dependencies, or workflows;
- coefficient fitting, hidden hand-tuning, residual-driven operator selection, or class-specific correction terms;
- changing scientific-admission, conservation-closure, or acceptance-oracle flags;
- declaring Issue #4 physically or numerically closed.

## 11. Completion definition

The tooling task is complete only when the reviewed observer:

1. implements the exact shared contract;
2. passes all fail-closed and zero-step tests;
3. produces reproducible manifests for both state classes, or precise fail-closed refusal records where an input class cannot satisfy the contract;
4. is documented with the inherited assumptions and limitations above.

Only after that review should the two frozen state-class runs be authorized. The three-class comparison and any scientific interpretation come afterward.
