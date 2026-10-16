# Native-module migration roadmap and geometry-core contract

**Status:** Draft engineering plan; no runtime switch or scientific-admission change.
**Coordination:** Issue #5; task inventory on Issue #135.
**Baseline checked for this branch:** main 0dcb86e57baaeaddd772130bbe96c7195f381f81.
**Current diagnostic production-source anchor:** 7f47669a1c0edfafef9fc5f27b2db80840d3f035.
**Vendor comparison pin:** RonBibb/bb-palatini-unified-r0, d6052d605673ce9d82cdc99b0df79855fa2ea215.

## 1. Decision and source-of-truth rule

Build the geometry, evolution, and constraint machinery as first-party modules in this repository. Independent equations and derivations define expected behavior. The pinned vendor implementation is a differential-verification aid, not the mathematical authority and not the long-term runtime owner.

Migrate by bounded layers. Keep the existing vendor-backed implementation intact while each native layer is derived, implemented, and tested in isolation. Do not redirect production call sites merely because the native and vendor arrays look similar or a short smoke test passes. Any runtime switch requires the gates below and an explicit owner/reviewer decision.

This roadmap concerns software implementation fidelity and numerical-operator verification. It does not reopen or change the project's existing scientific-admission decision: scientific_admission remains false.

## 2. Repository-specific implementation inventory

### Vendor boundary

- .gitmodules pins vendor/bb-palatini-unified-r0 to https://github.com/RonBibb/bb-palatini-unified-r0.git; docs/PRODUCTION_KERNEL.md records commit d6052d605673ce9d82cdc99b0df79855fa2ea215.
- engine/v55_pirk_adapter.py is the direct adapter boundary. Its vendor_modules() dynamically imports r0_operators, vacuum_bssn, and moving_puncture_pirk.
- geometry_stage_terms() wires moving-puncture explicit RHS, primary L2/L3 RHS blocks, and Lambda L2/L3 blocks. primary_l3_with_matter() and lambda_l3_with_matter() add local V5.5 matter projections; those local source terms must remain local during the native port.
- geometry_constraints() currently delegates to vacuum.constraints(). This is an explicit native-constraint target.
- engine/production_kernel.py owns ProductionState and synchronizes geometry, scalar fields, conservative matter, and observables. It is the orchestration seam, not the first module to rewrite.
- engine/production_contract.py defines runtime capability requirements. Do not advertise a new production capability before the corresponding native implementation and tests exist.

### Existing local contracts and dependencies

- engine/matter_system.py already declares BSSNMetricSlice and MetricDerivativeSet. geometry_metric_derivatives() records local ADM metric-derivative formulas from a,b,X,alpha,beta,Aa,K,Lambda,B. This is useful local convention evidence, but it is not a substitute for independently deriving the native curvature operators.
- The local metric reconstruction used by matter is g_rr = a/X^2 and g_theta_theta = b r^2/X^2. The shift belongs in the ADM spacetime line element; it must not be inserted again into the spatial 3-metric.
- Geometry-dependent consumers include engine/valencia.py, engine/v55_matter.py, engine/scalar_system.py, engine/cmc_gauge.py, engine/stress_energy.py, engine/worldtube.py, engine/invariant_diagnostics.py, and engine/radial_maxwell.py. Preserve their public data contracts while the core is being built; audit each explicitly before changing a consumer.
- Existing evidence to reuse rather than duplicate: Issue #4 and PR #275 (signed momentum/connection propagation tangent), PR #284 (analytic regular-centre Ricci calibration), PR #290 (shift-aware tensor/divergence helpers), and the Issue #292 / PR #296 HLL-to-total-stress reconciliation lane. These are adjacent tests and diagnostics, not proof that a native runtime replacement is complete.

## 3. Target module boundaries

The first native implementation should use pure functions and repository-owned dataclasses. It must not import the vendor package, modify sys.path, use importlib to resolve vendor modules, mutate a supplied state, or silently select a vendor fallback.

Proposed package boundary (names are provisional until the first focused PR):

- engine/native_geometry/fields.py — explicit native metric/state fields and shape/type validation, independent of vendor classes.
- engine/native_geometry/derivatives.py — the declared radial grid derivative contract: cell/face locations, stencil order, parity, centre treatment, and outer-boundary closure.
- engine/native_geometry/connection.py — metric-derived connection coefficients and regular-centre limits.
- engine/native_geometry/ricci.py — independently derived spatial Ricci components and scalar contractions, with signed component/profile outputs.
- engine/native_geometry/constraints.py — Hamiltonian, momentum, connection, determinant and trace-free witnesses assembled from native geometry and declared sources.
- engine/native_geometry/evolution.py — native explicit/L2/L3 and gauge stage RHS blocks only after the geometry primitives and signs are independently reviewed.
- engine/native_geometry/pirk.py — native sequential PIRK2 stage composition only after stage-level RHS checks pass.
- engine/v55_pirk_adapter.py remains the legacy comparison adapter during migration. No production call-site changes are part of this contract task.

The package can start with fewer files if that keeps each reviewed change small; the mathematical ownership boundaries must remain explicit. Native types may be converted to the existing consumer-facing shape in one compatibility adapter, rather than making native modules depend on vendor object types.

## 4. Conventions that must be written down before coding

The derivation note and tests must explicitly state:

1. Coordinates and units: radial coordinate versus areal radius, whether radius arrays denote cell centres or faces, units/normalization, angular factors, and where factors of 8 pi enter constraints or matter projections.
2. Metric variables: a,b,X, lapse alpha, radial shift beta, K, Aa, Lambda and gauge B; spatial metric components; determinant convention; trace-free and extrinsic-curvature sign conventions.
3. Index placement and shift: covariant, contravariant and mixed components must not be interchanged implicitly. Spatial metric operators are derived from the spatial 3-metric; shift enters spacetime and evolution formulas as specified by ADM kinematics.
4. Centre behavior: parity for every evolved/derived field, the even/odd extension used by each derivative, cell-centred half-cell treatment, and analytic limits for every expression with apparent 1/r or 1/r^2 factors. No arbitrary epsilon radius, floor, clamp, or centre projection is allowed to hide a singular expression.
5. Boundaries and stages: inner and outer closure, grid refinement convention, time-level used by each RHS term, and sequential PIRK2 ordering. These are part of the discrete operator, not implementation details to infer from vendor output.
6. Source ownership: geometry equations/operators live in the native geometry core; V5.5 matter projections/source terms remain local; scalar and conservative-fluid ownership remains with their existing local modules until a separately reviewed interface change.

## 5. Phased implementation roadmap

### Phase 0 — contract and inventory (this PR)

Record the vendor touchpoint map, local dependencies, field conventions, target APIs, verification matrix, and no-go boundaries. No production source or vendor files change.

### Phase 1 — native metric, connection and Ricci primitives

Derive the spherical spatial metric and connection from the declared metric ansatz independently. Derive each Ricci component by hand/symbolic algebra or an independently authored derivation script, then implement pure functions. Include the analytic regular-centre limits rather than evaluating singular-looking generic formulas at r=0. The first module must produce signed arrays and named intermediate terms so a reviewer can inspect each term.

### Phase 2 — independent constraint operators

Build H, radial M, Lambda/connection, determinant and trace-free checks on top of the native geometry output. Keep geometric-side and matter-subtracted total residuals separate. Preserve signed profiles and centre/fixed-annulus summaries; do not compare only maxima or unsigned norms.

### Phase 3 — native evolution and gauge RHS blocks

Port explicit, L2 and L3 geometry/gauge blocks one at a time from independent derivations. Keep matter additions as named local source projections. Test each block on frozen copies of a common input state, including zero/nonzero shift and parity-regular manufactured profiles. Do not change stage ordering until each block has passed in isolation.

### Phase 4 — native sequential integrator and state ownership

Implement the required sequential PIRK2 stage schedule using native state/RHS data structures. Check that each stage consumes the prescribed time-level and that the complete candidate state is finite, correctly shaped, and independent of the input object. Keep scalar and conservative-matter state ownership explicit; do not couple the migration to a new source law.

### Phase 5 — shadow-mode integration

Run native and legacy paths on identical copied states and pinned inputs without changing the accepted production state. Compare individual RHS blocks, stage intermediates, final one-step deltas, and signed constraints. Differences must be interpreted against independent discrete identities and truncation expectations; do not tune the native result to match vendor output.

### Phase 6 — explicit runtime switch, then dependency retirement

Only after all earlier gates pass, add a deliberate native-mode selection and require the existing standard tests plus the bounded production smoke on the exact candidate head. Keep the old path available for rollback during the agreed transition. Retire runtime vendor imports in a separate, reviewed change; preserve the pinned vendor as a test-only comparison oracle until the maintainers decide it no longer adds value. Removing the submodule is a separate decision, not a side effect of the first port.

## 6. Required verification gates

**G0 — provenance and interface.** Main/base SHA, production-source anchor, vendor pin, grid contract and conventions recorded. No ambiguity over which code owns a term.

**G1 — analytic/manufactured geometry.** Flat space in spherical coordinates; flat space represented in a nontrivial regular radial chart; constant/regular metric profiles; and smooth manufactured profiles with known derivatives/curvature. Verify tensor components and contractions separately.

**G2 — centre/parity/boundary.** Even/odd parity fixtures for every relevant field; centre half-cell and outer-boundary tests; no non-finite values or unreviewed 1/r singular evaluations. Report the centre and the fixed physical annulus separately.

**G3 — algebraic identities and negative controls.** Symmetry/trace identities, determinant normalization, connection definition and constraint identities must close at the declared discretization order. Include deliberate sign-flip, omitted-term, and index-placement negative controls so the tests prove they can catch plausible implementation mistakes.

**G4 — convergence and signed error accounting.** Use at least three resolutions for smooth manufactured data; report signed component profiles, pointwise absolute/relative error, L-infinity and L2 measures, and observed refinement rates. Fixed-annulus and centre behavior are separate witnesses. Tolerances must be justified from truncation/roundoff and fixed before seeing candidate results.

**G5 — independent vendor differential report.** On identical frozen metric/state inputs and exact pinned revisions, compare native and vendor components term by term. Keep a discrepancy ledger with sign, location, expected order, and interpretation. Vendor agreement is corroboration, not an acceptance oracle; a mismatch triggers equation/operator reconciliation, not automatic coefficient fitting.

**G6 — stage/integration checks.** Test each stage RHS independently, then one native sequential step and dt/dt-half response on copied, manufactured states. Include immutable input hashes and finite/shape checks. No long campaign is needed to close a local operator gate.

**G7 — integration and release.** Exact-head full standard workflows and the repository's bounded production smoke pass; interface consumers are covered; branch is current with main; independent derivation and implementation review are complete. Any source/runtime diff is inspected. The reviewer must explicitly confirm that admission-related fields and workflow thresholds were not altered.

A failed gate remains a failure with its signed report preserved. Do not relax tolerances after observing a mismatch; first identify its analytic/discrete source and add a negative-control or derivation explanation.

## 7. Scope exclusions and completion rule

This migration does not authorize changes to V5.5 physical equations, source laws, coefficients, gauge/initial-data choices, floors or clamps, damping, outer boundary physics, radiation admissibility/cone rules, campaign settings, archive content, or scientific-admission flags. It does not claim that a successful geometry port proves coupled trajectory stability, total stress-energy conservation, convergence of the full physical model, or scientific admission.

The Phase 0 deliverable is complete only when the owner/reviewer accepts this interface and gate contract. The next implementation PR should be limited to Phase 1: native spherical metric/connection/Ricci primitives plus analytic/manufactured tests, with no production call-site switch.
