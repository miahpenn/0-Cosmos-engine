# Native-module migration roadmap and geometry-core contract

**Status (checked 2026-10-17 CDT):** Active native-module migration contract; native CMC explicit-metric API and its CMC-stage observer are merged. No production runtime switch or scientific-admission change is merged.
**Coordination:** Issue #5; task inventory on Issue #135.
**Original roadmap baseline (before the native-module promotion PRs):** main `0dcb86e57baaeaddd772130bbe96c7195f381f81`.
**Current main after CMC API and observer increments:** `1a34f16333f3409cdaa41a44211cc7e322418bc2` (PR #358 API; PR #360 CMC observer).
**Pinned diagnostic production-source anchor (historical baseline retained intentionally):** `7f47669a1c0edfafef9fc5f27b2db80840d3f035`. This is not the current main SHA; keep the provenance distinction explicit.
**Vendor comparison pin:** RonBibb/bb-palatini-unified-r0, `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
**Production-migration gate:** Draft PR #359 (`9184df9b2928465804ab26578fdbc9c76e4bb75b`) passes the three core functional workflows but its source-pinned trace-free Ricci/proper-radial workflow fails closed on protected production-source changes. Keep it draft until an explicit trusted owner/reviewer authorizes a separate production/source-anchor transition. Do not weaken or bypass the existing pin/allowlist.

## 1. Decision and source-of-truth rule

Build the geometry, evolution, and constraint machinery as first-party modules in this repository. Independent equations and derivations define expected behavior. The pinned vendor implementation is a differential-verification aid, not the mathematical authority and not the long-term runtime owner.

Migrate by bounded layers. Keep the existing vendor-backed implementation intact while each native layer is derived, implemented, and tested in isolation. Do not redirect production call sites merely because the native and vendor arrays look similar or a short smoke test passes. Any runtime switch requires the gates below and an explicit owner/reviewer decision.

This roadmap concerns software implementation fidelity and numerical-operator verification. It does not reopen or change the project's existing scientific-admission decision: scientific_admission remains false.

## 2. Repository-specific implementation inventory

### Vendor boundary

- .gitmodules pins vendor/bb-palatini-unified-r0 to https://github.com/RonBibb/bb-palatini-unified-r0.git; docs/PRODUCTION_KERNEL.md records commit d6052d605673ce9d82cdc99b0df79855fa2ea215.
- engine/v55_pirk_adapter.py is the direct adapter boundary. Its vendor_modules() dynamically imports r0_operators, vacuum_bssn, and moving_puncture_pirk.
- `geometry_stage_terms()` remains the adapter boundary for the pinned explicit RHS and primary/Lambda blocks on current `main`; its local `primary_l3_with_matter()` and `lambda_l3_with_matter()` projections keep V5.5 matter-source ownership local.
- PR #358 adds `native_cmc_explicit_metric_rhs(grid, state)` in `engine/native_geometry/evolution.py`. It returns exactly `a/b/X`, requires the CMC convention `beta=0`, does not require the B-driver field, and shares its metric formula core with `native_explicit_rhs()`.
- PR #360 updates the existing CMC-stage observer to call that CMC-only API while comparing directly against the pinned vendor explicit RHS on copied stage inputs; exact array equality and accepted/copied-state immutability remain enforced.
- `geometry_constraints()` still delegates to `vacuum.constraints()`; this is a separate native target.
- `engine/production_kernel.py` owns `ProductionState` and orchestrates geometry, scalar fields, conservative matter, and observables. Its production call sites remain protected by the source-pinned diagnostic until a separate explicitly authorized migration path is established.
- `engine/production_contract.py` defines runtime capability requirements. Do not advertise a new production capability before the corresponding native implementation and tests exist.

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
- `engine/native_geometry/evolution.py` — native explicit metric RHS functions exist. The CMC-safe public `native_cmc_explicit_metric_rhs()` API is backed by exact-array, input-immutability, and nonzero-shift rejection tests. Primary/Lambda blocks, constraints, gauge evolution, and full PIRK2 stage composition remain separate candidates.
- `engine/native_geometry/pirk.py` — native sequential PIRK2 stage composition remains a separate candidate; the existing native moving-puncture step is not a drop-in CMC integrator.
- `engine/v55_pirk_adapter.py` remains the comparison/compatibility adapter for pinned primary/Lambda blocks and local sources. The proposed PR #359 override path is not merged; keep default direct adapter behavior and production stage behavior unchanged until the explicit owner/reviewer and source-provenance gates are satisfied.

The package can start with fewer files if that keeps each reviewed change small; the mathematical ownership boundaries must remain explicit. Native types may be converted to the existing consumer-facing shape in one compatibility adapter, rather than making native modules depend on vendor object types.

## 4. Conventions that must be written down before coding

The derivation note and tests must explicitly state:

1. Coordinates and units: radial coordinate versus areal radius, whether radius arrays denote cell centres or faces, units/normalization, angular factors, and where factors of 8 pi enter constraints or matter projections.
2. Metric variables: a,b,X, lapse alpha, radial shift beta, K, Aa, Lambda and gauge B; spatial metric components; determinant convention; trace-free and extrinsic-curvature sign conventions.
3. Index placement and shift: covariant, contravariant and mixed components must not be interchanged implicitly. Spatial metric operators are derived from the spatial 3-metric; shift enters spacetime and evolution formulas as specified by ADM kinematics.
4. Centre behavior: parity for every evolved/derived field, the even/odd extension used by each derivative, cell-centred half-cell treatment, and analytic limits for every expression with apparent 1/r or 1/r^2 factors. No arbitrary epsilon radius, floor, clamp, or centre projection is allowed to hide a singular expression.
5. Boundaries and stages: inner and outer closure, grid refinement convention, time-level used by each RHS term, and sequential PIRK2 ordering. These are part of the discrete operator, not implementation details to infer from vendor output.
6. Source ownership: geometry equations/operators live in the native geometry core; V5.5 matter projections/source terms remain local; scalar and conservative-fluid ownership remains with their existing local modules until a separately reviewed interface change.

## 4A. Concrete candidate conventions and API (Phase-0 review addendum)

**Review state (reconciled 2026-10-16 CDT):** This table is a candidate contract for the diagnostic interface, with equations sourced from the independently validated research library and documented local conventions. The absence of an APPROVED review on merged PR #298 is not itself a blocker: do not repeat scientific review merely to re-validate established derivations. Promotion gates are source/provenance traceability plus applicable implementation checks (analytic/manufactured tests, negative controls, exact-head CI, integration smoke, and scope checks). Reopen convention/scientific review only for a named concrete discrepancy, changed convention, failed check, provenance mismatch, or scope violation. Genuinely provisional source-law, index, and normalization choices remain explicitly provisional and cannot be treated as scientifically admitted by merging diagnostics. The code under `diagnostics/` remains a test prototype, not the canonical production API.

Pinned inputs for the table: Phase-1 PR #299 head `1665cc3c1c061e22589f1c7797f5ecc4ac22c798`; local convention source `engine/matter_system.py` blob `e52a117690253ddad7b24d116145a65474a351d2`; production diagnostic source anchor `7f47669a1c0edfafef9fc5f27b2db80840d3f035`; vendor comparator `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

### Grid, coordinate, and stencil contract

| Item | Concrete Phase-1 candidate convention | Evidence / boundary |
|---|---|---|
| Radial coordinate | `r` is the coordinate/isotropic radial coordinate used to sample metric profiles. It is **not** the areal radius. The derived areal radius is `R_areal = r*sqrt(b)/X`. Derivatives in the native grid act with respect to coordinate `r`, unless an expression explicitly converts to proper distance `dl = sqrt(a/X^2) dr`. | Local metric contract in `engine/matter_system.py`; nonlinear-flat-chart test in PR #301. |
| Locations and spacing | Uniform cell centres `r_i=(i+1/2)h`, `i=0..n-1`; faces `r_f=f h`, `f=0..n`; `h=dr=r_max/n`; `r_max` denotes the outer face. Cell profiles have shape `(n,)`; face profiles have shape `(n+1,)`. | `test_shell_geometry_matches_pinned_grid_contract`. |
| Coordinate shell weights | Face area weight `A_f=r_f^2`; per-unit-solid-angle shell weight `V_i=(r_{i+1}^3-r_i^3)/3`. A factor `4*pi` is added only by a caller that explicitly integrates over a full sphere. These are coordinate-grid FV weights, **not** the physical proper-volume element for a non-flat spatial metric. | Shell-volume sum and pinned differential check in `test_shell_geometry_matches_pinned_grid_contract`; divergence/SBP contract in `test_face_operators_match_pinned_grid_and_preserve_fv_identity`. |
| Parity | Parity parameter `p=+1` means even and `p=-1` means odd. At the origin, reflected cell ghosts are `u[-1]=p*u[0]` and `u[-2]=p*u[1]` under reflected cell indexing. Reject booleans and every parity value except `-1,+1`. Parity is specified per field/operator; never infer it from values. | `test_fourth_order_derivatives_reproduce_regular_polynomials` and `test_native_grid_rejects_invalid_resolution_shapes_and_parity`. |
| Cell derivatives | The compatibility methods `cell_derivative` and `cell_second_derivative` are second-order: parity-aware centre stencil, centred interior, one-sided outer closure. The `*_fourth` cell methods use five-point centred polynomial-moment stencils where available, parity-reflected origin ghosts, and one-sided polynomial-moment closures over the last two cells. First derivative requires at least 5 cells; fourth-order second derivative requires at least 6 cells. | Polynomial reconstruction and outer second-derivative observed-order witness in `tests/test_native_spherical_operators.py`. |
| Face gradient/interpolation | `gradient_faces` is second order in the interior and takes explicit inner/outer boundary gradients. `gradient_faces_fourth` and `interpolate_faces_fourth` have parity-aware inner handling and documented lower-order outer closures (first-difference / constant where currently implemented). Do **not** label these arrays globally fourth order; mask or separately establish boundary order for each use. | Current vendor differential checks exist. Independent smooth-profile observed-order tests for these face operators remain a promotion requirement. |
| Conservative divergence | `divergence_faces(F)_i=(r_{i+1}^2 F_{i+1}-r_i^2 F_i)/V_i`; diffusion is the divergence of a face coefficient times face gradient. Boundary gradients are explicit inputs, not hidden state extrapolation. | `test_face_operators_match_pinned_grid_and_preserve_fv_identity`. |
| Centre and outer closure | The grid is cell-centred; no field is evaluated at a centre cell located at `r=0` (the first centre is `h/2`). A term with apparent `1/r` or `1/r^2` is not regularized by an epsilon, floor, clamp, or projection: derive and test its analytic regular limit or use a parity-consistent expression. Report the first cells separately from a fixed physical annulus. The existing geometry differential fixture uses centre mask `2 <= r <= 6`; that mask must not be used to claim centre/outer closure validation. | Existing flat/nonlinear-chart controls are diagnostic evidence; geometry-core acceptance still needs explicit centre and outer witnesses. |

### Metric, tensor, connection, and units contract

| Item | Concrete local convention / candidate definition | Evidence / review boundary |
|---|---|---|
| Spatial metric | `gamma_rr=a/X^2`; `gamma_theta_theta=b*r^2/X^2`; `gamma_phi_phi=b*r^2*sin(theta)^2/X^2`. Require `a>0`, `b>0`, and `X>0` for the geometry diagnostic domain. The spatial 3-metric contains no separate shift term; shift enters the ADM spacetime metric and evolution kinematics. | `engine/matter_system.py` local reconstruction; native geometry input checks. |
| Determinant convention | The conformal determinant condition is `a*b^2=1`, with residual `a*b^2-1`. Do not confuse this with the physical determinant `det(gamma)=a*b^2*r^4*sin(theta)^2/X^6`. | Flat-chart and determinant controls in PR #302; source-linked convention traceability and passing tests are the default gate, with targeted review reopened only for a concrete discrepancy. |
| Proper radial derivative and Ricci sign | Set `A=gamma_rr=a/X^2`, `R_areal=r*sqrt(b)/X`, and `d/dl=A^(-1/2)d/dr`. With the candidate Ricci sign convention, `Ricci^r_r=-2*(d^2 R_areal/dl^2)/R_areal`; `Ricci^theta_theta=(1-(dR_areal/dl)^2-R_areal*d^2R_areal/dl^2)/R_areal^2`; scalar `Ricci=Ricci^r_r+2 Ricci^theta_theta`. The positive-curvature analytic witness is now included on main through PR #302. Its exact-head lightweight workflow passed on `668dfee2c394d5ddce19132500b12fffff0b666d` ([run 38985267332](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38985267332)); PR #302 is merged at `83bbdb5bb243e0d0abc05bf08d8e2b9098693709`. Flat-space zero curvature alone cannot determine this sign. | Independent warped-product derivation in PR #301; vendor comparison is supplementary only. |
| Extrinsic curvature and sign | The current local ADM kinematic convention is `partial_t gamma_ij=-2*alpha*K_ij + Lie_beta(gamma)_ij`. The local tensor reconstruction uses `K^r_r=K/3+Aa` and `K^theta_theta=K/3-Aa/2`; thus the trace-free angular mixed component is `Ab=-Aa/2` and `K=K^r_r+2*K^theta_theta`. `Aa` is the mixed radial trace-free component in this contract, not an implicitly lowered `A_rr`. | Derived from the explicit `Krr`, `Kthth`, and metric-time-derivative definitions in `engine/matter_system.py`; preserve the stated index placement and use the sign/parity tests as implementation gates, reopening convention review only if a concrete mismatch appears. |
| Connection variable | The candidate `Lambda` state profile represents the **contravariant radial** reference-metric contracted conformal connection `tilde_Lambda^r`, not a lowered `Lambda_r`. Its metric contraction is `Lambda_metric = a'/(2*a^2) - b'/(a*b) + 2*(1/b - 1/a)/r`; the connection residual is `Lambda-Lambda_metric`. Away from that constraint, keep metric-derived Ricci and BSSN connection-form Ricci as separate named outputs; do not enforce a discrete trace identity by construction. | Formula and wrong-sign offset control are present in PR #302, now merged; preserve the contravariant-radial index placement from the documented derivation and reopen review only on a concrete mismatch. |
| Hamiltonian diagnostic and source inputs | Candidate geometry-only expression is `H_geom = Ricci - Aa^2 - 2*Ab^2 + (2/3)*K^2`. If an explicit `rho` array is supplied, the current diagnostic candidate is `H=H_geom-2*kappa*rho`. Sources are explicit inputs; this API does not derive a source law or authorize production matter coupling. | Source-toggle sign test in PR #301 verifies the implementation contract only, not the physical provenance of `rho`. |
| Momentum diagnostic and source inputs | Candidate geometry residual is `M_geom = d_r(Aa) - (2/3)*d_r(K) + 2*(Aa-Ab)*(d_r R_areal)/R_areal`. The diagnostic currently returns `M=M_geom-kappa*S_r` for an explicitly supplied radial-momentum profile `S_r`. Treat its sign, lowered/raised index and normalization as **provisional until source-law/index provenance is established from the validated archive**; do not propagate it into a production source law from software equivalence alone. | API/negative controls are useful but do not settle the covariant-to-ADM momentum-source mapping. |
| `8*pi*G` / density units | Existing local matter code declares `KAPPA=8*pi`, i.e. the local normalized lane uses the `G=1` coefficient convention. The diagnostic default `kappa=8*pi` mirrors that local convention. This does **not** establish a conversion for archive/COSMOS density inputs: the cross-lane density/Hubble normalization question remains Issue #210 and must stay separate. | Current local source declaration; no new normalization law is authorized. |

### Callable prototype API and validation policy

These signatures describe the existing diagnostic candidates for review, not a promise that the `diagnostics/` modules become the production API.

| Callable | Input/output contract | Status / gap |
|---|---|---|
| `NativeSphericalCellGrid(n: int, r_max: float)` | Frozen grid descriptor; integer `n>=4`, finite positive `r_max`. Properties: `dr` scalar; `centers` shape `(n,)`; `faces`, `face_areas` shape `(n+1,)`; `volumes` shape `(n,)`. | Concrete Phase-1 diagnostic API. |
| `cell_derivative(values, *, parity)`, `cell_second_derivative(values, *, parity)`, `cell_derivative_fourth(values, *, parity)`, `cell_second_derivative_fourth(values, *, parity)` | Real numeric cell profile shape `(n,)` to new float profile shape `(n,)`; parity is `-1,+1`. The fourth-order methods enforce their minimum cell counts above. | PR #303 now rejects non-finite operator profiles and includes NaN/Inf negative controls; the exact-head regression workflow passed (run 38985509358). |
| `gradient_faces(values, *, inner_flux_gradient=0, outer_flux_gradient=0)`, `gradient_faces_fourth(values, *, parity, inner_flux_gradient=0, outer_flux_gradient=0)`, `interpolate_faces_fourth(values, *, parity)` | Real cell profile shape `(n,)` to face profile shape `(n+1,)`; parity is required only by the parity-aware variants; boundary-gradient values are explicit scalar inputs. | The face-array shape contract and finite boundary-scalar validation are tested (PR #303); independent smooth-profile observed-order tests for the fourth-order variants remain a concrete gap. |
| `divergence_faces(flux)`, `integration_by_parts_residual(values, flux)` | `flux` shape `(n+1,)`; `values` shape `(n,)`. Return shape `(n,)` float divergence or a finite scalar residual, respectively. | FV/SBP identity test exists; non-finite flux/profile values are now rejected by the shared validator (PR #303). |
| `diffusion(values, coefficient_faces=1, *, inner_flux_gradient=0, outer_flux_gradient=0)` | `values` shape `(n,)`; coefficient is a real scalar or shape `(n+1,)`; return shape `(n,)`. | PR #303 adds non-finite coefficient, profile, and boundary-gradient negative controls; exact-head regression passed (run 38985509358). |
| `metric_connection_lambda(grid, state)`, `native_bssn_geometry_terms(grid, state)` | Diagnostic geometry state fields are 1-D real finite profiles of shape `(n,)`. Geometry requires positive `a,b,X,alpha`; return named float profiles of shape `(n,)`, including metric/BSSN Ricci witnesses, connection residual, areal radius, radial lapse Hessian and lapse Laplacian. Inputs are read-only. | These remain diagnostic functions, not approved production APIs. The positive-curvature sign/refinement control is now on main via PR #302 and passed its exact-head checks. Direct finite-output checks for extreme-but-finite inputs remain a separate implementation gap. |
| `native_bssn_constraints(grid, state, *, rho=None, radial_momentum_source=None, kappa=8*pi)` | Requires `a,b,Aa,K` plus geometry fields. Each source may be omitted, scalar, or shape `(n,)`; source arrays are explicit and not derived. `kappa` must be finite and positive. Return named profiles `H,M,M_geometry,connection,determinant,R`, each shape `(n,)`. | PR #383 wires canonical local `rho/j` projections and checks the local consumer deltas; PR #384 checks rejection of a non-finite H witness from extreme-but-finite `Aa`. The momentum source index/sign/normalization remain provisional; these tests do not admit a cross-lane source law. |
| All native functions | Reject wrong shape, non-real numeric data, invalid parity, impossible resolution, and invalid positive metric parameters; return new arrays and never mutate caller state. Non-finite profiles/scalars must be rejected before promotion rather than silently allowed to propagate. | PR #303 implements finite-value checks for operator profiles/fluxes, boundary gradients, and diffusion coefficients; geometry helpers validate finite input profiles. PR #384 adds an extreme-finite-output rejection test for `native_bssn_constraints`; extreme-but-finite output checks for the native geometry witnesses remain an explicit gap. |

### Convention-to-test traceability and visible gaps

| Convention / failure mode | Existing or required discriminating control | Gate state |
|---|---|---|
| Coordinate `r` confused with areal radius `R_areal` | PR #301 nonlinear-coordinate flat metric with analytic lapse Hessian/Laplacian; a coordinate/areal-radius substitution should break its independent analytic result. | Existing; use the documented derivation and analytic witness as the default evidence, reopening review only on a concrete mismatch. |
| Incorrect shell volumes / omitted `r^2` face area | Shell geometry and summation-by-parts/FV identity tests. A deliberately unweighted divergence is the negative control. | Positive identity exists; add an explicit mutant control if a concrete discrimination gap remains, not merely to trigger another generic review. |
| Wrong origin parity or ghost sign | Even/odd regular-polynomial exact derivative tests; negative variant flips parity mapping at the first two cells. | Existing polynomial fixtures; confirm the parity-flip negative control is exercised before production-core promotion. |
| Wrong outer derivative order | Smooth `sin(r)` second-derivative witness over N=32/64/128 checks the last two cells for observed order `>3.8`; boundary cells remain explicit rather than hidden inside an interior norm. | Existing for second derivative; corresponding first-derivative and face-operator observed-order witnesses remain a gap. |
| Wrong metric determinant, Ricci sign or angular/radial component | Nonlinear-chart flat-space control and separate component checks; **include the positive-curvature warped-product witness** (`R_areal(l)=L*sin(l/L)`), which is now included on main by PR #302 and passed the exact-head lightweight check ([run 38985267332](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38985267332)), with three-resolution refinement; flat-space zero curvature alone cannot determine this sign. | Flat/nonlinear-chart controls and the positive-curvature sign/refinement witness pass on main-based PR #302. This is diagnostic evidence, not full-trajectory or production-runtime validation. |
| Wrong connection-correction sign/index placement | Nonzero evolved-`Lambda` offset compared to independently derived scalar/radial differences, plus a wrong-sign rejection control. | Present in PR #302 diagnostics, now merged with exact-head checks passing; targeted review is needed only if implementation evidence conflicts with the documented derivation. |
| Wrong source sign or hidden source law | Turn explicit `rho,S_r` inputs on/off and compare signed deltas against the declared diagnostic algebra; separately review the covariant-to-ADM source convention. | Software algebra test exists; the physical source-law and covariant-to-ADM mapping stay explicitly provisional until traceable archive evidence resolves the specific mapping. |
| Wrong API shape, non-mutation or NaN propagation | Wrong-shape/parity/resolution rejection tests and frozen-state input snapshots; add NaN/Inf profiles, non-finite face boundary values, and non-finite coefficient controls. | Shape/non-mutation coverage exists; PR #303 now exercises NaN/Inf rejection for operator inputs. Direct finite-output checks for extreme-but-finite inputs remain open. |

The table deliberately exposes remaining holes instead of implying that existing tests already cover them. Record each concrete convention with its source/provenance and discriminating test. If source support or an implementation check is absent, mark that specific convention unresolved and block only the dependent API. Green CI alone does not establish source-law acceptance, outer-face accuracy, or full production validity; however, do not require a blanket Phase-0 approval review merely to repeat already validated derivations.


## 5. Phased implementation roadmap

### Phase 0 — contract and inventory (this PR)

Record the vendor touchpoint map, local dependencies, field conventions, target APIs, verification matrix, and no-go boundaries. No production source or vendor files change.

### Phase 1 — native metric, connection and Ricci primitives

Derive the spherical spatial metric and connection from the declared metric ansatz independently. Derive each Ricci component by hand/symbolic algebra or an independently authored derivation script, then implement pure functions. Include the analytic regular-centre limits rather than evaluating singular-looking generic formulas at r=0. The first module must produce signed arrays and named intermediate terms so each term remains directly inspectable.

### Phase 2 — independent constraint operators

Build H, radial M, Lambda/connection, determinant and trace-free checks on top of the native geometry output. Keep geometric-side and matter-subtracted total residuals separate. Preserve signed profiles and centre/fixed-annulus summaries; do not compare only maxima or unsigned norms.

### Phase 3 — native evolution and gauge RHS blocks

The explicit metric `a/b/X` RHS implementation and CMC-specific interface are now present in `engine/native_geometry/evolution.py`; PRs #353/#356/#360 provide copied CMC-stage parity evidence and an exact-array regression contract. This is native-component readiness, not production caller promotion. Continue with bounded implicit/L2/L3 and gauge/source ownership tests only where they add coverage beyond existing evidence. Keep matter additions as named local source projections. Test blocks on frozen, same-input states including parity-regular profiles; do not alter CMC stage ordering or integrate moving-puncture alpha/beta RHS into the CMC path.

### Phase 4 — native sequential integrator and state ownership

Implement the required sequential PIRK2 stage schedule using native state/RHS data structures. Check that each stage consumes the prescribed time-level and that the complete candidate state is finite, correctly shaped, and independent of the input object. Keep scalar and conservative-matter state ownership explicit; do not couple the migration to a new source law.

### Phase 5 — shadow-mode integration

The CMC explicit-metric shadow observer is merged in #353 and uses the CMC-only API after #360. It compares `a/b/X` on source-reconstructed stage0, explicit-predictor, and primary-predictor slices against the direct pinned-vendor explicit RHS; the observed arrays were bitwise equal at the established N=80 fixture, and #356 makes native/adaptor equality a regression assertion. Reuse those records instead of rebuilding the same stage matrix. Any further shadow test should cover a named missing stage/input or interface defect and preserve accepted-state fingerprints. This evidence does not replace the protected source-pin gate for a production caller diff.

### Phase 6 — explicit runtime switch, then dependency retirement

A first explicit-metric production wiring proposal exists as draft PR #359, head `9184df9b2928465804ab26578fdbc9c76e4bb75b`, based before PR #358/#360. Its lightweight/one-step, strict regression/JSON, and integrated retry/full-machine workflows pass, but the trace-free Ricci/proper-radial source-pinned workflow fails closed because production kernel and adapter source blobs differ from the protected baseline and the allowlist currently accepts only additive native-core candidates. This is an intentional provenance gate, not a numerical mismatch. Do not merge the proposal, update the source pin, or weaken the workflow without the explicit trusted owner/reviewer decision required by the roadmap.

If that separate authorization is granted, start a new runtime candidate from the current main SHA, use `native_cmc_explicit_metric_rhs()` (three fields only), preserve the CMC lapse solve, `beta/B` resets, adapter-owned primary/Lambda and local matter-source blocks, and boundary ordering, and run every applicable exact-head gate. Keep the previous adapter available as an independent test comparator. Retire vendor runtime imports only in a later, separately reviewed change.

## 6. Required verification gates

**G0 — provenance and interface.** Main/base SHA, production-source anchor, vendor pin, grid contract and conventions recorded. No ambiguity over which code owns a term.

**G1 — analytic/manufactured geometry.** Flat space in spherical coordinates; flat space represented in a nontrivial regular radial chart; constant/regular metric profiles; and smooth manufactured profiles with known derivatives/curvature. Verify tensor components and contractions separately.

**G2 — centre/parity/boundary.** Even/odd parity fixtures for every relevant field; centre half-cell and outer-boundary tests; no non-finite values or unreviewed 1/r singular evaluations. Report the centre and the fixed physical annulus separately.

**G3 — algebraic identities and negative controls.** Symmetry/trace identities, determinant normalization, connection definition and constraint identities must close at the declared discretization order. Include deliberate sign-flip, omitted-term, and index-placement negative controls so the tests prove they can catch plausible implementation mistakes.

**G4 — convergence and signed error accounting.** Use at least three resolutions for smooth manufactured data; report signed component profiles, pointwise absolute/relative error, L-infinity and L2 measures, and observed refinement rates. Fixed-annulus and centre behavior are separate witnesses. Tolerances must be justified from truncation/roundoff and fixed before seeing candidate results.

**G5 — independent vendor differential report.** On identical frozen metric/state inputs and exact pinned revisions, compare native and vendor components term by term. Keep a discrepancy ledger with sign, location, expected order, and interpretation. Vendor agreement is corroboration, not an acceptance oracle; a mismatch triggers equation/operator reconciliation, not automatic coefficient fitting.

**G6 — stage/integration checks.** Test each stage RHS independently, then one native sequential step and dt/dt-half response on copied, manufactured states. Include immutable input hashes and finite/shape checks. No long campaign is needed to close a local operator gate.

**G7 — integration and release.** Exact-head full standard workflows and the repository's bounded production smoke pass; interface consumers are covered; branch is current with main; source/provenance and implementation evidence are traceable. Inspect any source/runtime diff and confirm admission-related fields and workflow thresholds remain unchanged. Do not require blanket independent scientific review for established derivations; reopen targeted review only for a named concrete discrepancy, changed convention, failed check, provenance mismatch, or scope violation.

A failed gate remains a failure with its signed report preserved. Do not relax tolerances after observing a mismatch; first identify its analytic/discrete source and add a negative-control or derivation explanation.

## 7. Scope exclusions and completion rule

This migration does not authorize changes to V5.5 physical equations, source laws, coefficients, gauge/initial-data choices, floors or clamps, damping, outer boundary physics, radiation admissibility/cone rules, campaign settings, archive content, or scientific-admission flags. It does not claim that a successful geometry port proves coupled trajectory stability, total stress-energy conservation, convergence of the full physical model, or scientific admission.

The native geometry core has progressed beyond the original Phase-0 plan: #358 adds the CMC-safe explicit-metric API, and #360 makes the existing copied CMC-stage observer use that API. Preserve the original scientific archive and vendor pin as provenance/comparison inputs; do not duplicate already-passed same-input matrices. The next production step is not authorized yet: PR #359 remains draft because the protected source-pinned diagnostic rejects production kernel/adapter changes. Resolve that blocker only through an explicit trusted owner/reviewer decision and a separately approved source-provenance transition; do not bypass the guard. Until then, prefer source-safe additive native-core work, exact same-input software checks, and documentation updates. `scientific_admission=false`; no production admission is implied.
