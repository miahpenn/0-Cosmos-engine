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

## 4A. Concrete candidate conventions and API (Phase-0 review addendum)

**Review state:** This table turns the current Phase-1 diagnostic prototype and existing local engine conventions into a concrete candidate interface for independent review. It is not acceptance by itself. PR #298 must remain draft until an owner/independent reviewer explicitly accepts or corrects the conventions below. The quarantined code under `diagnostics/` is a test prototype, not a canonical `engine/native_geometry/` API yet.

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
| Determinant convention | The conformal determinant condition is `a*b^2=1`, with residual `a*b^2-1`. Do not confuse this with the physical determinant `det(gamma)=a*b^2*r^4*sin(theta)^2/X^6`. | Flat-chart and determinant controls in PR #301; independent review still required. |
| Proper radial derivative and Ricci sign | Set `A=gamma_rr=a/X^2`, `R_areal=r*sqrt(b)/X`, and `d/dl=A^(-1/2)d/dr`. With the candidate Ricci sign convention, `Ricci^r_r=-2*(d^2 R_areal/dl^2)/R_areal`; `Ricci^theta_theta=(1-(dR_areal/dl)^2-R_areal*d^2R_areal/dl^2)/R_areal^2`; scalar `Ricci=Ricci^r_r+2 Ricci^theta_theta`. Positive-curvature analytic sign control is not yet present and is required before promotion; flat-space zero curvature alone cannot determine this sign. | Independent warped-product derivation in PR #301; vendor comparison is supplementary only. |
| Extrinsic curvature and sign | The current local ADM kinematic convention is `partial_t gamma_ij=-2*alpha*K_ij + Lie_beta(gamma)_ij`. The local tensor reconstruction uses `K^r_r=K/3+Aa` and `K^theta_theta=K/3-Aa/2`; thus the trace-free angular mixed component is `Ab=-Aa/2` and `K=K^r_r+2*K^theta_theta`. `Aa` is the mixed radial trace-free component in this contract, not an implicitly lowered `A_rr`. | Derived from the explicit `Krr`, `Kthth`, and metric-time-derivative definitions in `engine/matter_system.py`; independent sign/index review remains a G0 gate. |
| Connection variable | The candidate `Lambda` state profile represents the **contravariant radial** reference-metric contracted conformal connection `tilde_Lambda^r`, not a lowered `Lambda_r`. Its metric contraction is `Lambda_metric = a'/(2*a^2) - b'/(a*b) + 2*(1/b - 1/a)/r`; the connection residual is `Lambda-Lambda_metric`. Away from that constraint, keep metric-derived Ricci and BSSN connection-form Ricci as separate named outputs; do not enforce a discrete trace identity by construction. | Formula and wrong-sign offset control in PR #301; exact notation/index placement must be explicitly accepted by the independent reviewer. |
| Hamiltonian diagnostic and source inputs | Candidate geometry-only expression is `H_geom = Ricci - Aa^2 - 2*Ab^2 + (2/3)*K^2`. If an explicit `rho` array is supplied, the current diagnostic candidate is `H=H_geom-2*kappa*rho`. Sources are explicit inputs; this API does not derive a source law or authorize production matter coupling. | Source-toggle sign test in PR #301 verifies the implementation contract only, not the physical provenance of `rho`. |
| Momentum diagnostic and source inputs | Candidate geometry residual is `M_geom = d_r(Aa) - (2/3)*d_r(K) + 2*(Aa-Ab)*(d_r R_areal)/R_areal`. The diagnostic currently returns `M=M_geom-kappa*S_r` for an explicitly supplied radial-momentum profile `S_r`. Treat its sign, lowered/raised index and normalization as **provisional until independent source-equation review**; do not propagate this into a production source law from software equivalence alone. | API/negative controls are useful but do not settle the covariant-to-ADM momentum-source mapping. |
| `8*pi*G` / density units | Existing local matter code declares `KAPPA=8*pi`, i.e. the local normalized lane uses the `G=1` coefficient convention. The diagnostic default `kappa=8*pi` mirrors that local convention. This does **not** establish a conversion for archive/COSMOS density inputs: the cross-lane density/Hubble normalization question remains Issue #210 and must stay separate. | Current local source declaration; no new normalization law is authorized. |

### Callable prototype API and validation policy

These signatures describe the existing diagnostic candidates for review, not a promise that the `diagnostics/` modules become the production API.

| Callable | Input/output contract | Status / gap |
|---|---|---|
| `NativeSphericalCellGrid(n: int, r_max: float)` | Frozen grid descriptor; integer `n>=4`, finite positive `r_max`. Properties: `dr` scalar; `centers` shape `(n,)`; `faces`, `face_areas` shape `(n+1,)`; `volumes` shape `(n,)`. | Concrete Phase-1 diagnostic API. |
| `cell_derivative(values, *, parity)`, `cell_second_derivative(values, *, parity)`, `cell_derivative_fourth(values, *, parity)`, `cell_second_derivative_fourth(values, *, parity)` | Real numeric cell profile shape `(n,)` to new float profile shape `(n,)`; parity is `-1,+1`. The fourth-order methods enforce their minimum cell counts above. | Existing prototype checks shape/type/parity, but its shared input-profile helper does not currently reject NaN/Inf field values. Finite-input rejection and tests are required before any promotion. |
| `gradient_faces(values, *, inner_flux_gradient=0, outer_flux_gradient=0)`, `gradient_faces_fourth(values, *, parity, inner_flux_gradient=0, outer_flux_gradient=0)`, `interpolate_faces_fourth(values, *, parity)` | Real cell profile shape `(n,)` to face profile shape `(n+1,)`; parity is required only by the parity-aware variants; boundary-gradient values are explicit scalar inputs. | The face-array shape contract is tested. Finite boundary-scalar validation and independent observed-order tests for the fourth-order variants remain required. |
| `divergence_faces(flux)`, `integration_by_parts_residual(values, flux)` | `flux` shape `(n+1,)`; `values` shape `(n,)`. Return shape `(n,)` float divergence or a finite scalar residual, respectively. | FV/SBP identity test exists. Non-finite input policy still needs to be made explicit in implementation/tests. |
| `diffusion(values, coefficient_faces=1, *, inner_flux_gradient=0, outer_flux_gradient=0)` | `values` shape `(n,)`; coefficient is a real scalar or shape `(n+1,)`; return shape `(n,)`. | Current shape/type checks; add non-finite coefficient/profile/gradient negative controls before promotion. |
| `metric_connection_lambda(grid, state)`, `native_bssn_geometry_terms(grid, state)` | Diagnostic geometry state fields are 1-D real finite profiles of shape `(n,)`. Geometry requires positive `a,b,X,alpha`; return named float profiles of shape `(n,)`, including metric/BSSN Ricci witnesses, connection residual, areal radius, radial lapse Hessian and lapse Laplacian. Inputs are read-only. | These are PR #301 diagnostic functions, not approved production APIs. Positive-curvature sign control and direct finite-output assertions are still needed. |
| `native_bssn_constraints(grid, state, *, rho=None, radial_momentum_source=None, kappa=8*pi)` | Requires `a,b,Aa,K` plus geometry fields. Each source may be omitted, scalar, or shape `(n,)`; source arrays are explicit and not derived. `kappa` must be finite and positive. Return named profiles `H,M,M_geometry,connection,determinant,R`, each shape `(n,)`. | Input state is not mutated in the current witness test. Momentum source index/sign/normalization remain provisional; this is not an admitted source law. |
| All native functions | Reject wrong shape, non-real numeric data, invalid parity, impossible resolution, and invalid positive metric parameters; return new arrays and never mutate caller state. Non-finite profiles/scalars must be rejected before promotion rather than silently allowed to propagate. | The last finite-input requirement is **not fully implemented** by the Phase-1 operators at this snapshot; track it as an explicit follow-up, not an assumed property. |

### Convention-to-test traceability and visible gaps

| Convention / failure mode | Existing or required discriminating control | Gate state |
|---|---|---|
| Coordinate `r` confused with areal radius `R_areal` | PR #301 nonlinear-coordinate flat metric with analytic lapse Hessian/Laplacian; a coordinate/areal-radius substitution should break its independent analytic result. | Existing; independent equation review pending. |
| Incorrect shell volumes / omitted `r^2` face area | Shell geometry and summation-by-parts/FV identity tests. A deliberately unweighted divergence is the negative control. | Positive identity exists; explicit mutant control may be added if reviewer requires it. |
| Wrong origin parity or ghost sign | Even/odd regular-polynomial exact derivative tests; negative variant flips parity mapping at the first two cells. | Existing polynomial fixtures; ensure the negative case is exercised by the final independent review. |
| Wrong outer derivative order | Smooth `sin(r)` second-derivative witness over N=32/64/128 checks the last two cells for observed order `>3.8`; boundary cells remain explicit rather than hidden inside an interior norm. | Existing for second derivative; corresponding first-derivative and face-operator observed-order witnesses remain a gap. |
| Wrong metric determinant, Ricci sign or angular/radial component | Nonlinear-chart flat-space control and separate component checks; **add an analytic nonzero-curvature metric** (for example a constant-positive-curvature warped product with `R_areal(l)=L*sin(l/L)`) so a global Ricci-sign flip cannot pass the flat-space null test. | Flat null is existing; positive-curvature sign gate is open. |
| Wrong connection-correction sign/index placement | Nonzero evolved-`Lambda` offset compared to independently derived scalar/radial differences, plus a wrong-sign rejection control. | Present in PR #301 diagnostics; independent derivation review required. |
| Wrong source sign or hidden source law | Turn explicit `rho,S_r` inputs on/off and compare signed deltas against the declared diagnostic algebra; separately review the covariant-to-ADM source convention. | Software algebra test exists; physical source provenance stays open. |
| Wrong API shape, non-mutation or NaN propagation | Wrong-shape/parity/resolution rejection tests and frozen-state input snapshots; add NaN/Inf profiles, non-finite face boundary values, and non-finite coefficient controls. | Shape/non-mutation coverage exists in the relevant prototypes; NaN/Inf rejection is not yet fully implemented/tested. |

The table deliberately exposes remaining holes instead of implying that existing tests already cover them. A Phase-0 acceptance review should either accept the concrete local conventions, correct them with a cited equation/source, or mark them unresolved and block the dependent API. In particular, do not declare the prototype's source terms, outer-face accuracy, positive-curvature Ricci sign, or NaN/Inf behavior accepted merely because current CI is green.


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
