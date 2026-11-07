# Issue #4 coupled analytic matter-direction fixture — derivation-first specification

**Date:** 2026-11-07  
**Proposed class:** `coupled_analytic_matter_direction`  
**Status:** specification only; no fixture construction, archive, or predictor evaluation authorized by this document  
**Parent controls:** immutable `manufactured_regular_family`, immutable `equation_consistent_matter_direction` (frozen geometry), and immutable accepted-state archive  
**Scope:** research-only, instantaneous, zero production steps; no production initializer/RHS/stepper/vendor call; all admission and closure flags remain false.

## 1. Purpose and distinction between controls

The existing `manufactured_regular_family` is an arbitrary-family control. Its scalar-momentum and fluid-velocity variations are not required to satisfy the matter evolution equations; its residual must not be forced to zero.

The existing `equation_consistent_matter_direction` is a second matter-direction control, but it holds geometry, lapse, and connection fixed. It is therefore a frozen-geometry matter control, not the coupled control specified here.

This proposed third class is the missing **coupled geometric-and-matter direction**: scalar field/momentum rates satisfy the pinned scalar equations; fluid conservative rates satisfy the pinned Valencia system and their continuum continuity/acceleration identities; geometry, primary variables, connection, and CMC lapse directions are derived from their pinned instantaneous equations; and all six sector stress projections and their directional derivatives are evaluated with the same locked projection definitions as the observer.

It is not automatically a complete Einstein-matter solution. It is an analytic instantaneous control whose equation coverage and remaining constraints must be stated explicitly in the manifest.

## 2. Frozen provenance and immutable inputs

Use the provenance contract already established by the Issue #4 observer and coupled-control derivations:

- diagnostic source anchor: `4aec88d0e3ff028e0a65a8ac36b30948625df792`
- helper baseline: `7f47669a1c0edfafef9fc5f27b2db80840d3f035`
- vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- tangent helper blob: `3f00a4e0a2063b86d2e007c5a825f9edfff61a4c`
- total-stress helper blob: `c10f9e8b9e654d5195fd55c493c3fe06ba5c7082`

The manifest must pin the exact reviewed blob IDs for the independent constructor, scalar/Valencia/geometry transcription, six-sector projection derivative, CMC tangent implementation, fixture loader, observer, evaluator, and tests. Branch names are descriptive only. Any mismatch fails closed before construction or evaluation.

The three existing classes, archives, reports, and comparison outputs are read-only inputs. Do not regenerate, relabel, or modify them.

## 3. Base state and complete direction inventory

Define one regular base state (U_0) and one analytic instantaneous direction (dot U). Construct the symmetric members (U_0\pm\epsilon\dot U) independently for the locked values (epsilon=2\times10^{-5}) and (10^{-5}). No member may come from a production initializer, production RHS, kernel/vendor step, trajectory, or time integrator.

The state inventory must explicitly declare, field by field, both the base value and its direction:

- Geometry: (a,b,X,A_a,K,\Lambda); the dependent trace-free component (A_b=-A_a/2); metric determinant/volume factors; radial and angular metric coefficients.
- Gauge: (\alpha_{\rm in}), (\alpha_{\rm out}), and their distinct directions until an independently demonstrated identity permits equality. CMC gauge has (\beta=B=0); their zero directions are explicit gauge conditions, not omitted data.
- Scalars: (S,PS,D,PD,\phi,\Pi), including each field and conjugate-momentum direction.
- Fluids: dark matter, baryons, and radiation conservative arrays ((D,\tau,S_r)) (or exact canonical names at the pinned interface), their primitive recovery, and their directions.
- Six-sector projections: for each of dark matter, baryons, radiation, (S,D,\phi), retain ((\rho,p_r,p_t,j)) and ((\dot\rho,\dot p_r,\dot p_t,\dot j)).
- Derived source totals: raw six-sector sums and the separately time-density-accounted channels. Never replace the raw channels with accounted values.

The base must be regular, positive-metric, and parity-declared. Metric/scalar-even fields and the odd connection difference (C^r=\Lambda-\Lambda_{\rm metric}) receive explicit parity and centre regularity declarations. The base must satisfy the pinned discrete CMC lapse equation and the declared determinant/trace-free/connection identities before any directional witness is evaluated.

## 4. Scalar direction: pinned field equations

For each scalar (f\in\{S,D,\phi\}), use the pinned curved-slice spatial operator
[
\Delta_g f=\frac{X^2}{a}\left[D_2f-D_1f\left(\frac{D_1a}{2a}-\frac{D_1b}{b}+\frac{D_1X}{X}-\frac2r\right)\right].
]

The independent instantaneous transcription is
[
\dot f=\alpha p_f,
\qquad
\dot p_f=\alpha\Delta_g f+\alpha Kp_f
+\frac{X^2}{a}(D_1\alpha)(D_1f)-\alpha V_f'(f)+S_f.
]

Use the pinned potentials and signs, not a refit:
- (V_S=+S^2/2), hence (V'_S=S);
- (V_D=-D^2/2), hence (V'_D=-D);
- (V_\phi=V_0(e^{-\phi}-0.10)), with (V_0=8.242522415500654\times10^{-5}), hence (V'_\phi=-V_0e^{-\phi});
- the pinned scalar exchange source is (S_\phi=\beta_{DM}\rho_{DM}) in the conjugate-momentum equation, with (\beta_{DM}=-0.04). The corresponding dark-matter exchange covector has the opposite sector-exchange sign under the locked convention.

The direction of every term must be differentiated analytically, including metric, lapse, field, momentum, spatial-operator, and exchange-density responses. The scalar direction is not a freely assigned profile. Cross-check the independently coded rates against the pinned source equations and against symmetric plus/minus scalar reconstructions at both epsilons.

### 4.1 Detailed scalar/Valencia derivation note

The pinned source transcription and analytic directional formulas are now recorded in [`docs/ISSUE4_SCALAR_VALENCIA_DIRECTION_DERIVATION_2026-11-07.md`](ISSUE4_SCALAR_VALENCIA_DIRECTION_DERIVATION_2026-11-07.md), using the exact scalar, matter-system, matter-RHS, Valencia, matter-bundle, and production-kernel blob IDs from helper baseline `7f47669a1c0edfafef9fc5f27b2db80840d3f035`.

That note separates:
- the scalar field-rate map and its full analytic direction, including metric, lapse, spatial-operator, potential, and DM-density responses;
- the actual stored-array semi-discrete Valencia RHS (face interpolation, minmod reconstruction, HLL flux, geometric sources, and DM exchange);
- dust and radiation conservative-to-primitive directional recovery;
- continuum continuity/acceleration witnesses from the finite-volume discrete RHS contract.

The note also records a critical branch requirement: minmod, reconstruction guards, HLL extrema, and radiation recovery must remain on the same differentiable branch across each centered (+\epsilon/-\epsilon) pair. Branch ties or switches are inconclusive, not grounds to smooth or loosen the test. The actual `matter_system.py` stored-array arithmetic is authoritative for the RHS; similarly named alternate helper variables must not be substituted.

**Disposition:** this completes a source-pinned derivation/transcription milestone only. No independent implementation has yet been run, and no scalar/fluid direction gate is claimed as passed.

## 5. Fluid direction: continuity and acceleration, with pinned discrete realization

For each fluid, the continuum direction must satisfy its declared matter equations on the same base geometry:

[
\nabla_\mu(\rho_0u^\mu)=0,
\]
[
(\rho+p)u^\mu\nabla_\mu u_\nu
+ (\delta_\nu{}^\mu+u_\nu u^\mu)\nabla_\mu p
= (\delta_\nu{}^\mu+u_\nu u^\mu)Q_\mu.
]

Use the pinned equation of state and exchange assignment:
- baryons: pressureless, no dark-matter exchange source;
- dark matter: pressureless, with (Q_{DM,r}=+\beta_{DM}\rho_{DM}\partial_r\phi) in the locked covector convention and the scalar exchange opposite;
- radiation: (p=\rho/3), zero rest-current channel, and its pressure-gradient acceleration retained.

For the actual archived directional arrays, independently transcribe the pinned Valencia conservative RHS: primitive recovery, face metric interpolation, reconstruction, HLL fluxes, radial geometric sources, and sector-specific exchange/rest-current sources. Convert the conservative rates to primitive continuity and acceleration witnesses and retain both forms. A continuum Euler identity alone does not substitute for the pinned finite-volume directional contract.

The matter direction must include the metric/lapse response in conservative-to-primitive recovery. Do not set density derivatives to zero while assigning nonzero velocity derivatives unless the continuity equation independently gives that result.

## 6. Geometry and primary direction: pinned CMC equations

On the zero-shift CMC slice, with (A_b=-A_a/2), the explicit metric rates are
[
\dot a=-2\alpha aA_a,\qquad
\dot b=\alpha bA_a,\qquad
\dot X=\frac{\alpha XK}{3}.
]

The primary directions are the pinned L2 curvature/lapse-Hessian plus L3 matter-coupled blocks:
[
\dot A_a=-(D_rD_r\alpha-\tfrac13\Delta\alpha)
+\alpha(R^r{}_r-\tfrac13R)+\alpha KA_a
-\frac{16\pi}{3}\alpha(p_r-p_t),
]
[
\dot K=-\Delta\alpha
+\alpha(\tfrac32A_a^2+\tfrac13K^2)
+4\pi\alpha(\rho+p_r+2p_t).
]

The connection direction is the pinned Lambda L2 geometry block plus L3 connection and matter momentum source. Retain the exact parity-aware discrete derivatives and source normalization from the pinned implementation; do not replace the L2 block with a continuum simplification. The current independent transcription in `tools/issue4_cmc_lapse_tangent.py` is a source map and implementation starting point, not by itself proof that the new fixture has passed all coupled gates.

For the connection, define the pinned geometry-only momentum bracket
\[
M_{geom}=D_1A_a-\frac23D_1K+6A_a\chi'
 +(A_a-A_b)\left(\frac2r+\frac{D_1b}{b}\right),
\qquad \chi'=-\frac{D_1X}{2X}.
\]
With the locked \(\lambda_m=2\), retain the source-separated blocks
\[
\dot\Lambda_{L2}=
-\frac{2(A_aD_1\alpha+\alpha D_1A_a)}{a}
-\frac{4\alpha(A_a-A_b)}{rb}
+\frac{\lambda_m\alpha}{a}M_{geom},
\]
\[
\dot\Lambda_{L3}=
2\alpha A_a\Lambda
-\frac{8\pi\lambda_m\alpha}{a}j,
\qquad
\dot\Lambda=\dot\Lambda_{L2}+\dot\Lambda_{L3}.
\]
The matter current \(j\) is the raw sum of the six pinned sector projections; its direction is retained separately in the archive and is not substituted into the base \(\dot\Lambda\) expression.

All curvature, lapse-Hessian, connection, and matter operands must be reconstructed from the declared base arrays and locked operators. The geometry direction is not a vacuum or test-field direction: it uses the same six-sector source projections and their derived directions.

## 7. CMC lapse: preserve stage-aware input/output semantics

The stage-aware contract distinguishes the incoming lapse from the newly solved output lapse:
[
\lambda=P_w[K_{RHS}(\alpha_{in})],\qquad
L(a,b,X,\mathcal Q)\alpha_{out}=-\lambda,
]
where (w_i=r_i^2\sqrt{a_i}b_iX_i^{-3}\Delta r), (W=\sum_iw_i), and
[
\mathcal Q=\tfrac32A_a^2+\tfrac13K^2+4\pi(\rho+p_r+2p_t).
]

The incoming physical K-RHS direction must be derived from the pinned fieldwise stencil:
[
K_{RHS}(\alpha_{in})=
-\frac{X^2}{a}\left[D_2\alpha_{in}
-D_1\alpha_{in}\left(\frac{D_1a}{2a}-\frac{D_1b}{b}
+\frac{D_1X}{X}-\frac2r\right)\right]
+\alpha_{in}\mathcal Q.
]

Differentiate the complete expression, the proper-volume weights and normalization, and the output matrix/RHS. The output direction is obtained from
[
A\dot\alpha_{out}=\dot b-\dot A\alpha_{out}.
]
Keep (\dot\alpha_{in}) and (\dot\alpha_{out}) separate. The accepted-state archive's missing preceding incoming-lapse direction must not be replaced by zero. For this new analytic control, provide the incoming-lapse direction from an explicitly declared preceding analytic refresh/stage, or construct and independently validate a self-consistent fixed-point control. The latter is permitted only if both the base fixed-point identity and its two-epsilon directional cross-check pass; a small internal linear-system residual alone is insufficient.

### 7.1 Synthetic fixed-point lapse-control evidence (2026-11-07)

A refined research-only stage-history cross-check was supplied from the user's sandbox. Its two controls are explicitly separated:

- **History A — declared stage history:** smooth positive Gaussian incoming lapse and an independently declared nonzero incoming-lapse direction. The map residual is intentionally nonzero because this incoming lapse is not the fixed-point root.
- **History B — self-consistent fixed point:** the same base state and declared state rates, with the incoming lapse set to the affine fixed-point root and its direction set to the independently solved implicit fixed-point tangent.

The refined script reports the following results (user-reported; the sandbox script and run were not independently executable from this repository session):

| Gate | Check | Reported result |
|---|---|---:|
| 1 | Linear CMC solve residual at History A incoming lapse | 1.8e-15 |
| 1 | Map residual at History A incoming lapse | 3.8e-2 (expected nonzero; not a fixed point) |
| 2 | History A coupled state + incoming-lapse pure-map centered-difference error, two epsilons | 3.0e-10 |
| 2b | History A incoming-lapse-only pure-map centered-difference error, two epsilons | 9.0e-11 |
| 3 | History B fixed-point/map residual | 3.0e-14 |
| 4a | History B analytic implicit identity residual | 1.6e-17 |
| 4b | History B coupled state + implicit-lapse pure-map centered-difference error, two epsilons | 1.3e-10 / 6.5e-10 |

The reported separation is also explicit: the infinity-norm difference between declared and fixed-point incoming lapse is 2.6e-2; the difference between declared incoming-lapse direction and implicit fixed-point tangent is 1.0e-3. These are distinct controls and are not required to agree.

**Disposition:** this resolves the previously open *synthetic CMC fixed-point directional-control* gate, conditional on the reported run being reproducible against the pinned helper. The independent pure-map centered-difference check at the fixed point is the key added evidence; the implicit residual alone would not suffice. The earlier float64 finite-difference-of-solver-roots discrepancy is not used as the validation target here.

**Scope limit:** the declared geometric and matter rates in this synthetic history are not yet derived from the scalar field and Valencia fluid evolution equations. The script does not construct the physical coupled-matter fixture, accepted-state incoming-lapse history, six-sector archive, or conservation closure. It does not authorize production changes or admission-flag promotion. The script reports a nonzero centre gap for the Gaussian incoming lapse (about 1.2e-3); the base's parity/centre regularity requirements in this specification still need an explicit satisfactory construction for the eventual coupled fixture.

## 8. Six-sector projection and time-density accounting

Use the same pinned canonical projection helper for all six sectors. Independently differentiate its metric and primitive dependencies. At minimum, the volume response is
[
\frac{\dot{\sqrt\gamma}}{\sqrt\gamma}
=\frac12\frac{\dot a}{a}+\frac{\dot b}{b}-3\frac{\dot X}{X},
\qquad
\frac{\dot g^{rr}}{g^{rr}}
=2\frac{\dot X}{X}-\frac{\dot a}{a}.
]

For scalar sectors, retain the kinetic, gradient, potential, and current derivatives. For fluid sectors, retain conservative-to-primitive, Lorentz-factor, metric, pressure, and current responses, including the implicit radiation primitive relation. Require per-sector analytic directions to match independently reconstructed plus/minus projections at both epsilons.

Report separately:
1. raw six-sector reconstruction of total ((\rho,p_r,p_t,j));
2. time-density component accounting and its remainder;
3. pointwise divergence;
4. finite-volume divergence.

A raw reconstruction pass does not imply the time-density accounting pass. No remainder is absorbed into another sector or used to fit a coefficient.

## 9. Archive, zero-step, and evaluator contract

Use a new immutable archive class and dedicated fail-closed loader; do not add this class by weakening or overloading the existing arbitrary/frozen-flat schemas. Preserve the established archive members `report.json`, `tangent_arrays.npz`, and `SHA256SUMS.txt`, with exact array key/shape/dtype/hash records and canonical state digests for base and each plus/minus member.

The manifest must record:
- class and schema identity; analytic family definition and every fixed/varied field;
- all source/helper/vendor/reviewed-code pins;
- grid, resolution, parity, regularity, profile definitions, and epsilon inventory;
- scalar, fluid, geometry, primary, incoming/output lapse, and six-sector base/directional arrays;
- pointwise and FV divergence arrays and raw/accounted sector ledgers;
- per-equation residuals and plus/minus cross-checks;
- `production_steps_run=0`, `state_advanced=false`, `production_initializer_called=false`, `production_rhs_called=false`, `vendor_step_called=false`;
- `scientific_admission=false`, `acceptance_oracle_satisfied=false`, `physical_exchange_law_admitted=false`, and `conservation_closure_admitted=false`.

The evaluator must consume this class through a narrow adapter into the existing six-sector/zero-step machinery. Predictor terms, signs, masks, centre cells 0/1, cell-2 diagnostic, fixed-radius/annular summaries, epsilon handling, and zero-step tripwires remain unchanged. The input state digest must match before and after evaluation.

## 10. Ordered derivation and implementation gates

No fixture implementation begins until the derivation package is complete and reviewed against the pinned source blobs. Then proceed in this order:

1. **Equation transcription:** scalar and Valencia source equations and directional contracts are now documented in the linked derivation note; complete the independent code transcription and retain exact source/blob/operator pins, then document the CMC metric/primary/Lambda equations and stage-aware lapse map.
2. **Base construction:** define a regular analytic base and independently verify metric positivity, parity, centre regularity, determinant/trace-free/connection identities, six-sector base projections, and discrete CMC base residual.
3. **Direction construction:** derive every scalar, fluid, geometry, primary, connection, and incoming/output lapse direction. No free arrays, residual fitting, coefficient selection, or unexplained zero directions.
4. **Equation residuals:** evaluate scalar, continuity, acceleration, Valencia, geometric/primary, connection, and both CMC directional residuals on the base and direction.
5. **Two-epsilon cross-checks:** construct plus/minus members independently at both locked epsilons; compare every primitive, sector projection, geometry rate, K-RHS projection, lapse output, and six-sector directional projection against centered differences.
6. **Archive and loader:** implement the separate immutable schema and fail-closed provenance/digest/array validation.
7. **Observer:** only after all prior gates pass, run the unchanged pointwise and FV observer with zero-step tripwire armed.
8. **Comparison:** compare against the three immutable controls, keeping raw/accounted channels and pointwise/FV channels distinct.

A failure at any gate stops the sequence and is recorded without loosening tolerances or substituting a different operator.

## 11. Explicit non-goals

- No changes to production equations, initializer, solver, vendor, kernel, predictor, stencils, tolerances, fixtures, or existing archives.
- No use of the manufactured arbitrary-family residual as a target to force to zero.
- No use of the frozen-geometry matter control as a coupled tangent.
- No residual-based coefficient fitting or operator selection.
- No inferred missing lapse direction, silent fixed-point assumption, or unlabelled zero direction.
- No scientific-admission, acceptance-oracle, physical-exchange, or conservation-closure promotion based on this control.

**Readiness status:** derivation-first specification only. The synthetic self-consistent CMC fixed-point directional control is reported as passing, subject to reproduction against the pinned helper. The coupled fixture remains unbuilt. Equation-consistent scalar and Valencia fluid directions, the physical coupled geometry/connection direction, six-sector directional projections, regularity-complete base, immutable archive, and unchanged zero-step observer remain explicit gates.
