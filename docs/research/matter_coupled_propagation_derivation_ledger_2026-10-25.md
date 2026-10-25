# Matter-coupled momentum-propagation derivation ledger — 2026-10-25

**Status:** research-only derivation ledger; no identity admitted.  
**Parent convention record:** PR #421, head `6dfc53c48aaecc18b8389a75b6de7b14340747d1`.  
**Accepted-state production-source anchor recorded by #421:** `4aec88d0e3ff028e0a65a8ac36b30948625df792`.  
**Separate Brown-observer source pin:** the pinned `tools/audit_brown_vacuum_momentum_propagation.py` declares `SOURCE_COMMIT = 7f47669a1c0edfafef9fc5f27b2db80840d3f035`; do not conflate these identifiers. Reconcile the exact artifact-manifest source identity and observer source-root contract before executing any observer.  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`.  
**Accepted states:** existing full-continuous-profile states at t=0.60, N=60,120,240.  
**Scope:** organize the remaining source-backed derivation and identify which measurements already exist. No evolution, no production changes, no new convergence sweep, no scientific admission.

## 1. Decision and immediate work boundary

The next task is **derivation-first**, not another N=60/120/240 observer campaign. PR #420 already supplies accepted-state residual profiles, sector-source mapping, pointwise/FV divergence comparisons, and the Brown source bridge. PR #421 records the conventions that are resolved for the Brown candidate and the conventions that remain open.

This ledger does not declare the candidate identity correct because its numerical residual decreases or because the signed source bridge closes algebraically. It records the evidence needed to derive the full relation independently, then permits only targeted post-processing for gaps that remain.

## 2. Frozen candidate and notation

Use the PR #421 conventions without re-deciding them:

- Primary momentum residual: lower-index radial covector M_r.
- Brown connection residual: contravariant C^r = Lambda - Lambda_metric.
- Brown candidate Hamiltonian-gradient slots:
  `-(D_r alpha / 3) H + (alpha / 6) D_r H`, plus the existing C^r derivative terms.
- H in this Brown candidate is the mapped metric-only Hamiltonian expression; do not silently substitute an evolved-Lambda-dependent form.
- Accepted S-sector spatial operator: the already identified cell-centred, parity-aware operator. Do not replace it with a fluid HLL face-flux form.
- Geometry arm: full continuous-profile accepted states only.
- Existing resolutions: N=60,120,240, reuse verified artifacts.
- Production source and vendor pins: as stated above.

The candidate residual under audit is
`D_tot = dot(M_r)_measured - B_tot - S_Q_candidate`,
where `S_Q_candidate = -8 pi alpha Q_r` remains provisional. The map's Brown operator and candidate source must be evaluated with their signed operands and the exact accepted-state convention.

## 3. Evidence already available — do not rerun to reproduce

### 3.1 Full-profile Brown candidate residual

Existing #420 report values:

| Witness | N=60 | N=120 | N=240 |
|---|---:|---:|---:|
| Global max absolute residual | 8.975e-6 | 9.317e-7 | 4.133e-8 |
| Signed residual at r=2 | +8.975e-6 | +6.122e-7 | +4.030e-8 |
| Signed residual at r=8 | -1.194e-6 | -7.328e-8 | -4.591e-9 |

Reported fixed-radius orders are approximately 3.87, 3.93 at r=2 and 4.03, 4.00 at r=8. These support convergence of this *candidate residual* at those fixed radii. They do not establish the origin-cell continuum limit, the matter source law, or the full identity.

### 3.2 Brown source bridge

Existing #420 arrays define:
- B_tot: Brown construction on matter-inclusive H/M constraints.
- B_geom: same Brown construction using geometric-only constraints.
- S_act: existing geometric matter-injection plus current-tangent composite.
- S_Q_candidate: provisional stress-divergence source.
- Delta_S = S_act - S_Q_candidate.
- D_tot = dot(M_r)_measured - (B_tot + S_Q_candidate).
- D_geom = dot(M_r)_measured - (B_geom + S_act).

The recorded signed bookkeeping identity
`(B_tot - B_geom) - Delta_S = D_geom - D_tot`
closes to approximately 5.08e-21, 7.46e-21, and 1.06e-20 for N=60,120,240. The maximum gap in `(B_tot-B_geom)-Delta_S` is reported as 3.625e-10, 3.009e-11, and 2.015e-12.

Interpretation: this is an algebraic reconciliation of existing arrays. It is not an independent derivation of the continuum matter source, and it must not be used to fit or tune a source term.

### 3.3 Six-sector and divergence evidence

The existing sector observer covers S, D, COSMOS phi, dark matter, baryons, and radiation. Its sector additivity checks are tight, but additivity is not total stress-energy conservation. The pointwise and finite-volume divergences remain different discrete observables: the former differentiates cell-centred covariant flux with pointwise metric-source terms; the latter reconstructs spherical face fluxes and cell-integrates terms before dividing by cell measure.

Keep these representations separate unless a source-backed coordinate/measure/index derivation proves equivalence. In particular, the S-sector centre-cell discrepancy is not a fluid-HLL defect.

## 4. Exact-source code trace and pin discrepancy

Read-only inspection of the files at the recorded accepted-source anchor returned these blob SHAs:

| Path | Blob SHA at inspected ref | Directly established by the source |
|---|---|---|
| `tools/audit_brown_vacuum_momentum_propagation.py` | `d1e61abe2fd9d97760ebf2d7292e7fda9a373b08` | Declares a distinct `SOURCE_COMMIT=7f47669a1c0edfafef9fc5f27b2db80840d3f035` and vendor pin; builds `H_metric` from metric-only Ricci, `Aa), `Ab=-Aa/2), and `K`; explicitly avoids the vendor Hamiltonian to prevent overlap with explicit connection terms. |
| `engine/native_geometry/constraints.py` | `124d1f237e54d082b4d76213f9ebe144e02b02ba` | Defines `H=R-Aa^2-2Ab^2+2K^2/3-2*kappa*rho`; defines `M=M_geometry-kappa*radial_momentum_source`; its own docstring says source law, index placement, and normalization remain provisional. |
| `engine/scalar_system.py` | `e2c2dae9a016b4d6454cf74a9eae8b0167d5733d` | Scalar coordinate-time RHS includes shift advection; `Pi) RHS has `alpha*beta_dm*rho_dm` when the optional DM source is supplied. Scalar projections divide the COSMOS-`phi` sector by `KAPPA=8*pi`, while S/D use scale 1. |
| `engine/matter_rhs.py` | `0153c232ee15b5a9d3dbac37ea858e4d2d730011` | Fluid RHS is conservative face-flux differencing plus source terms; this does not define the scalar S spatial operator. |
| `engine/stress_energy.py` | `ac3e7518c62bb72e2218a61489fa486bb787809f` | Canonical total stress assembly delegates to the shared matter projection; optional EM projections must already be consumer-normalized. |
| `engine/v55_pirk_adapter.py` | `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3` | Matter source enters the primary/Lambda geometry RHS through total projections; the Lambda source shape is `-8*pi*lambda_m*alpha*total[j]/a`. This is code shape, not a derivation of the requested covector source law. |

**Fail-closed finding:** the Brown observer's declared source commit differs from the #421 accepted-state source anchor. This may be an intentional diagnostic/source-root separation, but the manifest and execution contract must prove it before running anything. No observer run is authorized until the identity is reconciled. The source code trace above is a code-level inventory only; it does not settle the continuum Bianchi/source derivation.

## 5. Derivation obligations

Each row must end with a derivation citation or an explicit unresolved statement. A small residual is not a substitute for the derivation.

| ID | Obligation | Required proof artifact | Current status |
|---|---|---|---|
| D1 | Exact accepted definitions and index/density convention for M_r, H, Q_r, alpha, shift, K and time derivative | Equation-to-code crosswalk with each raising/lowering and density factor explicit | Open |
| D2 | Translate Brown's vacuum structural equation to the accepted coordinate-time / shift / Lie split | Signed term-by-term derivation; identify which terms are already included in the measured tangent | Open |
| D3 | Derive matter-coupled Hamiltonian-gradient contribution | Derive from selected Einstein/Bianchi conventions and exact accepted H definition; do not infer from residual fitting | Open |
| D4 | Derive source coefficient/sign and map to accepted current tangent | Source law with index, lapse, density weight, unit normalization, and exchange sign explicit | Open |
| D5 | Derive reference-metric connection response in the full matter-coupled propagation relation | Coefficient/sign of each C^r derivative term, stencils/parity, and a check against double-counting Lambda-dependent H terms | Open |
| D6 | Relate pointwise, FV, and expanded stress divergences | Coordinate/measure transformation, boundary/origin treatment, and equivalence proof or explicit non-equivalence | Open |
| D7 | Close total six-sector stress-energy/exchange accounting | Sector equations and normalized total ledger, including COSMOS-phi normalization premise | Open |

### D2/D3 caution: do not transplant an unverified formula

PR #421 records Brown's spherical-BSSN Eq. (11b) Hamiltonian-gradient slots as
`-(D_r alpha / 3) H + (alpha / 6) D_r H`.
The different expression `-D_i(alpha) H - (alpha/2) D_i(H)` is not interchangeable as written; a global rescaling/sign flip of H does not fix the relative coefficient sign. A time-split translation must be derived from the accepted evolution equations and cannot be assumed to resolve this coefficient mismatch.

Brown Eq. (11b) is a vacuum structural formula. Matter-dependent terms must be derived independently from the selected source-aware Einstein/Bianchi conventions; they cannot be inferred solely from the Brown source-bridge residual.

### D4 caution: preserve candidate labels

Keep these as separate values:
- `S_Q_candidate = -8 pi alpha Q_r` — provisional candidate only.
- `S_act` — accepted-update geometric injection plus current tangent.
- `Delta_S = S_act - S_Q_candidate` — source-map difference.
- `B_tot - B_geom` — Brown matter-dependent operator difference.
- `bridge_gap = (B_tot - B_geom) - Delta_S` — algebraic comparison, not a physical closure criterion.

Do not adjust signs, coefficients, density weights, or normalizations to reduce `D_tot`. Derive them from the accepted-source equations and archive premises first.

### D7 caution: exchange is not automatically total conservation

The archive-locked DM–COSMOS-phi exchange pair and the stored-variable unit map apply to that sector pair. They do not by themselves establish the normalized divergence of the total six-sector stress tensor. Preserve the archive's documented conflicting earlier exchange-sign sentence as provenance debt; do not silently rewrite it. Resolve the total-sector statement against the canonical master checkpoint and accepted source equations with a cited premise.

## 6. Minimal next execution sequence

1. **Source extraction:** complete the archive master/checkpoint read at its exact valid path and resolve the accepted-artifact versus Brown-observer source-pin distinction. The initial code trace is recorded in §4; continue recording path, blob SHA, relevant equation/function, and premise for every source.
2. **Signed symbolic ledger:** write each candidate term with its sign, index position, density weight, units/normalization, derivative operator, and whether it belongs to the measured tangent or the Brown RHS. Unknown entries remain explicitly unknown.
3. **Derive before measuring:** complete D1–D5 as far as the sources permit. If a premise is genuinely missing or conflicting, report it as a blocker rather than choose by residual size.
4. **Reconcile existing arrays:** use only SHA-verified #420 artifacts and evaluate whether the derived term mapping is already measurable from saved operands. Do not rerun the existing resolution campaign.
5. **Targeted observer only if needed:** if a derived term cannot be audited from existing arrays, specify the exact missing observable and implement a diagnostic-only post-processor. It must not advance the kernel or write trajectories.
6. **Report outcome:** distinguish derived, algebraically checked, numerically observed, provisional, and unresolved statements. Include signed centre cells 0/1 separately from fixed-radius values and global maximum if any targeted residual is recomputed.

## 7. Fail-closed provenance and execution contract

- Branch base is PR #421 head `6dfc53c48aaecc18b8389a75b6de7b14340747d1`.
- Accepted-state production-source anchor recorded by #421 is `4aec88d0e3ff028e0a65a8ac36b30948625df792`; separately reconcile the Brown observer's declared source commit `7f47669a1c0edfafef9fc5f27b2db80840d3f035` against the accepted-state manifest before execution.
- Vendor pin must equal `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Verify each accepted artifact digest against its manifest before reading; preserve original files and record their SHA-256.
- Abort on any protected-path mismatch in `engine/`, `vendor/`, `.gitmodules`, or `requirements.txt`.
- No production-kernel call, no new trajectory, no production changes, no source/vendor pin changes, no runner-dispatch matrix.
- All scientific-admission flags remain false until independently authorized.

## 8. Machine-readable status seed

```json
{
  "schema": "matter-coupled-propagation-derivation-ledger/v1",
  "status": "research_only_open",
  "conventions_source_pr": 421,
  "conventions_source_head": "6dfc53c48aaecc18b8389a75b6de7b14340747d1",
  "accepted_state_source_pin": "4aec88d0e3ff028e0a65a8ac36b30948625df792",
  "vendor_pin": "d6052d605673ce9d82cdc99b0df79855fa2ea215",
  "accepted_states": {"time": 0.60, "resolutions": [60, 120, 240], "reuse_only": true},
  "existing_candidate_residual_campaign_reused": true,
  "new_evolution_performed": false,
  "production_changes": false,
  "derivation_obligations": {
    "accepted_variable_index_density_map": "open",
    "time_split_translation": "open",
    "hamiltonian_gradient_derivation": "open",
    "matter_source_sign_normalization": "open",
    "reference_metric_connection_response": "open",
    "pointwise_fv_expanded_equivalence": "open",
    "six_sector_total_stress_energy_conservation": "open"
  },
  "flags": {
    "full_matter_coupled_propagation_identity_closed": false,
    "matter_source_sign_and_normalization_reviewed": false,
    "reference_metric_connection_coupling_convention_frozen": false,
    "total_stress_energy_sector_conservation_closed": false,
    "scientific_admission": false
  }
}
```

## 9. Source records

- PR #421 convention map: `docs/research/full_matter_coupled_propagation_convention_map_2026-10-25.md`.
- Accepted-state sector matter-source mapping: `docs/research/accepted_state_sector_matter_source_mapping_2026-10-25.md`.
- Brown matter-source bridge and absolute residual: `docs/research/brown_matter_source_bridge_absolute_residual_2026-10-25.md`.
- Pointwise versus FV total-stress divergence: `docs/research/pointwise_vs_fv_total_stress_divergence_2026-10-25.md`.
- Brown structural reference cited by the map: J. David Brown, “BSSN in Spherical Symmetry,” arXiv:0705.3845, Eq. (11b).
- Archive master/checkpoint and the accepted source definitions must be read at the exact pinned revision before a convention can be promoted from open to resolved.

**Admission:** all five flags remain false.


## 10. Source-only derivation refinement — sector exchange is not total divergence (2026-10-25)

**Source pin inspected:** `4aec88d0e3ff028e0a65a8ac36b30948625df792`; vendor remains `d6052d605673ce9d82cdc99b0df79855fa2ea215`. This section records source semantics only. No RHS, initializer, observer, or evolution was executed.

### 10.1 Archive sign record and explicit provenance conflict

The Library's canonical `GEAR_V5_5_MASTER_CHECKPOINT_2026-10-05.md` locks the exchange pair as:

```text
nabla_mu T_DM^{mu nu} = + beta_DM * rho_DM * nabla^nu(phi)
nabla_mu T_phi^{mu nu} = - beta_DM * rho_DM * nabla^nu(phi)
```

These two sector divergences cancel at the continuum level. The older Library file `GEAR_V5_5_COVARIANT_GLOBAL_GATE_2026-10-05_v2.md`, §13, contains an internal sign contradiction: its first proposed covariant DM exchange sentence uses the opposite DM sign, while its later explicit archive-normalized sign statement agrees with the master checkpoint. Preserve the contradictory sentence as provenance debt; do not silently rewrite it or resolve the conflict by fitting a residual.

The pinned accepted-source code shape is consistent with the master checkpoint's later explicit pair:

- [`engine/valencia.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/valencia.py), blob `03ed86398f424c03d9d6a3e646da9222a34fac5c`, defines the lower-covector DM source as `q_t = beta_DM * rho_DM * partial_t(phi)`, `q_r = beta_DM * rho_DM * partial_r(phi)`.
- [`engine/matter_system.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/matter_system.py), blob `0969f4d996e4f8e0fc5b70a7dff147f024907302`, states the mixed-stress source equation and applies this covector to the conservative RHS. The accepted Eulerian-energy and radial-momentum source slots are `q_E = -sqrt(gamma) * (q_t - beta^r * q_r)` and `q_Sr = alpha * sqrt(gamma) * q_r`; the variable-mass rest-current source is a separate slot.
- [`engine/scalar_system.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/scalar_system.py), blob `e2c2dae9a016b4d6454cf74a9eae8b0167d5733d`, uses `partial_t(phi) = alpha * Pi + beta^r * partial_r(phi)` and adds `+alpha * beta_DM * rho_DM` to the coordinate-time Pi RHS.
- [`engine/v55_matter.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/v55_matter.py), blob `0ee3c51fac7dc46c63feba2747a750542bcee3d0`, identifies the RHS density as raw archive-normalized DM density before stress projection. `engine/matter_system.py` applies a `1/kappa` conversion to fluid Einstein-source projections, and `engine/scalar_system.py` applies the same `1/kappa` scale to the stored COSMOS-phi stress projection, where `kappa = 8*pi*G` and `G=1` in this source path.

Thus the normalized continuum sector-pair convention is:

```text
nabla_mu T_DM,phys^{mu nu}  = +(beta_DM * rho_DM,archive / kappa) * nabla^nu(phi_stored)
nabla_mu T_phi,phys^{mu nu} = -(beta_DM * rho_DM,archive / kappa) * nabla^nu(phi_stored)
```

This is a source-convention and unit-map result, not a measurement that the cell-centred scalar and HLL fluid updates cancel exactly at every discrete stage.

### 10.2 The key distinction for candidate S_Q

Use separate notation:

- `q_DM_nu = +beta_DM * rho_DM,archive * partial_nu(phi_stored)`: the **DM sector exchange covector** consumed by the Valencia RHS.
- `D_total_r = nabla_mu T_total^{mu}{}_{r}`: the **total stress-divergence residual** that enters a source-aware Bianchi/constraint-propagation derivation, if that derivation produces such a term in the selected convention.

They are not interchangeable. The scalar sector carries the opposite exchange source, while the total divergence also depends on all included sectors and on the chosen discrete representation. Therefore the sector value `q_DM_r` must not be substituted directly for `D_total_r`.

The candidate `S_Q = -8*pi*alpha*Q_r` remains unadmitted until the derivation defines `Q_r` explicitly as the appropriate **total** lower-index divergence in the same metric, density, lapse, time-split, and discrete-measure convention as `Mdot_r`. Do not infer that definition from a small residual.

The three existing representations remain distinct: pointwise covariant divergence, finite-volume/cell-integrated divergence, and expanded product-rule divergence. The source trace does not select one by residual magnitude or convergence order.

### 10.3 What this source trace resolves—and what it does not

| Sub-obligation | Result from pinned source/archive | Status for full propagation identity |
|---|---|---|
| DM–COSMOS-phi exchange sign at archive-sector level | Master checkpoint and accepted RHS source shape agree on DM positive / scalar negative; one contradictory earlier gate sentence remains preserved as provenance debt. | **Resolved at sector-pair convention level only.** |
| Stored density and stress projection scale | RHS consumes archive-normalized rho_DM; fluid Einstein-source projections and COSMOS-phi stress projection apply 1/kappa. | **Code map recorded; total Q_r normalization still open.** |
| Coordinate-time source operands | partial_t(phi) = alpha*Pi + beta^r*partial_r(phi); scalar source carries alpha; fluid source is mapped into its conservative energy/momentum variables. | **Source-placement map only; no propagation time-split derivation.** |
| Total stress divergence and S_Q | Sector exchange is not total divergence; pointwise/FV/expanded forms remain distinct. | **Open; do not promote S_Q = -8*pi*alpha*Q_r.** |
| Hamiltonian-gradient translation | Brown Eq. (11b) slots remain -(D_r alpha)*H/3 + alpha*D_r(H)/6; standard ADM subsidiary coefficients cannot be substituted without deriving the formulation/time split. | **Open.** |
| Reference-metric connection response in the full source-aware tangent | Section 10.5 derives the candidate matter-current pairing in one Lambda RHS channel. The exact time derivative of the reference-metric target and its discrete stencils/parity must still be carried into Mdot_r. | **Open.** |
| Total six-sector discrete conservation | Existing ledger verifies signed additivity/bookkeeping, not cancellation of the total discrete stress divergence across scalar and HLL fluid updates. | **Open.** |

### 10.4 Next derivation gate

No additional observer is justified by this source trace alone. Continue from the accepted source equations and frozen M_r/H/C^r definitions, and derive the Brown-vs-accepted time split while stating exactly which total-divergence representation is used. If the source record does not determine an identity or coefficient, leave it unresolved rather than infer it from the #420 residual.

### 10.5 Accepted Lambda-stage momentum-source combination (exact source algebra)

A further part of the reference-connection path can be resolved directly from the pinned vendor/adapter source. It narrows the open question but does not close the full propagation identity.

- In vendor `vacuum_bssn.py` at the pinned gitlink `d6052d605673ce9d82cdc99b0df79855fa2ea215` (file blob `099ce52ba9e6dcd38d50c2915248c76c4236f7f4`), `lambda_l2_rhs` includes `+lambda_m * alpha * momentum_constraint(grid, state) / a`. The vendor `momentum_constraint` is the geometry-only quantity; denote it `M_r,geom` for this same-state algebra.
- In accepted-source `engine/v55_pirk_adapter.py` (blob `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3`), `lambda_l3_with_matter` adds `-8*pi * lambda_m * alpha * total["j"] / a` to the Lambda RHS. Here `total["j"]` is the adapter's canonical combined scalar/fluid radial-momentum projection.
- Therefore the *paired terms as implemented*, with `G=1` and the candidate source convention, combine algebraically as:

```text
+ lambda_m * alpha / a * M_r,geom
- 8*pi * lambda_m * alpha / a * j_total
= + lambda_m * alpha / a * (M_r,geom - 8*pi*j_total)
= + lambda_m * alpha / a * M_r,total_candidate
```

This proves a source-code identity **within each same-stage RHS evaluation** when the vendor `lambda_l2` and adapter `lambda_l3` blocks are evaluated on the same state. For the existing #420 *semi-discrete* `Mdot` tangent, that same-state RHS pairing is the relevant source-level result: it rules out a simple omission of the matter momentum projection from this Lambda-stage channel. It does not prove the full `C^r` tangent or the complete `Mdot_r` propagation identity. Separately, the finite accepted step in `production_kernel.py` uses a PIRK combination of start-state `lambda_l2`, final-primary `lambda_l2`, start-state `lambda_l3`, and primary-predictor `lambda_l3`; finite-step conservation therefore cannot be inferred from the same-stage RHS pairing alone.

### 10.6 Chain rule for the reference-metric connection target under the accepted CMC metric RHS

The pinned native CMC explicit-metric RHS implementation is `engine/native_geometry/evolution.py`, blob `7b47d700680b12ab3c82cb6a556d40580b103ce9`. Its shared metric core requires `beta=0` for the CMC override and returns:

```text
a_dot = -2*alpha*a*Aa
b_dot = -2*alpha*b*Ab = +alpha*b*Aa,  where Ab = -Aa/2
X_dot = alpha*X*K/3
```

The frozen reference-metric target uses only a and b; X is absent:

```text
Lambda_metric = D_r(a)/(2*a^2) - D_r(b)/(a*b) + 2/r*(1/b - 1/a)
C^r = Lambda - Lambda_metric
```

Substituting the accepted CMC `a_dot` and `b_dot` into the exact chain rule gives the operator expression:

```text
Lambda_metric_dot
  = -D_r(alpha*a*Aa)/a^2
    + 2*alpha*Aa*D_r(a)/a^2
    - D_r(alpha*b*Aa)/(a*b)
    - alpha*Aa*D_r(b)/(a*b)
    - 4*alpha*Aa/(r*a)
    - 2*alpha*Aa/(r*b)

C_dot^r = Lambda_dot - Lambda_metric_dot
```

Here `D_r` acting on products means the same pinned linear derivative operator applied to the product array; do not replace it with a finite-difference product-rule expansion. This is the exact time derivative of the target's discrete expression at a frozen CMC slice if the same fixed-grid operator is used for both the target and RHS operands. It is a derived operator formula, **not a measured closure result**.

I inspected the existing #420 candidate NPZ inventories: they carry the saved `Mdot`, Brown/source terms and six-sector arrays, but do not expose `a_dot`, `b_dot`, or `Lambda_metric_dot` profiles. So the connection tangent cannot currently be evaluated from those bundles alone. The next step is to finish mapping `Lambda_dot` through the accepted PIRK split (including which stage supplies each L2/L3 term) and then decide whether one targeted, RHS-only postprocessor is necessary to emit these missing operands. No evolution or campaign is warranted.

No production/source/vendor/dependency pin or accepted-state artifact changed. All five scientific-admission flags remain false.


## 10.7 Accepted PIRK split for Lambda and the reference-metric response (source-only)

**Finding:** the accepted step's Lambda time split can now be written explicitly from the pinned production kernel. This closes the *stage-assignment bookkeeping question*, not the full \(M_r\) propagation derivation.

**Exact source anchors inspected**

- engine/production_kernel.py — blob 1720615628a98a8d7cfc9fcadcc984fe28c88be1; accepted source anchor 4aec88d0e3ff028e0a65a8ac36b30948625df792.
- engine/v55_pirk_adapter.py — blob be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3.
- engine/native_geometry/evolution.py — blob 7b47d700680b12ab3c82cb6a556d40580b103ce9.
- Pinned vendor vacuum_bssn.py — gitlink d6052d605673ce9d82cdc99b0df79855fa2ea215, file blob 099ce52ba9e6dcd38d50c2915248c76c4236f7f4.

### Stage names and accepted update

Use these labels; they refer to different stage states, not interchangeable samples:

- \(0\): start state after the step refreshes the CMC lapse and sets \(\beta=B=0\).
- \(p\): primary predictor \(g_1\), formed after the first explicit \(a,b,X\) prediction, scalar/fluid prediction, CMC-lapse solve, and the primary \(A_a,K,\Lambda\) update.
- \(f\): final-primary operand for the last implicit \(L_2\) evaluation. It uses the final explicit \(a,b,X\), the predictor-resolved lapse, and the final \(A_a,K\); it is a mixed split-stage object, not itself the accepted endpoint.
- \(n+1\): the completed accepted state, after the final Lambda update and outer-boundary/CMC handling.

The final accepted Lambda update assembled in ProductionKernel.step is

\[
\Lambda_{n+1}^{\mathrm{pre-boundary}}-\Lambda_0
=\frac{\Delta t}{2}\left(
L_{2,0}^{\Lambda}+L_{2,f}^{\Lambda}
+L_{3,0}^{\Lambda}+L_{3,p}^{\Lambda}
\right).
\]

The predictor computation of \(g_1^\Lambda\) separately uses \(L_{2,0}^{\Lambda}\), \(L_{2,p}^{\Lambda}\), and \(L_{3,0}^{\Lambda}\). Its \(L_{2,p}^{\Lambda}\) is **not** the second \(L_2\) operand in the final accepted update. That operand is \(L_{2,f}^{\Lambda}\), evaluated after the final implicit-primary \(A_a,K\) update. This distinction is important when reconstructing the tangent or connection budget.

The CMC metric fields use a different split: \(a,b,X\) are updated by the explicit trapezoid from the start and primary-predictor explicit RHS values. For the pinned native CMC RHS with \(\beta=0\), \(A_b=-A_a/2\),

\[
E_a=-2\alpha a A_a,\qquad
E_b=-2\alpha b A_b=+\alpha b A_a,\qquad
E_X=\alpha X K/3.
\]

Thus the exact accepted metric increments use \(E_{a,0},E_{a,p}\) and \(E_{b,0},E_{b,p}\), not the \(f\) state used by final \(L_2^\Lambda\).

### Lambda source content and the finite-step mismatch that must not be hidden

With zero shift and the pinned vendor's \(KO\_EPSILON=0\), its Lambda terms reduce to

\[
L_{2}^{\Lambda}
=-\frac{2}{a}\left(A_a D_r\alpha+\alpha D_r A_a\right)
-\frac{4\alpha(A_a-A_b)}{rb}
+\lambda_m\frac{\alpha}{a}M_{r,\mathrm{geom}},
\]

\[
L_{3,\mathrm{vendor}}^{\Lambda}=2\alpha A_a\Lambda.
\]

The first-party adapter adds to each stage's \(L_3^\Lambda\) the matter-current term

\[
L_{3,\mathrm{matter}}^\Lambda=-8\pi\lambda_m\frac{\alpha}{a}j_{\mathrm{total}}.
\]

At one and the same RHS stage, the two momentum-source terms therefore pair algebraically as
\[
\lambda_m\frac{\alpha}{a}
\left(M_{r,\mathrm{geom}}-8\pi j_{\mathrm{total}}\right).
\]
This is a source-code identity for that stage only. In the *accepted finite-step* update the geometry-only momentum term is quadratured at \(0\) and \(f\), whereas \(j_{\mathrm{total}}\) is quadratured at \(0\) and \(p\). Its signed source contribution is therefore

\[
\frac{\lambda_m\Delta t}{2}
\left[
\frac{\alpha_0}{a_0}M_{r,\mathrm{geom},0}
+\frac{\alpha_f}{a_f}M_{r,\mathrm{geom},f}
-8\pi\left(
\frac{\alpha_0}{a_0}j_0+
\frac{\alpha_p}{a_p}j_p
\right)
\right].
\]

Because the second-stage operands are evaluated on different states, the same-stage cancellation cannot be promoted to a finite-step conservation or \(M_r\)-propagation identity. The index, normalization, and Bianchi source map in the earlier obligations remain open.

### Reference-metric response: exact discrete comparison versus instantaneous chain rule

The existing chain-rule expression in §10.6 remains the source-level instantaneous derivative of the sampled target \(\Lambda_{\mathrm{metric}}(a,b)\) when the same fixed \(D_r\) is applied to the product arrays. For an *accepted finite step*, the non-ambiguous discrete quantity is instead

\[
\Delta\Lambda_{\mathrm{metric}}
=\Lambda_{\mathrm{metric}}(a_{n+1},b_{n+1})
-\Lambda_{\mathrm{metric}}(a_0,b_0),
\qquad
\Delta C^r=\Delta\Lambda-\Delta\Lambda_{\mathrm{metric}}.
\]

Evaluate that difference from the actual start and completed accepted arrays, not by replacing it with a single-stage \(\Delta t\,\dot\Lambda_{\mathrm{metric}}\). The target is nonlinear in \(a,b\), and discrete product/chain/composition defects are separately observed; they must not be silently folded into the physical source term.

### Existing-data availability check — completed read-only schema audit

The read-only schema inventory is complete. It reused downloadable archived artifacts only; no observer, RHS, or evolution was run.

- **#420 full-profile states:** the N=60/120/240 artifacts from run 40156226086 (artifact IDs 12584064724, 12584164579, 12584980156) contain endpoint profiles at \(t=0.60\). Their profiles.npz arrays include the completed accepted \(a,b,X,\alpha,A_a,K,\Lambda\), scalars, matter, and endpoint constraints; they do not contain the full start/predictor/final-implicit stage profiles.
- **#415 saved-state stage reports:** artifacts 12577357052 (N=120) and 12577886263 (N=240) are report.json-only, schema **saved_state_cmc_operator_contract_one_step_v1**. They include summaries for start-lapse-refreshed, explicit-predictor, and primary-predictor stages, plus final-implicit and Lambda/connection budget summaries at selected cells. They do not include a full radial NPZ of stage-resolved scalar/fluid profiles or stage-resolved \(j_{\mathrm{total}}\). The matter-\(L_3\) budget remains increment-inferred in that report.
- **#408 raw witness:** artifact 12160420638 contains 131 arrays and a stored \(J_{\mathrm{total}}\)/Lambda RHS snapshot, but it is the separate N=60 state at \(t=0.015\) with state digest c6acb371…, and declares the older production-source baseline 7f47669a1c0edfafef9fc5f27b2db80840d3f035, not the #420 accepted-state source anchor 4aec88d0e3ff028e0a65a8ac36b30948625df792. It is useful source/array-schema evidence but is not the target \(t=0.60\) start/predictor witness. Its companion shared_state_signed_arrays.npz has 20 signed residual/divergence arrays, not the missing stage profiles.

**Exact missing operands for the target accepted step:** full radial \(j_0\) and \(j_p\) arrays, reconstructed from the same-stage scalar/fluid fields and geometry; the matching \(\alpha/a\) arrays; and full radial \(M_{r,\mathrm{geom}}\) values for the \(0\) and \(f\) L2 operands. The complete \(a,b,\Lambda\) stage arrays are also required to compute the actual discrete \(\Delta\Lambda_{\mathrm{metric}}\) and \(\Delta C^r\), rather than rely on selected-cell summaries or a single-stage chain-rate approximation.

### Stage-capture implementation and result — completed diagnostic-only copied step

The exact-source capture is now implemented in tools/audit_pirk_stage_resolved_momentum_source.py, with its guarded matrix workflow in .github/workflows/issue4-pirk-stage-resolved-momentum-source.yml and helper tests in tests/test_pirk_stage_resolved_momentum_source.py.

**Successful workflow:** [run 40235790388](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388). All three lanes passed the helper tests, accepted-source/vendor guards, copied-step capture, raw-operand closure assertions, and artifact upload. Inputs were the three checksum-verified #420 endpoint-profile artifacts from [run 40156226086](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40156226086).

Output artifacts:
- N=60: [artifact 12614314612](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12614314612). report.json SHA-256: bf777499cc9e20b664c7cf906e3ddfb297fba53f44f312d65d62b80a23fae94b; stage_arrays.npz SHA-256: 942325e3bfd9cfcd95ad7b2f0b5511c6303045f2a8eb36c4ec6f5144a0cc3306.
- N=120: [artifact 12615190935](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12615190935). report.json SHA-256: fd33c9fe07847df7eb102244ca344021e70530a918ac5d5d0e65f70ccd65a343; stage_arrays.npz SHA-256: 8d449fbfd224f5608292c375e2c47e625b1e9a70f8562686cbdafd996503fd90.
- N=240: [artifact 12614533930](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40235790388/artifacts/12614533930). report.json SHA-256: e52f852694c5a5c36d208a2439af813e699b2e0f7c53e67bfecef08bd21f1658; stage_arrays.npz SHA-256: 8b85b445c4b1213cb6904734ed6fcf7cd9344d857a112685e718b6bef922d141.

The input profile and report member hashes were checked against each artifact's SHA256SUMS.txt. Every lane confirmed production source anchor 4aec88d0e3ff028e0a65a8ac36b30948625df792, vendor pin d6052d605673ce9d82cdc99b0df79855fa2ea215, and all 11 protected source blobs matching the report's pin table. The saved input state digest was unchanged by the copied-step run. No production, initializer, dependency, vendor, or source-pin edit was made.

**Scope precision:** each lane restores the saved accepted endpoint at t approximately 0.60, evaluates one copied in-memory PIRK step, and emits arrays only. The copied output is at t=0.64, 0.62, and 0.61 for N=60, 120, and 240, respectively. This does not replay the historical accepted interval ending at t=0.60 and writes no accepted trajectory. Captured arrays include geometry, scalar and fluid states and RHS terms for the start-lapse-refreshed, explicit-predictor, and primary-predictor stages; the final-primary mixed-stage L2 geometry and result; independently reconstructed total j at start and primary predictor; geometric momentum constraint, alpha/a, and separated vacuum/matter L3 terms. The NPZ contains full radial profiles rather than only centre summaries.

#### Numerical and source-accounting results

| Resolution | Maximum absolute Lambda update residual, r≤10 | Maximum raw-operand bound, r≤10 | Maximum L3 reconstruction residual |
|---:|---:|---:|---:|
| N=60 | 1.36e-19 | 7.62e-16 | 2.71e-20 |
| N=120 | 1.81e-19 | 7.75e-16 | 2.71e-20 |
| N=240 | 2.01e-19 | 7.78e-16 | 5.42e-20 |

The copied-step Lambda update reconstructs from the exact captured L2-at-start, L2-at-final-primary, L3-at-start, and L3-at-primary-predictor operands within the raw-operand envelope. Independently assembled canonical total j reproduces the adapter's L3 matter term. This verifies stage-source bookkeeping and implementation consistency only; it does not admit the physical normalization of total j, prove S_Q=-8π α Q_r, or close momentum-constraint propagation.

The difference between the accepted-split momentum-source combination and a counterfactual that substitutes the primary-predictor geometric momentum operand for the exact final-primary L2 operand is nonzero. This is expected evidence that the second-stage source operands are not interchangeable, not a residual to fit away:

| N | Difference at cell 0 | Cell 1 | Cell 2 | Cell 6 |
|---:|---:|---:|---:|---:|
| 60 | +1.3483e-8 | +1.4570e-8 | +1.5435e-8 | +2.6263e-9 |
| 120 | +2.4748e-9 | +1.7133e-9 | +1.6554e-9 | +1.7734e-9 |
| 240 | +5.2274e-10 | +3.0602e-10 | +1.8436e-10 | +2.2799e-10 |

**Critical interpretation limit:** the report's finite-step connection residual is zero by definition because it reassembles C-out minus C-in from those same endpoints using delta Lambda minus delta Lambda-metric. It is not an independent propagation test and must not be counted as one. The nontrivial checks here are the stage-resolved Lambda update reconstruction, direct canonical total-current/L3 decomposition, provenance checks, and input immutability. This capture did not evaluate an independent full time derivative of M-r or a full matter-coupled propagation RHS.

### Next derivation gate

The stage-current/operand availability gap is now filled for this copied diagnostic step. Next, finish the signed equation-to-code crosswalk for the independent full partial-time M-r identity and total-stress source, then join that derivation to an exactly state-matched existing tangent only after the stage and lapse-state mapping is verified. Do not infer or fit S_Q from the captured values and do not launch another resolution or trajectory campaign merely to reproduce this result.

**Status remains unchanged:** the time-split assignments above are source-derived; the full matter-coupled \(M_r\) identity, physical \(S_Q\), reference-metric coupling admission, and total six-sector discrete stress-energy conservation remain open. Pointwise/FV/expanded divergence measures remain separate, and all five scientific-admission flags remain false.
