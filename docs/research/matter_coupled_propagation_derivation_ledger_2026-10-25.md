# Matter-coupled momentum-propagation derivation ledger — 2026-10-25

**Status:** research-only derivation ledger; no identity admitted.  
**Parent convention record:** PR #421, head `6dfc53c48aaecc18b8389a75b6de7b14340747d1`.  
**Accepted-state production-source anchor:** `4aec88d0e3ff028e0a65a8ac36b30948625df792`.  
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

## 4. Derivation obligations

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

## 5. Minimal next execution sequence

1. **Source extraction:** read the exact pinned accepted update, native constraint helper, Brown observer, sector source mapper, divergence comparison, and archive master/checkpoint documents. Record path, blob SHA, relevant equation/function, and premise for every source.
2. **Signed symbolic ledger:** write each candidate term with its sign, index position, density weight, units/normalization, derivative operator, and whether it belongs to the measured tangent or the Brown RHS. Unknown entries remain explicitly unknown.
3. **Derive before measuring:** complete D1–D5 as far as the sources permit. If a premise is genuinely missing or conflicting, report it as a blocker rather than choose by residual size.
4. **Reconcile existing arrays:** use only SHA-verified #420 artifacts and evaluate whether the derived term mapping is already measurable from saved operands. Do not rerun the existing resolution campaign.
5. **Targeted observer only if needed:** if a derived term cannot be audited from existing arrays, specify the exact missing observable and implement a diagnostic-only post-processor. It must not advance the kernel or write trajectories.
6. **Report outcome:** distinguish derived, algebraically checked, numerically observed, provisional, and unresolved statements. Include signed centre cells 0/1 separately from fixed-radius values and global maximum if any targeted residual is recomputed.

## 6. Fail-closed provenance and execution contract

- Branch base is PR #421 head `6dfc53c48aaecc18b8389a75b6de7b14340747d1`.
- Accepted-state production-source anchor must equal `4aec88d0e3ff028e0a65a8ac36b30948625df792`.
- Vendor pin must equal `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Verify each accepted artifact digest against its manifest before reading; preserve original files and record their SHA-256.
- Abort on any protected-path mismatch in `engine/`, `vendor/`, `.gitmodules`, or `requirements.txt`.
- No production-kernel call, no new trajectory, no production changes, no source/vendor pin changes, no runner-dispatch matrix.
- All scientific-admission flags remain false until independently authorized.

## 7. Machine-readable status seed

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

## 8. Source records

- PR #421 convention map: `docs/research/full_matter_coupled_propagation_convention_map_2026-10-25.md`.
- Accepted-state sector matter-source mapping: `docs/research/accepted_state_sector_matter_source_mapping_2026-10-25.md`.
- Brown matter-source bridge and absolute residual: `docs/research/brown_matter_source_bridge_absolute_residual_2026-10-25.md`.
- Pointwise versus FV total-stress divergence: `docs/research/pointwise_vs_fv_total_stress_divergence_2026-10-25.md`.
- Brown structural reference cited by the map: J. David Brown, “BSSN in Spherical Symmetry,” arXiv:0705.3845, Eq. (11b).
- Archive master/checkpoint and the accepted source definitions must be read at the exact pinned revision before a convention can be promoted from open to resolved.

**Admission:** all five flags remain false.
