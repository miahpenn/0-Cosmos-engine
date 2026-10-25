# Full matter-coupled momentum-propagation convention map — 2026-10-25

**Status:** diagnostic convention inventory only. This document does not admit the continuum identity or authorize a production change.

## 1. Provenance and scope

- Documentation branch: `diag/full-matter-coupled-convention-freeze-2026-10-25`.
- Branch point / inspected PR #420 diagnostic head: `e98b8fa5bd67e9cf023e001e81a6792b5ce469d2`.
- PR #420 base SHA: `19fbc814a041a7aaf2678c0337a02f02905f40a4` (not the accepted-state production-source anchor).
- Accepted-state production-source anchor recorded by the existing #420 observers: `4aec88d0e3ff028e0a65a8ac36b30948625df792`.
- Pinned vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Accepted states reused: (t=0.60), (N=60,120,240), with existing artifact manifests verified by the #420 workflows.
- Scope: source-backed convention inventory and prior evidence crosswalk only. No evolution steps are advanced here.

## 2. Residual object and index position

**Resolved within the Brown-BSSN candidate only:** use the lower-index radial covector component (M_r) for the Brown structural evaluator. The repository's [vacuum Brown propagation observer](../tools/audit_brown_vacuum_momentum_propagation.py) explicitly describes its target as lower-index (M_r) and its connection constraint as contravariant (C^r). Its formula is sourced to J. David Brown, [“BSSN in Spherical Symmetry,” arXiv:0705.3845, Eq. (11b)](https://arxiv.org/abs/0705.3845).

This does **not** establish that all other residual arrays named `M`, (M^r), or a densitized (M_r) are identical. In particular, the first-party [native constraint helper](../engine/native_geometry/constraints.py) uses an explicit lower-radial source API but documents that the source law, index placement, and normalization are provisional. Any crosswalk between these objects must show the metric raising/lowering, any density factor, and source normalization explicitly. Do not silently equate them.

## 3. Accepted S-sector spatial RHS operator

**Resolved from the pinned accepted update; do not repeat this identification.** The production S-field spatial operator is the cell-centred, parity-aware derivative form

```text
L_S[S] = (X^2/a) * (
    D2(S)
    + (D1(b)/b - D1(a)/(2*a) + 2/r - D1(X)/X) * D1(S)
)
```

The corresponding isolated spatial term in the S momentum RHS is `alpha * L_S[S]`; lapse-gradient, curvature, potential, shift, and any conditional exchange operands are separate terms. This is not a fluid-HLL or finite-volume face-flux update. The reconstruction of the pinned S-current RHS spatial contribution closes to roundoff in the existing crosswalk.

Evidence: [S-sector RHS/operator crosswalk](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/docs/research/s_sector_divergence_operator_crosswalk_2026-10-25.md), [full-profile scalar-operator addendum](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/docs/research/full_profile_scalar_operator_addendum_2026-10-25.md).

## 4. Stress-divergence representation and source mapping

**Unresolved as a physical propagation convention.** Keep these as separate named quantities:
- pointwise covariant (
abla_mu T^mu{}_r), using the pointwise flux derivative and its metric-source terms;
- the cell-integrated / finite-volume divergence reconstructed from face fluxes and integrated sources;
- expanded covariant expressions assembled from differentiated factors.

The existing diagnostics identify a discrete radial-flux product-rule gap between some pointwise and expanded forms. The S-sector pointwise and FV values also differ substantially in the origin cell because they are not the same discrete observable. Neither the smaller magnitude nor a good observed order alone selects the physical authority. Any future equality claim must derive the measure, index, density, coordinate, and boundary/origin treatment from the accepted formulation.

The provisional source (S_Q=-8pialpha Q_r) used by the Brown crosswalk must remain labeled **candidate**, not accepted oracle. The measured matter/current composite (S_{m act}) is not to be interpreted in isolation: the existing Brown bridge shows that its difference from (S_Q) is nearly balanced by the Brown matter-dependent difference between total-constraint and geometry-only evaluations. The signed reassembly closes at approximately (5.08	imes10^{-21}), (7.46	imes10^{-21}), and (1.06	imes10^{-20}) for (N=60,120,240); this is algebraic reconciliation of the existing arrays, not independent proof of the continuum source law.

Evidence: [pointwise versus FV total-stress divergence](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/docs/research/pointwise_vs_fv_total_stress_divergence_2026-10-25.md), [sector source mapping](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/docs/research/accepted_state_sector_matter_source_mapping_2026-10-25.md), [Brown matter-source bridge](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/docs/research/brown_matter_source_bridge_absolute_residual_2026-10-25.md).

## 5. Brown references and Hamiltonian-gradient normalization

**Reference disambiguation resolved; candidate H definition mapped; matter-coupled gradient relation unresolved.** Two distinct publications must not be combined under one citation:
- J. David Brown, [“BSSN in Spherical Symmetry,” arXiv:0705.3845, Eq. (11b)](https://arxiv.org/abs/0705.3845), which is the structural formula implemented by the repository's Brown observer.
- J. David Brown and Lisa L. Lowe, [“Modifying the Einstein Equations off the Constraint Hypersurface,” gr-qc/0606008](https://arxiv.org/abs/gr-qc/0606008), a separate off-constraint evolution-equation construction.

Brown defines its vacuum Hamiltonian constraint as
`H = K^2 - K_ij K^ij + R`. With `K_ij = A_ij + (1/3) gamma_ij K` and trace-free `A_ij`, this becomes `H = R - A_ij A^ij + (2/3)K^2`. The repository's metric-only candidate `H_metric = R_metric - Aa^2 - 2*Ab^2 + (2/3)*K^2` has the same algebraic normalization when `A_ij A^ij = Aa^2 + 2*Ab^2`, with the spatial metric/Ricci and mixed-component conventions mapped consistently.

In the repository's Brown-BSSN formula, the Hamiltonian-gradient slots are
```text
-(D_r(alpha)/3) * H + (alpha/6) * D_r(H)
```
with `H` reconstructed as metric-only, independently of the evolved-Λ connection terms. The original proposal's `-D_i(alpha)*H -(alpha/2)*D_i(H)` is **not interchangeable as written**. More strongly, a single overall rescaling or sign flip of `H` cannot reconcile the formulas by itself: the relative coefficient sign between the lapse-gradient and Hamiltonian-gradient terms differs. Coordinate-time versus normal/Lie-split derivative conventions affect shift/advection bookkeeping, but do not alone establish a translation of these two slots. Do not use the proposed formula as the Brown Eq. (11b) expression without a full derivation from the selected field/evolution conventions.

**Remaining derivation before freezing the matter-coupled identity:** map the exact accepted-source definitions of `H`, `M_r`, `Q_r`, `K`, lapse, shift and time-derivative split; then derive the complete matter-coupled propagation relation from the selected Einstein/Bianchi conventions. Brown Eq. (11b) is a vacuum structural formula, so the matter terms must be independently derived and may not be inferred from the source-map residual. Keep the candidate's algebraic `H` normalization separate from whether the full matter-coupled gradient term, source term, and reference-metric response are accepted.

## 6. Reference-metric connection coupling

The existing Brown observer uses (chi=X^2) and a contravariant radial connection-constraint residual (C^r), with explicit (D_r^2 C^r) and (D_r C^r) terms. The metric-derived target for the connection state is documented as
```text
Lambda_metric = D1(a)/(2*a^2)
             - D1(b)/(a*b)
             + 2/r * (1/b - 1/a)
C^r = Lambda - Lambda_metric
```

The presence of this structural candidate is known; its complete mapping into the accepted, matter-coupled reference-metric propagation identity is **not** yet frozen. Keep metric-only Ricci/H and evolved-Λ-dependent Hamiltonian representations separate, as the existing observer does. Any future derivation must identify derivative parity/stencils and the coefficient/sign of every (C^r) derivative term. Do not infer this coupling from residual fitting or from the geometric closure alone.

## 7. Six-sector exchange and conservation

The existing per-sector observer covers (S), (D), COSMOS φ, dark matter, baryons, and radiation. Its sector additivity checks close to tight numerical tolerances. Additivity is bookkeeping, not total stress-energy conservation. The full total-(T)/sector-exchange gate remains open; the archived COSMOS-φ normalization/comparator needs its source-premise review. Keep this gate distinct from the Brown candidate residual and from the origin-cell regularity witness.

## 8. Numerical evidence correctly interpreted

For the full-continuous-profile accepted states, the existing matter-inclusive Brown candidate has global maximum absolute residuals (8.975	imes10^{-6}), (9.317	imes10^{-7}), and (4.133	imes10^{-8}) at (N=60,120,240). At fixed (r=2), its observed orders are approximately 3.87 and 3.93; at fixed (r=8), approximately 4.03 and 3.99. This supports a decreasing absolute residual for this candidate on these states. It does not alone prove the continuum identity, settle the regular-origin limit, or close total sector conservation. Relative RMS is secondary when the reference derivative becomes very small.

## 9. Explicit unresolved items / next gate

1. Derive the Hamiltonian-gradient normalization and derivative-splitting translation against the exact accepted-source equations and Brown Eq. (11b).
2. Establish the source-backed reference-metric connection coupling and show no double-counting of Λ-dependent H terms.
3. Derive the allowed (Q_r)/source representation under the exact index/density convention; keep pointwise, FV, and expanded diagnostics separate unless equivalence is demonstrated.
4. Complete the archive-backed six-sector total stress-energy/exchange identity, including the COSMOS-φ normalization premise.
5. Only after the above are source-pinned, decide whether any additional observer is needed. Do not rerun the existing (N=60/120/240) campaign or a runner-dispatch matrix just to reproduce the measured evidence.

## 10. Guardrails / admission flags

No production equation, initializer, gauge, boundary, coupling, floor/clamp, damping, dependency, vendor, or source pin is changed or authorized by this document. No evolution is advanced. The source/vendor pins above must be checked by any future run. All gates remain false:

- `full_matter_coupled_propagation_identity_closed = false`
- `matter_source_sign_and_normalization_reviewed = false`
- `reference_metric_connection_coupling_convention_frozen = false`
- `total_stress_energy_sector_conservation_closed = false`
- `scientific_admission = false`


## 11. Reference-metric connection-variable translation — candidate-level mapping

The pinned vendor and first-party native geometry helpers use the reference-metric target

```text
Lambda_metric = a'/(2*a^2) - b'/(a*b) + 2/r*(1/b - 1/a)
C_code = Lambda - Lambda_metric
```

Brown's spherical-BSSN connection constraint (Brown, arXiv:0705.3845, Eq. (10c)) is

```text
G_B^r = -a'/(2*a^2) + Gamma_B^r + b'/(a*b) + 2/(a*r)
Gamma_B,metric^r = a'/(2*a^2) - b'/(a*b) - 2/(a*r)
```

The algebraic translation

```text
Gamma_B^r := Lambda - 2/(r*b)
``

implies `Gamma_B,metric^r = Lambda_metric - 2/(r*b)` and therefore

```text
G_B^r = Gamma_B^r - Gamma_B,metric^r
      = Lambda - Lambda_metric
      = C_code
```

Thus, for the same accepted state and the same differentiation operators, the constraint residual entering Brown's `(C^r)'` and `(C^r)''` terms maps algebraically to the code's `Lambda - Lambda_metric`; the radial metric factor in Brown's coefficient is `g_theta_theta = b*r^2`, so `(g_theta_theta)'/g_theta_theta = b'/b + 2/r`. This removes the apparent raw-variable naming mismatch for the **Brown candidate's constraint field**. It is not a proof that the complete production matter-coupled tangent obeys Brown Eq. (11b), nor does it derive the matter-source law, Hamiltonian-gradient translation, or total stress-energy conservation.

**Status refinement:** `reference_metric_connection_variable_mapping = resolved_for_Brown_candidate_only`; `reference_metric_connection_coupling_convention_frozen = false` remains unchanged. The latter stays open until the complete source-aware propagation relation and exact accepted-update derivative/time-split conventions are independently derived. No source code, evolution, or admission state is changed by this algebraic note.


## 12. Archive-locked DM–COSMOS-phi exchange: signs versus normalization

The canonical `GEAR_V5_5_MASTER_CHECKPOINT_2026-10-05.md` locks the covariant exchange pair as

```text
nabla_mu T_DM^{mu nu}  = + beta_DM * rho_DM * nabla^nu(phi)
nabla_mu T_phi^{mu nu} = - beta_DM * rho_DM * nabla^nu(phi)
```

The same pair appears in the later explicit sign statement of `GEAR_V5_5_COVARIANT_GLOBAL_GATE_2026-10-05_v2.md` §13. That file's earlier sentence gives the opposite signs, so the archive retains a real provenance conflict; it must not be silently edited away. The implementation follows the master checkpoint / later explicit pair in its DM covector source and scalar `Pi)-RHS source shape: DM uses covector components built from `beta_DM*rho_DM*(partial_t phi, partial_r phi)`, while the scalar `Pi) RHS receives `alpha*beta_DM*rho_DM`. This code-shape match is not a substitute for deriving the exact total-`T) divergence with the accepted normal/time split.

**Normalization seam remains unresolved:** `rho_DM` supplied to the coupled RHS is the archive-normalized conservative density; the scalar and fluid Einstein stress projections apply the stored `1/(8*pi*G)` conversion on their stress arrays. A complete source map must derive the corresponding conversion and index/lapse factors between the archived exchange equation and the physical `Q_r) used in the matter-coupled constraint propagation identity. Do not infer the conversion by fitting the Brown residual. Until the exact archive-normalized coupling and accepted projection conventions are reconciled, keep `matter_source_sign_and_normalization_reviewed = false` and `scientific_admission = false`.
