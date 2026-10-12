# Radiation source and boundary diagnostic inventory

Scope: diagnostic-only radiation four-force counterfactual in PR #9. This note inventories existing production-kernel witnesses to reuse; it does not add a new diagnostic layer or promote a physical model.

## Existing radiation and boundary witnesses

### 1. Outer-cell radiation characteristic state

The production diagnostics expose:

- `radiation_outer_vhat`: orthonormal radial velocity, √(γ_rr) v_r, at the outer cell.
- `radiation_outer_s_minus` and `radiation_outer_s_plus`: radiation characteristic speeds using c_s = 1/√3.
- `radiation_outer_ratio`: |S|/|E| at the outer cell.

The counterfactual adds `radiation_margin_i0` and `radiation_margin_admitted`, using the same ratio to report 1 - |S|/|E| at the centre and over the admitted interior. A positive margin means the local |S|/|E| admissibility inequality is satisfied.

### 2. Radiation transport flux witnesses

The existing diagnostics path computes the outer-face Valencia flux and HLL fluxes at selected interior faces, including a nearby interior face used to help distinguish boundary-closure effects from an outgoing transport feature. Reuse these measurements rather than introducing a parallel flux definition.

### 3. Geometric source terms

The conservative Valencia RHS includes the metric-source contributions for each species. In the radiation counterfactual, the optional four-force is an additional source term and is entered only when `beta_rad != 0.0`; the null coefficient leaves the existing source path unchanged.

When interpreting results, compare the existing transport/flux and metric-source witnesses wherever the saved stage-level trace makes that split available. Do not infer a term-by-term attribution from aggregate endpoint diagnostics alone.

### 4. Outer boundary condition

The incoming-light characteristic boundary reconstruction (`apply_light_constraint_boundary`) is applied during the predictor Lambda update and again on the completed slice. Candidate/null differences must be read alongside the outer-face flux and characteristic-speed witnesses so a boundary response is not mistaken for a direct four-force effect.

### 5. Constraint and geometry witnesses

The existing diagnostics include Hamiltonian residuals (including admitted-region reporting), momentum residual, connection residual C_Λ, determinant residual, and CMC residual / target-K-dot witnesses. These are the principal decision-table entries. Any influence of the radiation source on these geometry witnesses is mediated through the total stress-energy and the geometry evolution stages.

### 6. Explicitly unavailable: Misner–Sharp current/balance residual

The Misner–Sharp mass helper exists, but a consistent current and balance residual has not been assembled. The runner must retain explicit NaN values for current fields; do not substitute a mass, finite difference, or other proxy and label it a current.

## First-run interpretation protocol

1. **Null arm:** `beta_rad = 0.0` must preserve the current radiation RHS. The focused regression test checks omitted/default versus explicit-zero RHS arrays and preserves the null call pattern with `dphi_t=None`, `dphi_r=None`. After that test passes, compare the outer-face/boundary witnesses between null executions if those values are included in the saved artifact.
2. **Candidate arm:** `beta_rad = -0.04` may change radiation admissibility margins, admitted H/M residuals, and the total source seen by geometry. Such changes are observations, not evidence by themselves that the coefficient or coupling is physically correct.
3. **Attribution:** use a saved stage-level flux-versus-metric-source split when available. If the artifact contains only endpoint aggregates, describe only endpoint differences and do not claim a term-level causal decomposition.
4. **Admission:** preserve the four one-step comparisons (N=80/160, null/candidate), JSON artifact, SHA-256, and explicit decision table. No long campaign, parameter fitting, production promotion, or merge follows automatically from this diagnostic.

## Execution status

This note and the focused regression test are source-level additions only. Unless a subsequent run artifact or CI result is linked here, tests and the four one-step comparisons must be treated as **not executed / not admitted**.
