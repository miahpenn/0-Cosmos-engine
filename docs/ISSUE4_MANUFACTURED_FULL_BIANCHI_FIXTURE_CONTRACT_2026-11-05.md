# Issue #4 manufactured full-Bianchi fixture contract

**Date:** 2026-11-05  
**Class:** `manufactured_regular_family`  
**Status:** manufactured fixture constructed and fail-closed loader validated; evaluator integration pending  
**Scope:** research-only, immutable, zero production steps; all scientific-admission and conservation-closure flags remain false.

## 1. Purpose and boundary

Supply a frozen, analytically prescribed regular-profile state and a declared one-parameter family around it to the existing Issue #4 full-Bianchi predictor. The predictor's signed equations, normalizations, centre/fixed-radius masks, pointwise/FV channels, and zero-step guards remain unchanged.

The manufactured class is a controlled mathematical discriminator. It is not a production-initialized state, a trajectory, or evidence of physical realism. It does not establish scientific admission or total-stress conservation.

The accepted-state fixture loader, accepted archive hashes, accepted-state regression values, and accepted-state report path remain unchanged.

## 2. Exact provenance contract

The manufactured archive must fail closed unless all of the following identities match:

| Identity | Required value |
|---|---|
| Diagnostic/source anchor | `4aec88d0e3ff028e0a65a8ac36b30948625df792` |
| Protected helper baseline | `7f47669a1c0edfafef9fc5f27b2db80840d3f035` |
| Vendor | `d6052d605673ce9d82cdc99b0df79855fa2ea215` |
| Tangent helper blob | `3f00a4e0a2063b86d2e007c5a825f9edfff61a4c` |
| Total-stress helper blob | `c10f9e8b9e654d5195fd55c493c3fe06ba5c7082` |

The fixture manifest must additionally record exact Git blob identities for the reviewed evaluator, observer, accepted loader, manufactured constructor, manufactured loader, and schema tests. Branch names are descriptive only and never substitute for commit/blob identities. Any required identity mismatch must stop before state construction, predictor evaluation, or report generation.

## 3. State and array interface

All cell-centred arrays are one-dimensional, finite `float64`, shape `(N,)`, and use the pinned vendor grid with centres `r_i=(i+1/2) r_max/N`. Resolution, `r_max`, `dr`, parity, and profile coefficients are explicit manifest fields.

The frozen base state must provide the accepted observer's canonical state keys:

- `grid.centers`
- `geometry.a`, `geometry.b`, `geometry.X`, `geometry.Aa`, `geometry.K`, `geometry.Lambda`, `geometry.alpha`, `geometry.beta`, `geometry.B`
- `scalars.S`, `scalars.PS`, `scalars.D`, `scalars.PD`, `scalars.phi`, `scalars.Pi`
- `matter.dark_matter.{rest,energy_t,momentum_r}`
- `matter.baryons.{rest,energy_t,momentum_r}`
- `matter.radiation.{rest,energy_t,momentum_r}`

The pinned total-stress projection must reconstruct the six named sectors, with their existing conventions and signs: dark matter, baryons, radiation, S, D, and COSMOS-`phi`. The constructor must not substitute a parallel hand-coded sector normalization. The report retains each sector's `rho`, `p_r`, `p_t`, and `j` projections (or the exact pinned helper's equivalent primitives), raw six-sector sums, and separately reported time-density component accounting.

The primary fixture has exactly zero shift and no electromagnetic carrier or EM-on content. Metric parity is declared per field: regular metric/scalar-even fields are even; the connection difference `C^r = Lambda - Lambda_metric` is odd. Centre regularity declarations are explicit, not inferred from sampled values.

## 4. Tangent semantics — essential distinction

The accepted-state tangent is a directional derivative along the pinned production RHS. Its helper obtains production rates, forms `U +/- epsilon * Udot_RHS`, and differences the measured total momentum constraint. The manufactured constructor is prohibited from calling that RHS.

Therefore the manufactured archive must **not** label its tangent a production-RHS or semi-discrete-production tangent.

Instead, define an explicit analytic one-parameter manufactured family `U(s)`, with `U(0)=U_0`. For each locked epsilon (`2e-5` and `1e-5`), independently construct and freeze `U(+epsilon)` and `U(-epsilon)`. The family must specify every changing geometry, lapse, scalar, and matter primitive; fields held fixed must be declared. No member may be obtained by calling the production initializer, production RHS, kernel step, vendor evolution step, or a time integrator.

Measure the manufactured directional witness using the same pinned discrete momentum-constraint and six-sector projection definitions:

```
Mdot_manufactured(epsilon)
  = [M_total(U(+epsilon)) - M_total(U(-epsilon))] / (2*epsilon)
```

Store the two signed arrays under the evaluator's existing paths:
- `tangent_eps2e-05.Mdot`
- `tangent_eps1e-05.Mdot`

The manifest must identify this as `analytic_family_central_difference`, record the family parameterization and plus/minus state digests, and state explicitly that it is **not** the accepted-state production-RHS tangent. This preserves the predictor's input interface without conflating the physical meaning of the two tangent classes.

## 5. Total-stress divergence witnesses

For each epsilon, pass the exact frozen triplet `(U_0,U(+epsilon),U(-epsilon))` through the pinned total-stress helper's existing directional projection/divergence routines. Preserve the helper's pointwise and finite-volume channels independently:

- `total_T_eps2e-05.total.div_T_r`
- `total_T_eps1e-05.total.div_T_r`
- `total_T_eps2e-05.total.div_T_r_fv`
- `total_T_eps1e-05.total.div_T_r_fv`

Do not synthesize either divergence from the Bianchi predictor, from the measured `Mdot`, or from the other divergence channel. Preserve the sector-resolved decomposition and its raw and time-density-accounted forms. Keep the two epsilon channels separate.

## 6. Archive, digest, and immutability contract

Use the same three-member archive layout as the accepted fixture:

- `report.json`
- `tangent_arrays.npz`
- `SHA256SUMS.txt`

The report declares a distinct manufactured schema and class label; includes all pins, exact code blob identities, array records (path/key/shape/dtype/SHA-256), resolution/grid/profile parameters, parity controls, sector names/sign conventions, tangent semantics, and construction method; and records:

- `production_steps_run = 0`
- `state_advanced = false`
- `production_initializer_called = false`
- `production_rhs_called = false`
- `vendor_step_called = false`
- `scientific_admission = false`
- all admission/closure flags false unless explicitly outside this fixture's scope and still reported false.

The loader verifies the archive/member hashes, every array record, required keys, shape, dtype, finiteness, exact pins, and canonical base-state digest before constructing an evaluator state. It also verifies the stored plus/minus family digests and confirms the input digest is unchanged after evaluation. Missing or extra members, arrays, or required metadata fail closed.

## 7. Evaluator integration boundary

Add a separate manufactured loader and a narrow fixture-class dispatch/adapter. The adapter returns the evaluator's existing internal state/array interface; it does not alter predictor formulas, signs, masks, epsilon handling, zero-step tripwires, or digest checks.

The accepted loader and accepted-state regression remain unchanged and execute only for the accepted class. The manufactured class receives its own schema and structural checks. Reports retain identical signed-term and summary structure, with an explicit `state_class` and `tangent_semantics` so the two tangent meanings cannot be mistaken for one another.

## 8. Construction prohibitions and negative controls

The constructor must not call:

- `V55ProductionKernel.initialize`
- `V55ProductionKernel._rhs` or any production RHS entry point
- `V55ProductionKernel.step`
- vendor PIRK/evolution steps
- any trajectory/campaign/capture path

The existing analytic metric/scalar profile definitions may be reused only after their formulas and provenance are identified. Any fluid profile must be specified analytically and independently; do not carry fluid arrays forward from an initialized seed. The old Ricci-oracle tool is not a Bianchi witness generator.

A parity-breaking negative-control arm, if added, is a separate artifact with its own manifest and digest. It never modifies or enters the primary fixture. No residual fitting, coefficient selection by residual size, EM-on content, production modification, or admission-flag change is permitted.

## 9. Completion gates

The manufactured fixture is ready for observer evaluation only when:

1. A complete initializer-free six-sector analytic base state and explicit one-parameter family are documented.
2. Both epsilon tangent arrays are measured from that family using the pinned discrete momentum constraint and projection definitions, and are honestly labelled as manufactured-family tangents.
3. Pointwise and FV total-stress divergence arrays are independently generated by the pinned helper from the same base/plus/minus triplets.
4. The archive, member, array, state, family, source, helper, vendor, and reviewed-code identities all pass fail-closed validation.
5. The evaluator consumes the class through the adapter with no predictor or zero-step-guard changes.
6. The accepted-state regression remains unchanged and green.
7. A manufactured report is produced with the same signed structure, separate channel summaries, centre cells 0/1 (cell 2 diagnostic), fixed-radius/annular masks, and six-sector raw/accounted summaries.

Only after these gates pass may the manufactured class be evaluated. The t=0 fixture remains a separate subsequent task.


## 10. Construction and loader validation record — 2026-11-05

The contract's constructor, projection, tangent, divergence, and archive/loader gates have now been exercised in bounded GitHub Actions workflows.

- Constructor workflow: [run #42072872069](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42072872069), six tests passed in the pinned Python 3.12.15 environment.
- Five-state six-sector projection workflow: [run #42073309019](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42073309019), all independent/canonical projection cross-checks passed.
- Final manufactured fixture and loader workflow: [run #42074621280](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42074621280), constructor and schema tests passed; archive generation and fail-closed round-trip loader passed.

Final archive:
- Class: `manufactured_regular_family`
- Resolution/domain: `N=60, r_max=40`
- State members: `U(0), U(±2e-5), U(±1e-5)`
- Tangent semantics: `analytic_family_central_difference`; explicitly not a production-RHS tangent.
- Array count: 593
- ZIP SHA-256: `ad8b376881dfeab692af5804870f1e3743cd4fe1fffbc05ca5ec17c9d98415c3`
- Artifact: [run #42074621280](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42074621280)

The runner reports `production_initializer_called=false`, `production_rhs_called=false`, `vendor_step_called=false`, and `production_steps_run=0`. The accepted-state loader, accepted archive, accepted-state regression, production source, and vendor files were not modified.

The raw six-sector divergence reconstruction passed with maximum raw gap approximately `5.56e-18`. The helper's separate time-density-component-accounting flag remains false at both epsilons (the unexplained remainder is approximately `5.56e-18), and is retained as a diagnostic finding rather than promoted to a pass. This is not a physical conservation-closure result.

**Next gate:** add a narrow manufactured-class adapter into the existing full-Bianchi predictor, keeping the accepted loader and accepted regression path unchanged; then run the manufactured report under the same locked predictor and zero-step guards.
