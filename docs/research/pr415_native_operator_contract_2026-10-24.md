# PR #415 research note: native operator contract exposed by the CMC centre-row counterfactual

Date: 2026-10-24
Status: diagnostic design only; not a production proposal.

## Finding

The current-source saved-state counterfactual at N=60/120/240 shows that the CMC centre row and the production geometry RHS do not use the same discrete operator contract.

The production-local `engine/cmc_gauge.py::solve_cmc_lapse` currently:
- computes first derivatives through `_d1`;
- assembles a second-order centred three-point interior lapse operator;
- imposes the centre row `alpha[0] = alpha[1]` as its discrete representation of `alpha'(0)=0`;
- normalizes the outer lapse to one.

The production geometry RHS obtains its L2 terms through the pinned vendor geometry operator. The saved-state counterfactual instead constructs the lapse row with the same parity-aware fourth-order D1/D2 stencil used by the production K L2 operator, while holding the existing CMC target fixed. This copied solve reduced the centre K_L2 centre-minus-even-fit mismatch by 93.2%, 97.8%, and 97.4% at N=60/120/240. It does not prove the replacement is a correct global gauge solve or an evolved fix.

## Native-kernel blueprint hypothesis

The result supports investigating a native, engine-owned *operator contract*, rather than adding more local patches around borrowed modules. It does not yet establish that the vendor kernel as a whole is inadequate. The concrete incompatibility currently demonstrated is between the local CMC gauge discretization and the production geometry RHS stencil at the origin.

A future native operator interface should make these items explicit and shared:
1. Grid placement and radial support (cell-centred positive radii; no implicit point at r=0).
2. Parity-aware D1 and D2, including the origin closure and outer-boundary rows.
3. One documented spherical/reference-metric geometry operator convention for the lapse PDE and the K/Aa RHS.
4. Stage-local geometry, scalar, matter and electromagnetic projections, all evaluated from the same named PIRK stage.
5. CMC target projection and its discrete volume weights, with the exact K RHS used to define the target.
6. Explicit update-stage ordering and signed operand capture so the accepted update can be reconstructed independently.
7. Fail-closed provenance: source/vendor/dependency pins and saved-state digests, with diagnostic counters kept separate from scientific admission.

This is a specification emerging from observed interfaces, not a request to replace the vendor kernel or rewrite production now. Keep the vendor provenance and production equations intact until the cross-stage experiment establishes which interfaces actually fail.

## Next falsifiable experiment

Artifact availability correction (2026-10-24): the verified baseline/candidate bundles include complete coupled state checkpoints only at t=0.6. The histories near the predictor crossover (N=120: t≈0.22; N=240: t≈0.10) contain stage/profile diagnostics but not the full scalar and conserved-matter fields needed to restart a faithful coupled step. Therefore the first executable one-step comparison is deliberately limited to the complete t=0.6 checkpoints. We will not synthesize missing state from plotted profiles or rerun the long campaign merely to recreate these snapshots. A full checkpoint at a crossover stage remains a separate evidence gap.

Using the existing baseline and half-cell complete t=0.6 checkpoints at N=120 and N=240, run a *copied-state one-step* comparison:

- Arm A: reproduce the current production step exactly.
- Arm B: in a copy only, replace the lapse row with the parity-aware D1/D2 row at every CMC refresh within that one PIRK step; leave the CMC target projection and all other equations unchanged.
- Capture the lapse, signed K_L2/K_L3, Aa lapse-Hessian and trace-free-Ricci terms, and H/M/connection residuals at the start, explicit predictor, primary predictor, final implicit operand, and accepted output.
- Reconstruct the actual Aa/K updates from captured operands and require closure to the same raw-operand-scaled roundoff envelope used by the existing audit.
- Verify bitwise that Arm A's saved input and all original artifacts remain unchanged. Do not chain counterfactual steps; this is a local sensitivity test, not a new trajectory.


## First copied one-step result at t=0.6

The stage observer and update reconstruction passed on both lineages at N=120 and N=240 (workflow run 40136649511). All four arms per resolution reported that the saved input remained unchanged; each captured four CMC refreshes, three geometry RHS stages, and raw-operand-scaled A_a/K update closure.

Maximum absolute stage-local residual |K_rhs - CMC target|, across baseline and origin-half-cell lineages:

| Stage | N=120 production | N=120 operator-compatible | N=240 production | N=240 operator-compatible |
|---|---:|---:|---:|---:|
| Start-lapse refresh | 2.374e-4 | 1.78e-15–2.48e-15 | 2.406e-4 | 1.09e-14–1.40e-14 |
| Explicit predictor | 3.382e-4 | 1.55e-15–1.69e-15 | 2.913e-4 | 6.89e-15–2.01e-14 |
| Primary predictor | 3.382e-4 | 1.95e-15–2.96e-15 | 2.913e-4 | 1.04e-14–1.37e-14 |
| Accepted-output refresh | 3.380e-4 | 1.94e-15–2.75e-15 | 2.913e-4 | 9.33e-15–1.24e-14 |

The improvement is consistent across both initial-data lineages and all captured stages. The copied solve shares the derivative operators and target definition with the K RHS, so this is a **discrete operator-compatibility result**, not an independent physical validation of the CMC target.

The other source channels did not close as a side effect. On the baseline lineage at the start-lapse refresh, the cell-0 A_a lapse-Hessian term changed from -5.904e-6 to +4.824e-7 at N=120, and from -6.175e-6 to +2.511e-7 at N=240. The trace-free-Ricci term remained approximately -1.2001e-3 at N=120 and +1.0312e-3 at N=240. At the accepted output, H/M/connection maxima changed only slightly (roughly 1e-7, 1e-6, and 1e-8 respectively), while remained materially nonzero. Thus the lapse-row change fixes the sampled K-target compatibility mismatch but does not remove the separate A_a/Ricci or constraint residuals.

Artifact records from run 40136649511:
- N=120: artifact 12577505582; uploaded ZIP SHA-256 `5170ac602b21fe4c2f0aaf78d905a27839b9445bde9dacbed086473e4e205a67`.
- N=240: artifact 12577005887; uploaded ZIP SHA-256 `269a97a691a56ce47a1b8a87be49064f231bfe1ea92190aceac9771582836c89`.


## Exact t=0.6 evolved-state Ricci channel decomposition

The next observer reuses the existing pinned methods in `audit_evolved_centre_ricci_stage_budget.py` and `audit_current_main_centre_proper_radial_ricci.py`, but applies them to the **same four saved-state copied steps above** rather than initializing a new trajectory. The report now records, at the start-lapse refresh, explicit predictor, primary predictor, final mixed-state L2 operand, and accepted output:

- pinned-vendor TF Ricci versus the independent proper-radial TF Ricci formula;
- the derivative-only `D4(Lambda)` chain-rule counterfactual gap (the algebraic `Lambda D1(a)` term is deliberately retained);
- separate signed `D4(Lambda)`, algebraic-`Lambda`, remaining metric-bracket, and A_a lapse-Hessian contributions;
- full-array reconstruction closure for both the TF Ricci identity and the A_a L2 source split.

The identity closures are at floating-point roundoff: maximum TF-Ricci channel closure is below `1.7e-19`, and maximum A_a L2 reconstruction closure is below `5.5e-19`. This proves the ledger reconstructs the implemented source; it does not prove either Ricci route is the correct continuum equation.

### Accepted-output centre discriminator

Signed values below are at cell 0. “Raw gap” is pinned-vendor TF Ricci minus the independently reconstructed proper-radial TF Ricci; “chain-rule gap” applies only the saved-state derivative representation counterfactual.

| Resolution / saved lineage | Raw gap at cell 0 | Chain-rule gap at cell 0 | Effect of counterfactual |
|---|---:|---:|---|
| N=120 baseline | -1.415477e-3 | +4.377803e-4 | magnitude reduced 69.1% (3.23x) |
| N=120 origin-half-cell candidate | +2.668611e-4 | +1.682555e-4 | magnitude reduced 36.9% (1.59x) |
| N=240 baseline | +3.598385e-3 | +1.149145e-3 | magnitude reduced 68.1% (3.13x) |
| N=240 origin-half-cell candidate | +7.151231e-4 | +9.893254e-4 | magnitude **increased 38.3%** |

The response is emphatically not a universal “replace D4 by the chain rule” fix. It reduces the baseline lineage's cell-0 mismatch, does much less for the N=120 half-cell candidate, and worsens cell 0 for the N=240 half-cell candidate. Yet after the counterfactual the fixed physical annulus maximum is only about `1.4e-9` at N=120 and `6.9e-11` at N=240. Cell 6 is likewise about `1.2e-9` / `2.4e-9`; the material residual remains in the first few centre cells. The centre cells do not track the smooth-annulus response.

The signed pattern across cells 0/1/2 also changes with resolution and initialization lineage. For example, after the counterfactual the N=120 baseline residuals are approximately `(+4.38e-4, -7.52e-5, +1.08e-5)`; N=240 baseline is `(+1.149e-3, -2.097e-4, +3.94e-5)`. For the N=240 half-cell lineage the corresponding values are `(+9.893e-4, -1.695e-4, +2.821e-5)`. This supports retaining the centre as a separate regularity/operator-compatibility problem rather than fitting a single global correction.

### What the A_a source budget adds

The A_a L2 channel split remains closed to roundoff, but its signed components are strongly state- and resolution-dependent. On the accepted N=120 baseline, the cell-0 A_a contribution from the `D4(Lambda)` TF-Ricci channel is about `+1.013e-3`, the remaining metric bracket is about `-2.185e-3`, and the lapse-Hessian contribution is about `-8.51e-6`. On the N=240 baseline, the D4 channel is about `-7.03e-5` while the remaining metric bracket is `+1.061e-3`. The algebraic-`Lambda` term is only O(1e-7) or smaller in these examples.

The operator-compatible lapse row changes the accepted cell-0 lapse-Hessian term to roughly `+5.8e-7` at N=120 and `+2.6e-7` at N=240, but leaves the large signed TF-Ricci components and Ricci gap largely unchanged. This is direct evidence that the CMC/K row mismatch and the A_a/Ricci centre residual are distinct channels. The half-cell lineage also changes the A_a component balance substantially, but neither counterfactual yields a uniformly smaller evolved Ricci gap at both resolutions.

### Provenance and reproducibility

The channel decomposition is from workflow run [40136928892](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40136928892), reusing the verified complete t=0.6 artifacts. Both N=120/N=240 jobs passed, including stage capture and accepted-update closure.

- N=120 artifact [12577046128](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40136928892/artifacts/12577046128); ZIP SHA-256: `44139de32bdce0c37bd4bf990d88958fcf0554acef4cdc757b95d30cbcd3caee`.
- N=240 artifact [12577411251](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40136928892/artifacts/12577411251); ZIP SHA-256: `36e5a7b7d8a4905da629a3f34aa89019030a3d243f81490b8711d6b5f9730b98`.

The input checkpoint SHA-256 values are reproduced in each report; baseline/candidate checkpoints were loaded from the archived artifact bundles, were not re-evolved to t=0.6, and each A/B arm's source input remained unchanged. The report leaves `scientific_admission=false`.

**Interpretation:** current evolved-state data do not support treating the derivative representation counterfactual as a universal cure, despite its strong initialized-slice effect in the archive. The defensible next step is to keep the `D4(Lambda)`/metric-bracket mismatch, the origin half-cell/regularity behavior, the CMC lapse/K compatibility, and the coupled H/M/connection propagation budgets independently observable under one shared native operator contract. Do not promote any one counterfactual into production on this evidence.


### Separate the evolved-connection and metric-stencil responses

The latest saved-state observer further splits the change to the TF-Ricci derivative channel into two signed responses:

`connection response = (2/3) X^2 [D4(Lambda_evolved) - D4(Lambda_metric)]`

`metric-stencil/chain response = (2/3) X^2 [D4(Lambda_metric) - d_chain(Lambda_metric)]`

Their sum is the total `D4(Lambda_evolved) - d_chain(Lambda_metric)` response. Here `Lambda_metric` is the established discrete metric reconstruction, while `d_chain` is the existing analytic chain-rule derivative built from the same D1/D2 metric operands. This is a decomposition, not a claim that either operator should replace the evolved connection in production.

Accepted-output cell-0 signed TF-Ricci responses:

| Resolution / saved lineage | Evolved-connection response | Metric-stencil/chain response | Sum (evolved D4 minus chain rule) |
|---|---:|---:|---:|
| N=120 baseline | +2.41318e-4 | -2.09458e-3 | -1.85326e-3 |
| N=120 origin-half-cell candidate | -1.84870e-5 | +1.17093e-4 | +9.86057e-5 |
| N=240 baseline | -1.53235e-3 | +3.98159e-3 | +2.44924e-3 |
| N=240 origin-half-cell candidate | -6.35024e-4 | +3.60822e-4 | -2.74202e-4 |

The two response pieces sum to the total within at most `8.7e-19`; the independent identity `raw vendor-minus-direct gap - chain-rule gap = total D4/chain response` closes within `4.4e-19`. The algebra is therefore accounted for; the varying signs and amplitudes are features of the captured discrete states, not ledger non-closure.

The decomposition does not identify one universally dominant term. The metric-stencil/chain response is larger in magnitude on both N=120 lineages and on the N=240 baseline, but the evolved-connection response is larger on the N=240 half-cell lineage. At cell 1 the terms can oppose each other differently again. This is consistent with the prior observed centre-localized derivative representation issue, while also demonstrating that the evolved connection response cannot be assumed negligible on every accepted state.

**Next diagnostic gate:** preserve both signed channels in the stage contract and compare them with `Lambda_evolved - Lambda_metric`, H/M/connection propagation, and A_a source terms on the same named PIRK operands. Any candidate native kernel should expose these operands and their identities rather than silently setting `Lambda_evolved = Lambda_metric` or substituting the chain rule. Only after an independently derived propagation identity closes on the exact accepted state should we infer which source term or closure needs a physical change.


## Accepted-step Lambda update and connection-residual budget

The next experiment asks whether a one-step change in the connection residual
`C = Lambda_evolved - Lambda_metric(a,b)`
is carried by the evolved-`Lambda` source, the metric update, or both. This uses the same four copied t=0.6 checkpoint arms at N=120/240; it does not chain a counterfactual step.

The observer independently reconstructs the actual accepted Lambda update on the interior `r <= 10`:
[
Lambda_{n+1}=Lambda_n+rac{Delta t}{2}
left(L_{2,mathrm{start}}+L_{2,mathrm{final}}+
L_{3,mathrm{start}}+L_{3,mathrm{primary}}ight).
]
The final L2 source is reevaluated on the exact final mixed state, with accepted A_a and K operands, matching the production split. The outer characteristic boundary is excluded from this arithmetic comparison because it separately updates the finite-radius boundary mode.

**All four arms pass.** The maximum absolute Lambda-update reconstruction residual over `r <= 10` is between `1.38e-19` and `2.16e-19`, against raw-operand-scaled bounds of approximately `7.73e-16`–`7.78e-16`. The independent identity
`delta(C) = delta(Lambda_evolved) - delta(Lambda_metric)`
closes exactly or within `4.1e-20`. All original saved inputs remain unchanged. This is arithmetic closure of the implemented update, not proof that its source terms are physically correct.

### Signed cell-0 movement across the accepted step

| Resolution / lineage | C at input | delta Lambda evolved | delta Lambda metric from a,b | delta C | C at output |
|---|---:|---:|---:|---:|---:|
| N=120 baseline, production | +7.11697e-4 | -1.44802e-5 | +3.05539e-5 | -4.50342e-5 | +6.66663e-4 |
| N=120 baseline, operator-compatible lapse | +7.11697e-4 | -1.45586e-5 | +3.04940e-5 | -4.50526e-5 | +6.66645e-4 |
| N=120 half-cell, production | -9.77600e-5 | -1.44771e-5 | -1.91146e-5 | +4.63751e-6 | -9.31225e-5 |
| N=120 half-cell, operator-compatible lapse | -9.77600e-5 | -1.45557e-5 | -1.91724e-5 | +4.61671e-6 | -9.31433e-5 |
| N=240 baseline, production | -8.58895e-4 | -1.50865e-5 | -3.77644e-6 | -1.13100e-5 | -8.70205e-4 |
| N=240 baseline, operator-compatible lapse | -8.58895e-4 | -1.51227e-5 | -3.80326e-6 | -1.13194e-5 | -8.70214e-4 |
| N=240 half-cell, production | -3.60770e-4 | -1.50859e-5 | +1.83253e-5 | -3.34111e-5 | -3.94181e-4 |
| N=240 half-cell, operator-compatible lapse | -3.60770e-4 | -1.51221e-5 | +1.82982e-5 | -3.34203e-5 | -3.94190e-4 |

The `Lambda` source movement at cell 0 is consistently negative (about `-1.45e-5` at N=120 and `-1.51e-5` at N=240), while the independently calculated metric-derived-`Lambda` movement changes sign across lineage/resolution. Thus neither source lane alone explains the residual change:

- On the N=120 baseline, the metric-derived connection moves positive by about `3.06e-5`, opposing the negative evolved-`Lambda` change; the residual decreases by `4.50e-5`.
- On the N=120 half-cell lineage, both changes are negative, with the metric-derived change more negative, so the residual moves toward zero by about `4.64e-6`.
- On the N=240 baseline, the evolved-`Lambda` change dominates, producing a `-1.13e-5` residual change.
- On the N=240 half-cell lineage, the metric-derived connection moves positive by `1.83e-5` against the negative evolved-`Lambda` change; their difference pushes the residual more negative by `3.34e-5`.

The copied operator-compatible lapse changes these `delta C` values only marginally (at most about `9.4e-9` at cell 0 in these comparisons); it repairs the K-RHS/CMC-target row consistency, but is not the cause or cure of the Lambda-connection residual. Cell 1/2/6 source and residual increments are retained in the signed per-arm artifact JSON.

This evidence points to a stage-coupled response between the metric update and the independently evolved connection variable, with a balance that depends on the saved state and resolution. It argues for an operator interface that captures `Lambda_evolved`, metric-derived `Lambda`, their separate source/update operands, and the accepted `delta C` in the same ledger. Do not enforce the constraint by setting the two Lambda profiles equal; expose the propagation defect and derive its full matter-coupled identity.

At cell 0, splitting the Lambda update into its two integration lanes shows that L2 drives the negative evolved-`Lambda` movement, while L3 offsets it positively:

| Resolution / lineage | delta Lambda from L2 | delta Lambda from L3 | net delta Lambda |
|---|---:|---:|---:|
| N=120 baseline, production | -1.51561e-5 | +6.75858e-7 | -1.44802e-5 |
| N=120 baseline, operator-compatible lapse | -1.52345e-5 | +6.75842e-7 | -1.45586e-5 |
| N=120 half-cell, production | -1.51559e-5 | +6.78784e-7 | -1.44771e-5 |
| N=120 half-cell, operator-compatible lapse | -1.52344e-5 | +6.78766e-7 | -1.45557e-5 |
| N=240 baseline, production | -1.52563e-5 | +1.69815e-7 | -1.50865e-5 |
| N=240 baseline, operator-compatible lapse | -1.52925e-5 | +1.69815e-7 | -1.51227e-5 |
| N=240 half-cell, production | -1.52561e-5 | +1.70265e-7 | -1.50859e-5 |
| N=240 half-cell, operator-compatible lapse | -1.52923e-5 | +1.70265e-7 | -1.51221e-5 |

Here each contribution is formed from the actual accepted trapezoidal split, `0.5*dt*(L2_start+L2_final)` and `0.5*dt*(L3_start+L3_primary)`. L3 offsets the L2 decrease by approximately 4.5% at N=120 and 1.1% at N=240. The lapse-row counterfactual slightly changes the L2 contribution but not the L3 term meaningfully; the resulting difference in `delta Lambda` is still small compared with the lineage-dependent metric-derived Lambda movement.


### Artifact provenance

Workflow [40137296859](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40137296859) completed successfully for N=120 and N=240 from the saved t=0.6 artifacts.

- N=120 artifact [12577357052](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40137296859/artifacts/12577357052), ZIP SHA-256: `d1f52473efd61b90459b221068070954b7dd6ff586ede0b5021024a04a5f5e9b`.
- N=240 artifact [12577886263](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40137296859/artifacts/12577886263), ZIP SHA-256: `13d9c493127d7fbc14adb7479def9f981b8021fb2768aad9dbbea063f7e92b5f`.

The report has `scientific_admission=false`; no production, initializer, vendor, or dependency path changed.


## Same-accepted-state momentum tangent plus independent total-T divergence

**Completed:** workflow [40138324805](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40138324805) reuses the checksummed baseline and origin-half-cell candidate checkpoints at N=120 and N=240, all at t=0.6. No campaign states were re-evolved. The observer uses the exact current `kernel._geometry_stage_terms` interface (including the native explicit a/b/X RHS), then evaluates the two-sided momentum tangent at epsilon = 2e-5 and 1e-5. It joins the existing independent total stress-energy divergence and six-sector decomposition on those **same directional state pairs**.

### 1. The Lambda-connection rate is strongly coupled to M, but not yet a closed full propagation identity

For the residual `C = Lambda_evolved - Lambda_metric(a,b)`, the candidate relation `Cdot = (2 alpha/a) M` is much closer than its opposite sign `Cdot = -(2 alpha/a) M` on all four accepted states:

| Resolution / lineage | M at cell 0 | Max abs [Cdot - (2 alpha/a) M], full grid | Max abs opposite-sign residual |
|---|---:|---:|---:|
| N=120 baseline | -1.031329e-3 | 1.081e-7 | 4.121e-3 |
| N=120 origin-half-cell candidate | +1.107020e-4 | 1.006e-8 | 4.423e-4 |
| N=240 baseline | -6.508156e-4 | 3.697e-8 | 2.600e-3 |
| N=240 origin-half-cell candidate | -1.639339e-3 | 4.166e-8 | 6.549e-3 |

The same-sign relation is consistent at both directional epsilons; geometric five-term tangent closure and six-sector current projection closure pass. However, the same-sign residual remains finite and centre-dominated (cell 0/1/2), so this is a **strong leading-order relation, not exact full matter-coupled propagation closure**. Keep the remainder explicit; do not assume it is roundoff or drop it.

### 2. The independent total-T projection agrees with the canonical projection

For both epsilons and both lineages at N=120/240, the independently reconstructed `rho, pr, pt, j` projections pass their comparison with the canonical total stress-energy projection. The total-T divergence itself is nonzero and centre-localized in the pointwise evaluator. Representative cell-0 values are:

| Resolution / lineage | \\(\\nabla_\\mu T^\\mu{}_{t}\\) (pointwise) | \\(\\nabla_\\mu T^\\mu{}_{r}\\) (pointwise) | \\(\\nabla_\\mu T^\\mu{}_{t}\\) (finite-volume) | \\(\\nabla_\\mu T^\\mu{}_{r}\\) (finite-volume) |
|---|---:|---:|---:|---:|
| N=120 baseline | +5.41630e-7 | +1.55747e-7 | +3.83121e-9 | -6.29328e-9 |
| N=120 origin-half-cell candidate | +5.43384e-7 | +1.67217e-7 | +3.97157e-9 | +7.72126e-11 |
| N=240 baseline | +1.46405e-7 | +5.32677e-8 | +8.20604e-10 | -9.43553e-9 |
| N=240 origin-half-cell candidate | +1.46403e-7 | +5.49858e-8 | +8.34340e-10 | -7.79914e-9 |

The pointwise fixed annulus (5<=r<=10) is smaller than the first-cell response: the largest pointwise time-component annulus maxima are about 7.52e-10 at N=120 and 4.88e-11 at N=240. The finite-volume and pointwise forms differ materially in the first few cells; neither is silently selected as the authoritative propagation operator here. Both signed forms remain in the artifact.

These are directional derivatives of the independently projected total tensor along the current semi-discrete RHS, not a zero-divergence claim. The fact that the two methods yield different centre amplitudes is a reason to keep both discretization contracts visible.

### 3. A previous sector-sum discrepancy was observer accounting, not an unexplained source

Summing the six per-sector directional divergences initially differs from the independently calculated total divergence by about `2.8e-16`–`7.6e-16`. The report separately measures the difference between the *sum of per-sector time-density rates* and the total-projection time-density rate. Propagating that difference through the same divergence operator predicts the observed sector-sum gap.

After this explicit time-density accounting, the maximum unexplained remainder is only about `2.8e-20`–`6.5e-20` across the N=120/N=240 lineages and epsilons, and all four divergence channels report `component_corrected_passed=true`. Thus this particular sector-additivity gap is explained by the nonlinear finite-difference time-density projection; it should not be mislabeled as evidence of physical nonconservation. Conversely, the corrected **sector accounting** does not make the total divergence itself vanish—the residuals above remain.

The total-T observer is six-sector (S, D, COSMOS phi, dark matter, baryons, radiation) and zero-shift, matching these saved states; it does not cover optional electromagnetic/charged-dust branches and does not admit the conditional dark-matter exchange law.

### Provenance / exact artifacts

The run passed both N=120 and N=240 jobs, with input digests unchanged, exact current-engine source pins checked, the historical provenance reference for the total-T helper checked independently, and `scientific_admission=false`.

- N=120 artifact [12577267713](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40138324805/artifacts/12577267713), ZIP SHA-256: `80e4c5d94929458487a8dd05db4f118a9a4bd61acfd7ca3922b4b27f3c31049e`; embedded `report.json` SHA-256: `432206b08347c8110031765b4e2846afbb77eadec79a0206d62b1bc39f46143d`.
- N=240 artifact [12577152799](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40138324805/artifacts/12577152799), ZIP SHA-256: `8dfe4f894727ef856ee57d580b99619aa9165bb4b61aa4b58721cd0bb45a3543`; embedded `report.json` SHA-256: `76057a6ec5330c3418cf1ebe0a48314fe9a7e73b4542b912ccb07f7a56c2d8b8`.

**Disposition:** the `Cdot ≈ (2 alpha/a) M` relation is corroborated but not exact at the centre; independent total-T projection is cross-checked, and its sector-sum discrepancy is accounted for, but the total-T divergence remains nonzero. The next gate remains a frozen, convention-pinned full matter-coupled momentum-constraint propagation RHS, including the normalized Hamiltonian-gradient term, reference-metric connection response, and complete stress-energy exchange/conservation residuals. No production change is justified by this diagnostic.


Decision: if the K-row benefit disappears at a specific stage, that identifies a stage interaction to isolate. If it persists but Aa/Ricci and constraints remain unchanged or worsen, retain the two defect channels as separate. If the copied step cannot be reconstructed exactly, stop and fix the observer rather than interpreting the result.

## Guardrails

- No production equation, initializer, vendor source, or dependency changes.
- No new scientific-admission claims; admission remains false.
- PR #415 remains draft.
- Reuse existing saved profiles/artifacts; do not repeat the completed full evolution or dispatch matrix.
