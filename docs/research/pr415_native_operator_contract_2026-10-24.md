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


Decision: if the K-row benefit disappears at a specific stage, that identifies a stage interaction to isolate. If it persists but Aa/Ricci and constraints remain unchanged or worsen, retain the two defect channels as separate. If the copied step cannot be reconstructed exactly, stop and fix the observer rather than interpreting the result.

## Guardrails

- No production equation, initializer, vendor source, or dependency changes.
- No new scientific-admission claims; admission remains false.
- PR #415 remains draft.
- Reuse existing saved profiles/artifacts; do not repeat the completed full evolution or dispatch matrix.
