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

Using the existing baseline and half-cell saved histories at N=120 and N=240, run a *copied-state one-step* comparison at a small set of captured times bracketing the known predictor crossover (N=120: t≈0.22; N=240: t≈0.10), plus the common t=0.6 checkpoint:

- Arm A: reproduce the current production step exactly.
- Arm B: in a copy only, replace the lapse row with the parity-aware D1/D2 row at every CMC refresh within that one PIRK step; leave the CMC target projection and all other equations unchanged.
- Capture the lapse, signed K_L2/K_L3, Aa lapse-Hessian and trace-free-Ricci terms, and H/M/connection residuals at the start, explicit predictor, primary predictor, final implicit operand, and accepted output.
- Reconstruct the actual Aa/K updates from captured operands and require closure to the same raw-operand-scaled roundoff envelope used by the existing audit.
- Verify bitwise that Arm A's saved input and all original artifacts remain unchanged. Do not chain counterfactual steps; this is a local sensitivity test, not a new trajectory.

Decision: if the K-row benefit disappears at a specific stage, that identifies a stage interaction to isolate. If it persists but Aa/Ricci and constraints remain unchanged or worsen, retain the two defect channels as separate. If the copied step cannot be reconstructed exactly, stop and fix the observer rather than interpreting the result.

## Guardrails

- No production equation, initializer, vendor source, or dependency changes.
- No new scientific-admission claims; admission remains false.
- PR #415 remains draft.
- Reuse existing saved profiles/artifacts; do not repeat the completed full evolution or dispatch matrix.
