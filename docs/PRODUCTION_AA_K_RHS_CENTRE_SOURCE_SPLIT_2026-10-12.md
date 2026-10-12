# Production \(A_a\)/\(K\) RHS and matter-current tangent audit — 2026-10-12

**Run:** [38294736014](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294736014)  
**Artifact:** \`production-Aa-K-centre-source-tangent\`, ID \`11768515284\`  
**Artifact SHA-256:** \`2128028ad4fcd80c4d603107d481a78a1add233d18b35eaa3d1d19a43cae2a91\`  
**Diagnostic commit:** \`58e35a0b64d0f17fc053dfaa2aff737365424d69\`  
**Frozen source:** \`6d8b3f2790644591326df57c55908863f994cb0c\`  
**Pinned vendor:** \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`  
**Production main:** unchanged at \`809228c380d80af83858412557c78b6329ffa5c8\`.

## Executive finding

The new signed component audit narrows the centre-series problem further.

1. **The pinned initial profile's first-cell \(\dot A_a\) is dominated by the trace-free Ricci term.** At N=60/80/120, the first-cell Ricci contribution is roughly \(+1.05\times10^{-3}\), while the trace-free lapse-Hessian term is only \(3.6\)–\(4.1\times10^{-5}\).
2. **The uniform fourth-order radial-quadrature initial profile largely removes the first-cell Ricci contribution, but the lapse-Hessian term remains.** At N=80–160 it stays around \(+3.87\)–\(4.15\times10^{-5}\); the total \(\dot A_a\) remains near \(+3.8\)–\(3.9\times10^{-5}\). Because \(r_0=\Delta r/2\), \(\dot A_a/r_0^2\) rises from \(1.53\times10^{-4}\) at N=80 to \(6.21\times10^{-4}\) at N=160. This does not meet the required bounded/convergent centre quotient.
3. **The first-cell \(K\) RHS has a separate, even-scalar centre-series incompatibility.** On the uniform arm, \(\dot K(r_0)\) is approximately \(-1.38\) to \(-1.52\times10^{-3}\), whereas the next cell is about \(+1.0\)–\(1.45\times10^{-4}\). The lapse-Laplacian component changes from negative at cell 0 to about \(+0.96\)–\(1.29\times10^{-3}\) at cell 1, while the matter source is much more continuous. The quotient \(D_r\dot K/r\) at cell 0 increases from \(5.11\times10^{-4}\) (N=40) to \(8.68\times10^{-3}\) (N=160).
4. **This points to centre compatibility of the CMC lapse solution / its central-row treatment as the next discriminator, not to a generic Ricci derivative sign error.** Manufactured smooth-profile tests show the pinned Ricci and trace-free lapse-Hessian operators converge on their fixed-annulus control. Earlier CMC centre-row counterfactuals were mixed and worsened an admitted-inner momentum witness, so this is a diagnostic lead—not permission to patch the row.
5. Signed RHS reconstruction closes to roundoff, and epsilon-direction differences are small. The production source, gauge, coefficients, boundary conditions and centre rules were not edited. All 17 states are frozen; the five uniform-quadrature states are t=0 only.

**Decision:** audit the exact CMC lapse equation/central row and its produced \(\alpha\) series alongside the existing \(A_a/K\) identities. Do not apply a centre projector, tune the row, or resume a long production campaign until a single discrete identity is demonstrably inconsistent.

## Method and coverage

The workflow restored 12 saved baseline/discrete-consistent-initial-metric checkpoints (baseline and candidate at \(t=0\) and saved \(t=2.16\), N=60/80/120), plus five newly reconstructed uniform-fourth-order quadrature initial states (N=40/60/80/120/160). Each is a frozen state: no trajectory was continued by this audit.

On each state the tool:
- decomposed the actual production \(A_a\) RHS into trace-free lapse Hessian, trace-free Ricci, advection/extrinsic, KO-dissipation and anisotropic-stress terms;
- decomposed the actual production \(K\) RHS into lapse Laplacian, advection/extrinsic, KO-dissipation and matter-source terms;
- retained the signed geometric momentum and per-sector scalar/fluid-current terms;
- evaluated the production semi-discrete RHS in two directions (\(\epsilon=2\times10^{-5}\) and \(10^{-5}\)), re-solving the CMC lapse for each perturbed state, then compared directional tangents;
- recorded cells \(i=0,\ldots,5\), fixed radii \(r=2,3,4,5\), and the quotients \(\dot A_a/r^2\) and \(D_r\dot K/r\).

The source checkpoint and vendor are pinned exactly; the uniform initial profile uses the same analytic cubic-Lagrange interval integral described in the quadrature handoff. No \(t=2.16\) evolution is performed for the uniform arm.

## Signed first-cell \(A_a\) budget at \(t=0\)

The relevant identity is the production sum
\[
\dot A_a =
-\left(D_rD_r\alpha-\frac13\Delta\alpha\right)
+\alpha\left(R^r{}_r-\frac13R\right)
+\text{advection/extrinsic}
+\text{KO}
+\text{anisotropic-stress source}.
\]
The reported terms are from the pinned production RHS and vendor operators, with the source term assembled from total stress-energy.

| N | Initial arm | \(\dot A_a\) | TF lapse Hessian | TF Ricci | Matter anisotropic stress | \(\dot A_a/r_0^2\) |
|---:|---|---:|---:|---:|---:|---:|
| 60 | pinned baseline | \(+1.09159e{-3}\) | \(+3.60758e{-5}\) | \(+1.05674e{-3}\) | \(-1.228e{-6}\) | \(+2.45608e{-3}\) |
| 80 | pinned baseline | \(+1.08913e{-3}\) | \(+3.87376e{-5}\) | \(+1.05109e{-3}\) | \(-6.989e{-7}\) | \(+4.35652e{-3}\) |
| 120 | pinned baseline | \(+1.08795e{-3}\) | \(+4.07553e{-5}\) | \(+1.04751e{-3}\) | \(-3.129e{-7}\) | \(+9.79154e{-3}\) |
| 60 | discrete-consistent initial metric | \(+3.30274e{-5}\) | \(+3.60804e{-5}\) | \(-1.825e{-6}\) | \(-1.228e{-6}\) | \(+7.43117e{-5}\) |
| 80 | discrete-consistent initial metric | \(+3.65671e{-5}\) | \(+3.87392e{-5}\) | \(-1.473e{-6}\) | \(-6.990e{-7}\) | \(+1.46269e{-4}\) |
| 120 | discrete-consistent initial metric | \(+3.97270e{-5}\) | \(+4.07556e{-5}\) | \(-7.157e{-7}\) | \(-3.130e{-7}\) | \(+3.57543e{-4}\) |
| 40 | uniform fourth-order quadrature | \(+1.45701e{-4}\) | \(+2.94348e{-5}\) | \(+1.18916e{-4}\) | \(-2.649e{-6}\) | \(+1.45701e{-4}\) |
| 60 | uniform fourth-order quadrature | \(+5.18390e{-5}\) | \(+3.60813e{-5}\) | \(+1.69859e{-5}\) | \(-1.228e{-6}\) | \(+1.16638e{-4}\) |
| 80 | uniform fourth-order quadrature | \(+3.83681e{-5}\) | \(+3.87397e{-5}\) | \(+3.274e{-7}\) | \(-6.990e{-7}\) | \(+1.53472e{-4}\) |
| 120 | uniform fourth-order quadrature | \(+3.72922e{-5}\) | \(+4.07558e{-5}\) | \(-3.151e{-6}\) | \(-3.130e{-7}\) | \(+3.35630e{-4}\) |
| 160 | uniform fourth-order quadrature | \(+3.88433e{-5}\) | \(+4.14856e{-5}\) | \(-2.466e{-6}\) | \(-1.765e{-7}\) | \(+6.21493e{-4}\) |

The algebraic reconstruction is not the issue: maximum \(\dot A_a\)-term reconstruction error over the 17 frozen states is \(2.17\times10^{-19}\); the maximum directional-tangent reconstruction discrepancy is \(1.09\times10^{-14}\). The source and geometry term definitions therefore agree with the implemented RHS to stored precision.

### Interpretation

The uniform quadrature and discrete-consistent initial-metric candidates both greatly reduce the central Ricci term, while the first-cell trace-free lapse-Hessian contribution persists at approximately \(4\times10^{-5}\). The production trace-free lapse-Hessian operator is not thereby proved wrong: on manufactured smooth regular profiles the independent analytic control converges at approximately fourth order on a fixed annulus. The unresolved question is whether the **actual CMC-solved lapse profile and its central discrete row produce a regular even scalar**, so that its trace-free Hessian is \(O(r^2)\) at the first cell.

The baseline/candidate \(t=2.16\) states also show that evolved profiles can regain sizeable first-cell Ricci contributions; the uniform initialization was not evolved to \(t=2.16\) in this run. Thus initial quadrature alone is not an evolution repair.

## Signed first-cell \(K\) RHS and its radial regularity witness

For an even, regular scalar \(K(r)\), the spatial derivative of its RHS should be compatible with the even centre series; \(D_r\dot K/r\) should remain finite and converge rather than grow strongly with refinement.

| N | Arm | \(\dot K_0\) | Lapse term \(-\Delta\alpha\) at \(i=0\) | Matter source at \(i=0\) | \(\dot K_1\) | Lapse term at \(i=1\) | \(D_r\dot K/r\) at \(i=0\) |
|---:|---|---:|---:|---:|---:|---:|---:|
| 40 | uniform quadrature | \(-1.37837e{-3}\) | \(-2.20562e{-4}\) | \(-1.25538e{-3}\) | \(+1.45008e{-4}\) | \(+9.62686e{-4}\) | \(+5.10738e{-4}\) |
| 60 | uniform quadrature | \(-1.45770e{-3}\) | \(-2.70491e{-4}\) | \(-1.28489e{-3}\) | \(+1.29676e{-4}\) | \(+1.14797e{-3}\) | \(+1.19503e{-3}\) |
| 80 | uniform quadrature | \(-1.48816e{-3}\) | \(-2.90475e{-4}\) | \(-1.29541e{-3}\) | \(+1.18402e{-4}\) | \(+1.21724e{-3}\) | \(+2.14878e{-3}\) |
| 120 | uniform quadrature | \(-1.51086e{-3}\) | \(-3.05634e{-4}\) | \(-1.30299e{-3}\) | \(+1.07738e{-4}\) | \(+1.26777e{-3}\) | \(+4.86945e{-3}\) |
| 160 | uniform quadrature | \(-1.51900e{-3}\) | \(-3.11122e{-4}\) | \(-1.30565e{-3}\) | \(+1.03413e{-4}\) | \(+1.28561e{-3}\) | \(+8.67666e{-3}\) |

The centre-to-first-neighbour jump in the lapse-Laplacian term is large and does not shrink under refinement; the matter source changes much less between those cells. This makes the CMC lapse solve/central-row compatibility a high-priority discriminator. It is not yet a proved centre-row defect: a proper audit must compare the exact solved lapse to the intended discrete CMC equation and regular radial series, identify how the centre row is constructed, and distinguish equation residual from derivative-evaluation effects.

The earlier centre-row/lapse counterfactuals on separate branches reduced some first-cell \(K\)-RHS mismatch measures but worsened an admitted-inner momentum-constraint witness after one normal step by approximately \(3\times10^{-6}\). That negative result must remain part of the record; it forbids treating a local \(K\)-row improvement as an adequate repair.

## Tangent stability and closures

Across all 17 cases:
- Base \(\dot A_a\) reconstruction maximum error: \(2.17\times10^{-19}\).
- Base \(\dot K\) reconstruction maximum error: \(2.17\times10^{-19}\).
- Momentum geometric reconstruction: exactly 0 in stored arrays.
- Per-sector momentum source reconstruction maximum error: \(4.07\times10^{-20}\).
- Directional tangent reconstruction maximum error: \(1.09\times10^{-14}\) for \(A_a\) and \(K\).
- Largest difference between the two epsilon tangents is \(2.84\times10^{-10}\) for the full \(\dot A_a\) term, \(8.54\times10^{-10}\) for the full \(\dot K\) term, and \(1.10\times10^{-12}\) for the momentum residual tangent.

These are algebraic and directional-consistency checks, not continuum conservation proofs. The production scalar/fluid current was kept as signed, separately reported sectors; this report does not claim that their arithmetic reconstruction proves covariant stress-energy conservation.

## Next bounded audit

1. **Read-only CMC lapse centre-row audit:** inspect the exact \`_solve_lapse\` discrete equation and boundary/centre row for the pinned source. On frozen baseline, discrete-consistent, and uniform-quadrature profiles, report \(\alpha_i\), even-series residuals, centre-row equation residual, and \(-\Delta\alpha\) at \(i=0,1,2\), fixed \(r=2,3,4,5\), and across N=40–160. Compare against independent regular even profiles. No row is changed in this test.
2. **Coupled centre identities:** reconstruct the expected regular-series limits of trace-free \(D_iD_j\alpha\), trace-free Ricci, \(\dot A_a\), \(\dot K\), and \(D_r\dot K/r\) with consistent metric/conformal conventions. Show which identity first fails.
3. **Radiation admissibility remains a separate blocker:** matched one-step N=240 and N=320 cases fail at \(i=0\) for all initializer arms (\(|S|/E\approx1.0203\) and \(1.3603\), respectively). Do not clip primitives or loosen the cone.
4. Only after a specific identity fails should a one-change numerical counterfactual be proposed. Recheck H, M, \(C_\Lambda\), determinant, radiation admissibility and Misner–Sharp balance. No long campaign, especially no \(t=160\), until the admission gates pass.

## Guardrails

No production source, physics equation, scalar/matter source, coefficient, gauge, boundary condition, centre projection, damping, floor, clamp or radiation admissibility rule was changed. \`main\` remains untouched. The quadrature arms and all CMC/centre counterfactuals remain diagnostic. Workflow success means that this frozen-state audit executed and reconstructed its signed terms—not that the production machine has passed scientific admission.
