# Issue #4 coupled geometric-and-matter tangent control

**Date:** 2026-11-05  
**Branch:** `research/issue4-coupled-geometric-matter-control-2026-11-05`  
**Status:** source-equation map and derivation preflight only; no coupled fixture has been constructed or evaluated  
**Scope:** research-only, instantaneous, zero-step; no production RHS, initializer, vendor step, evolution, predictor modification, or admission change.

## Purpose and immutable inputs

This is a new control class. The arbitrary manufactured family, equation-consistent frozen-flat matter direction, accepted-state archive, their reports, and the three-class comparison remain immutable. This branch must not regenerate or rewrite any of those inputs.

The control target is a single analytic base state and a first-order direction that includes both (i) the six-sector matter tangent and (ii) the geometry tangent dictated by the pinned CMC metric/primary/Lambda equations. It is not a claim of a complete Einstein-matter solution.

Locked provenance remains:
- diagnostic source: `4aec88d0e3ff028e0a65a8ac36b30948625df792`
- helper baseline: `7f47669a1c0edfafef9fc5f27b2db80840d3f035`
- vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- tangent helper blob: `3f00a4e0a2063b86d2e007c5a825f9edfff61a4c`
- total-stress helper blob: `c10f9e8b9e654d5195fd55c493c3fe06ba5c7082`

## Pinned instantaneous geometry equation map

The CMC branch is not the vacuum moving-puncture evolution. The source path is:
`V55ProductionKernel._geometry_stage_terms` ->
`native_cmc_explicit_metric_rhs` for (a,b,X) ->
`v55_pirk_adapter.geometry_stage_terms` for the primary and Lambda blocks.
The CMC gauge sets (eta=B=0) on the slice. The lapse is a positive, elliptically determined CMC field, not an independently evolved 1+log variable.

With (eta=0), (A_b=-A_a/2), and the pinned fourth-order parity-aware operators, the explicit metric rates reduce to
[
dot a=-2alpha a A_a,qquad
dot b=+alpha b A_a,qquad
dot X=rac{alpha X K}{3}.
]
The primary rates are the sum of the pinned L2 curvature/lapse-Hessian block and the matter-coupled L3 block:
[
dot A_a=
-left(D_rD_ralpha-rac13Deltaalphaight)
+alphaleft(R^r{}_r-rac13 Right)
+alpha K A_a
-rac{16pi}{3}alpha(p_r-p_t),
]
[
dot K=-Deltaalpha+
alphaleft(rac32 A_a^2+rac13K^2ight)
+4pialpha(ho+p_r+2p_t).
]
The connection rate is the sum of the pinned (Lambda_{L2}) derivative/singular block, the (Lambda_{L3}) advection/connection block, and the matter momentum source:
[
dotLambda=
Lambda_{L2}
+left(2alpha A_aLambdaight)
-rac{8pilambda_malpha j}{a},
qquad lambda_m=2.
]
Here (Lambda_{L2}) retains the pinned (D_r^2eta), shift-divergence, (A_a)/lapse-gradient, and geometry-only momentum-constraint terms. In the CMC slice (eta=0), but its parity-aware discrete derivatives must still be evaluated under the same operator contract; do not replace the finite-difference expression with a continuum simplification without an explicit identity.

These expressions are an equation map, not yet an executable implementation. Every curvature, lapse-Hessian, connection, and source operand must be reconstructed independently from the declared base arrays and pinned operators. The research constructor must not import or call `V55ProductionKernel`, `_rhs`, `_geometry_stage_terms`, `solve_cmc_lapse`, vendor `full_rhs`, or any stepper.

## CMC lapse is a separate derivation gate

The CMC lapse is not a free profile. The pinned gauge condition determines it from
[
(D^2-mathcal Q)alpha=-dot K_{m CMC},
]
with the proper-volume projected (K)-rate defining (dot K_{m CMC}), and
[
mathcal Q=rac32 A_a^2+rac13K^2+
4pi(ho+p_r+2p_t).
]
The finite-grid operator, centre regularity row, outer normalization, and proper-volume weights are part of the locked convention. A coupled direction therefore also requires the directional derivative of this elliptic CMC solve, including the variations of its coefficients, source, and proper-volume projection. Simply holding (alpha) fixed or differentiating only the pointwise equation is not sufficient.

The independent implementation must solve the pinned discrete linear system and its differentiated system from reconstructed arrays, without importing the production CMC solver. It must retain a separate residual witness for the original and differentiated CMC equations. If this derivative cannot be independently derived and validated, the coupled fixture is blocked; do not substitute an arbitrary lapse tangent.

## Matter coupling and source conventions

The matter direction must be reconstructed with the same metric-dependent six-sector projection used by the locked observer. The scalar and fluid instantaneous rates must be consistent with their pinned scalar and Valencia equations on the same base metric. The dark-matter exchange covector retains the archived sign:
[
Q_{DM,r}=+eta_{DM}ho_{DM},partial_rphi,
]
with the scalar exchange opposite. The total source operands (ho,p_r,p_t,j) used by the geometry rates must be generated independently from the six sectors and checked against the existing projection contract. Do not infer them from residuals or tune them.

## Required implementation gates

1. **Independent equation transcription:** implement pure research-only functions for the CMC lapse solve and its directional derivative, the (a,b,X,A_a,K,Lambda) geometry rates, and matter rates. Record source file/blob provenance and operator parity for every derivative.
2. **Analytic base compatibility:** choose a regular, positive, metric-compatible base satisfying the declared discrete determinant/trace-free/connection conventions and the discrete CMC lapse equation. Demonstrate the geometric constraints and CMC residuals before any tangent evaluation.
3. **Directional compatibility:** construct plus/minus analytic states from the derived coupled direction. Verify metric positivity, centre parity/regularity, determinant and trace-free identities, and both original/differentiated CMC residuals. No coefficient fitting.
4. **Six-sector and tangent witnesses:** run the same pinned projection helper; retain all six sector pieces, signed semi-discrete momentum tangent at both epsilons, pointwise total-stress divergence, and independent finite-volume divergence.
5. **Fail-closed archive and evaluation:** use a new class label, immutable archive, dedicated loader, and the same full-Bianchi predictor. Verify exact provenance and zero-step tripwire; production step count must remain zero and state digest unchanged.
6. **Comparison:** only after all preceding gates pass, compare the new coupled class with the three existing immutable classes. Preserve centre cells 0/1, cell 2, fixed-radius and annular witnesses, both epsilons, and raw versus time-density-corrected sector accounting.

## Explicit non-goals

- No production equation, initializer, vendor, kernel, or accepted archive changes.
- No production RHS, initializer, CMC production solver, vendor RHS, vendor step, PIRK step, or time integrator call.
- No predictor retuning, coefficient fitting, residual-based term selection, or operator substitution.
- No scientific admission, acceptance-oracle, physical-exchange, or conservation-closure flag changes.
- No claim that a passing software or analytic-control test establishes a complete Einstein-matter solution.

## Independent CMC discrete solve and tangent gate

A pure research implementation now exists at `tools/issue4_cmc_lapse_tangent.py`. It reconstructs the discrete tridiagonal CMC operator from supplied arrays, solves the base lapse with the centre Neumann row and fixed unit outer normalization, and independently differentiates:
- the proper-volume weights and their total;
- the proper-volume projected CMC target `kdot_cmc`, including the supplied instantaneous `K RHS dot`;
- the radial first-derivative coefficient (c);
- the potential (mathcal Q), inverse radial metric coefficient, and (m=mathcal Q/(X^2/a));
- the complete discrete matrix and right-hand side;
- the lapse tangent from (A\,\dot\alpha=\dot b-\dot A\,\alpha).

The differentiated residual is retained directly. The plus/minus witness independently resolves the perturbed discrete systems and compares their central difference with the analytic lapse tangent. The module imports no engine, production CMC, vendor, RHS, or stepper code. Its cell-centred fourth-order parity-ghost derivative is independently transcribed from vendor commit `d6052d605673ce9d82cdc99b0df79855fa2ea215`, file `r0_operators.py`, blob `7f92fd9857981561dc2dff568addae87790c989f`; the implementation preserves the two parity ghosts, centred five-point stencil, and polynomial one-sided outer rows. The CMC tests now use this independent transcription rather than a generic test derivative.

CI: [run #42087257240](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42087257240) — **5 tests passed** with the independent pinned-operator transcription, compile passed, and the static import guard confirmed no production/vendor imports. The earlier algebra-only runs are retained in workflow history. The tests cover all discrete rows, positive lapse, differentiated-system residual, central-difference agreement, proper-volume projection response, and fail-closed invalid inputs.

Important scope limitation: this is the isolated discrete CMC algebra gate. Its current tests inject a fixed linear test derivative operator; the production-pinned fourth-order parity-aware operator and independently reconstructed physical (K)-RHS directional array have not yet been wired into this module. In particular, `k_rhs_dot` is an explicit required input and must be derived from the full instantaneous geometry/matter tangent—not estimated from the lapse residual. Thus this passing CI does **not** yet establish a pinned physical CMC lapse tangent or a coupled control.

## Current preflight conclusion

The source map and isolated CMC algebra gate establish the correct differentiated system and its testable contract. The coupled fixture is **not yet built**. The next milestone is to bind the independent solver to the exact pinned fourth-order parity-aware derivative operator, independently reconstruct the base (K)-RHS and its directional derivative from the CMC geometry plus six-sector source equations, and verify the differentiated CMC identity on a regular metric-compatible base. Do not label a frozen-lapse or vacuum-geometry shortcut as the coupled control.
