# Issue #4 physical CMC K-RHS directional derivation

**Date:** 2026-11-05  
**Branch:** `research/issue4-coupled-geometric-matter-control-2026-11-05`  
**Status:** analytic derivation; implementation and physical fixture not yet complete  
**Scope:** research-only, instantaneous, zero-step. No production RHS, production CMC solver, vendor RHS, or stepper calls.

## 1. Physical K equation

Define
[
\mathcal Q=\frac32 A_a^2+\frac13K^2+
4\pi(\rho+p_r+2p_t),\qquad \nu=\frac{X^2}{a},
]
and
[
c=-\frac{D_1a}{2a}+\frac{D_1b}{b}
-\frac{D_1X}{X}+\frac2r.
]
On the zero-shift CMC slice,
[
K_{\rm RHS}(\alpha)=-\Delta_g\alpha+\alpha\mathcal Q
=-\nu(D_2\alpha+cD_1\alpha)+\alpha\mathcal Q.
]
This is the trace of the pinned primary L2 plus L3 equation: L2 supplies (-\Delta_g\alpha), while L3 supplies (alpha(3A_a^2/2+K^2/3)) and the six-sector matter source (4\pi\alpha(\rho+p_r+2p_t)).

## 2. Stage-aware lapse contract: keep input and output distinct

The inspected CMC implementation is stage-aware. Its target routine projects the actual K RHS using the **incoming** `geometry.alpha`; the lapse solver then solves for a new **output** lapse. Thus write these as (\alpha_{\rm in}) and (\alpha_{\rm out}), respectively:
[
\lambda=P_w[K_{\rm RHS}(\alpha_{\rm in})],
\qquad
L(a,b,X,\mathcal Q)\alpha_{\rm out}=-\lambda,
]
with the fixed centre regularity row and outer normalization. Here (P_w[f]=\sum_iw_if_i/W).

The source contract was inspected at repository commit `0c40172f32b72a8149549adadf63d668a2f24f3e`, file `engine/cmc_gauge.py`; the same file is present on this research branch. The source confirms:
- `target_kdot` forms L2 plus L3 K RHS from the incoming geometry lapse and projects it over all proper-volume cells;
- `solve_cmc_lapse` then builds the lapse operator and solves for the output lapse;
- (w_i=r_i^2\sqrt{a_i}b_iX_i^{-3}\Delta r);
- (c) is assembled using the grid's fourth-order parity-aware first derivative;
- the interior lapse matrix uses the three-point second-order radial Laplacian row with the centred (cD_1) coefficient;
- the centre row is (\alpha_0-\alpha_1=0), and the outer row fixes (\alpha_{N-1}=1).

Therefore the three-point interior matrix in `tools/issue4_cmc_lapse_tangent.py` matches the inspected CMC solver's matrix structure; it is not a guessed fourth-order replacement. Its coefficient derivative must nevertheless use the same fourth-order parity-aware (D_1), and all row placements must match the source.

**Do not silently identify** (\alpha_{\rm in}) and (\alpha_{\rm out}). They are separate arguments/outputs of the stage-aware map. A self-consistent analytic fixed-point control may deliberately impose (\alpha_{\rm in}=\alpha_{\rm out}), but that is an additional control condition, not an automatic property of the stage-aware implementation. The direction must carry both (\dot\alpha_{\rm in}) and (\dot\alpha_{\rm out}) until that condition is explicitly established.

## 3. Directional derivative of the physical K RHS

For a fixed grid and fixed linear discrete derivative operators,
[
\frac{\dot\nu}{\nu}=2\frac{\dot X}{X}-\frac{\dot a}{a},
]
[
\dot c=-\frac12\left(\frac{D_1\dot a}{a}
-\frac{(D_1a)\dot a}{a^2}\right)
+\frac{D_1\dot b}{b}-\frac{(D_1b)\dot b}{b^2}
-\frac{D_1\dot X}{X}+\frac{(D_1X)\dot X}{X^2},
]
[
\dot{\mathcal Q}=3A_a\dot A_a+\frac23K\dot K
+4\pi(\dot\rho+\dot p_r+2\dot p_t).
]
The physical projected K-RHS direction is
[
\boxed{
\begin{aligned}
\dot K_{\rm RHS}^{\rm in}={}&
-\dot\nu(D_2\alpha_{\rm in}+cD_1\alpha_{\rm in})\\
&-\nu(D_2\dot\alpha_{\rm in}
+cD_1\dot\alpha_{\rm in}+\dot cD_1\alpha_{\rm in})\\
&+\dot\alpha_{\rm in}\,\mathcal Q
+\alpha_{\rm in}\dot{\mathcal Q}.
\end{aligned}}
]
Here (dot K) inside (dot{\mathcal Q}) is the physical base (K_{\rm RHS}^{\rm in}). The six-sector source directions and all other dotted geometry/primary fields must be derived independently from the pinned instantaneous equations. No free `k_rhs_dot` profile is permitted.

## 4. Differentiate the proper-volume projection

With
[
w_i=r_i^2\frac{\sqrt{a_i}b_i}{X_i^3}\Delta r,
\qquad W=\sum_iw_i,
\qquad
\frac{\dot w}{w}=\frac12\frac{\dot a}{a}
+\frac{\dot b}{b}-3\frac{\dot X}{X},
]
the target derivative is
[
\boxed{
\dot\lambda=\frac1W\sum_i
\left[\dot w_i((K_{\rm RHS}^{\rm in})_i-\lambda)
+w_i(\dot K_{\rm RHS}^{\rm in})_i\right].
}
]
The projection includes the complete proper-volume domain, including boundary cells, exactly as the source implementation does. It is differentiated from the physical K RHS built with (\alpha_{\rm in}), not from the solved output lapse unless the fixed-point identity is part of the control definition.

## 5. Differentiate the output CMC lapse solve

The discrete output equation is
[
(D_2+cD_1-\mathcal Q/\nu)\alpha_{\rm out}
=-\lambda/\nu
]
on interior rows, with the fixed centre and outer rows above. Let its assembled matrix and RHS be (A) and (b). Its directional system is
[
\boxed{A\dot\alpha_{\rm out}
=\dot b-\dot A\alpha_{\rm out}.}
]
The interior RHS derivative is
[
\dot b_i=-\frac{\dot\lambda}{\nu_i}
+\frac{\lambda\dot\nu_i}{\nu_i^2}.
]
The matrix derivative contains (\dot c), (\dot{\mathcal Q}), and (\dot\nu), with the exact three-point row placement. The centre and fixed outer-normalization rows have zero directional RHS.

This is the stage-aware coupling: the incoming lapse direction enters (dot K_{\rm RHS}^{\rm in}), then (dot\lambda) enters the output lapse tangent equation. If the analytic control additionally imposes (\alpha_{\rm in}=\alpha_{\rm out}) and (\dot\alpha_{\rm in}=\dot\alpha_{\rm out}), these equations become one coupled linear system in the common lapse tangent and (dot\lambda). Otherwise retain the two lapse directions and do not force a fixed point.

## 6. Implementation and validation gates

1. Reconstruct the incoming physical K RHS from the independent L2 trace, L3 trace, and six-sector projections. Cross-check its algebraic form against the formula in §1 without importing or calling production/vendor routines.
2. Reconstruct the proper-volume projection and its derivative from the incoming lapse, including (\dot w) and (\dot K_{\rm RHS}^{\rm in}).
3. Bind the output lapse matrix to the inspected source stencil and exact row convention. Retain base and differentiated row-wise residuals.
4. Solve for (\dot\alpha_{\rm out}) from the differentiated output system. Preserve the distinction between input and output lapse directions.
5. If and only if the chosen analytic base is a fixed point, impose and test the input/output lapse equality and solve the corresponding coupled block system.
6. Independently construct plus/minus states at both required epsilons and compare the incoming K-RHS projection, target rate, output lapse, and output lapse tangent. Check positivity, parity, determinant/trace-free identities, six-sector projections, and all CMC residuals.
7. Only after these pass, construct the new immutable coupled archive and run the existing zero-step two-channel observer.

No production code, accepted archive, existing comparison class, predictor, admission flag, or conservation flag is changed by this derivation.


## 7. Independent stage-aware implementation checkpoint (2026-11-05)

Implemented in `tools/issue4_cmc_lapse_tangent.py`:
- `physical_k_rhs_and_direction` independently transcribes the native warped-product lapse-Laplacian identity using (A=a/X^2), (R=r\sqrt b/X), and the parity-aware (D_1,D_2) operators. It differentiates (A), (R), the lapse Hessian, and the matter/primary potential analytically.
- `solve_stage_aware_cmc_tangent` constructs incoming (K_{\rm RHS}) and its direction internally, projects the incoming source over the full proper-volume measure, solves the distinct output lapse, and differentiates that output solve. Its API has no `k_rhs_dot` argument.
- Incoming `alpha_in_dot) and output `alpha_dot) remain distinct. The existing output operator keeps the source-matched three-point interior row and fixed centre/outer rows.
- No engine, production, or vendor imports are used.

The bounded workflow run for commit `a21363a3557f16006bacedb31691b89ac85856c2` completed the compile, all seven tests, and the static no-production/no-vendor import guard successfully:
[CI run #42092298042](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42092298042).

**Scope limitation:** these new plus/minus tests use a manufactured smooth control and combined ((\rho,p_r,p_t)) input profiles. They verify the independent directional algebra and stage-aware solve, not yet the actual pinned accepted-state six-sector arrays. The next binding gate is to feed the locked six-sector projection outputs and independently derived coupled geometry/primary directions into this API, compare the reconstructed incoming (K_{\rm RHS}) against the pinned L2+L3 component ledger, and run the required two-epsilon plus/minus checks on the accepted (t=0.60) state. Do not claim that physical-state binding or coupled-fixture admission has passed yet.


### Discrete-operator equivalence boundary

The compact (-\nu(D_2+cD_1)\alpha) expression is the radial scalar-Laplacian form of the continuum identity. The new independent code currently differentiates the native warped-product transcription ((A=a/X^2), (R=r\sqrt b/X)) using parity-aware fourth-order (D_1,D_2). Because finite-difference product/quotient operations need not commute, this is **not yet asserted to be bitwise or tolerance-equivalent** to the pinned vendor `primary_l2_rhs()["K"]) used by the stage-aware production projection. That crosswalk must be measured on the exact pinned states before the new independent reconstruction is accepted as the production-discrete K-RHS transcription.

Likewise, the present manufactured regression passes combined total ((\rho,p_r,p_t)) arrays. It does not yet independently reconstruct and bind each of the six named sectors (dark matter, baryons, radiation, S, D, phi) and their directional projections from the accepted archive. These two source-fidelity checks are explicit prerequisites for the accepted-state gate; the seven-test synthetic pass does not waive them.
