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
and the radial scalar-Laplacian coefficient
[
c=-\frac{D_1a}{2a}+\frac{D_1b}{b}
-\frac{D_1X}{X}+\frac2r.
]
For the locked zero-shift CMC slice, the instantaneous trace equation is
[
K_{\rm RHS}=-\Delta_g\alpha+\alpha\mathcal Q
=-\nu(D_2\alpha+cD_1\alpha)+\alpha\mathcal Q.
]
The same declared discrete operators must be used for the base equation, its direction, and the plus/minus witnesses. The base (K_{\rm RHS}) is therefore reconstructed from the lapse and geometry/matter arrays; it is not a free profile.

## 2. Directional derivative

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
Thus
[
\boxed{
\begin{aligned}
\dot K_{\rm RHS}={}&
-\dot\nu(D_2\alpha+cD_1\alpha)\\
&-\nu(D_2\dot\alpha+cD_1\dot\alpha+\dot cD_1\alpha)\\
&+\dot\alpha\,\mathcal Q+\alpha\dot{\mathcal Q}.
\end{aligned}}
]
Here (dot K) inside (dot{\mathcal Q}) is the physical base (K_{\rm RHS}). It is not an independently supplied tangent. All other dotted geometry and six-sector source fields must come from the independently derived coupled instantaneous direction.

This expression is affine-linear in the unknown (dot\alpha). Its known terms contain the metric, primary, and six-sector directions; its (dot\alpha)-dependent terms are the linear lapse-feedback operator.

## 3. Proper-volume projection and coupled base solve

The locked proper-volume weights are
[
w_i=r_i^2\frac{\sqrt{a_i}b_i}{X_i^3}\,\Delta r,
\qquad W=\sum_iw_i.
]
The CMC target is
[
\lambda=\dot K_{\rm CMC}
=\frac{\sum_iw_i(K_{\rm RHS})_i}{W}.
]
Because (K_{\rm RHS}) itself depends on (alpha), a physical base solve must solve for ((\alpha,\lambda)) together, not compute (lambda) from an externally prescribed (K_{\rm RHS}) and then solve for (alpha).

The discrete augmented system has (N+1) unknowns:
- the centre regularity row for (alpha);
- the (N-2) interior CMC equations (K_{\rm RHS}=\lambda);
- the fixed outer lapse normalization;
- the proper-volume projection equation (\lambda-\sum_iw_i(K_{\rm RHS})_i/W=0).

The last row uses the same declared discrete (K_{\rm RHS}) operator at every cell, including boundary cells. This makes the boundary closure and projection part of the system rather than silently dropping boundary contributions. The exact boundary stencil and centre/outer rows must be pinned to the source contract before implementation.

## 4. Coupled directional solve

Differentiate the projection:
[
\boxed{
\dot\lambda=
\frac{1}{W}\sum_i\left[
\dot w_i((K_{\rm RHS})_i-\lambda)
+w_i(\dot K_{\rm RHS})_i
\right],
}
]
where
[
\frac{\dot w}{w}=\frac12\frac{\dot a}{a}
+\frac{\dot b}{b}-3\frac{\dot X}{X}.
]
Together with the differentiated interior equations
[
\dot K_{\rm RHS}=\dot\lambda,
]
the differentiated centre row and fixed-normalization row, this is an augmented linear system for ((\dot\alpha,\dot\lambda)). The (dot K_{\rm RHS}) expression in §2 supplies the lapse-feedback coefficients. This replaces the current free `k_rhs_dot` input with a derived coupled unknown.

The resulting (dot\lambda) must agree with the direct differentiated proper-volume projection above. Retain both the original and differentiated equation residuals, the projection residuals, and the plus/minus independently re-solved witness.

## 5. Important implementation audit before coding

The current `tools/issue4_cmc_lapse_tangent.py` correctly differentiates its own supplied-array algebra, but its assembled interior radial matrix is a three-point, second-order (D_2+cD_1) stencil. Its coefficient-derivative helper is a separate fourth-order parity-aware transcription. Those facts do **not** establish that the matrix is the exact pinned CMC operator.

Before binding the augmented solve, inspect and transcribe the actual pinned CMC operator, including:
- cell-centred first- and second-derivative stencils and parity ghosts;
- centre and outer boundary rows;
- radial coefficient placement and any face/cell averaging;
- the precise proper-volume projection domain and weights;
- whether the outer normalization is fixed or otherwise specified by the source.

Do not call the current three-point matrix the pinned physical CMC operator unless that audit establishes identity. If the pinned operator is wider-band, the augmented system must preserve that exact structure.

## 6. Validation gates

1. Independently reconstruct the physical base (K_{\rm RHS}) from (alpha), metric, primary fields, and six-sector sources.
2. Solve the augmented base ((\alpha,\lambda)) system and retain row-wise residuals.
3. Independently derive the full geometry and matter direction; insert the physical base (K_{\rm RHS}) as (dot K) in (dot{\mathcal Q}).
4. Assemble and solve the augmented ((\dot\alpha,\dot\lambda)) system. No free `k_rhs_dot`.
5. Compare against independently reconstructed plus/minus augmented base solves at both required epsilons; check metric positivity, parity, determinant/trace-free identities, six-sector projections, and both CMC residuals.
6. Only after these pass, construct the new immutable coupled archive and run the existing zero-step two-channel observer.

No production code, accepted archive, existing comparison class, predictor, admission flag, or conservation flag is changed by this derivation.
