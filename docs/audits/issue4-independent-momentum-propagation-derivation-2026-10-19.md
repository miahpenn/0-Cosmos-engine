# Issue #4 — independent momentum-propagation RHS derivation

Date: 2026-10-19  
Status: derivation in progress; no closure or scientific-admission claim.

## Scope and provenance

This note extends the diagnostic-only investigation on [PR #408](https://github.com/miahpenn/0-Cosmos-engine/pull/408). The accepted-state replay remains pinned to historical production source `7f47669a1c0edfafef9fc5f27b2db80840d3f035` and vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. The recorded initial/accepted digests and the raw witness artifact remain the comparison state. PR #406 is untouched. No production equations, pins, or admission rules are changed.

This is an algebraic derivation plan, not a report of a completed numerical closure test. The measured finite-difference tangent remains a comparison measurement only.

## 1. Constraint decomposition

The pinned observer uses

\[
M_r=M_{\rm geom}-8\pi J_r.
\]

The independently predicted tangent must therefore be assembled as

\[
\dot M_r^{\rm pred}=\dot M_{\rm geom}^{\rm RHS}-8\pi\dot J_r^{\rm RHS}.
\]

Both terms must be derived from the same accepted state and the recorded production RHS rates. Do not compute the predicted tangent by finite-differencing the constraint observer.

## 2. Geometric chain rule

With the production conventions \(A_b=-A_a/2\), \(\chi_r=-\tfrac12 D_1X/X\), and the existing parity-aware discrete derivative \(D_1\), the observed geometric constraint is

\[
M_{\rm geom}=D_1A_a-\frac23D_1K+6A_a\chi_r+
\frac32A_a\left(\frac2r+\frac{D_1b}{b}\right).
\]

At fixed grid coordinates, let \(\dot A_a,\dot K,\dot X,\dot b\) be the accepted RHS rates. The independent chain-rule prediction is

\[
\dot\chi_r=-\frac12\left(\frac{D_1\dot X}{X}
-\frac{(D_1X)\dot X}{X^2}\right),
\]

\[
\begin{aligned}
\dot M_{\rm geom}^{\rm RHS}={}&D_1\dot A_a-\frac23D_1\dot K
+6(\dot A_a\chi_r+A_a\dot\chi_r)\\
&+\frac32\left[
\dot A_a\left(\frac2r+\frac{D_1b}{b}\right)
+A_a\left(\frac{D_1\dot b}{b}
-\frac{(D_1b)\dot b}{b^2}\right)
\right].
\end{aligned}
\]

The grid coordinate \(r\) and derivative operator are fixed in this directional derivative. The same derivative stencil and parity convention as the production observer must be used; do not substitute a continuum derivative or a different boundary stencil.

## 3. Scalar-current chain rule

The audited scalar current pieces are \(J_S=-P_S D_1S\), \(J_D=-P_D D_1D\), and \(J_\phi=-\Pi D_1\phi/\kappa\), with \(\kappa=8\pi\). Their RHS derivatives are

\[
\dot J_S=-\dot P_S D_1S-P_S D_1\dot S,
\qquad
\dot J_D=-\dot P_D D_1D-P_D D_1\dot D,
\]

\[
\dot J_\phi=-\frac{\dot\Pi D_1\phi+\Pi D_1\dot\phi}{\kappa}.
\]

These expressions use the accepted scalar rates and the exact production discrete derivative. The \(S\), \(D\), and COSMOS-\(\phi\) terms must remain separately visible in the signed report.

## 4. Conservative-fluid current: required independent derivative

The fluid current is obtained from the production `project_species` mapping, which depends on both the conservative species state and the metric slice. Its time derivative must account for all of the following:

- each species' accepted conservative rates (rest, coordinate energy, and radial momentum);
- evolving metric fields and metric-derived quantities used by primitive recovery/projection;
- the exact species-specific projection and sign/normalization conventions;
- all enabled baseline species (dark matter, baryons, radiation), with any unsupported opt-in source causing a fail-closed error.

Do not approximate this term by differentiating the total constraint observer. The implementation should expose the per-species signed \(\dot J\) arrays and a separately computed total. If a derivative of the projection map is evaluated numerically, it must be a separately specified, converged projection-Jacobian check with declared step-size/error bounds—not the constraint finite-difference tangent relabeled as independent.

## 5. Closure protocol

On the exact accepted state and rates:

1. Independently assemble \(\dot M_{\rm geom}^{\rm RHS}\), each signed scalar \(\dot J\), each signed fluid \(\dot J\), and \(\dot M_r^{\rm pred}\).
2. Compare predicted and measured signed \(\dot M_r\) at cells 0 and 1; also report each cell separately across fixed-coordinate annuli. Do not conflate fixed-coordinate annuli with fixed-areal-radius worldtubes.
3. Preserve raw arrays and their hashes; record the source and vendor pins, state digests, derivative/stencil conventions, and numerical bounds.
4. Include manufactured controls for the geometric chain rule, scalar-current chain rule, and fluid projection derivative. Establish tolerances before examining the closure residual.
5. Fail closed for missing sectors, non-finite arrays, source-pin or digest mismatch, failed controls, or failed signed closure. Do not tune coefficients to match the measured tangent.

## 6. Interpretation gates

Keep three outcomes distinct:

- **Implementation consistency:** the existing measured tangent is internally stable.
- **Independent propagation closure:** the independently assembled RHS reproduces the measured signed tangent within justified, preregistered bounds.
- **Scientific admission:** existing project gates are satisfied; a one-state closure result alone does not establish this.

Until the independent calculation and controls run successfully, `independent_full_propagation_identity_closed=false`, `scientific_admission=false`, and `acceptance_oracle_satisfied=false`. PR #406 remains untouched.
