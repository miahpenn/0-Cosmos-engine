# Continuous origin half-cell integral oracle — 15 October 2026

## Status and scope

**Research-only. Scientific admission remains NOT ADMITTED.** This work supplies a source-specific numerical reference for the origin half-cell moment missing from the sampled production profiles. It does not alter the production initializer, engine, vendor, dependencies, physical equations, campaign configuration, or saved state. The audit takes zero evolution steps.

## Motivation

The existing origin-fit studies establish that source-fit family and fit window affect the inferred half-cell correction. They do not, by themselves, establish the true integral of the continuous source. This audit adds a conditional continuous extension of the initializer's declared Gaussian initial fields, integrates that source on the interval from the origin to the first cell centre, and compares it with the production point-sample moment and even-polynomial fit candidates.

## Continuous source derivation

For the default uncharged initializer, the analytic fields are

\[
S(r)=A_S e^{-(r/w)^2},\qquad D(r)=A_D e^{-(r/w)^2},
\]

with constant \(\phi=\phi_0\), \(\Pi=\Pi_0\), and \(K_t=-H_0\). Their derivatives are \(S'=-2rS/w^2\) and \(D'=-2rD/w^2\). The initializer source is

\[
s(r)=1-8\pi r^2\rho+2r^2K_rK_t+r^2K_t^2.
\]

Because the scalar-gradient contributions to \(\rho\) are linear in the metric variable \(B\), and \(P_D=D\) cancels the two D-field terms independent of \(B\), this can be written as

\[
s(r)=A(r)-q(r)B(r),\qquad
q(r)=4\pi r^2\left[(S')^2+(D')^2\right].
\]

The radial metric reconstruction satisfies \(I'(r)=s(r)\) and \(B(r)=I(r)/r\), giving

\[
I'(r)+\frac{q(r)}{r}I(r)=A(r),\qquad I(0)=0,\quad B(0)=1,\quad s(0)=1.
\]

The audit integrates this regular initial-value problem with an implicit trapezoid update and successive halvings of the substep size. Richardson extrapolation estimates the limiting integral. The change between the final extrapolants is recorded as a *numerical convergence estimate*, not as a rigorous analytic remainder bound.

## Comparison layers

For each production resolution, the report separately records:

1. The continuously integrated half-cell moment \(I_{\mathrm{cont}}\).
2. The continuous-source midpoint moment, evaluated at the continuum fixed point.
3. The actual production point-sampled half-cell moment from a six-pass source replay, including the pinned grid derivative stencils.
4. The continuous quadrature delta, the continuum-versus-sampled-source delta, and their algebraic sum.
5. Fit-derived moments for the unanchored and origin-anchored degree-3 fit on the reference window, plus higher-degree/window controls where there are enough samples to keep the fit overdetermined.
6. The source-replay agreement with \(B=X^6\), state digests before/after observation, all source/vendor pins, and a SHA-256 manifest.

The comparison is deliberately split because the continuum source and the sampled discrete source use different derivative representations. The total difference cannot be labelled purely a quadrature error.

## Review gates and interpretation limits

The result should be used to assess how well each fit family approximates the specified continuous source integral. It does not prove this analytic continuation is the uniquely intended physical model, nor does a small fit residual imply a rigorous integral-error bound.

Before a production change is proposed, the team still needs independent numerical-method review, a stated acceptable error criterion, source-specific uncertainty propagation through the self-consistent initializer, signed Ricci/constraint witnesses, and an independent review of any production patch. The existing offline Lambda derivative counterfactual remains a separate evidence layer and is not merged into the quadrature claim.

## Reproduction

Run the specialized workflow on branch research/continuous-origin-half-cell-oracle-2026-10-15. It pins engine source 285a2ee0c1a3b97cddeb4519d8e749109b115aab and vendor d6052d605673ce9d82cdc99b0df79855fa2ea215, executes the unit tests, verifies that the production six-pass source replay matches the initialized state, emits deterministic JSON/CSV reports, and records checksums.
