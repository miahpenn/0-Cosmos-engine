# GEAR-63 explicit-\u03be equation regression — conditional model only

**Date:** 15 October 2026  
**Base used for the research branch:** `614167904cc8cb80cf0d55521375bc52ca1362ff`  
**Scope:** analytic/numeric regression only; no production equation, vendor, dependency, or archive modification.

## Declared assumptions

This audit tests the action

\[
S=\int d^4x\sqrt{-g}\left[\frac{R}{2\kappa}
-\frac12(\nabla d)^2-V(d)-\frac12\xi R d^2\right],
\qquad
V(d)=-\frac12 m^2d^2,
\]

with signature \((-+++ )\), the standard Ricci convention, \(\kappa=8\pi G>0\), a canonical scalar, and no other matter or constant-potential offset. Here \(m^2\) is the positive inverted-quadratic parameter (the \(\alpha^2\) used in the GEAR-63 proposal). This is a *conditional test model*, not a recovered authoritative GEAR/COSMOS action.

The independent reference for the field-equation form is Hrycyna, “The non-minimal coupling constant and the primordial de Sitter state,” *European Physical Journal C* **80**, 817 (2020), especially Eqs. (5)–(9), which state the nonminimal stress tensor, scalar equation, and flat-FLRW constraint/acceleration equations: https://link.springer.com/article/10.1007/s10052-020-8397-5. The paper uses its own \(\kappa^2=8\pi G\) notation; this note writes \(\kappa=8\pi G\). Its discussion of \(\xi=3/16\) is explicitly the five-dimensional conformal value; the four-dimensional conformal value tested below is \(\xi=1/6\).

## Algebra under those assumptions

Define two distinct factors:

\[
D_{\rm grav}=1-\kappa\xi d^2,
\qquad
D_{\rm trace}=1-\kappa\xi(1-6\xi)d^2.
\]

\(D_{\rm grav}\) is the effective Einstein prefactor. \(D_{\rm trace}\) is the denominator produced when the Einstein trace and scalar equation are combined to eliminate \(\Box d\). They are **not interchangeable**. The regression reports them separately and rejects points at either pole. A negative but nonzero prefactor remains explicit in the returned state; nonsingularity alone must not be mistaken for a physically acceptable positive-effective-gravity regime.

For field velocity \(v=\dot d\), elimination gives

\[
R=-\kappa\,
\frac{(1-6\xi)v^2+(2-6\xi)m^2d^2}
{D_{\rm trace}}.
\]

The curvature contribution to the scalar acceleration is \(F_\xi=-\xi R d\). Including the inverted-potential term, the net force in \(\ddot d+3H\dot d=F_{\rm net}\) is

\[
F_{\rm net}=m^2d-\xi R d
=\frac{m^2d(1+\kappa\xi d^2)+\kappa\xi(1-6\xi)v^2d}
{D_{\rm trace}}.
\]

The code checks the original trace residual after inserting \(\Box(d^2)=-2m^2d^2+2\xi R d^2-2v^2\), compares the direct force with the separately written eliminated expression, and checks the sign of the homogeneous scalar equation.

The complete flat-FLRW constraint is

\[
3H^2D_{\rm grav}
=\kappa\left(\frac12v^2+V(d)+6\xi Hdv\right).
\]

The linear-in-\(H\) term is retained when solving its quadratic, rather than checking only \(H^2\) from a minimally coupled approximation.

## Analytic limits and the negative-\(\xi\) root

- **\(\xi=0\):** \(R=-\kappa(v^2+2m^2d^2)\), \(F_\xi=0\), and \(F_{\rm net}=m^2d\).
- **\(\xi=1/6\) in four dimensions:** \(D_{\rm trace}=1\) and \(R=-\kappa m^2d^2\), independent of \(v\). This does not remove the separate possible pole in \(D_{\rm grav}\).
- **Reduced rest-force root for \(\xi<0\):** \(d_*^2=-1/(\kappa\xi)\), with \(D_{\rm grav,*}=2\), \(D_{\rm trace,*}=2-6\xi\), \(R_*=m^2/\xi<0\), and \(F'_{\rm net}(d_*)=m^2/(3\xi-1)\). The derivative is a one-dimensional reduced-force criterion only.

For that root, at \(v=0\), the full flat-FLRW constraint requires

\[
H^2=\frac{\kappa V(d_*)}{3D_{\rm grav,*}}
=-\frac{\kappa m^2d_*^2}{12}<0.
\]

Thus the nonzero \(\xi<0\) root is **a reduced-force root but not a real spatially flat FLRW state** for this negative zero-offset potential. An additional source/offset or a different background changes the problem and must be derived consistently; it may not be inserted just to preserve the root. Full coupled metric-scalar stability remains **unknown**.

## Regression coverage and limits

`tests/test_gear63_explicit_xi_regression.py` exercises a nonsingular parameter grid, original-trace closure, force-form equivalence, scalar-equation sign, \(\xi=0\) and \(\xi=1/6\) limits, each denominator pole independently, the full FLRW quadratic, and the negative-\(\xi\) root counterexample. The implementation lives in `tools/gear63_explicit_xi_regression.py` and is isolated from the production runtime.

These checks can expose algebra/sign/domain errors *within this declared action*. They do not establish that the archive intended this action or any particular \(\xi\); do not infer a localized \(Rd/6\) force from a curvature trace; do not choose an S/D normalization or source coefficient; and do not claim stability, a cosmological bounce, or scientific admission from an algebraic regression. Issue #210 remains open until primary-action/projection provenance or an explicit owner scope decision resolves that separate gate.
