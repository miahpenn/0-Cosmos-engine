# Origin half-cell fit-estimator order audit — 15 October 2026

## Status

Research-only manufactured-source test on a branch cut from live main at
`6398620b85b7f5a027de0532ded80bdf04a9b53b`. No production
engine, vendor kernel, dependency, initialized state, or campaign setting is
changed. No evolution is performed.

This audit follows the fit-family sensitivity identified in draft PR #257.
It checks a methodological question in isolation: does fixed-physical-radius
polynomial fitting preserve the expected order for the integral over the
origin half-cell?

## Derivation

For a regular even source, write

\[
s(r)=s_0+s_2r^2+s_4r^4+\cdots,\qquad h=\Delta r,\quad r_0=h/2.
\]

The exact origin-cell contribution is

\[
I(h)=\int_0^{h/2}s(r)\,dr
=s_0\frac{h}{2}+s_2\frac{h^3}{24}+O(h^5).
\]

The point-sampled half-cell used by the baseline is

\[
I_{\rm point}=r_0s(r_0)
=s_0\frac{h}{2}+s_2\frac{h^3}{8}+O(h^5),
\]

so its leading half-cell error is \(O(h^3)\).

A degree-three fit in \(r^2\), performed on a fixed physical window
\(0<r\le L\), converges under refinement toward a least-squares projection
over that non-shrinking window. For a generic smooth source that is not
exactly cubic in \(r^2\) throughout the window, the fitted intercept need not
converge to \(s(0)\). If its limiting bias is \(\epsilon_0\), the estimated
half-cell moment contains an error \(\epsilon_0 h/2+O(h^3)\), i.e. first-order
moment error instead of the expected third-order behaviour.

This can be more serious in the centre-cell metric than the moment order
alone suggests. The initializer forms the radial metric source integral
average as \(B_i=I(r_i)/r_i\). At the first cell, \(r_0=h/2\), so a moment
error \(\delta I=O(h)\) yields \(\delta B_0=\delta I/r_0=O(1)\). The error may
still decay at any fixed positive radius while failing to converge uniformly
as the first cell approaches the origin. This is a structural consequence of
the formula under a persistent fitted-intercept bias, not proof that the
production Ricci discrepancy has this cause.

A local window whose physical radius shrinks with \(h\) is a different
estimator. This audit uses six cell-centre samples for the local window. It
does not claim that this local fit is automatically the right production
prescription; the point is to measure the convergence order and intercept
bias against a known answer.

## Reproducible manufactured source

The source is

\[
s(r)=1+0.3r^2-0.04r^4+0.004r^6+0.001r^8+0.0002e^{-r^2}.
\]

Its exact moment is evaluated analytically, including the error-function term
from the exponential. A separate polynomial control omits the last two terms,
making it exactly cubic in \(r^2\). The fitting basis is scaled to
\(z=(r/L)^2\), which spans the same polynomial space as raw powers of \(r^2\)
while avoiding unnecessary conditioning issues for shrinking windows.

## Results

Signed errors below are estimated minus exact half-cell moment for the
non-polynomial smooth source.

| Grid spacing \(h\) | Fixed window \(L=2.5\): signed error | Fixed-window order | Local window \(L=6h\): signed error | Local-window order |
|---:|---:|---:|---:|---:|
| 0.125 | \(-4.8859530\times10^{-4}\) | — | \(-1.5675682\times10^{-8}\) | — |
| 0.0625 | \(-2.5571562\times10^{-4}\) | 0.934 | \(-3.0649490\times10^{-11}\) | 8.998 |
| 0.03125 | \(-1.2930933\times10^{-4}\) | 0.984 | \(-5.9869\times10^{-14}\) | 9.000 |

The last local-window value is close to floating-point cancellation; the
0.125-to-0.0625 comparison is the stronger observed-order measurement.

The exact degree-three polynomial control integrates to roundoff for both
the fixed and local windows across \(h=0.5,0.25,0.125,0.0625\) (largest
absolute error approximately \(4.5\times10^{-16}\)). This shows the quadrature
formula and implementation reproduce their declared polynomial oracle.

The fixed-window fitted-intercept bias approaches approximately
\(-8.28\times10^{-3}\) relative to the source value at the origin for this
manufactured source. Because the first-cell metric average is the half-cell
moment divided by \(r_0\), the corresponding first-cell metric bias likewise
stays finite rather than vanishing in the refinement sequence. The local
six-cell fit's intercept and first-cell errors fall rapidly with refinement.
These values are specific to this source and selected fit family, not
estimates of the production Ricci discrepancy.

## Interpretation and next gate

**Established:** a smooth even source exists for which fixed-physical-window,
fixed-degree fitting loses the expected half-cell moment order, while a
shrinking local window recovers much higher order in this manufactured
control. Finite iteration defect, fit residual, and convergence of the
source/metric fixed point are not substitutes for a quadrature-order test.

**Not established:** that the production source has the same intercept bias,
that the local-window method is uniquely preferred, or that changing the
initializer would fix the entire Ricci/constraint discrepancy.

Next, compare the two estimators on current source-pinned initializer
profiles at several resolutions, report intercept bias and signed half-cell
moment, and check held-out regular-source controls before proposing any
initializer change. Keep the origin-moment and Lambda/Ricci derivative
channels separate. Scientific admission remains **NOT ADMITTED**.
