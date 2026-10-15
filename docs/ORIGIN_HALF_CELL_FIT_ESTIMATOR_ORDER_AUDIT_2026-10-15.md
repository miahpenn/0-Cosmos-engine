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


## Additional control — exact origin value available

For the production initializer's integrand,

\[
s(r)=1-8\pi r^2\rho+2r^2K_rK_t+r^2K_t^2,
\]

regularity and finite source fields give \(s(0)=1\) exactly. This does not make
the rest of the source known, but it means the constant coefficient need not be
inferred from a finite-radius fit. The standalone fitter now offers an
origin-anchored option that fixes \(c_0=s(0)\) and fits only the higher powers.

For the same manufactured smooth source used above (whose exact origin value is
\(1.0002\)), the anchored results are:

| Grid spacing \(h\) | Fixed window: signed error | Fixed-window order | Local window: signed error | Local-window order |
|---:|---:|---:|---:|---:|
| 0.125 | \(+2.8073478\times10^{-6}\) | — | \(+1.3056263\times10^{-9}\) | — |
| 0.0625 | \(+3.6155624\times10^{-7}\) | 2.957 | \(+2.5526248\times10^{-12}\) | 8.999 |
| 0.03125 | \(+4.5530052\times10^{-8}\) | 2.989 | \(+4.9925342\times10^{-15}\) | 8.998 |

The exact degree-three polynomial control remains at roundoff for both windows
with and without anchoring. **Anchoring the origin changes the fixed-window
asymptotic result from first order to third order on this manufactured source**:
the remaining leading error comes from the fitted quadratic coefficient, not
from a fitted constant term. The local six-cell fit remains approximately ninth
order for this source.

### Offline reconciliation against the existing production-source artifact

This is a post-processed check of the existing source-pinned PR #257 artifact,
not a new engine run. I downloaded artifact 11993174037 from
[workflow 38776498208](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38776498208)
and verified its ZIP SHA-256 as
`ff667892633d20c17fda62ceb9752de420f79276d0681ed9e8195dd1cb65f5c2`.
Its production-source anchor is
`285a2ee0c1a3b97cddeb4519d8e749109b115aab`; the initializer blob
`7db2e1a10cf16fb40e4221b0df504fdf39aae611` is identical at that anchor and
current main. Using the artifact's actual last-iteration source samples at
N=160/240/320 and fitting degree 3 over \(r\le2.5\), the known \(s(0)=1\)
gives:

| Resolution | Unconstrained fit \(c_0-1\) | Intercept term \((c_0-1)r_0\) | Fit-minus-point moment correction | Intercept term / correction |
|---:|---:|---:|---:|---:|
| N=160 | \(-3.233\times10^{-8}\) | \(-8.081\times10^{-9}\) | \(+1.30448\times10^{-5}\) | \(-0.062\%\) |
| N=240 | \(-3.352\times10^{-8}\) | \(-5.587\times10^{-9}\) | \(+3.86825\times10^{-6}\) | \(-0.144\%\) |
| N=320 | \(-7.352\times10^{-8}\) | \(-9.190\times10^{-9}\) | \(+1.62615\times10^{-6}\) | \(-0.565\%\) |

These are a decomposition of the **intercept component**, not an exact error in
the unknown production half-cell integral. The current source-profile intercept
bias is many orders smaller than the adversarial manufactured example and is
below 0.6% of the fitted correction even at N=320. This weakens any claim that
the generic fixed-window O(h) mechanism is the leading explanation of the
present fit correction. It does not rule out other fit-coefficient or
operator/geometry errors, and does not quantify the Ricci response to an
intercept-only perturbation.

## Interpretation and next gate

**Established:** a smooth even source exists for which fixed-physical-window,
fixed-degree fitting loses the expected half-cell moment order, while a
shrinking local window recovers much higher order in this manufactured
control. Finite iteration defect, fit residual, and convergence of the
source/metric fixed point are not substitutes for a quadrature-order test.

**Not established:** that the production source has the same intercept bias,
that the local-window method is uniquely preferred, or that changing the
initializer would fix the entire Ricci/constraint discrepancy.

The offline source-profile check now bounds the fitted-intercept contribution
for the available N=160/240/320 profiles, but it does not provide an exact
production half-cell integral or a Ricci response bound. Next, a source-pinned
anchored-intercept counterfactual would need to report both the moment and the
resulting signed centre-Ricci gap on matching profiles, while preserving the
origin-moment and Lambda/Ricci derivative channels as separate diagnostics.
Do not select a production fit/window from the manufactured controls alone.
Scientific admission remains **NOT ADMITTED**.
