# Continuous origin half-cell integral oracle — 15 October 2026

## Status and scope

**Research-only. Scientific admission remains NOT ADMITTED.** This work supplies a source-specific numerical reference for the origin half-cell moment missing from sampled production profiles. It does not alter the production initializer, engine, vendor, dependencies, physical equations, campaign configuration, or saved state. The audit takes zero evolution steps.

## Motivation

The existing origin-fit studies establish that source-fit family and fit window affect the inferred half-cell correction. They do not, by themselves, establish the continuous-source integral. This audit extends the initializer's declared Gaussian initial fields continuously, integrates the resulting source from the origin to the first cell centre, and compares it with the production point-sample moment and even-polynomial fit candidates.

## Continuous source derivation

For the default uncharged initializer, the analytic fields are

\[
S(r)=A_S e^{-(r/w)^2},\qquad D(r)=A_D e^{-(r/w)^2},
\]

with constant \(\phi=\phi_0\), \(\Pi=\Pi_0\), and \(K_t=-H_0\). Their derivatives are \(S'=-2rS/w^2\) and \(D'=-2rD/w^2\). The initializer source is

\[
s(r)=1-8\pi r^2\rho+2r^2K_rK_t+r^2K_t^2.
\]

Because the scalar-gradient contributions to \(\rho\) are linear in \(B\), and \(P_D=D\) cancels the two D-field terms independent of \(B\), this source can be written as

\[
s(r)=A(r)-q(r)B(r),\qquad
q(r)=4\pi r^2\left[(S')^2+(D')^2\right].
\]

The radial metric reconstruction satisfies \(I'(r)=s(r)\) and \(B(r)=I(r)/r\), giving

\[
I'(r)+\frac{q(r)}{r}I(r)=A(r),\qquad I(0)=0,\quad B(0)=1,\quad s(0)=1.
\]

The audit solves this regular initial-value problem with an implicit trapezoid update and successive halvings of the substep size. Richardson extrapolation gives a convergence estimate; that estimate alone is not an analytic remainder bound.

## Analytic integral and estimator enclosure

The source term independent of the metric feedback is

\[
A(r)=1+c_{\rm bg}r^2-4\pi A_S^2 r^2e^{-ar^2}
+\frac{16\pi K_tA_D^2}{w^2}r^4e^{-ar^2},\qquad a=2/w^2,
\]

where \(c_{\rm bg}=3K_t^2-8\pi(\rho_\phi+\rho_{\rm matter})\) cancels algebraically for the initializer's declared Friedmann value of \(H_0\). Its primitive is evaluated in closed form as

\[
I_A(r)=r+\frac{c_{\rm bg}r^3}{3}
-4\pi A_S^2F_2(r)
+\frac{16\pi K_tA_D^2}{w^2}F_4(r),
\]

with \(F_p(r)=\int_0^r t^p e^{-at^2}dt\) for \(p=2,4\), evaluated using a stable small-argument series or the error function. Also,
\[
q(r)=C r^4e^{-ar^2}\le Cr^4,\qquad
C=\frac{16\pi(A_S^2+A_D^2)}{w^4}.
\]

On each tested half-cell the code bounds \(A\) below by a positive value, and above by \(A_{\max}\). The regular integrating-factor solution therefore has \(0<B(r)\le A_{\max}\). Writing \(Q=\int_0^{r_0}q(r)B(r)dr\) gives the analytic inequalities

\[
0\le Q\le Q_{\max}=\frac{A_{\max} C r_0^5}{5},\qquad
I_A-Q_{\max}\le I(r_0)\le I_A.
\]

Thus any fitted half-cell estimator \(F\) has a **signed estimator-minus-integral error enclosure**

\[
F-I(r_0)\in[F-I_A,\;F-I_A+Q_{\max}].
\]

The continuum midpoint defect \(\Delta=I-r_0s(r_0)\) also has the enclosure
\[
\Delta_A-Q_{\max}\le\Delta\le\Delta_A+r_0q(r_0)A_{\max},
\qquad \Delta_A=I_A-r_0A(r_0).
\]

These are analytic inequalities conditional on the declared uncharged Gaussian source and \(K_t<0\); optional charged-dust/Maxwell terms are excluded. The reported endpoints are evaluated in float64, **not** with directed-rounding interval arithmetic, so this is not yet a formally machine-certified floating-point proof. It does, however, bound the source-specific estimator error independently of the Richardson estimate.

### Analytic enclosure results

| Resolution | Integral lower bound | Integral upper bound | \(Q_{\max}\) | Midpoint-defect lower bound | Midpoint-defect upper bound |
|---:|---:|---:|---:|---:|---:|
| 160 | 0.24999346461512656 | 0.24999346502401773 | \(4.089\times10^{-10}\) | \(1.304954388\times10^{-5}\) | \(1.305199202\times10^{-5}\) |
| 240 | 0.16666472867678334 | 0.16666472873062910 | \(5.385\times10^{-11}\) | \(3.873181923\times10^{-6}\) | \(3.873504693\times10^{-6}\) |
| 320 | 0.12499918217712129 | 0.12499918218989914 | \(1.278\times10^{-11}\) | \(1.634981594\times10^{-6}\) | \(1.635058221\times10^{-6}\) |

For the selected estimator family, the maximum absolute error allowed by the analytic enclosure is:

| Resolution | Degree 3, unanchored, window 2.5 | Degree 3, origin-anchored, window 2.5 | Degree 4, window 3.0 | Degree 5, window 3.5 |
|---:|---:|---:|---:|---:|
| 160 | \(7.206\times10^{-9}\) | \(9.905\times10^{-10}\) | \(7.956\times10^{-10}\) | \(4.331\times10^{-10}\) |
| 240 | \(1.689\times10^{-8}\) | \(6.911\times10^{-10}\) | \(7.865\times10^{-10}\) | \(1.727\times10^{-10}\) |
| 320 | \(8.911\times10^{-9}\) | \(2.209\times10^{-10}\) | \(7.421\times10^{-10}\) | \(8.171\times10^{-11}\) |

The signed lower/upper fit-error intervals are written to `fit_family.csv`. This is more informative than treating a tiny in-sample residual as an error bound: for this declared source and these windows, it yields a model-conditional upper bound on the actual integral-estimator error. It does not establish that this continuous extension is the unique physical model, nor bound evolved-trajectory errors.

## Verified run and provenance

The successful source-pinned audit run is [GitHub Actions run 38792318422](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38792318422); its [artifact 12000553218](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38792318422/artifacts/12000553218) contains the JSON report, moment CSV, fit-family CSV (including leave-one-out diagnostics and analytic error intervals), and SHA-256 manifest. All eight dedicated analytic-oracle/regression tests passed on that head.

- Engine source anchor: 285a2ee0c1a3b97cddeb4519d8e749109b115aab.
- Vendor commit: d6052d605673ce9d82cdc99b0df79855fa2ea215.
- Audited research head: 397f920fc1b3238515717d1067bede5fc2cc48a2 (audit checkout head: 52836b9b4235ee61dc0eca75e6d6e63c916802be).
- Artifact ZIP SHA-256: 151fa36698c50333a810cf3a3620aa51faa654ec4cb98287db09cdd50d347c2c.
- Report JSON SHA-256: c9f916205fde3297f0cc44bc57b8dc59d634cc89868877e3fef03e3d8c1bf204.
- Moment CSV SHA-256: 9c6e07948c0e7c5e2d14c73a229d6b437dc862dfe5746ead5ad06dce0ade22a3.
- Fit-family CSV SHA-256: f6cd7910331b03626c38e4203d6e2a0e6912f4a39d4255e8f8ea327036331ee1.

The run verified all declared engine-source blobs and the vendor pin, required an empty engine/vendor/dependency diff, reproduced the initializer's six-pass B profile against \(X^6\), checked the Richardson reference against the analytic integral enclosure and every eligible fit-error interval, confirmed finite output and unchanged initialized-state digests, and took zero evolution steps.

## Numerical results

The continuous quadrature delta is \(I_{\mathrm{cont}}-r_0s_{\mathrm{cont}}(r_0)\), using the continuum fixed point at the first cell centre \(r_0=\Delta r/2\). The sampled-source delta compares this continuous-source midpoint moment with the point-sampled moment used by the production initializer.

| Resolution | \(\Delta r\) | Continuous half-cell integral | Production point-sample moment | Continuous quadrature delta | Continuum-vs-sampled-source delta | Max source replay error |
|---:|---:|---:|---:|---:|---:|---:|
| 160 | 0.5000 | 0.24999346461587835 | 0.24998041303237156 | \(+1.3051583823\times10^{-5}\) | \(-3.1664\times10^{-13}\) | \(5.55\times10^{-16}\) |
| 240 | 0.3333… | 0.16666472867682744 | 0.16666085522595117 | \(+3.8734508880\times10^{-6}\) | \(-1.1741\times10^{-14}\) | \(1.22\times10^{-15}\) |
| 320 | 0.2500 | 0.12499918217712812 | 0.12499754713168003 | \(+1.6350454493\times10^{-6}\) | \(-1.2490\times10^{-15}\) | \(6.66\times10^{-16}\) |

The decomposition closure residual is zero at reported precision for all three resolutions. The sampled-source delta is tiny compared with the quadrature delta, isolating the dominant term under this declared continuous extension. The correction values scale with observed orders 2.996 and 2.998 between successive resolutions, consistent with the \(O(\Delta r^3)\) midpoint-rule error expected from a regular even source on the origin half-cell. This is a local quadrature scaling result, not a proof of full-solution convergence.

## Fit-family comparison

The signed fit error below is the fitted integral minus the continuous numerical reference. Positive and negative signs are retained.

| Resolution | Degree 3, window 2.5, unanchored | Degree 3, window 2.5, origin-anchored | Degree 4, window 3.0 | Degree 5, window 3.5 |
|---:|---:|---:|---:|---:|
| 160 | \(-6.7974\times10^{-9}\) | \(+9.8975\times10^{-10}\) | \(-3.8740\times10^{-10}\) | \(-2.4922\times10^{-11}\) |
| 240 | \(-1.6835\times10^{-8}\) | \(+6.9098\times10^{-10}\) | \(-7.3264\times10^{-10}\) | \(-1.1887\times10^{-10}\) |
| 320 | \(-8.8977\times10^{-9}\) | \(+2.2086\times10^{-10}\) | \(-7.2923\times10^{-10}\) | \(-6.8936\times10^{-11}\) |

For this source and these sample windows, anchoring the known \(s(0)=1\) improves the degree-3 moment estimate. Higher-degree fits are closer still in the recorded family. This is evidence about the estimator under a specified continuous source; it does not make the highest-degree/window choice a proven universal prescription. In-sample residuals and Richardson changes are diagnostics; unlike those quantities, the analytic inequality above bounds the fitted estimator's integral error under the stated continuous-source assumptions.

### Leave-one-out source prediction

To avoid using in-sample residuals as the only fit-family discriminator, the audit also refits each candidate after withholding one source sample at a time and predicts the held-out value. The table reports the largest absolute held-out source residual; the CSV additionally preserves RMS error, worst-fold condition number, and the minimum residual degrees of freedom in a training fold.

| Resolution | Degree 3, unanchored, window 2.5 | Degree 3, origin-anchored, window 2.5 | Degree 4, window 3.0 | Degree 5, window 3.5 |
|---:|---:|---:|---:|---:|
| 160 | \\(1.775\\times10^{-6}\\) | \\(1.808\\times10^{-6}\\) | \\(3.710\\times10^{-7}\\) | \\(8.906\\times10^{-8}\\) |
| 240 | \\(1.642\\times10^{-6}\\) | \\(1.745\\times10^{-6}\\) | \\(1.784\\times10^{-7}\\) | \\(5.336\\times10^{-8}\\) |
| 320 | \\(8.110\\times10^{-7}\\) | \\(8.722\\times10^{-7}\\) | \\(1.203\\times10^{-7}\\) | \\(2.041\\times10^{-8}\\) |

The held-out errors substantially exceed the in-sample residuals. At N=160, several LOO folds have zero residual degrees of freedom, and the degree-5 fit's worst-fold condition number is approximately \\(4.67\\times10^4\\); its LOO prediction error is still much lower than the degree-3 controls, but is not negligible. This means fit-window policy remains a model-selection question. LOO values are prediction diagnostics for sampled source values—not confidence intervals and not substitutes for the analytic half-cell moment-error enclosure. No fit family is promoted to production on this evidence alone.

## Comparison layers

The report records:
1. The continuously integrated half-cell moment.
2. The continuous-source midpoint moment at the continuum fixed point.
3. The actual production point-sampled moment from a six-pass source replay using the pinned grid derivative stencils.
4. The continuous quadrature delta, continuum-versus-sampled-source delta, and their algebraic sum.
5. Fit-derived moments for an unanchored and origin-anchored degree-3 fit, plus higher-degree/window controls when the data overdetermine the fit.
6. Source replay against \(B=X^6\), state digests before/after observation, source/vendor pins, and SHA-256 checksums.

The comparison is split because the continuum source and sampled discrete source use different derivative representations. The total difference cannot be labelled purely a quadrature error without this decomposition.

## Review gates and interpretation limits

The evidence supports the origin half-cell quadrature hypothesis as the leading local discrepancy under the stated continuous extension. It does not prove the extension is the uniquely intended physical model, establish trajectory convergence, or authorize a production equation change.

Before a production patch is proposed, require independent numerical-method review, an explicit acceptable error criterion, source-specific uncertainty propagation through the self-consistent initializer, signed Ricci/constraint witnesses, and independent review of the isolated patch. The offline Lambda derivative counterfactual remains a separate evidence layer and is not merged into the quadrature claim.

## Reproduction

The specialized workflow runs the continuous-source estimator tests, validates source/vendor pins, compiles the audit, and runs the zero-evolution audit at resolutions 160, 240, and 320 on \(r_{\max}=80\). The report and CSV files are written beneath results/continuous-origin-half-cell-oracle-2026-10-15/ and uploaded as the linked artifact.
