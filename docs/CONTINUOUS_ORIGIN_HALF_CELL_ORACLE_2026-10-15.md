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

The audit solves this regular initial-value problem with an implicit trapezoid update and successive halvings of the substep size. Richardson extrapolation estimates the limiting integral. The reported numerical uncertainty is the greater of the final extrapolant-change estimate and an eight-epsilon floating-point floor. Neither is a rigorous analytic remainder bound.

## Verified run and provenance

The successful source-pinned audit run is [GitHub Actions run 38790339132](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38790339132); its [artifact 12000470732](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38790339132/artifacts/12000470732) contains the JSON report, moment CSV, fit-family CSV, and SHA-256 manifest.

- Engine source anchor: 285a2ee0c1a3b97cddeb4519d8e749109b115aab.
- Vendor commit: d6052d605673ce9d82cdc99b0df79855fa2ea215.
- Audited code head: 7f4c32ac9718e985593d166fa2f1556708528adb.
- Artifact ZIP SHA-256: 3c261e50c02d753d026ef17517d27fe69452e75ca37d1370186c351cf25efcaa.
- Report JSON SHA-256: ce9bc124f6401fa7d75c7ffabf30bb9ac48fe2d8c4b0e91bf0a619adf4712ab6.
- Moment CSV SHA-256: fea5243f4e28642c0fc935fd772e9155c4c7be21cd11f8e696c2aa4ea9c327f5.
- Fit-family CSV SHA-256: 636d69f2f1756cb89e0463b45ed51907bee59df067c014cef2a0ab6a43c73841.

The run verified all declared engine-source blobs and the vendor pin, required an empty engine/vendor/dependency diff, reproduced the initializer's six-pass B profile against \(X^6\), confirmed finite output and unchanged initialized-state digests, and took zero evolution steps.

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

For this source and these sample windows, anchoring the known \(s(0)=1\) improves the degree-3 moment estimate. Higher-degree fits are closer still in the recorded family. This is evidence about the estimator under a specified continuous source; it does not make the highest-degree/window choice a proven universal prescription. In-sample residuals and Richardson changes are reported for diagnostics, not promoted to rigorous error bounds.

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
