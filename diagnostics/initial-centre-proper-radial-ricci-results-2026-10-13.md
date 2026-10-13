# Initial-centre proper-radial Ricci comparison — 13 October 2026

**Status:** completed diagnostic-only; scientific admission remains false.  
**PR:** [#74](https://github.com/miahpenn/0-Cosmos-engine/pull/74)  
**Exact-head workflow:** [38342420905](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38342420905)  
**Artifact:** [11787675751](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38342420905/artifacts/11787675751)

## Question and scope

The signed current-source Hamiltonian/momentum-constraint audits and PIRK L2 stage split localized a persistent first-cell \(\dot A_a\) to the trace-free Ricci term, but did not determine whether the underlying cause is the initialized centre profile, the discrete centre operators, or the BSSN Ricci representation. A manufactured smooth-metric control showed that the vendor Ricci expression can reproduce regular \(O(r^2)\) trace-free curvature on a compatible profile. This comparison tests the actual initialized profiles against a separate spherical-geometric reconstruction.

This is a **no-evolution** comparison at \(N=60,80,120\), \(r_{\max}=80\), with default scalar amplitude/width, \(D_{\rm amplitude}=10^{-10}\), and radiation enabled. Initialization includes the ordinary CMC lapse solve, but no production time step is taken. The observer verifies pre/post state digests and reports signed profiles and first-six-cell values.

## Reproducibility and guardrails

- Production source baseline: \`29d828b4d72b017c12a20b6dc90549f2511f025c\`.
- Vendor submodule: \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.
- Six production-file blobs (kernel, initializer, PIRK adapter, CMC gauge, matter and stress-energy assembly) are checked exactly against the source baseline.
- The workflow fails closed on engine/vendor/submodule/dependency drift; no production file, vendor file, physical source, equation, coefficient, gauge, boundary, initializer, or admission rule changes in this PR.
- Exact-head software result: **233 tests passed**; provenance checks passed; all three observer states stayed unchanged and finite.
- Artifact ZIP SHA-256: \`bf6cf98c99b3a9f3eb434c1ff684302ca5a897f0ece405bbf98cbc7058a16316\`.
- Verified internal artifact hashes:
  - \`report.json\`: \`b8b4a85e227f756ce9e04cd3271306bf49af90bcc30d7e7732b6856c8c7ae226\`
  - \`profiles.npz\`: \`5e0d8451deaabed1445ec769f0252f2c94c607b66bebac5bee85af11f75e8b3b\`
  - \`profiles.csv\`: \`9d0bdbac7769442159cd77a703035218e023190ec0e22258b95bfc8712637c20\`

## Independent geometric reconstruction

The physical spatial metric can be written
\[
d\ell^2 = f(r)^2\,dr^2 + R(r)^2 d\Omega^2,\quad
R(r)=\frac{r\sqrt{b}}{X},\quad f(r)=\frac{\sqrt a}{X}.
\]
With derivatives with respect to proper radius \(s\), \(d s=f\,dr\),
\[
R_s=\frac{R_r}{f},\quad
R_{ss}=\frac{R_{rr}-R_r f_r/f}{f^2},
\]
\[
{}^{(3)}R^s{}_s=-\frac{2R_{ss}}{R},\quad
{}^{(3)}R^\theta{}_\theta=\frac{1-R_s^2-RR_{ss}}{R^2}.
\]
The corresponding trace-free radial Ricci is
\[
{}^{(3)}R^s{}_s-\frac{{}^{(3)}R}{3}
=\frac23\left({}^{(3)}R^s{}_s-{}^{(3)}R^\theta{}_\theta\right).
\]
This evaluates the geometry through areal radius and proper radial distance, without using evolved \(\Lambda\) or the vendor's BSSN Ricci brackets. It **does share the pinned radial D1/D2 grid operators**, so it is an independent representation, not an independent discretization.

The regression suite now includes an analytic constant-positive-curvature areal slice, \(R=r\), \(f=(1-kr^2)^{-1/2}\), with \(k=10^{-3}\). On \(r\le2\), both the vendor and proper-radial trace-free Ricci, and their difference, are bounded below \(5\times10^{-6}\). This checks the geometric identities/sign convention on an analytic compatible metric; it does not establish correctness on the production initial profile.

## Signed result at the first cell

The values below are the unscaled spatial trace-free Ricci \(R^r{}_r-R/3\). The lapse factor is near unity and does not explain the sign difference.

| \(N\) | \(r_0\) | Vendor BSSN TF Ricci | Proper-radial TF Ricci | Vendor minus direct |
|---:|---:|---:|---:|---:|
| 60 | 0.666667 | +1.04290558e-3 | -5.71662047e-4 | +1.61456763e-3 |
| 80 | 0.500000 | +1.03688498e-3 | -5.65396364e-4 | +1.60228135e-3 |
| 120 | 0.333333 | +1.03302627e-3 | -5.61366204e-4 | +1.59439247e-3 |

Across the first six cells, the maximum absolute difference divided by the larger of the two representations' maximum absolute values is 1.5481, 1.5453, and 1.5434 at \(N=60,80,120\). The first-cell mismatch stays near \(1.6\times10^{-3}\) while \(r_0\) halves from \(2/3\) to \(1/3\), rather than displaying the \(O(r^2)\) scaling expected of a smooth trace-free-centre series. The direct result itself also remains nonzero in this first-cell sequence; this is **not** evidence to patch the vendor formula alone.

## Localization and identities

Linear interpolation of the full CSV profiles to fixed coordinate radii gives the following signed vendor-minus-direct differences:

| Radius | \(N=60\) | \(N=80\) | \(N=120\) |
|---:|---:|---:|---:|
| \(r=2\) | -2.24234e-4 | -1.04684e-4 | +8.51923e-6 |
| \(r=5\) | +5.48024e-6 | +1.69491e-6 | +3.40731e-7 |
| \(r=8\) | -1.76295e-6 | -6.91733e-7 | -1.56592e-7 |

Thus the disagreement is strongly localized to the inner-centre layer; at fixed \(r=5\) and \(r=8\), the two representations approach one another rapidly under refinement. The CSV contains the full signed profiles for independent reconstruction of these interpolated samples.

The initialized metric identities close to roundoff:
- \(\max|ab^2-1|\le3.34\times10^{-16}\);
- \(\max|b-X^2|\le2.23\times10^{-16}\);
- \(\max|a-X^{-4}|\le4.45\times10^{-16}\);
- the discrete initialized \(\Lambda\)-from-metric identity residual is exactly zero under the same formula;
- \(\max|R-r|\le1.43\times10^{-14}\).

These checks rule out a simple failure of the algebraic conformal metric identities but do not prove centre differentiability/regularity of the initialized sequence, nor do they force two different discrete curvature constructions to obey a continuum product rule.

## Source-level origin half-cell check

The current production initializer's first cumulative-integral entry is \`I_0=(dr/2)\,s_0^{\mathrm{sample}}\), while the first cell centre is \(r_0=dr/2\), and it then forms \(B_0=I_0/r_0\). Thus, for a fixed smooth even source profile
\[
s(r)=s_0+s_2 r^2+s_4 r^4+O(r^6),
\]
the current first-cell rule gives \(B_0=s(r_0)=s_0+s_2 r_0^2+s_4 r_0^4+\cdots\). The actual origin-to-centre average would be
\[
\frac1{r_0}\int_0^{r_0}s(r)\,dr
=s_0+\frac{s_2}{3}r_0^2+\frac{s_4}{5}r_0^4+O(r_0^6).
\]
For a fixed source, the leading first-cell differences are therefore \(+(2/3)s_2 r_0^2\) and \(+(4/5)s_4 r_0^4\). This is a source-derived explanation for why a first-cell profile can carry incorrect centre-series coefficients even while \(a,b,X\) satisfy their algebraic identities. In production, \(s\) is itself coupled to the iterated \(B\), so this identity localizes a plausible numerical mechanism but is not a full nonlinear error budget.

This is **not** approval to switch to the candidate quadrature. Draft PR #18 already has analytic even-series origin weights and four-arm matched-time evidence, but its quadrature arm worsens the production-defined Misner–Sharp current residual at the tested resolutions. The right next experiment is to apply this no-evolution, same-profile Ricci observer to baseline and the existing quadrature-only candidate source profiles (with explicit allowed source diffs) and read that against PR #18's signed centre/flux/work balance findings. Avoid another long trajectory and do not fit quadrature weights to a residual.

## Interpretation and next discriminator

The supported conclusion is an **inner-centre discrete compatibility failure** on the actual initialized profiles. The analytic constant-curvature control argues against a universal sign or formula mistake in the independent reconstruction. The fixed-radius comparison argues against a general, domain-wide disagreement between the Ricci representations. What remains unresolved is which centre-local layer is responsible: the initializer's origin-to-first-cell radial-source quadrature, the parity/near-origin stencil closure acting on that profile, or their interaction with the BSSN connection formulation.

Do not patch production physics, replace the Ricci equation, tune a coefficient, or claim scientific admission from this result. Next, independently audit the first half-cell integration used to construct \(B(r)=r^{-1}\int_0^r s(x)\,dx\) against the regular even Taylor series of the actual source \(s(r)\), then compare the implied leading coefficients of \(f(r)\) and the signed radial/tangential Ricci on the same fixed source profiles. Any candidate change should first be evaluated as an offline counterfactual and checked against the already-open initial-quadrature/worldtube balance audits, not promoted from this centre-only result.
