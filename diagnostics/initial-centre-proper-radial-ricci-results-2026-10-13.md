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

## Cross-check against the team's current-main origin-only A/B

This section uses two separate source-pinned team artifacts rather than a duplicate run:

1. Current-main initialized-profile and origin-only A/B: workflow [38342627215](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38342627215), source capture branch head \`ca81723220bae953d0f281ac3e88b1dfafb2e6f5\`, artifact [11787378026](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38342627215/artifacts/11787378026), ZIP SHA-256 \`c9914bca21722c5060343c7a986de1060e3145047c0c13b8c2be38e8c3cad637\`. The A/B report SHA-256 is \`d15762d1aabbd9ecc7b0e31ef4181757e2553978b86df906d16e06794be6e384\`; profile archive SHA-256 is \`ddf77f43d9a0949ae70c6388bb380e64413e5811ee02cab6283bfcac2a2d49d2\`.
2. Current-main PIRK L2 independent metric cross-check: workflow [38342662323](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38342662323), source commit \`c01f2467347e8aebd57d220374fd7bc4e93045bc\`, kernel blob \`2954e7a79656849b5f088053ea0f4ac97d4bb443\`, vendor \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`, artifact [11787118604](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38342662323/artifacts/11787118604), ZIP SHA-256 \`38262677ba5bc501d615b3329782de2e88f48754432e24ee65dfad074b8f1c84\`. Its report SHA-256 is \`7ed52f357de4290827d5c2f3e83ff52a1e30804fc8ff1b2a2bc40b67b67e146c\`.

The first artifact's A/B is deliberately narrower than the broader PR #18 radial-integration change: only the origin half-cell treatment is counterfactually changed; inter-centre trapezoid intervals remain intact. The independent geometry comparison was recomputed from its stored \`a,b,X,r\` profiles with the same pinned grid derivative stencil. A second independent metric-route implementation in the L2 artifact applies \(qD_r(qD_rR)\), \(q=X/\sqrt a\), rather than the \(D_{rr}R-D_rR\,D_rf/f\) sequence used above. Both recover the same localization.

| \(N\) | Origin treatment | Vendor TF Ricci, cell 0 | Proper-radial TF Ricci, cell 0 | Vendor minus direct |
|---:|---|---:|---:|---:|
| 60 | Production half-cell trapezoid | +1.04291e-3 | -5.71662e-4 | +1.61457e-3 |
| 60 | Even-series half-cell only | +1.67635e-5 | -1.72906e-5 | +3.40541e-5 |
| 80 | Production half-cell trapezoid | +1.03688e-3 | -5.65396e-4 | +1.60228e-3 |
| 80 | Even-series half-cell only | +3.22967e-7 | -4.50507e-6 | +4.82804e-6 |
| 120 | Production half-cell trapezoid | +1.03303e-3 | -5.61366e-4 | +1.59439e-3 |
| 120 | Even-series half-cell only | -3.10710e-6 | -9.84612e-8 | -3.00864e-6 |

For the alternate metric-route sequence the post-counterfactual signed gaps are \(+3.41290\times10^{-5}\), \(+4.85652\times10^{-6}\), and \(-3.00225\times10^{-6}\) at \(N=60,80,120\). The close agreement between two derivative sequences is a useful robustness check: the origin-only change lowers the first-cell discrepancy by roughly factors \(47\), \(330\), and \(530\), respectively. This is substantially stronger evidence for a half-cell/initial-data mechanism than the baseline comparison alone, but still a finite-resolution result, not proof of continuum regularity.

The origin-only A/B also reduces the first-cell Hamiltonian residual from \(3.93161\times10^{-3}\) to \(1.00308\times10^{-4}\) at \(N=60\), from \(3.98859\times10^{-3}\) to \(3.67656\times10^{-5}\) at \(N=80\), and from \(4.02811\times10^{-3}\) to \(9.09006\times10^{-6}\) at \(N=120\). The central full \(\dot A_a\) witness falls from \(1.09159\times10^{-3}\), \(1.08913\times10^{-3}\), \(1.08795\times10^{-3}\) to \(5.18390\times10^{-5}\), \(3.83681\times10^{-5}\), \(3.72922\times10^{-5}\). Momentum, connection, and determinant witnesses remain at roundoff in both branches of this *initialization-only* A/B.

These findings narrow, but do not close, attribution: the production origin half-cell treatment is enough to reproduce most of the centre-layer anomaly on the initialized slice. The separate PR #18 four-arm trajectory still warns that the broader quadrature change (including inter-centre interval integration) worsens the Misner–Sharp current residual. The origin-only A/B has not advanced time, so it cannot establish conservation behavior under evolution. The next discriminating candidate is therefore an explicitly isolated **origin-half-cell-only** matched-time run, compared against baseline with the existing fixed-coordinate/worldtube flux-work ledger; do not conflate that result with the broader all-interval quadrature. Check the concurrent centre-series branch before launching, since its workflow may be rerunning after a harness-only change.

## Interpretation and next discriminator

The supported conclusion is an **inner-centre discrete compatibility failure** on the actual initialized profiles. The analytic constant-curvature control argues against a universal sign or formula mistake in the independent reconstruction. The fixed-radius comparison argues against a general, domain-wide disagreement between the Ricci representations. What remains unresolved is which centre-local layer is responsible: the initializer's origin-to-first-cell radial-source quadrature, the parity/near-origin stencil closure acting on that profile, or their interaction with the BSSN connection formulation.

Do not patch production physics, replace the Ricci equation, tune a coefficient, or claim scientific admission from this result. Next, independently audit the first half-cell integration used to construct \(B(r)=r^{-1}\int_0^r s(x)\,dx\) against the regular even Taylor series of the actual source \(s(r)\), then compare the implied leading coefficients of \(f(r)\) and the signed radial/tangential Ricci on the same fixed source profiles. Any candidate change should first be evaluated as an offline counterfactual and checked against the already-open initial-quadrature/worldtube balance audits, not promoted from this centre-only result.
