# Offline centre-Ricci derivative-channel decomposition

**Status:** research-only; diagnostic decomposition, not a production fix or scientific admission.

## Provenance

- Engine source anchor: 285a2ee0c1a3b97cddeb4519d8e749109b115aab
- Pinned vendor commit: d6052d605673ce9d82cdc99b0df79855fa2ea215
- Vendor operator Git blob SHA-1: `r0_operators.py` — `7f92fd9857981561dc2dff568addae87790c989f` (Git object identifier, not a SHA-256 file digest).
- Source-pinned PR #249 stage-artifact ZIP hashes: N=40 1cd1ea6a07ec59477729b3724df12cab2414e17eb9777efb3f9fdb264878429c; N=80 7f64d7fd15347c25924f83381133a6df4347c692239299ac1f58c2827f43ca0c.
- Source-pinned PR #250 accepted-step ZIP hash: 2b6486b6159af3a272cb11f2351c9079e835946440726de983ba0cf67011e333.
- All three ZIP hashes and all embedded report/checkpoint/NPZ/CSV SHA-256 values were checked against the extracted files. No new production step was run.

## Method

On each saved N=40/N=80 stage, apply the pinned parity-odd discrete derivative to the saved Lambda_evolved and Lambda_metric_reconstructed profiles. Independently reconstruct Lambda_metric from the same D1(a), D1(b), then evaluate its analytic product/chain-rule derivative using the same-grid D2(a), D2(b). Split the derivative gap as

\[
D_4\Lambda_{\rm evolved}-\partial_r\Lambda_{\rm metric}
=(D_4\Lambda_{\rm evolved}-D_4\Lambda_{\rm metric})
+(D_4\Lambda_{\rm metric}-\partial_r\Lambda_{\rm metric}|_{\rm chain}).
\]

The pinned vendor TF-Ricci dependence is linear in the D4-Lambda operand with coefficient (2/3) X^2; the algebraic term 0.5 X^2 (D1(a)/a) Lambda is independently retained. The full vendor call on a metric-Lambda copy changes both terms, so it is reported separately from the derivative-only counterfactual.

## Result: centre and fixed-annulus examples

| Resolution / stage | TF response from evolved-vs-metric D4 (cell 0) | TF response from metric-D4-vs-chain-rule (cell 0) | Annulus max, connection channel | Annulus max, stencil/chain channel |
|---|---:|---:|---:|---:|
| N=40 primary_predictor | 0.000000e+00 | 1.655012e-03 | 0.000000e+00 | 3.497183e-05 |
| N=40 final_u_old_v | 2.084627e-06 | 1.641369e-03 | 5.608782e-07 | 3.450637e-05 |
| N=40 final_state | 1.841446e-06 | 1.641369e-03 | 6.532852e-07 | 3.450637e-05 |
| N=80 primary_predictor | 0.000000e+00 | 1.601732e-03 | 0.000000e+00 | 1.110025e-06 |
| N=80 final_u_old_v | 2.191107e-06 | 1.588448e-03 | 3.052304e-08 | 1.110108e-06 |
| N=80 final_state | 1.748457e-06 | 1.588448e-03 | 7.659308e-10 | 1.110108e-06 |

Here “connection channel” means the TF-Ricci response to D4(Lambda_evolved)-D4(Lambda_metric), while “stencil/chain channel” means the response to D4(Lambda_metric)-the analytic chain-rule derivative. The two signed channels add to the total evolved-D4 minus analytic-chain-rule derivative response.

## Three-term signed budget for the metric-Lambda channel

A follow-up pass decomposes the second (stencil/chain-rule) channel into the three terms used to reconstruct

\[
\Lambda_{\rm metric}=
\underbrace{\frac{D_1 a}{2a^2}}_{\lambda_a}
\;\underbrace{-\frac{D_1 b}{ab}}_{\lambda_b}
\;+\underbrace{\frac{2}{r}\left(\frac1b-\frac1a\right)}_{\lambda_r}.
\]

For each term the report records the signed derivative-only TF-Ricci response
\((2/3)X^2[D_4(\lambda_i)-\partial_r\lambda_i|_{D_1/D_2}]\).
The three pieces are summed before taking any norm; their sum is checked against the unsplit
metric-Lambda channel on every saved stage.

### Primary-predictor centre values

| Resolution | Total metric-Lambda channel | \(\lambda_a\) response | \(\lambda_b\) response | \(\lambda_r\) response | Recombination max error |
|---:|---:|---:|---:|---:|---:|
| N=40 | \(+1.655011946\times10^{-3}\) | \(-8.043721740\times10^{-6}\) | \(-8.055575911\times10^{-6}\) | \(+1.671111244\times10^{-3}\) | \(2.17\times10^{-19}\) |
| N=80 | \(+1.601732028\times10^{-3}\) | \(-9.708243950\times10^{-7}\) | \(-9.792177506\times10^{-7}\) | \(+1.603682070\times10^{-3}\) | \(2.17\times10^{-19}\) |

### Fixed-annulus maxima

| Resolution | \(\max_{5\le r\le10}|\lambda_a\text{ response}|\) | \(\max_{5\le r\le10}|\lambda_b\text{ response}|\) | \(\max_{5\le r\le10}|\lambda_r\text{ response}|\) |
|---:|---:|---:|---:|
| N=40 | \(7.1530\times10^{-6}\) | \(7.1635\times10^{-6}\) | \(2.0655\times10^{-5}\) |
| N=80 | \(4.7444\times10^{-7}\) | \(4.6991\times10^{-7}\) | \(2.2262\times10^{-7}\) |

On the first cell of these two primary-predictor states, the radial-geometry reconstruction term dominates the signed chain-rule difference; the two derivative-product terms contribute much smaller negative offsets. The fixed-annulus maxima show that this is not a claim of radial-term dominance at every radius or resolution. The diagnostic now emits all three signed profiles and rejects a term-sum closure error above the same \(2\times10^{-12}\) absolute gate.

**Interpretation limit:** this locates the largest first-cell component in the discrete representation comparison; it does not prove that the radial term or either derivative representation is defective. The next numerical-method review should compare the parity/regular-origin expansion of \(\delta=1/b-1/a\) with the sampled profile and the pinned stencil, and assess whether the stage profile resolves that expansion. No production equation or operator is changed.

## Accepted-step connection-derivative history (PR #250 saved artifact)

The accepted-step artifact contains Lambda_evolved, Lambda_metric, and the sampled connection residual, but does not include a, b, or X. It therefore supports the exact parity-odd D4 derivative split, but not the analytic chain-rule or X^2 Ricci-response split at each later step.

| Accepted step | Time | D4(Lambda_evolved-Lambda_metric) cell 0 | Cell 1 | Max abs, first two cells | Max abs, fixed 5<=r<=10 |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 |
| 1 | 0.06 | 2.765190e-06 | -3.589692e-06 | 3.589692e-06 | 9.825768e-07 |
| 2 | 0.12 | 1.098457e-05 | -1.426672e-05 | 1.426672e-05 | 3.914362e-06 |
| 3 | 0.18 | 2.451200e-05 | -3.185266e-05 | 3.185266e-05 | 8.763505e-06 |
| 4 | 0.24 | 4.311432e-05 | -5.606336e-05 | 5.606336e-05 | 1.548014e-05 |
| 5 | 0.30 | 6.647599e-05 | -8.651450e-05 | 8.651450e-05 | 2.399717e-05 |
| 6 | 0.36 | 9.420502e-05 | -1.227287e-04 | 1.227287e-04 | 3.423123e-05 |
| 7 | 0.42 | 1.258402e-04 | -1.641441e-04 | 1.641441e-04 | 4.608390e-05 |
| 8 | 0.48 | 1.608590e-04 | -2.101241e-04 | 2.101241e-04 | 5.944301e-05 |

At each saved step, applying the same linear D4 operator to the stored residual agrees with D4(Lambda_evolved)-D4(Lambda_metric) to at most 3.44e-19 (roundoff). This is a consistency identity, not a causal attribution of the full Ricci gap.

## Interpretation / limits

- On these saved states, the derivative-only centre TF-Ricci response of the analytic-chain-rule versus discrete-D4 channel is about 1.59–1.64e-3, roughly three orders of magnitude larger than the evolved-versus-metric connection-derivative response at the captured completed state (about 1.75–1.84e-6). This localizes strong sensitivity to the discrete-stencil/product-chain-rule representation difference, but does not uniquely explain the total vendor-vs-direct Ricci gap or establish which representation should be used.
- At the captured predictor, evolved and metric Lambda are equal on the sampled data, so the connection channel is exactly zero there; it becomes small but nonzero at the final mixed/final state.
- In the fixed coordinate annulus 5<=r<=10, the D4-chain channel is much smaller at N=80 than N=40, and the connection-residual derivative channel is also smaller. Two resolutions are not enough to infer a formal convergence order.
- Reconstructed-Lambda full-vendor response includes the changed algebraic Lambda*D1(a) term and is not the same as the derivative-only response. The full pinned vendor TF-Ricci response identity was checked directly; its maximum closure residual over all six stages was below 5.2e-19. The reconstructed metric Lambda exactly matches the saved stage profile, the recalculated evolved-Lambda TF-Ricci matches the stored profile, and the old chain-rule response matches the prior artifact to below 4.4e-19.
- PR #250 accepted-step history verifies D4(Lambda_evolved)-D4(Lambda_metric) = D4(Lambda_evolved-Lambda_metric) from saved values through eight steps; its artifact lacks a, b, X, so it does not independently determine the analytic-chain-rule channel or X-weighted response at later accepted steps. No new evolution was run.
- The discrepancy remains diagnostic-only. Do not change production equations or claim scientific admission without independent methodological review and wider controlled convergence / invariant evidence.
