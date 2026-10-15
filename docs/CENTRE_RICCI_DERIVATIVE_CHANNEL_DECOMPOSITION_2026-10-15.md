# Offline centre-Ricci derivative-channel decomposition

**Status:** research-only; diagnostic decomposition, not a production fix or scientific admission.

## Provenance

- Engine source anchor: 285a2ee0c1a3b97cddeb4519d8e749109b115aab
- Pinned vendor commit: d6052d605673ce9d82cdc99b0df79855fa2ea215
- Vendor operator source blob: r0_operators.py SHA 7f92fd9857981561dc2dff568addae87790c989f
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

- In both resolutions, the derivative-only centre response associated with the analytic-chain-rule vs discrete-D4 representation is about 1.6e-3, while the evolved-vs-metric connection-derivative channel at the captured completed stage is about 2e-6. This strongly suggests the leading observed centre discrepancy on these saved states is dominated by the sampled-derivative/product-chain-rule representation gap, with a smaller connection-residual derivative contribution.
- At the captured predictor, evolved and metric Lambda are equal on the sampled data, so the connection channel is exactly zero there; it becomes small but nonzero at the final mixed/final state.
- In the fixed coordinate annulus 5<=r<=10, the D4-chain channel is much smaller at N=80 than N=40, and the connection-residual derivative channel is also smaller. Two resolutions are not enough to infer a formal convergence order.
- Reconstructed-Lambda full-vendor response includes the changed algebraic Lambda*D1(a) term and is not the same as the derivative-only response. The full pinned vendor TF-Ricci response identity was checked directly; its maximum closure residual over all six stages was below 5.2e-19. The reconstructed metric Lambda exactly matches the saved stage profile, the recalculated evolved-Lambda TF-Ricci matches the stored profile, and the old chain-rule response matches the prior artifact to below 4.4e-19.
- PR #250 accepted-step history verifies D4(Lambda_evolved)-D4(Lambda_metric) = D4(Lambda_evolved-Lambda_metric) from saved values through eight steps; its artifact lacks a, b, X, so it does not independently determine the analytic-chain-rule channel or X-weighted response at later accepted steps. No new evolution was run.
- The discrepancy remains diagnostic-only. Do not change production equations or claim scientific admission without independent methodological review and wider controlled convergence / invariant evidence.
