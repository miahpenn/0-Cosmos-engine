# Continuous-oracle half-cell correction: Ricci response — 15 October 2026

## Status and provenance

**Research-only. Scientific admission remains NOT ADMITTED.** This is a zero-evolution frozen-source metric/Ricci counterfactual driven by the continuous half-cell moment, not a production-initializer change.

The dedicated audit passed on code head \`6d364655a0f9e7f0edd38c6f9e5d1c3e9a243a8f\`: [workflow run 38791262954](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38791262954), artifact [12000266806](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38791262954/artifacts/12000266806). Artifact ZIP SHA-256 is \`99f8315d020af1c70cc306beebb96ce32392b2c47d502fa0660bd9a198d872c9\`; report SHA-256 is \`a73541328be62805bdcc8913ea68e4b62b9bcb41a5e520c9c18c6088db7e4df4\`.

All eleven required engine-source blob pins match the source anchor \`285a2ee0c1a3b97cddeb4519d8e749109b115aab\`; vendor pin \`d6052d605673ce9d82cdc99b0df79855fa2ea215\` matches and the engine/vendor/dependency diff is empty. The initialized states are unchanged, profiles finite, and no evolution is performed.

The continuous half-cell integral was independently solved by DOP853 plus adaptive quadrature and cross-checked against the implicit-trapezoid/Richardson oracle on [PR #259](https://github.com/miahpenn/0-Cosmos-engine/pull/259). Integrals agree to \(\le9.5\times10^{-16}\). The ODE/quadrature residual is at most \(2.8\times10^{-17}\), and the regular \(s_2/s_4/s_6\) series moment residual at most \(5.3\times10^{-15}\). This is a high-accuracy numerical cross-check, not a rigorous interval enclosure.

## Counterfactual

For each resolution, define the continuous-oracle half-cell correction relative to the production point sample as

\[
\Delta I_{\rm cont}=I_{\rm cont}(r_0)-r_0\,s_{\rm production}(r_0),\qquad r_0=\Delta r/2.
\]

Apply it once to the same pinned production metric as \(\Delta B(r)=\Delta I_{\rm cont}/r\). The source samples and later inter-centre trapezoid intervals are held fixed. The copied geometry's metric identities, Lambda, pinned vendor trace-free Ricci and independent proper-radial Ricci are recomputed. This isolates the direct geometric effect of the best available numerical moment; it is not a self-consistent discrete initializer re-solve.

## Results

| Resolution | Continuous integral \(I_{\rm cont}\) | \(\Delta I_{\rm cont}\) | Production cell-0 signed Ricci gap | Oracle-corrected cell-0 signed gap | Gap magnitude reduction from production baseline |
|---:|---:|---:|---:|---:|---:|
| \(N=160\) | 0.24999346461587850 | \(+1.30515835\times10^{-5}\) | \(+1.59181308\times10^{-3}\) | \(-1.88494083\times10^{-4}\) | \(8.44\times\) |
| \(N=240\) | 0.16666472867682738 | \(+3.87345088\times10^{-6}\) | \(+1.59002921\times10^{-3}\) | \(-1.93165363\times10^{-4}\) | \(8.23\times\) |
| \(N=320\) | 0.12499918217712717 | \(+1.63504545\times10^{-6}\) | \(+1.58941644\times10^{-3}\) | \(-1.94789646\times10^{-4}\) | \(8.16\times\) |

Thus the continuous correction removes most of the original positive centre discrepancy, but leaves a sign-reversed residual of order \(10^{-4}\). It does not close the signed Ricci mismatch.

### Exact continuous moment vs fitted corrections

| Resolution | Fixed unanchored fit cell-0 gap | Continuous-oracle moment cell-0 gap | Oracle correction's increase in gap magnitude vs unanchored fit |
|---:|---:|---:|---:|
| \(N=160\) | \(-1.87566879\times10^{-4}\) | \(-1.88494083\times10^{-4}\) | \(+0.494\%\) |
| \(N=240\) | \(-1.85415084\times10^{-4}\) | \(-1.93165363\times10^{-4}\) | \(+4.180\%\) |
| \(N=320\) | \(-1.85080224\times10^{-4}\) | \(-1.94789646\times10^{-4}\) | \(+5.246\%\) |

This is the key discriminator. The source-specific continuous integral agrees with the independent Richardson reference, and the fit errors are already small in absolute moment units. Yet applying the *continuous-oracle* moment to the production geometry still worsens the remaining centre-Ricci residual relative to the unconstrained fit, especially at finer resolutions. The geometric mismatch is therefore not explained by half-cell quadrature error alone.

The fixed physical annulus \(5\le r\le10\) is essentially unchanged: the maximum-gap difference between the oracle-corrected arm and the unconstrained fit arm is approximately \(3.8\times10^{-14}\), \(1.7\times10^{-14}\), and \(1.2\times10^{-14}\) at \(N=160/240/320\). This response is tightly centre-localized.

## Interpretation and next gate

**Supported:** the production point-sampled half-cell moment is missing the expected positive midpoint-rule correction of \(O(\Delta r^3)\); a high-accuracy continuous estimate agrees with an independent Richardson reference; applying that correction reduces the original centre Ricci discrepancy substantially.

**Not supported:** claiming that the half-cell quadrature correction solves the Ricci discrepancy or that the exact continuous integral gives the best discrete Ricci/geometry closure. The oracle-corrected geometry retains a signed residual near \(-1.9\times10^{-4}\), despite an accurate moment estimate.

The next target should be the **discrete geometry/operator compatibility**: evaluate the independently reconstructed warped-product Ricci and the pinned vendor trace-free Ricci on the same corrected metric, then split the remaining cell-0/cell-1 residual by the Lambda derivative/chain-rule channel and the remaining geometric derivative/algebraic terms. Keep the integral estimate fixed during that decomposition, so quadrature is no longer a moving variable.

No production correction is proposed. Scientific admission remains **NOT ADMITTED**.
