# Analytic half-cell injection into the self-consistent copied initializer — 15 October 2026

**Disposition: research-only; scientific admission NOT ADMITTED.** This experiment tests whether replacing the polynomial-fitted origin half-cell moment with the independent continuous-source integral changes the initialized centre-Ricci mismatch, while reporting the signed constraints and a separate Lambda derivative counterfactual.

## Exact-head validation and provenance

The source-pinned counterfactual and manufactured Lambda symbol sweep both passed in [GitHub Actions run 38796154469](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38796154469) on code head afe359e9322f516666214b5817ca325c76431560 (audit checkout head 68f48423ec646fc3a50db767b97c4a40b0c934eb). All four focused tests passed; the lower/nominal/upper self-consistent half-cell arms, derivative-mode sweep, input-state hash checks and finite/response-closure gates succeeded; no evolution steps were taken.

The [analytic half-cell/constraint artifact 12001559466](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38796154469/artifacts/12001559466) includes signed Ricci/constraint witnesses and all profiles. ZIP SHA-256: 4aa79998226efbec9b90582212c0fe948b6ae06224287fb1681cec33d183f7e5; report.json SHA-256: 109069a2703e92af0277747191d85a00008da51f59709eeedc4aa59b35325cf8; signed_witnesses.csv SHA-256: d26712f1a8b1e115babf965f3eaeba02fb19b58017ddb803c95ce7d7a633af14; profiles.npz SHA-256: 9db4febab5da0c680a3fd4db3f8cb88335a1803f1ec9d9f43055eeb103c7f64a.

The separate [manufactured Lambda symbol artifact 12001124815](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38796154469/artifacts/12001124815) contains the wavelength sweep, full response profiles and checksums. ZIP SHA-256: 7b82033ad349e23db15678a6fb15ecfab679f5e2a4fd37f1c09c1a0f37c32ee6; report.json SHA-256: 771ec4ac1bfa6c665c415f25c625ebb7c8637059217aa03ec9b9164a2fbf8a70; mode_sweep.csv SHA-256: 41d0a489048621746f15f8058f3c1f37334e7a60307574e3a65dc666518de4bb; profiles.npz SHA-256: 0854fbb0d0ea24d3eba8fd97e1917d168272b22debe5500504c5ea37b9449b8d.

- Engine source anchor: 285a2ee0c1a3b97cddeb4519d8e749109b115aab.
- Vendor commit: d6052d605673ce9d82cdc99b0df79855fa2ea215.
- No production engine, vendor, or dependency file changed; the pin check confirmed the production/vendor/dependency diff is empty.

## Method

For the default uncharged Gaussian fields, the continuous source is \(s=A-qB\), \(I=rB\), and \(I'+(q/r)I=A\), with \(I(0)=0\). The independent continuous-source audit provides a Richardson nominal value and the model-conditional analytic interval \(I_A-Q_{\max}\le I(r_0)\le I_A\), where \(r_0=\Delta r/2\).

At each \(N=160,240,320\), the new copied-state experiment injects three half-cell moments: analytic lower endpoint, Richardson nominal, and analytic upper endpoint. It changes only the first half-cell moment. Source and metric are recomputed together on the copied state; all inter-centre trapezoids and sampled outer-grid derivative stencils remain unchanged. It reports both the six-update initializer analogue and the converged fixed point. The polynomial fit is still calculated as a comparator but is not applied to these three arms.

Each candidate metric is tested against an independent proper-radial reconstruction of trace-free Ricci. The report preserves signed cell-0/cell-1 Ricci gaps and native-grid Hamiltonian, momentum, connection and determinant witnesses. A separate Lambda arm changes only the Ricci term induced by comparing \(D_4(\Lambda)\) with the metric-to-Lambda chain-rule derivative on the same candidate metric. That arm is not folded into the primary Ricci gap or the primary constraint witnesses.

## Result 1 — the analytic half-cell change is real, but not the complete correction

The primary Ricci-gap sign is vendor trace-free Ricci minus independently reconstructed proper-radial trace-free Ricci. Values below are signed, not absolute values.

| Resolution | Production cell-0 gap | Analytic nominal cell-0 gap | Absolute-gap reduction | Corrected cell-1 gap | Gap span across analytic integral endpoints |
|---:|---:|---:|---:|---:|---:|
| \(N=160\) | \(+1.5918131\times10^{-3}\) | \(-1.8849409\times10^{-4}\) | \(8.45\times\) | \(+8.4574090\times10^{-5}\) | \(5.578\times10^{-8}\) |
| \(N=240\) | \(+1.5900292\times10^{-3}\) | \(-1.9316537\times10^{-4}\) | \(8.23\times\) | \(+8.5262824\times10^{-5}\) | \(2.479\times10^{-8}\) |
| \(N=320\) | \(+1.5894164\times10^{-3}\) | \(-1.9478965\times10^{-4}\) | \(8.16\times\) | \(+8.5516744\times10^{-5}\) | \(1.394\times10^{-8}\) |

**Interpretation:** replacing the point-sampled/fitted origin contribution with the analytic continuous integral removes most of the original positive cell-0 Ricci gap, but reverses the sign and leaves a robust residual near \(-1.9\times10^{-4}\). The next cell's signed gap becomes positive. This is evidence that the half-cell contribution is important, but it is not a complete solution to the Ricci mismatch.

The remaining Ricci gap is not explained by uncertainty in the continuous integral interval. The interval widths are only \(4.089\times10^{-10}\), \(5.385\times10^{-11}\), and \(1.278\times10^{-11}\); the corresponding cell-0 Ricci-gap spans are below \(5.6\times10^{-8}\), far smaller than the \(\sim1.9\times10^{-4}\) residual. This sensitivity is conditional on this copied discrete map and does not bound errors in evolved trajectories.

## Result 2 — the Hamiltonian maximum improves, but the residual relocates

| Resolution | Production cell-0 H residual | Analytic nominal cell-0 H residual | Analytic nominal cell-1 H residual | Production full-grid \(\max|H|\) | Analytic nominal full-grid \(\max|H|\) |
|---:|---:|---:|---:|---:|---:|
| \(N=160\) | \(+4.04156\times10^{-3}\) | \(+6.25313\times10^{-5}\) | \(+3.55752\times10^{-4}\) | \(4.04156\times10^{-3}\) | \(3.55752\times10^{-4}\) |
| \(N=240\) | \(+4.05101\times10^{-3}\) | \(+6.54767\times10^{-5}\) | \(+3.70172\times10^{-4}\) | \(4.05101\times10^{-3}\) | \(3.70172\times10^{-4}\) |
| \(N=320\) | \(+4.05429\times10^{-3}\) | \(+6.64722\times10^{-5}\) | \(+3.75276\times10^{-4}\) | \(4.05429\times10^{-3}\) | \(3.75276\times10^{-4}\) |

The full-grid Hamiltonian maximum drops by about \(10.8\)–\(11.4\times\), but the maximum moves from cell zero to cell one. Cell one is worse than its original cell-1 value in absolute magnitude and its sign reverses. Consequently, it would be misleading to report this as “constraints fixed”: the peak is reduced, yet local residual structure remains and must be explained.

The native-grid momentum witness remains at float64 roundoff (\(\max|M|\approx4.7\times10^{-18}\), \(2.0\times10^{-17}\), and \(7.6\times10^{-18}\)); the connection identity closes to zero in this metric-derived Lambda construction, and the determinant witness is at about \(2.22\times10^{-16}\). These are useful consistency checks, but the momentum field is expected to remain nearly unchanged in this particular copy because the experiment does not re-solve or evolve the extrinsic-curvature fields.

The fixed outer annulus \(5\le r\le10\) changes negligibly in the primary Ricci-gap maximum: approximately \(8.87\times10^{-8}\), \(1.89\times10^{-8}\), and \(6.19\times10^{-9}\) across the three resolutions. The Hamiltonian annulus witness is likewise unchanged within the reported digits.

## Result 3 — the separate Lambda derivative channel almost cancels the remaining gap

For the exact same analytic-half-cell candidate geometry, the offline chain-rule substitution yields:

| Resolution | Primary analytic nominal cell-0 gap | Chain-rule counterfactual cell-0 gap | Chain-rule counterfactual cell-1 gap |
|---:|---:|---:|---:|
| \(N=160\) | \(-1.88494\times10^{-4}\) | \(-1.492\times10^{-9}\) | \(+3.224\times10^{-10}\) |
| \(N=240\) | \(-1.93165\times10^{-4}\) | \(-6.089\times10^{-10}\) | \(+1.960\times10^{-10}\) |
| \(N=320\) | \(-1.94790\times10^{-4}\) | \(-3.311\times10^{-10}\) | \(+1.224\times10^{-10}\) |

This is a powerful diagnostic discriminator: on these initialized copied metrics, the Lambda derivative representation accounts for nearly all of the *remaining Ricci gap* after analytic half-cell injection. It does **not** authorize replacing the production \(D_4(\Lambda)\) term. Lambda is an independently evolved BSSN variable, and the chain-rule substitution may alter the off-constraint principal structure. Constraint witnesses in this report belong to the primary metric-derived Lambda candidate, not to an evolved run with that operator replacement.

## Result 4 — manufactured odd-Lambda stencil-symbol control

To determine whether the near-zero same-slice chain-rule gap justifies replacing the production D4(Lambda) operator, a separate manufactured test holds the initialized spatial metric fixed and injects smooth odd off-constraint modes delta-Lambda = epsilon sin(k r). It compares the pinned parity-odd fourth-order derivative with the known analytic derivative, and independently checks the fixed-metric vendor Ricci-response identity. It also verifies that the derivative-only chain-rule arm removes the D4(delta-Lambda) response while retaining the algebraic Lambda*D1(a) response.

The measured norm ratio matches the standard fourth-order modified-wavenumber symbol to floating-point precision; the full vendor gap-response identity closes to about 10^-18. As expected for a fourth-order centred derivative, the derivative accuracy depends strongly on the number of cells per mode:

| Resolution | dr | Wavelength 8: cells / relative derivative error | Wavelength 4: cells / relative derivative error | Wavelength 2: cells / relative derivative error |
|---:|---:|---:|---:|---:|
| N=160 | 0.5000 | 16 / 0.07783% | 8 / 1.1785% | 4 / 15.1174% |
| N=240 | 0.3333… | 24 / 0.01553% | 12 / 0.24249% | 6 / 3.51744% |
| N=320 | 0.2500 | 32 / 0.004932% | 16 / 0.07783% | 8 / 1.17848% |

The fixed physical wavelength-2 derivative error drops from 15.12% at N=160 to 3.52% at N=240 and 1.18% at N=320. At eight cells per wavelength the tested stencil's response is within about 1.18%; four cells per wavelength is plainly under-resolved in this control. Longer wavelength modes have smaller error, consistent with the stencil symbol.

**What this does establish:** the Lambda-response algebra is internally consistent, and D4(Lambda) has the expected, resolution-dependent response to a manufactured off-constraint mode. **What it does not establish:** the cause of the production Ricci residual, stability of the coupled evolved BSSN principal symbol, or that a chain-rule replacement is mathematically or physically admissible. The chain-rule arm removes the imposed off-constraint derivative response by construction; a near-zero residual on this frozen metric is not evidence that this is a valid evolution operator.




## Result 5 — frozen-stage RHS bridge closes in the Ricci-bearing L2 block

A first attempt required the lower-order L3 RHS to equal the Ricci derivative response; that test correctly failed because it attributed a Ricci/principal contribution to the wrong PIRK split. The test was corrected rather than relaxing its tolerance: the Ricci-bearing block is L2, and L3 remains a separately recorded lower-order channel.

On code head 3818206110a95c25276f73bb3fee62256bb2e105, the dedicated exact-head workflow [38797535339](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38797535339) passed source/vendor boundary checks, the four focused regression tests, the lower/nominal/upper analytic half-cell counterfactual, and the manufactured odd-Lambda mode sweep. The signed half-cell witnesses and the expanded manufactured-mode profiles are saved in [artifact 12001974762](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38797535339/artifacts/12001974762) and [artifact 12002229340](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38797535339/artifacts/12002229340), respectively.

For all 15 resolution/wavelength combinations (N=160, 240, 320; wavelengths 32, 16, 8, 4, 2), the measured change in the frozen primary Aa L2 RHS matches alpha times the independently predicted Ricci response to about \(2\times10^{-18}\) or better. After removing the imposed D4(delta-Lambda) derivative contribution, the L2 response matches the remaining algebraic Lambda*D1(a) term to about \(2\times10^{-18}\) or better. This closes the **frozen-stage algebraic accounting identity** and shows the manufactured off-constraint mode enters the Ricci-bearing Aa L2 block exactly where expected.

By contrast, the Aa L3 RHS perturbation is zero for the manufactured vacuum test and the Lambda L3/L2 deltas are at numerical roundoff in the test setup. The growing Aa L2 response as wavelength is shortened is consistent with the measured derivative-symbol error; for wavelength 2 and N=160 the maximum Aa L2 RHS delta is \(3.142\times10^{-4}\), with a 15.12% relative error in D4 compared with the analytic derivative.

This is a useful bridge between a spatial Ricci response and the stage RHS, but it remains a frozen-metric test. It does not test predictor/corrector coupling over time, the evolving Lambda equation's constraint propagation, stability of the coupled off-constraint principal symbol, or conservation over a full horizon. The chain-rule substitution still removes D4(delta-Lambda) response by construction and is not a production fix.

## Limits and next gate

The source/vendor pins pass, the observer leaves the initialized state unchanged, all profiles and endpoint sensitivity fields are finite, and the audit takes zero evolution steps. The analytic enclosure is conditional on the default uncharged Gaussian initializer with constant \(\phi,\Pi,K_t\), radiation enabled and \(K_t<0\). Charged-dust/Maxwell alternatives are not covered. The outer grid still uses the pinned discrete derivative/source reconstruction; this is a deliberately hybrid continuum-half-cell/discrete-outer-grid experiment, not a full continuum solution. Enclosure endpoint calculations use float64 rather than directed-rounding interval arithmetic.

The evidence now separates four matters: (1) the large origin-half-cell sampling contribution; (2) a persistent, sign-reversed Ricci/centre-Hamiltonian residual after that contribution is replaced; (3) a strong same-slice Lambda derivative-channel response; and (4) the expected, resolution-dependent D4 symbol error for manufactured off-constraint modes. The next valid step is independent numerical-method review of the signed residual relocation and chain-rule identity, followed by a manufactured test of the *coupled evolved* off-constraint principal symbol and constraint propagation—not just the now-closed frozen-stage L2 RHS identity. Then run high-resolution and full-horizon constraint/conservation diagnostics under an explicitly reviewed hypothesis. Do not combine the half-cell and Lambda counterfactuals into a production patch based only on the near-zero composed Ricci gap.

No production initializer, engine equation, vendor, dependency or campaign setting was modified. Scientific admission remains **NOT ADMITTED**.
