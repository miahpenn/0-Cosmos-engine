# Analytic half-cell injection into the self-consistent copied initializer — 15 October 2026

**Disposition: research-only; scientific admission NOT ADMITTED.** This experiment tests whether replacing the polynomial-fitted origin half-cell moment with the independent continuous-source integral changes the initialized centre-Ricci mismatch, while reporting the signed constraints and a separate Lambda derivative counterfactual.

## Exact-head validation and provenance

The latest dedicated run, [GitHub Actions 38795244429](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38795244429), succeeded on PR #265 head \`0f18fabe72c92dcc3b0e9076ae86eb449641bad8\). Both focused tests passed, all declared source blobs and the vendor pin matched, the recurrence report was finite, the input state hash was unchanged, and zero evolution steps were taken. [Artifact 12001317982](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38795244429/artifacts/12001317982) includes the report, signed witnesses, full profiles, and checksum manifest.

- Engine source anchor: \`285a2ee0c1a3b97cddeb4519d8e749109b115aab\`.
- Vendor commit: \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.
- Audit checkout head: \`596568e3d38d7ba10f12d1850eaa180a5a01a16d\`.
- Artifact ZIP SHA-256: \`83021bd8e06ca40b8a95c329e7779cd354b1abdf5dfef2521a00edfca973ec51\`.
- Report JSON SHA-256: \`a41dbfa04293db3eaa99d85f2bbb377ff20f3253f7e18e9d0e82ed5f59066ccf\`.
- Signed witness CSV SHA-256: \`568e40562454b8c5342a100295c3a5e120846cb794cb326fa791cd14df297285\`.
- Full profile NPZ SHA-256: \`3db720972c5b480b24b17a9a60d571312425f27b5b9270bc86673da8afd78477\`.

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

## Limits and next gate

The source/vendor pins pass, the observer leaves the initialized state unchanged, all profiles and endpoint sensitivity fields are finite, and the audit takes zero evolution steps. The analytic enclosure is conditional on the default uncharged Gaussian initializer with constant \(\phi,\Pi,K_t\), radiation enabled and \(K_t<0\). Charged-dust/Maxwell alternatives are not covered. The outer grid still uses the pinned discrete derivative/source reconstruction; this is a deliberately hybrid continuum-half-cell/discrete-outer-grid experiment, not a full continuum solution. Enclosure endpoint calculations use float64 rather than directed-rounding interval arithmetic.

The evidence now separates three effects more sharply: (1) the large origin-half-cell sampling contribution, (2) a persistent, sign-reversed Ricci/centre-Hamiltonian residual after that contribution is replaced, and (3) a strong same-slice Lambda derivative-channel response. The next valid step is independent mathematical/numerical review of the chain-rule identity and the signed residual relocations, plus a manufactured off-constraint stability/principal-symbol control before any operator modification is contemplated. Then re-run high-resolution and full-horizon constraint/conservation diagnostics on an explicitly declared hypothesis. Do not combine the two counterfactuals into a production patch based only on the near-zero composed Ricci gap.

No production initializer, engine equation, vendor, dependency or campaign setting was modified. Scientific admission remains **NOT ADMITTED**.
