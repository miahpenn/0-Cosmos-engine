# Centre Ricci investigation — evidence and next gates (15 October 2026)

## Current disposition

**Scientific admission: NOT ADMITTED.** No production or vendor equation has been changed. The evidence below establishes a discrepancy between Ricci representations and the behaviour of bounded offline counterfactuals; it does not identify an approved production correction.

PR #244 was merged at commit \`6a8f03d2fb7619b7e331f043091e6593ee8d5b59\`. GitHub reports \`merged_by=miahpenn\`, and the PR's formal review list is empty. The merge contains six diagnostic/test/documentation/workflow files; the exact-head pinned audit confirms no changes to production engine files, vendor code, dependencies, or physics. The repository's qualified-merge policy requires a trusted independent exact-head approval for non-documentation paths, so the observed merge path still needs reconciliation. Tracking issue: [#246](https://github.com/miahpenn/0-Cosmos-engine/issues/246). A merge and green CI are not scientific review.

Source anchor: \`285a2ee0c1a3b97cddeb4519d8e749109b115aab\`.
Vendor submodule: \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.

## 1. Initialized production-profile discrepancy

The source-pinned no-evolution observer compares vendor BSSN trace-free Ricci against an independent warped-product reconstruction from the proper-radial geometry. The initial centre discrepancy remains roughly constant under refinement, while the fixed outer annulus gap decreases.

| N | Baseline cell-0 gap | Cell-0 gap after offline chain-rule substitution | Max full-grid gap after substitution | Max gap on \(5\le r\le10\) |
|---:|---:|---:|---:|---:|
| 160 | \(+1.591813077\times10^{-3}\) | \(-2.449681397\times10^{-10}\) | \(2.449681397\times10^{-10}\) | \(7.761541355\times10^{-12}\) |
| 240 | \(+1.590029207\times10^{-3}\) | \(-5.056151277\times10^{-11}\) | \(5.056151277\times10^{-11}\) | \(1.629645868\times10^{-12}\) |
| 320 | \(+1.589416445\times10^{-3}\) | \(-1.625895577\times10^{-11}\) | \(1.625895577\times10^{-11}\) | \(5.214383883\times10^{-13}\) |

This is strong localization evidence for a discrete product/chain-rule compatibility mismatch on the initialized slice. The substitution is offline-only: it changes the \(D_4(\Lambda)\) operand in the diagnostic expression, not the production BSSN equation.

The refreshed exact-head run is [38690211345](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38690211345), artifact [11951678234](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38690211345/artifacts/11951678234), on code head \`607c07de47024601aad28704d937decd3a0204e9\). ZIP SHA-256: \`7dd17311fb83b6fac849e2161ec3e532954452237ee334d6982a909e9583e3a9\`. Internal hashes verified:
- \`report.json\`: \`dda18bb07e45f7c0121fe3819e95efec12a91816355771c638af211ed206642b\`
- \`profiles.npz\`: \`78de519b85050b82ba4ab8966ffba4e2f9debad78ccebef164e439aa50a0fb8d\`
- \`profiles.csv\`: \`a0208a3ca3d63906473ca8782d56810f9bfcce9413e98a84bda87aa9427826d0\`

The report confirms all 11 engine-source blob pins and the vendor pin, empty engine/vendor/dependency diff, finite profiles, unchanged initialized states, zero evolution steps, and no errors.

## 2. Origin half-cell hypothesis (secondary)

The pinned initializer's first accumulated radial integral entry uses \`0.5*dr*source[0]\` at \(r_0=\Delta r/2\). For a smooth even source \(s(r)=s_0+s_2r^2+s_4r^4+\cdots\), the exact moment terms are \(s_0r_0+s_2r_0^3/3+s_4r_0^5/5+\cdots\); point sampling overweights the \(r^2\) moment.

The frozen-source correction on a copied metric reduces the cell-0 gap magnitude by about 8.5×, but leaves a residual near \(-1.85\times10^{-4}\); cell 1 remains near \(+8.4\times10^{-5}\), and the fixed annulus barely changes. The combined frozen-source half-cell and offline chain-rule counterfactual reduces the first-cell gap to \(O(10^{-9}\!-\!10^{-10})\), with annulus gaps \(O(10^{-12}\!-\!10^{-13})\). This is **not** a self-consistent initializer re-solve, a convergence proof, or a production recommendation. The combined result is preserved in [run 38687830131 / artifact 11950013978](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38687830131).

## 3. Corrected manufactured Lambda-parity discriminator

**Correction:** the first manufactured-profile artifact preserved determinant one but used \(X=q^{-2/3}\), inconsistent with the production identities \(b=X^2\) and \(a=X^{-4}\). Those manufactured-profile values are superseded and must not be cited as evidence.

The corrected family uses
\[
q(r)=1+c_3r^2+c_5r^4,\qquad
a=q^{-4/3},\quad b=q^{2/3},\quad X=q^{1/3}.
\]
It satisfies \(ab^2=1\), \(b=X^2\), \(a=X^{-4}\), areal radius \(R=r\), and proper-radial factor \(f=1/q\). At N=160/240/320, identity residuals are at floating-point roundoff scale (typically \(\le 6\times10^{-16}\)); the maximum \(|R-r|\) is below \(3\times10^{-14}\).

The corrected artifact is part of [run 38690211345 / artifact 11951678234](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38690211345/artifacts/11951678234). It uses four fixed-metric connection arms: analytic continuum metric-derived \(\Lambda\), initializer-style discrete metric-derived \(\Lambda\), a smooth-odd off-constraint perturbation, and a separate even-parity negative control.

Representative N=320 first-two-cell signed gaps:

| Manufactured profile | Continuum arm | Discrete arm | Smooth-odd off-constraint arm | Even-parity negative control |
|---|---:|---:|---:|---:|
| q1 | \(+1.297\times10^{-13},+9.535\times10^{-14}\) | \(+1.078\times10^{-13},+7.368\times10^{-14}\) | \(+1.329\times10^{-5},+1.298\times10^{-5}\) | \(+6.222\times10^{-5},-8.889\times10^{-6}\) |
| q2 | \(-3.559\times10^{-14},-3.411\times10^{-14}\) | \(-2.669\times10^{-14},-2.475\times10^{-14}\) | \(-9.961\times10^{-6},-9.658\times10^{-6}\) | \(-4.667\times10^{-5},+6.667\times10^{-6}\) |

The direct warped-product Ricci is identical across arms because the metric is held fixed. Regular continuum/discrete arms agree with that independent geometric witness to roundoff; the smooth-odd off-constraint arm produces a much larger response while preserving the expected odd centre parity. The even control deliberately breaks regular parity.

**Limit of this artifact:** it is a parity/gap discriminator, not evolved stability evidence. The signed derivative-plus-algebraic-\(\Lambda a'\) response-closure test discussed on an earlier intermediate head is not part of the corrected exact-head artifact and is not claimed as verified for these corrected metric profiles. If that operator-budget regression is reintroduced, it must be rerun against this corrected metric family.

## 4. Provenance-matched evolved-state comparison

The older \(t=2.16\) sensitivity artifact [workflow 38289657512](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289657512) comes from engine source \`6d8b3f2790644591326df57c55908863f994cb0c\), which differs from the current pin in production files. It remains historical methodology context only; it is **not** the current-source evolved-state result.

The new bounded diagnostic [run 38690211348](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38690211348), [artifact 11951463175](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38690211348/artifacts/11951463175), uses the exact source/vendor pins above. ZIP SHA-256: \`016e3841218282a553631e5d39e3987b2f58f5051179bc02affb012ca65767d6\`. All embedded hashes match:
- \`report.json\`: \`a73c455f3d06482a4a59169c06ffe24da0197963322f5025f0ef8470618f33df\`
- \`evolved_state.npz\`: \`760320e4f18551198db0e9f807b9b16aab12ba484d8153da35e478d8e0364cb2\`
- \`profiles.npz\`: \`9b01f16d3d7cc199233045ed2e6e7235a6d4b206115c063c609f064ea001af43\`
- \`profiles.csv\`: \`1d030342b69e3d8f48e992b097cf78c883037a886470b0c0989bf268622cb813\`

Configuration: N=40, \(r_{\max}=80\), \(dt=0.06\), eight steps, final coordinate time \(t=0.48\). The current connection residual reaches \(\max|\Lambda_{\rm evolved}-\Lambda_{\rm metric}|=5.79987\times10^{-4}\) in the first two cells, and \(2.18072\times10^{-5}\) on the fixed annulus. The saved checkpoint and profile digests were verified, the observer left the evolved state unchanged, all profiles are finite, and no production/vendor/dependency paths changed.

| Signed TF-Ricci gap versus direct warped-product Ricci | Cell 0 | Cell 1 | Max on \(5\le r\le10\) |
|---|---:|---:|---:|
| Vendor BSSN using evolved \(\Lambda\) | \(+9.98026\times10^{-4}\) | \(-1.96837\times10^{-4}\) | \(4.46499\times10^{-5}\) |
| Same BSSN expression using metric-reconstructed \(\Lambda\) | \(+8.91279\times10^{-4}\) | \(-5.76272\times10^{-5}\) | \(8.46876\times10^{-6}\) |
| Offline chain-rule derivative for \(D_4(\Lambda)\) only, retaining evolved \(\Lambda\) algebraically | \(+4.84792\times10^{-5}\) | \(-7.89407\times10^{-6}\) | \(9.20251\times10^{-7}\) |

The derivative-only counterfactual reduces the centre discrepancy substantially but does **not** close it to roundoff on the evolved state. Reconstructing \(\Lambda\) in the full BSSN expression changes both derivative and algebraic terms, so it is a distinct counterfactual. This eight-step endpoint confirms that the mismatch persists away from the initial connection-constraint surface; it does not locate the first step at which the residual grows, diagnose a specific production RHS term, prove long-time stability, or authorize a production substitution.

## 5. Next research gates

1. Resolve the PR #244 merge/review-gate discrepancy in [issue #246](https://github.com/miahpenn/0-Cosmos-engine/issues/246); keep engineering/merge qualification separate from scientific review.
2. Obtain a distinct independent scientific review of the signed operator decomposition, parity test, and interpretation limits. No formal review of PR #244 is recorded.
3. Continue with a separate, source-pinned bounded diagnostic that records every intermediate step, connection residual growth, and signed \(A_a\) RHS components by stage; only then assess which term first drives the evolved centre discrepancy.

No production equation, vendor source, coefficient, gauge, boundary, threshold, dependency, or campaign setting is changed by the diagnostics. Scientific admission remains **NOT ADMITTED**.
