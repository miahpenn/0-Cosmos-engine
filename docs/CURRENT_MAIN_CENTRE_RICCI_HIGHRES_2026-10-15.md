# High-resolution centre Ricci replay — 15 October 2026

## Question

Do the existing first-cell vendor-minus-warped-product trace-free Ricci differences at N=60/80/120 enter a stable refinement regime when the same production initializer is sampled more finely?

This is a bounded no-evolution diagnostic. It is not a production physics correction or a scientific-admission claim.

## Pins and configuration

- PR base / live main snapshot at branch creation: `3971d9a7e959e4128c030319ee7305d390df105b`.
- Production-source comparison anchor: `285a2ee0c1a3b97cddeb4519d8e749109b115aab`.
- Vendor submodule: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Resolutions: N=160, 240, 320; outer radius r_max=80.
- The audit validates the 11 expected engine blob hashes, vendor pin, and empty engine/vendor/dependency diff before it evaluates any profiles.
- The observer reports state digests and zero evolution steps. Its result bundle contains signed CSV/NPZ profiles, strict JSON, source-pin metadata, and checksums.

## Why this is distinct

The merged #218 artifact already records N=60/80/120 on the same declared engine-source anchor and vendor. PR #242 tests a synthetic c5 profile and an exactly predicted direct-metric stencil truncation term; it does not sample the initialized production geometry. This run extends only the initialized production-profile resolutions and writes to a fresh output path, preserving prior artifacts.

## Reported comparisons

For each N, inspect the signed vendor and direct warped-product TF-Ricci values for cells i=0 and i=1; their difference and difference divided by r²; the lapse-weighted Ricci contribution to the A_a RHS; and the maximum signed/absolute gap in the fixed coordinate annulus 5≤r≤10. Keep these center and annulus results separate and do not assume a convergence order ahead of the measurements.


## Primary discriminator: offline chain-rule closure for the BSSN connection `Lambda`

A newer independent offline decomposition of the original high-resolution artifact (team handoff: [Issue #5 comment 6185063664](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6185063664)) suggests the strongest current lead is the discrete derivative of the composite metric-to-`Lambda` relation, not the origin half-cell quadrature alone.

The pinned initializer sets
[
\Lambda = \frac{D_1(a)}{2a^2}-\frac{D_1(b)}{ab}+\frac{2}{r}\left(\frac1b-\frac1a\right).
]
The vendor BSSN Ricci calculation uses a separate (D_4(\Lambda)) term because (\Lambda) is an independently evolved connection variable. For a controlled **offline-only counterfactual**, this observer computes the ordinary chain-rule derivative of the displayed metric definition using the same sampled (a',b',a'',b'') stencil operands:
[
\Lambda'_{\rm chain}=\frac{a''}{2a^2}-\frac{(a')^2}{a^3}-\frac{b''}{ab}
+\frac{b'(a'b+ab')}{a^2b^2}
-\frac{2(1/b-1/a)}{r^2}
+\frac{2(-b'/b^2+a'/a^2)}{r}.
]
The sole Ricci-term substitution in this counterfactual is (D_4(\Lambda)\mapsto\Lambda'_{\rm chain}); the induced radial trace-free Ricci shift is (	frac23 X^2(\Lambda'_{\rm chain}-D_4(\Lambda))). The new report records both derivatives, their signed difference, the induced shift, and baseline-versus-counterfactual gaps for cells 0/1, the full grid and fixed (5\le r\le10) annulus, with all original source/vendor/state/checksum gates retained.

The refreshed observer reproduced the counterfactual on exact PR head `d99cf6109e3ba76ef8da940bbbef064a25a3a309`: dedicated no-evolution workflow [38687362269](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38687362269), artifact [11950443446](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38687362269/artifacts/11950443446). The ZIP SHA-256 is `0d6be6b33c0d07c7585c22ac273d3f06c0cb7552662f4097e96d8bd9bef5751c`; embedded report, NPZ and CSV checksums all verify. Source anchor `285a2ee0c1a3b97cddeb4519d8e749109b115aab` and vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215` were pinned; source/vendor/dependency diff is empty, the repository regression suite passed, all profiles are finite, initialized-state digests are unchanged, there are zero evolution steps, and the report has no errors.

| N | baseline cell-0 TF-Ricci gap | chain-rule cell-0 gap | maximum full-grid gap after substitution | fixed-annulus gap after substitution |
|---:|---:|---:|---:|---:|
| 160 | +1.591813077e-3 | -2.449681397e-10 | 2.449681397e-10 | 7.761541355e-12 |
| 240 | +1.590029207e-3 | -5.056151277e-11 | 5.056151277e-11 | 1.629645868e-12 |
| 320 | +1.589416445e-3 | -1.625895577e-11 | 1.625895577e-11 | 5.214383883e-13 |

The cell-1 residuals are (-2.1383\times10^{-10}, -4.7629\times10^{-11}, -1.5717\times10^{-11}), respectively. The run is diagnostic-only; a reduced metric-representation gap does not show that substituting the chain-rule derivative into the evolved production BSSN operator is safe.

**Critical limit:** this is strong evidence of a discrete product/chain-rule closure mismatch between two representations on the initialized slice. It does not authorize replacing the evolved BSSN connection derivative in production. The vendor explicitly relies on independent `Lambda` as part of the BSSN principal structure, and direct metric Ricci need not equal the evolved-connection expression away from the connection-constraint surface. No vendor or production change is part of this PR.

## Added frozen-source origin half-cell counterfactual

The pinned production initializer computes its first cumulative integral entry as `0.5*dr*source[0]` at the first cell centre (r_0=\Delta r/2). For a smooth even source
[
s(r)=s_0+s_2r^2+s_4r^4+O(r^6),
]
the regular half-cell integral is
[
I_0=s_0r_0+\frac{s_2r_0^3}{3}+\frac{s_4r_0^5}{5}+O(r_0^7),
]
but the point-sampled implementation yields (r_0s(r_0)=s_0r_0+s_2r_0^3+s_4r_0^5+O(r_0^7)). The quadratic moment is therefore overweighted by a factor of three. Because (X=B^{1/6}), the leading first-cell coefficient in (X) inherits the same factor-three tendency.

The expanded observer reconstructs the final source of the pinned six-iteration initializer and verifies the reconstructed radial integral (B) against the actual geometry value (X^6). It fits an even (s_0+s_2r^2+s_4r^4+s_6r^6) model over fixed physical radii r≤2.5 and r≤3.0, with the difference reported as a fit-width sensitivity check. It then corrects only the origin half-cell integral and propagates its (delta I/r) effect to a **copied metric**. All inter-centre trapezoids remain unchanged, the source is frozen rather than iterated again, and no production state or source file is edited. The report records baseline and counterfactual signed Ricci gaps for the first two cells and fixed annulus.

The origin-half-cell mechanism remains a secondary hypothesis. This is deliberately a local mechanism discriminator, not a self-consistent initializer candidate. Earlier origin-half-cell experiments improved initial constraints but did not restore evolved centre convergence; their source blobs differ from this current source anchor. A reduced frozen-source gap therefore supports the quadrature hypothesis but does not authorize a production patch, prove which curvature representation is wrong, or establish scientific admission.

## Combined frozen-source half-cell + Lambda chain-rule counterfactual

The same audit now composes the two candidate mechanisms on the copied origin-half-cell metric: it applies the even-source half-cell moment correction with the source frozen at the final iteration, then recomputes the offline derivative of the metric-to-Lambda relation using the same D1/D2 operands. This new arm does not iterate the initializer again, alter production data, or evolve the state.

The first combined report is [workflow run 38687830131](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38687830131), artifact [11950013978](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38687830131/artifacts/11950013978), artifact ZIP SHA-256 `9c13084da4b8aef3da40a37e2a6c91c5044c7a385d3031362bec5ad1612cb158`. Its internal report/NPZ/CSV hashes match the embedded manifest. All 11 engine-source blobs, vendor pin, state-digest checks and zero-step gates passed; the artifact reports no errors and no production/vendor/dependency changes.

| N | original cell-0 gap | half-cell-only cell-0 gap | Lambda-chain-only cell-0 gap | combined cell-0 gap | combined max gap on fixed 5≤r≤10 |
|---:|---:|---:|---:|---:|---:|
| 160 | +1.591813e-3 | −1.875669e-4 | −2.449713e-10 | −1.490857e-9 | 7.762128e-12 |
| 240 | +1.590029e-3 | −1.854151e-4 | −5.055663e-11 | −6.062155e-10 | 1.630068e-12 |
| 320 | +1.589416e-3 | −1.850802e-4 | −1.624953e-11 | −3.292417e-10 | 5.205571e-13 |

This combined counterfactual reduces the initialized first-two-cell gap from O(10⁻³) to O(10⁻⁹–10⁻¹⁰). The half-cell-only arm reduces that gap about 8.5× but leaves a residual; the chain-rule arm applied on its own nearly closes the original gap, and the combination leaves a comparably small residual on the changed metric. All combined full-grid profiles are finite, the initialized state digest remains unchanged, and the reported evolution-step count remains zero.

**Historical-data guardrail:** the local Library's prior [discrete-consistent initial Lambda/Ricci sensitivity run](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289657512) found that the stored initial Lambda matches the discrete metric-reconstructed Lambda to roundoff at t=0, while on saved t=2.16 states the independently evolved Lambda can depart from that identity and its frozen-RHS substitution can materially change the first-cell A_a RHS. The present chain-rule arm is therefore only an initialized-slice representation diagnostic. It does not authorize changing the evolved BSSN Lambda derivative or promoting the origin quadrature rule. The joint collapse supports a numerical compatibility hypothesis, not a unique cause, self-consistent initializer solution, evolved stability, or scientific admission.


## Manufactured Lambda-parity discriminator (new bounded test)

The observer now includes a manufactured, no-evolution operator-sensitivity study using two smooth-even metric profiles
\[
q(r)=1+c_3r^2+c_5r^4,\quad a=q^{-4/3},\quad b=q^{2/3},\quad X=q^{-2/3}.
\]
These preserve \(ab^2=1\), proper factor \(\sqrt a/X=1\), and areal radius \(R=rq\). For each fixed metric the audit evaluates four connection-field arms:

- **Continuum:** analytic metric-defined contracted connection sampled at cell centres.
- **Discrete:** initializer-style metric-derived connection formed from the pinned fourth-order derivative operators.
- **Off-constraint odd:** the discrete connection plus \(\epsilon r\exp[-(r/w)^2]\), which has the expected odd centre parity while keeping the metric fixed.
- **Parity-breaking negative control:** the discrete connection plus a constant even offset.

Every arm reports signed vendor/direct TF-Ricci differences in cells 0 and 1, their \(r^2\)-scaled values, the fixed \(5\le r\le10\) annulus maximum, and the gap after the same offline metric-chain-rule comparison. The report also exports the individual vendor radial and scalar Ricci components, plus signed annulus minima/maxima for the response budget.

For fixed metric samples, the pinned vendor expression contains both a derivative dependence through \(D_4(\Lambda)\) and an algebraic dependence through \(\Lambda a'\). The signed TF-Ricci response relative to the discrete-metric-\(\Lambda\) arm is therefore audited as
\[
\Delta(R^r{}_r-R/3)=\frac{2}{3}X^2D_4(\Delta\Lambda)+\frac{1}{2}X^2\frac{D_1(a)}{a}\Delta\Lambda.
\]
The regression checks that the observed gap change closes against the separate derivative and algebraic terms to roundoff-scaled tolerance. This prevents attributing the full off-constraint response to the chain-rule derivative alone. The direct metric Ricci must be identical across the four arms because the metric is held fixed; only the supplied connection field changes.

This distinguishes three effects without an evolved trajectory: representation sampling for a regular connection, a regular off-constraint connection perturbation, and a deliberate parity violation. It still cannot validate a production substitution; the operator’s off-constraint behaviour and later-state connection evolution remain review gates. The new signed response decomposition and closure checks are included in this diagnostic-only branch; exact-head workflows remain the gate for the amended head.

## Acceptance and limits

A useful result requires all source/vendor pins to match, finite profiles, unchanged state digests, zero evolution steps, and a reproducible checksum manifest. Passing those integrity gates means the diagnostic completed, not that either formula is correct. No production equation, vendor code, coefficient, gauge, boundary, threshold, dependency, or campaign setting is changed. Scientific admission remains **NOT ADMITTED**.


## Manufactured metric consistency correction — 15 October 2026 UTC

A review of the first manufactured-profile artifact found that its metric preserved determinant one but used a conformal-factor sample inconsistent with the production identities \(b=X^2\) and \(a=X^{-4}\). Its manufactured-profile results are superseded and are not used as scientific evidence. The corrected family keeps \(q(r)=1+c_3r^2+c_5r^4\) but sets \(a=q^{-4/3}\), \(b=q^{2/3}\), and \(X=q^{1/3}\), which satisfies \(ab^2=1\), \(b=X^2\), and \(a=X^{-4}\). The independently reconstructed geometry then has areal radius \(R=r\) and proper-radial factor \(f=1/q\). The revised audit exports the manufactured metric samples and checks these identities explicitly. All conclusions await exact-head validation of the corrected artifact; production and vendor code remain unchanged.


## Current-main freshness rerun — 15 October 2026 UTC

The diagnostic branch was refreshed onto live `main` at `eb09400e3bdc3fb3f3129f8d80883356c9f444f9`; the merge commit is `da30e32699e8eb29bbf403bcaee5441ca610779d`. On that refreshed commit, all three standard PR workflows passed: Lightweight validation [38689407069](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38689407069), strict JSON [38689407091](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38689407091), and retry/integration [38689406991](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38689406991).

The branch synchronization itself changed only `docs/STATUS.md` at its tip, so the path-filtered no-evolution workflow was not guaranteed to run for the refreshed commit. This documentation-only freshness note deliberately touches a watched path so that the pinned audit can be rerun against the refreshed branch head. It does not change the frozen production source anchor `285a2ee0c1a3b97cddeb4519d8e749109b115aab`, the vendor pin `d6052d605673ce9d82cdc99b0df79855fa2ea215`, or any production/vendor/dependency file. The refreshed-head specialized artifact must be inspected independently; prior artifacts remain tied to their recorded exact heads.


## Provenance-matched evolved-state gate — 15 October 2026 UTC

The historical t=2.16 Lambda-sensitivity artifact (workflow 38289657512) is retained as methodology context only: its recorded engine source commit is `6d8b3f2790644591326df57c55908863f994cb0c`, which differs from the current pinned source anchor `285a2ee0c1a3b97cddeb4519d8e749109b115aab` in production engine files. It is therefore not admissible as the provenance-matched follow-up state.

A separate source-pinned diagnostic workflow now creates and saves a short evolved state using the exact anchored production kernel and vendor. It advances only eight steps (`N=40`, `r_max=80`, `CFL=0.03)), then compares the BSSN Ricci with independently evolved `Lambda`, the same BSSN formula with `Lambda` reconstructed from the current metric, direct warped-product Ricci, and the offline derivative-only chain-rule counterfactual. The run preserves signed cell-0/cell-1 and fixed-annulus summaries, current connection residuals, a checkpoint, raw profiles, source pins, and checksums. It is diagnostic-only and does not change any production/vendor equation or launch a long campaign.

The corrected manufactured-profile arm also enforces all production metric identities (`b=X^2`, `a=X^-4`, `ab^2=1)); the earlier inconsistent manufactured metrics are explicitly superseded by the corrected artifact. The provenance-matched evolved-state interpretation remains pending its exact-head artifact and source/checksum validation.
