# Centre-Ricci origin even-series audit — 15 October 2026

**Disposition: research-only; scientific admission NOT ADMITTED.** This audit analyses saved profiles only. It does not change production equations, the initializer, vendor operators, dependencies, campaign settings, or any evolved state.

## Provenance and fail-closed input

The input is the existing source-pinned, no-evolution initialized-state diagnostic from [workflow run 38690211345](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38690211345), [artifact 11951678234](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38690211345/artifacts/11951678234).

- Artifact ZIP SHA-256: `7dd17311fb83b6fac849e2161ec3e532954452237ee334d6982a909e9583e3a9`.
- `report.json` SHA-256: `dda18bb07e45f7c0121fe3819e95efec12a91816355771c638af211ed206642b`.
- `profiles.npz` SHA-256: `78de519b85050b82ba4ab8966ffba4e2f9debad78ccebef164e439aa50a0fb8d`.
- Source anchor: `285a2ee0c1a3b97cddeb4519d8e749109b115aab`; vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Audit head: `607c07de47024601aad28704d937decd3a0204e9`. All eleven declared engine-source blob pins match; production/vendor/dependency diff is empty.
- Resolutions are (N=160,240,320), (r_{\max}=80), with zero evolution steps. The report records finite profiles, unchanged input state, and roundoff-scale geometric identities.

The new reader `tools/audit_centre_ricci_origin_series.py` refuses to run unless the report/profile hashes, source/vendor pins, eleven source blobs, resolutions, no-evolution marker, and no-source-change marker all match the recorded artifact.

## Method

For each initialized profile, set

\[
\delta(r)=\frac1{b(r)}-\frac1{a(r)}.
\]

Smooth even metric fields at a regular spherical centre imply \(\delta=O(r^2)\). Thus a local regular expansion has

\[
\frac{\delta(r)}{r^2}=c_2+c_4r^2+c_6r^4+\cdots.
\]

The diagnostic performs a least-squares fit to \(c_2+c_4r^2+c_6r^4\) on the **fixed physical window** \(0<r\le3\) for (N=160,240,320). It reports the signed residual at every fitted sample in two ways: including the first cell, and excluding cell zero then evaluating the fitted curve at the actual first-cell radius. The latter is an extrapolation discriminator, not assumed truth.

The comparison is against the existing `origin_cf_a/b/Lambda` fields from the source-pinned artifact. These are a *frozen-source half-cell-integral counterfactual*: the initializer source is held at the saved final iteration, the half-cell moment correction is propagated as \(\Delta I/r\), and the metric / Lambda / Ricci are recomputed on a copied geometry. It is not a self-consistent re-solve.

## Results: \(\delta/r^2\) even-series fit

| (N) | \(\Delta r\) | Production cell-0 \(\delta/r^2\) | Counterfactual cell-0 \(\delta/r^2\) | Production fit \(c_2\), including cell 0 | Counterfactual fit \(c_2\) |
|---:|---:|---:|---:|---:|---:|
| 160 | 0.5 | \(1.253599\times10^{-3}\) | \(4.187033\times10^{-4}\) | \(9.495792\times10^{-4}\) | \(4.463301\times10^{-4}\) |
| 240 | 0.3333… | \(1.255286\times10^{-3}\) | \(4.222438\times10^{-4}\) | \(7.872296\times10^{-4}\) | \(4.491359\times10^{-4}\) |
| 320 | 0.25 | \(1.255877\times10^{-3}\) | \(4.232819\times10^{-4}\) | \(7.011154\times10^{-4}\) | \(4.465756\times10^{-4}\) |

The fit spans the same physical radius on all resolutions, but has 6, 9, and 12 samples respectively because the grid is cell-centred. Fit residuals show a stark difference:

| (N) | Production max signed-residual magnitude | Production residual / profile peak | Counterfactual max residual | Counterfactual residual / profile peak | Max-residual reduction |
|---:|---:|---:|---:|---:|---:|
| 160 | \(3.213554\times10^{-4}\) | 25.63% | \(3.283499\times10^{-5}\) | 6.96% | 9.79× |
| 240 | \(4.734917\times10^{-4}\) | 37.72% | \(3.270865\times10^{-5}\) | 6.86% | 14.48× |
| 320 | \(5.571477\times10^{-4}\) | 44.36% | \(3.508373\times10^{-5}\) | 7.33% | 15.88× |

The inclusive-fit quadratic coefficient has a 30.57% relative span across resolutions for the production profiles, compared with 0.627% for the frozen-source counterfactual. The production residual becomes a larger fraction of the profile as resolution increases; the counterfactual residual remains near 7%.

## First-cell comparison using only cells 1 and outward

The second fit excludes cell zero and evaluates the remaining-cell fit at the true first-cell radius. This tests whether the first sample is compatible with the local continuation inferred from its neighbours.

| (N) | Production fitted continuation at (r_0) | Production signed first-cell residual | Counterfactual fitted continuation at (r_0) | Counterfactual signed first-cell residual | Residual-magnitude reduction |
|---:|---:|---:|---:|---:|---:|
| 160 | \(5.124079\times10^{-4}\) | \(+7.411908\times10^{-4}\) | \(4.805099\times10^{-4}\) | \(-6.180662\times10^{-5}\) | 11.99× |
| 240 | \(4.854646\times10^{-4}\) | \(+7.698212\times10^{-4}\) | \(4.650719\times10^{-4}\) | \(-4.282809\times10^{-5}\) | 17.97× |
| 320 | \(4.705205\times10^{-4}\) | \(+7.853563\times10^{-4}\) | \(4.556629\times10^{-4}\) | \(-3.238103\times10^{-5}\) | 24.25× |

The production first-cell value remains close to \(1.256\times10^{-3}\) while its radius shrinks from 0.25 to 0.125. Its signed disagreement with the neighbour-based continuation grows across these resolutions. The frozen-source correction brings that first-cell value much closer to the fitted continuation, by factors of roughly 12–24 in residual magnitude.

The odd-connection companion \(\Lambda/r\), which should also have an even continuation for an odd regular \(\Lambda\), shows the same direction: its first-cell value changes from about \(3.06\times10^{-3}\) in the production profiles to about \(1.51\text{–}1.52\times10^{-3}\) in the counterfactual. On the same (r\le3) fit, the maximum residual is lower by about 15–18× across these three resolutions. This is corroborating evidence, not an independent metric observable, because \(\Lambda\) is reconstructed from the same metric/stencil.

## Mechanistic consistency with the origin half-cell moment

The existing source report fits an even expansion of the initializer source and compares its analytic half-cell integral against the production point-sample approximation. The frozen-source correction \(\Delta I\) is:

| (N) | \(\Delta I\) |
|---:|---:|
| 160 | \(1.304479\times10^{-5}\) |
| 240 | \(3.856616\times10^{-6}\) |
| 320 | \(1.626148\times10^{-6}\) |

Those corrections decrease by factors 3.38 and 2.37 as \(\Delta r\) decreases by factors 1.5 and 1.333, respectively—consistent with an \(O(\Delta r^3)\) half-cell moment error from the \(r^2\) source term. That matches the expected leading scaling when replacing \(r_0s(r_0)\) with \(s_0r_0+s_2r_0^3/3+\cdots\). The two source-fit windows in the saved report change \(\Delta I\) by only \(3.7\times10^{-8}\), \(2.9\times10^{-8}\), and \(3.0\times10^{-8}\) at the three resolutions.

## Interpretation and next gate

**Result:** the multi-resolution fit, signed residuals, coefficient stability, and corroborating \(\Lambda/r\) profile all support the origin half-cell quadrature as a leading explanation for the first-cell regularity discrepancy in these initialized profiles. The frozen-source counterfactual improves even-series consistency on the same saved states.

**Not established:** this does not prove that the full initialized continuum solution is wrong, establish evolved trajectory convergence, or authorize promotion of the counterfactual to production. The counterfactual holds the initializer source fixed after correcting the half-cell moment; it does not re-solve the initializer source/metric fixed point.

**Next bounded action:** independently review the fit/term convention, then build a *separate diagnostic-only copied initializer* that iterates the corrected source and metric together at (N=160,240,320), checks all metric and matter invariants, and reruns the same Ricci/constraint discriminator. Only after that should the team consider an isolated production-initializer correction. No production files should be changed for that test; scientific admission stays **NOT ADMITTED**.

## Reproduction

```bash
python tools/audit_centre_ricci_origin_series.py \
  --profiles /path/to/profiles.npz \
  --report /path/to/report.json \
  --output results/centre_ricci_origin_series_audit.json
python -m pytest -q tests/test_centre_ricci_origin_series.py
```

The audit accepts only the exact SHA-256-pinned saved profile/report pair named above, refuses source/vendor pin drift, reports signed fit residuals, and records `scientific_admission=false`.


## Follow-up: self-consistent copied-initializer test — 15 October 2026

The next bounded test was run in the isolated diagnostic branch, using the same pinned initialization source and no-evolution setup at (N=160,240,320): [workflow run 38774311507](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38774311507), [artifact 11992384031](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38774311507/artifacts/11992384031).

Artifact ZIP SHA-256: `21a734250db76b2e7f8583661a53066702acc2674d33f64c388bbcbbd3d386d5`.
- `report.json`: `7b496123539b101ebefa96b151dd649e1be22c696037f22085dbb1e5faca6219`
- `profiles.npz`: `89fbd6513805988ab13eff5b22cea9b542e0d13ac78507f4696411d3de9dc84d`
- `profiles.csv`: `99f8d8c092b4e97c7466a16f3bd3158970486b78ea605695ccbf9e672fbe2574`

The run passed the complete repository regression suite, source-pin checks, the copied-initializer audit, and artifact upload. Its audit head is `6c0db3e82dbba707b9d2d7fb23c75fb125444048`; engine source anchor `285a2ee0c1a3b97cddeb4519d8e749109b115aab); vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215). All eleven engine-source blobs and the vendor pin match; production/vendor/dependency diff is empty. All three runs have zero evolution steps and confirm the input state hashes are unchanged.

### Self-consistent source–metric loop

The diagnostic first recreated the ordinary six-pass point-sampled source recurrence. Its (B) profile matches the initialized (X^6) to at most (6.7\times10^{-16}). It then replaced *only* the first half-cell integral with the analytically integrated even-source fit, recomputing the source from each candidate metric on every pass. The fixed point was reached within six passes; recomputed map defects are zero at reported floating-point precision. The metric identities (ab^2=1), (b=X^2), (a=X^{-4}), (X^6=B), and (R=r) close at roundoff, with the largest reported radius identity residual (1.43\times10^{-14}). Scalar/matter fields are held fixed at their initialized samples; their source contributions are reevaluated from the candidate metric. No field is evolved.

| Resolution | Production \(\delta/r^2\), cell 0 | Corrected self-consistent \(\delta/r^2\), cell 0 | Production signed centre TF-Ricci gap | Corrected signed centre gap | Corrected / baseline gap magnitude |
|---:|---:|---:|---:|---:|---:|
| 160 | \(+1.253599\times10^{-3}\) | \(+4.187033\times10^{-4}\) | \(+1.591813\times10^{-3}\) | \(-1.875669\times10^{-4}\) | 11.78% |
| 240 | \(+1.255286\times10^{-3}\) | \(+4.222438\times10^{-4}\) | \(+1.590029\times10^{-3}\) | \(-1.854151\times10^{-4}\) | 11.66% |
| 320 | \(+1.255877\times10^{-3}\) | \(+4.232819\times10^{-4}\) | \(+1.589416\times10^{-3}\) | \(-1.850802\times10^{-4}\) | 11.64% |

Thus the self-consistent half-cell correction changes the first-cell \(\delta/r^2\) by about a factor of three and reduces the *magnitude* of the first-cell TF-Ricci discrepancy by factors of approximately 8.49, 8.58, and 8.59. The residual changes sign and remains O(\(10^{-4}\)); the half-cell term alone does not close the Ricci gap. In the fixed \(5\le r\le10\) annulus, the maximum gap changes only slightly: (8.8675\times10^{-8}\to8.8603\times10^{-8}) at N=160, (1.88824\times10^{-8}\to1.88780\times10^{-8}) at N=240, and (6.19092\times10^{-9}\to6.19030\times10^{-9}) at N=320.

The corrected fixed-point half-cell integral differs from the frozen-source correction by only (2.91\times10^{-13}), (8.24\times10^{-14}), and (2.24\times10^{-14}), respectively. On these profiles, recomputing source feedback therefore barely shifts the integral relative to the earlier frozen-source counterfactual.

### Combined derivative-representation check on the corrected metric

The previously separated Lambda chain-rule counterfactual was then recomputed on the **self-consistent corrected metric**, not inferred from the frozen-source copy. The first-cell combined vendor-minus-direct TF-Ricci gaps are (-1.491\times10^{-9}), (-6.062\times10^{-10}), and (-3.293\times10^{-10}) for N=160/240/320; the second-cell residuals are (+3.221\times10^{-10}), (+1.949\times10^{-10}), and (+1.216\times10^{-10}). This is a strong algebraic localization result: the origin half-cell counterfactual removes much of the central discrepancy, and the derivative-only Lambda representation shift nearly cancels the remaining central gap on the same copied geometry.

It remains a **two-counterfactual decomposition**, not validation of a production equation change. The Lambda chain-rule derivative is not the production BSSN evolution operator; the tiny combined residual cannot establish that it should replace the evolved-connection representation, nor does it prove evolved-trajectory convergence.

### Important source-fit sensitivity qualification

A separate sensitivity scan on the checksum-pinned corrected-source profile compared fits of degree 1, 2, and 3 in (r^2) over several fixed windows. For the chosen degree-3 even series ((s=s_0+s_2r^2+s_4r^4+s_6r^6)), the half-cell correction varies by at most about 0.28%, 0.75%, and 1.86% across admissible windows ending at (r=3.0) for N=160/240/320 (N=160 has too few samples for a degree-3 fit at (r\\le1.5)). If the fitting window is widened to (r\le3.5), the difference from the (r\le2.5) reference rises to about 1.31%, 5.56%, and 7.58%.

Lower-order fits over the wider windows are much less stable: a degree-2 source fit ending at (r=2.5) shifts the correction by (-5.35\%), (-23.22\%), and (-30.34\%) from the degree-3 reference at N=160/240/320; at (r\le3.0) those differences are (-17.82\%), (-47.45\%), and (-89.59\%). The source-fit degree and support are therefore substantive modeling assumptions, not harmless numerical details.

**Sensitivity scan verification:** The degree/window scan is now included in the regenerated source-pinned report from [workflow run 38775067240](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38775067240), [artifact 11992219877](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38775067240/artifacts/11992219877). Artifact ZIP SHA-256: `690f47afac9ce67ec59909a28fea4863e8947801c7d9c41d816cc0c3183fc600`; report JSON SHA-256: `ba06e662eb15104089bbc82b4f5011e5d83a88a2e8e3e3c648d3cdc09e6e553c`. The exact diagnostic head `1e37fbc762de42579b50129d75fcd784551c6e07` passed compilation, source/vendor/dependency pin validation, the full repository suite (**673 passed**), and the no-evolution profile audit. The new report has `errors=[]`, `scientific_admission=false`, zero evolution steps, and confirms the input initialized state was unchanged.




### Higher-order fits and leave-one-cell-out check

The next scan expands the family to degrees 1–5 in \(r^2\) and reports leave-one-cell-out (LOO) prediction errors wherever the remaining points determine the fit. Exact/saturated polynomials (no residual degrees of freedom) are labelled and are not used in the overdetermined spread bound.

Run [38775616302](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38775616302), [artifact 11993781423](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38775616302/artifacts/11993781423), regenerated the report on diagnostic head `548003538a0946be2a232d3de8ca7fdcef83cde9`. Artifact ZIP SHA-256: `400b3c17fed510f31e6c15144aa9b7d7b043004fec592e813238ef4a6e04b7c9`; report JSON SHA-256: `8e79fe9c6a2c2a8429face4cf0380f488cb1bbb16ac0136048b9f6319d5afda3`. The exact-head full test suite passed (**673 passed**), source/vendor/dependency pin checks passed, and the audit completed diagnostic-only with `errors=[]`, zero evolution steps, and unchanged input state.

Across **overdetermined** degree-3/4/5 fits and windows ending at \(r\le3.0\), the total span of the fitted half-cell correction, relative to the degree-3 \(r\le2.5\) reference, is:

| Resolution | Correction span across eligible degree-3/4/5 fits, windows \(r\le3\) | Reference correction |
|---:|---:|---:|
| N=160 | 0.331% | \(1.304479\times10^{-5}\) |
| N=240 | 1.181% | \(3.856616\times10^{-6}\) |
| N=320 | 2.411% | \(1.626148\times10^{-6}\) |

The degree-4/5 estimates on eligible windows are about +0.05% from the reference at N=160, about +0.42–0.44% at N=240, and about +0.50–0.55% at N=320. The larger total spread at N=320 comes from the degree-3 \(r\le3\) fit itself (−1.864%), while the degree-4/5 estimates stay grouped near +0.5%. This is a useful truncation diagnostic; it is not proof that degree 4 or 5 is the correct physical prescription.

LOO source prediction adds a separate check. For the degree-3 \(r\le2.5\) reference fit, the maximum absolute LOO source residual is \(1.775\times10^{-6}\), \(1.642\times10^{-6}\), and \(8.110\times10^{-7}\) at N=160/240/320. For the degree-4 \(r\le3\) fits it is \(3.709\times10^{-7}\), \(1.784\times10^{-7}\), and \(1.203\times10^{-7}\); for degree 5 at \(r\le3\), no overdetermined LOO estimate exists at N=160, while the maxima are \(7.88\times10^{-9}\) and \(4.43\times10^{-9}\) for N=240/320. These errors favour the higher-order local source representation on the sampled profiles, but small in-sample/LOO residuals alone do not establish a continuum error bound.

**Updated numerical conclusion:** degree-1/2 fits are not defensible candidates for this correction on the tested windows; degree-3 to degree-5 fits are much more consistent, with a measured spread around 0.3–2.4% over the admitted window family. That bounds, but does not eliminate, source-fit model uncertainty. The half-cell quadrature hypothesis remains strongly supported; the exact even-series truncation and fit support still need independent numerical-method review before any production change.

**Next gate:** independently justify the even-source truncation/window against the saved source’s signed residuals and smoothness, then bound the half-cell moment over defensible fit families before considering any production patch. No production/vendor/dependency file has changed; scientific admission remains **NOT ADMITTED**.
