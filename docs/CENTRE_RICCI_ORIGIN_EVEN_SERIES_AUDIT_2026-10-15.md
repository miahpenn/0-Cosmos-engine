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
