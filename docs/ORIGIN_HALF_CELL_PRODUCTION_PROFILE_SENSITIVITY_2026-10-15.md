# Origin half-cell estimator sensitivity on pinned production profiles — 15 October 2026

## Status and scope

**Research-only; scientific admission remains NOT ADMITTED.** This check reads an existing checksum-pinned, no-evolution initialized-state artifact. It performs no new engine initialization, no evolution, and no production, vendor, dependency, or campaign-setting edits.

The input is artifact `11993174037` from [workflow run 38776498208](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38776498208). The ZIP SHA-256 is `ff667892633d20c17fda62ceb9752de420f79276d0681ed9e8195dd1cb65f5c2`; `report.json` SHA-256 is `f2907932df5d76a0fb3377da7dcdb4e9bc11cf516b83b0769e3143178b4fe64f`; `profiles.npz` SHA-256 is `b9d2a849a3f32867dead913ca86eb61e8a6210b7db25629a9590315e67553587`.

The artifact reports source anchor `285a2ee0c1a3b97cddeb4519d8e749109b115aab` and vendor pin `d6052d605673ce9d82cdc99b0df79855fa2ea215`. All eleven recorded production-source blobs were independently compared with PR #258's declared base `main` SHA `6398620b85b7f5a027de0532ded80bdf04a9b53b`; every blob hash matches. The source profiles are therefore from the same pinned production implementations as the declared base, although the saved artifact retains its original source-anchor commit.

## Estimator comparison

For each saved final-iteration initializer source sample, two degree-three polynomials in \(r^2\) were fit with a scaled basis:

- **Fixed physical support:** every cell-centred source sample with \(r \le 2.5\).
- **Local support:** the first six cell-centred source samples, using support radius \(6\Delta r\).

Each fitted polynomial was analytically integrated from zero to \(r_0=\Delta r/2\). The signed quantity below is *fitted moment minus production point-sampled moment* (\(r_0\,s(r_0)\)). It is a method correction on the saved discrete profile, **not a true signed quadrature error**, because no exact continuous integral is available for these production samples.

| Resolution | \(\Delta r\) | Fixed-support sample count | Fixed-fit correction to point moment | Local-minus-fixed moment estimate | Local-minus-fixed \(B_0\) estimate |
|---:|---:|---:|---:|---:|---:|
| 160 | 0.5 | 5 | \(+1.30447861\times10^{-5}\) | \(-3.68044705\times10^{-8}\) | \(-1.47217882\times10^{-7}\) |
| 240 | 0.333333… | 8 | \(+3.85661575\times10^{-6}\) | \(+1.56010206\times10^{-8}\) | \(+9.36061239\times10^{-8}\) |
| 320 | 0.25 | 10 | \(+1.62614774\times10^{-6}\) | \(+8.80225266\times10^{-9}\) | \(+7.04180213\times10^{-8}\) |

Here \(B_0=I(r_0)/r_0\) is the first-cell initializer source average. The fixed-window fitted-moment correction decreases with observed order approximately three across both refinement pairs, while the associated \(B_0\) shift decreases at approximately second order. That is consistent with the leading \(O(\Delta r^3)\) point-sample half-cell moment discrepancy and the division by \(r_0\).

The fixed-vs-local estimator difference is about \(-0.282\%\), \(+0.405\%\), and \(+0.541\%\) of the fixed-window correction at \(N=160/240/320\). Its sign changes across the first refinement pair, so it should not be mistaken for a convergent signed error estimate. The four-fit degree/window-family spread (degree 3 at \(r\le2.5\), degree 3 at \(r\le3\), degree 4 at \(r\le3\), degree 5 at \(r\le3.5\)) is approximately \(0.334\%\), \(1.178\%\), and \(2.407\%\) of the fixed degree-three correction at those resolutions. This broader family remains more informative than comparing only two nearby degree-three windows.

## First-cell geometry sensitivity

The estimator difference propagates to the first-cell \(B_0\) as \(\delta I/r_0\). Applying the declared algebraic identities \(X^6=B\), \(a=X^{-4}\), and \(b=X^2\) to the two alternative \(B_0\) values gives an algebraic estimate of the resulting \(X/a/b\) sensitivity. This is **not** a self-consistent copied-initializer re-solve: it holds the saved production \(B_0\) fixed except for the half-cell moment delta and does not recompute the source/metric fixed point.

## Interpretation and next gate

The production samples do not show the large fixed-vs-local spread seen in the deliberately selected manufactured source. In the production profiles, the two degree-three estimators are close at these resolutions, while the broader degree/window family shows a larger but still modest spread. This is evidence against declaring the manufactured example's large persistent-intercept effect the measured production failure mode. It does not prove the fixed-window estimator is exact; low sample counts, profile sampling, and unmeasured continuous-source behaviour still matter.

The next decisive step is to establish an independent continuous/source oracle or a rigorous interpolation remainder bound for the actual initializer source, and to propagate candidate moment uncertainty through the full copied initializer at \(N=160/240/320\). Until then, use the moment differences as estimator corrections/sensitivity measurements, not a claimed error bound. No production correction is proposed.

## Reproduction

Download artifact `11993174037` from the linked workflow run, then run:

```bash
python tools/audit_origin_half_cell_fit_production_profiles.py \\
  --artifact path/to/current-main-centre-ricci-profiles.zip \\
  --output results/origin_half_cell_production_profile_sensitivity.json
```

The tool rejects any different ZIP, inner report/profile hashes, source/vendor pins, resolution set, evolution marker, or production-change marker. It uses NumPy only and records `scientific_admission=false`.
