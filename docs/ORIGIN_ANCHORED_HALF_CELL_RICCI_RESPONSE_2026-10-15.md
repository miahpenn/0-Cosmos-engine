# Origin-anchored half-cell Ricci response — 15 October 2026

## Status

**Research-only. Scientific admission remains NOT ADMITTED.** The diagnostic initialized and inspected source-pinned production profiles, then ran three copied source/metric recurrences per resolution. It took no evolution steps, modified no production/vendor/dependency files, and verified unchanged initialized-state hashes.

The exact-head audit completed successfully on `0a4a9f1839d853b874a333d36d2a5025c5bed16b` in [workflow run 38787838951](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38787838951), artifact [11997834796](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38787838951). Artifact ZIP SHA-256: `62d3e2102c8e759009ad1c8db8e42502e7666a2633c1199bdd5818535afbabe6`. `report.json` SHA-256: `cee8972c3fe3e180f0d60588b5e630766d43e51b13c8deda9edb79089db94ce9`; `profiles.npz` SHA-256: `38af797f7da93c434bcb47b3a6e0a3b680b428c22d80066de6303a8f9cbdb703`.

All eleven declared engine-source blob pins matched the current-main base `6398620b85b7f5a027de0532ded80bdf04a9b53b`; vendor pin `d6052d605673ce9d82cdc99b0df79855fa2ea215` matched and the engine/vendor/dependency diff was empty.

## Method

For a regular initializer source

[
s(r)=1-8pi r^2ho+2r^2K_rK_t+r^2K_t^2,
]

the exact regular-centre value is (s(0)=1). Three copied recurrences were compared:

1. **Fixed, unanchored:** degree-three fit in (z=(r/2.5)^2), using cell centres (rle2.5) inclusive; infer the intercept from samples.
2. **Fixed, origin anchored:** the same support and polynomial degree, but impose (c_0=1) and refit all higher coefficients under that constraint.
3. **Local, origin anchored:** degree-three fit with (c_0=1), using the first six cell-centre samples ((rle6Delta r)).

Each iteration recomputes the initializer source for its candidate (B), changes only the origin half-cell moment and preserves all later trapezoid intervals. Both the sixth iteration and converged copied fixed point are checked. The Ricci witness is the signed difference between the pinned vendor trace-free Ricci and the independent proper-radial warped-product reconstruction.

## Main results

The table compares the fixed, origin-anchored fit with the unconstrained fit at the converged copied fixed point.

| Resolution | Unanchored fitted intercept bias (c_0-1) | Moment correction, unanchored | Moment correction, anchored | Correction change | Cell-0 signed Ricci gap, unanchored | Cell-0 signed Ricci gap, anchored | Increase in gap magnitude |
|---:|---:|---:|---:|---:|---:|---:|---:|
| (N=160) | (-3.2325	imes10^{-8}) | (1.30447864	imes10^{-5}) | (1.30525733	imes10^{-5}) | (+0.0597%) | (-1.87566917	imes10^{-4}) | (-1.88629091	imes10^{-4}) | (+0.5663%) |
| (N=240) | (-1.06489	imes10^{-7}) | (3.85661583	imes10^{-6}) | (3.87414185	imes10^{-6}) | (+0.4544%) | (-1.85415121	imes10^{-4}) | (-1.93483464	imes10^{-4}) | (+4.3515%) |
| (N=320) | (-7.3523	imes10^{-8}) | (1.62614776	imes10^{-6}) | (1.63526630	imes10^{-6}) | (+0.5607%) | (-1.85080249	imes10^{-4}) | (-1.95030651	imes10^{-4}) | (+5.3763%) |

The signed second-cell gap also increases in positive magnitude under anchoring: (8.44117	imes10^{-5}	o8.45977	imes10^{-5}), (8.39051	imes10^{-5}	o8.53185	imes10^{-5}), and (8.38158	imes10^{-5}	o8.55590	imes10^{-5}) at (N=160/240/320).

All corrected recurrences converged in six updates at the reported precision, with a recomputed map defect of zero. The initialized state hash remained unchanged, and all checked copied metric profiles were finite. The fixed (5le rle10) annulus Ricci gap changed only at about (10^{-14}) or less in absolute value between the two fixed-window variants.

## Interpretation

**Anchoring is a valid regularity constraint, but it is not a Ricci-gap fix by itself.** It removes a small fitted-intercept bias, changes the estimated moment correction by only (0.06%)–(0.56%), and nevertheless increases the magnitude of the residual centre-Ricci gap by (0.57%)–(5.38%). The effect is strongly centre-localized; the fixed outer annulus is essentially unaffected.

This result does not invalidate the broader origin half-cell mechanism. Relative to the uncorrected production first-cell gap (about (+1.59	imes10^{-3})), the anchored copied recurrences still reduce the gap magnitude to roughly (1.9	imes10^{-4}), about an eightfold reduction. They do, however, show that imposing (s(0)=1) does not close the residual or improve it relative to the unconstrained fit.

Crucially, the constrained fit re-estimates higher coefficients as well as setting (c_0=1). Therefore the Ricci response in this experiment cannot be attributed to the constant coefficient alone. The exact continuous production-source half-cell integral is still unknown; fitted-minus-point moments remain corrections/sensitivity measures, not true quadrature errors.

## Follow-up: intercept-only and frozen-source controls — 15 October 2026

The next experiment ran on exact audit head `f8b957da02463ceac2f84f80ab2c4e83f4881bb4` in [workflow run 38788699294](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38788699294), artifact [11998384076](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38788699294). ZIP SHA-256: `32378a2ad273622a58870b195304306ea4734eff358c74675daab569b93fcf3e`; report SHA-256: `4e6dd42a9db967be1a22461bba95ec12258360b4f43b81ad7158158e5d60f410`. Status is `completed_diagnostic_only`, `errors=[]`, all eleven source blob pins match, vendor pin matches, production/vendor/dependency diff is empty, zero evolution steps, input state unchanged, all profiles finite, and `scientific_admission=false`. Seven manufactured-control tests passed on this head.

The new intercept-only arm changes only the fitted constant coefficient from its unanchored value to (c_0=1), while retaining the higher coefficients from the unanchored fit. It is compared with both the unconstrained and constrained least-squares fits. Two layers are recorded separately: (1) all three fits evaluated on exactly the same saved source samples, and (2) the fully iterated copied source/metric recurrences.

| Resolution | Unanchored fixed-point moment correction | Constrained-anchored correction | Intercept-only correction | Cell-0 gap magnitude increase: constrained | Cell-0 gap magnitude increase: intercept-only |
|---:|---:|---:|---:|---:|---:|
| (N=160) | (1.30447864	imes10^{-5}) | (1.30525733	imes10^{-5}) | (1.30528676	imes10^{-5}) | (0.5663%) | (0.5877%) |
| (N=240) | (3.85661583	imes10^{-6}) | (3.87414185	imes10^{-6}) | (3.87436401	imes10^{-6}) | (4.3515%) | (4.4067%) |
| (N=320) | (1.62614776	imes10^{-6}) | (1.63526630	imes10^{-6}) | (1.63533812	imes10^{-6}) | (5.3763%) | (5.4186%) |

The centre gap changes are near-identical for constrained and intercept-only anchoring: the constrained/intercept-only gap difference is below (0.06%) of the unanchored residual at each resolution. Matching-source comparisons likewise show the constrained moment minus the intercept-only moment is only (-2.94	imes10^{-10}), (-2.22	imes10^{-10}), and (-7.18	imes10^{-11}) at (N=160/240/320). These matched-source differences are much smaller than the change from imposing the constant term. This supports the conclusion that, in this fit family and source, the direct (c_0) replacement dominates the changed moment; it is not a universal claim about production physics.

### Frozen-source geometry/Ricci response

The audit also took the four moment corrections and applied each once to the **same unmodified production metric**, as (Delta B(r)=Delta I/r), with the source samples held identical between arms. This eliminates nonlinear source/metric feedback from the immediate comparison. The resulting cell-0 signed gaps are:

| Resolution | Free fixed fit | Constrained anchor | Intercept-only anchor | Local six-cell anchor |
|---:|---:|---:|---:|---:|
| (N=160) | (-1.87566879	imes10^{-4}) | (-1.88629090	imes10^{-4}) | (-1.88669235	imes10^{-4}) | (-1.88994682	imes10^{-4}) |
| (N=240) | (-1.85415084	imes10^{-4}) | (-1.93483465	imes10^{-4}) | (-1.93585739	imes10^{-4}) | (-1.93212900	imes10^{-4}) |
| (N=320) | (-1.85080224	imes10^{-4}) | (-1.95030651	imes10^{-4}) | (-1.95109018	imes10^{-4}) | (-1.94798344	imes10^{-4}) |

These frozen-source values match the fully iterated copied-recurrence cell-0 gaps to roughly (10^{-12}) at (N=160), (10^{-11}) at (N=240), and (10^{-12}) or better at (N=320). Thus source/metric feedback barely changes the Ricci witness for this specific counterfactual; the direction of the response is not caused by the self-consistent re-fit. The (5le rle10) annulus remains essentially unchanged (absolute changes in its maximum gap are about (10^{-14}) or smaller).

### Updated interpretation and remaining limitation

1. **Origin anchoring is not sufficient to repair the residual Ricci gap.** It enforces a mathematically known source value, but it changes the inferred half-cell correction and pushes the remaining central residual modestly farther from zero.
2. **The immediate effect is dominated by replacing (c_0), not by the constrained refit of higher coefficients.** The intercept-only and constrained-anchored Ricci results are nearly coincident.
3. **The response is strongly centre-localized.** The fixed physical annulus is practically unaffected.
4. **This still does not give the exact continuous production half-cell integral.** Fit-minus-point quantities remain estimator corrections, not true signed quadrature errors. Scalar/matter samples stay fixed while their sources are reevaluated at candidate metrics; no fields evolve.

The present result weakens the claim that a small fitted-intercept bias is the leading cause of the observed residual Ricci gap. The half-cell counterfactual still reduces the much larger baseline centre discrepancy by about eightfold, but the remaining sign-reversed (O(10^{-4})) residual requires a separate explanation—likely another source-fit coefficient, discrete derivative/geometry representation, or both. The next scientific gate is an independent continuous-source/remainder bound and a separate decomposition of the first-two-cell Ricci response by the even-series coefficients while keeping the Lambda derivative channel separate. Do not promote a quadrature rule from this evidence alone.

No production correction is proposed. Scientific admission remains **NOT ADMITTED**.
