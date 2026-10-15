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

## Next discriminator

Add an explicitly labelled **intercept-only replacement** control: for the same source samples, first perform the unconstrained fit, then replace only its fitted constant with (1) while leaving the fitted higher coefficients unchanged. Compare it against (a) the unconstrained fit and (b) the constrained least-squares anchored fit. Report both same-source moment differences and the fully iterated copied-initializer signed Ricci gaps. This separates the direct constant-term effect from the higher-coefficient refit and the nonlinear source/metric feedback.

No production correction is proposed. Scientific admission remains **NOT ADMITTED**.
