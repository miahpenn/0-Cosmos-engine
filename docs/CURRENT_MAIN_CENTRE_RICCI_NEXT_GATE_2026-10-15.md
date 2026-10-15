# Current-main centre Ricci discriminator — refreshed 15 October 2026

## Purpose

Refreshed after PR #237's signed momentum-tangent budget merged. Current-source evidence shows a geometry-dominated centre momentum tendency and a Ricci-dominated primary L2 Aa term. This observer compares the independent warped-product geometry formula with the vendor Ricci on the same initialized profiles; any discrepancy is localization evidence, not proof that either continuum equation is defective.

Diagnostic-only, no-evolution comparison of the current production initializer against an independent spherical warped-product Ricci reconstruction. This is a refreshed continuation of the archived #74/#78/#84/#140 evidence, not a repeat of their frozen source runs.

## Exact source

- Current `main` source anchor: `285a2ee0c1a3b97cddeb4519d8e749109b115aab`.
- Vendor gitlink: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- The tool pins the production kernel, initializer, BSSN adapter, CMC gauge, matter, stress, scalar, and fluid-source blobs. It refuses source/vendor/dependency drift and fails if it mutates any initialized state.
- Resolutions: (N=60,80,120), (r_{max}=80); no evolution step is taken.

## Measured terms

For the BSSN radial line element (d\ell=(\sqrt{a}/X)dr) and areal radius (R=r\sqrt{b}/X), the independent reference uses

[
R_{\ell}=\frac{R_r}{f},\qquad
R_{\ell\ell}=\frac{R_{rr}-R_r f_r/f}{f^2},\qquad f=\sqrt a/X,
]

[
\mathrm{Ric}^{\ell}{}_{\ell}=-2R_{\ell\ell}/R,\qquad
\mathrm{Ric}^{\theta}{}_{\theta}=\frac{1-R_{\ell}^2-RR_{\ell\ell}}{R^2}.
]

Its trace-free radial Ricci is compared with the pinned BSSN source (R^r{}_r-R/3), which uses the evolved conformal connection (Lambda). The report keeps the first two centre cells distinct from fixed coordinate-radius bands (2\le r<5) and (5\le r\le10). It also records the signed (\dot A_a) source budget: trace-free lapse Hessian, trace-free Ricci, vacuum L3, matter L3 increment, full initialized-slice RHS, and (\dot A_a/r^2), including an operand-scaled L2 reconstruction closure witness.

## Regression controls and limits

Tests retain the constant-curvature oracle, require a smooth manufactured regular profile to have bounded trace-free Ricci divided by (r^2), and include an intentionally parity-violating radius profile as a negative control.

The direct warped-product reference is independent of the evolved (Lambda) brackets, but shares the same pinned cell-centred radial differentiation operator and sampled metric arrays. This is an independent formula route, not an independent discretization or a continuum result. An initialized-profile comparison does not establish trajectory convergence, identify the responsible equation by itself, or authorize a production kernel patch.

No physical law, gauge, boundary, initializer, production/vendor source, campaign setting or scientific-admission rule is changed. Scientific admission remains **NOT ADMITTED**.
