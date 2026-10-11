# Domain-size audit addendum: sub-cell ridge and local lapse check

**Date:** 2026-10-11  
**Parent findings:** [strong-D domain-size residual audit findings](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/independent-metric-ricci-audit/diagnostics/curvature_consistency/strong-D-domain-size-residual-audit-findings-2026-10-11.md)  
**Source artifact:** workflow [#38108719482](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38108719482), artifact `strong-D-domain-size-residual-audit`, ID `11692468908`, SHA-256 `98395c9934e950dda46a70e5a011bcc27778641b8ca02b670880bc367e4204b4`.

## Sub-cell outer-peak fit at t=12

To test whether the one-cell sampled shift is merely grid quantization, I fitted a quadratic to the three neighboring samples of (|H_{vendor}|) around the outer maximum. This is a local diagnostic interpolation, not a change to the evolution.

| Initializer | Domain | Grid maximum radius | Quadratic-fit radius | Fit peak amplitude |
|---|---:|---:|---:|---:|
| Production | (r_{max}=80) | 9.75 | 9.76551 | (2.829204\times10^{-4}) |
| Production | (r_{max}=160) | 9.25 | 9.34210 | (2.598720\times10^{-4}) |
| Regular-(F) shadow | (r_{max}=80) | 9.75 | 9.76573 | (2.828140\times10^{-4}) |
| Regular-(F) shadow | (r_{max}=160) | 9.25 | 9.34091 | (2.597733\times10^{-4}) |

The fitted displacement is about (-0.423) in radius for both initializer variants, rather than a full (-0.5) grid step. Therefore the shift is **not explained entirely by the grid maximum snapping to a neighboring cell**. It remains a modest, grid-scale shift; do not overinterpret the sub-cell interpolation as independent resolution convergence.

The fitted outer-peak amplitude is lower on the larger domain by about (8.15\%) for both initializer variants, essentially the same domain effect as the sampled-maximum comparison.

## Lapse at the ridge

The lapse values requested at the ridge are now extracted from the same full profiles. At each domain's fitted peak:

| Initializer | (alpha) at fitted peak, (r_{max}=80) | (alpha) at fitted peak, (r_{max}=160) | Difference |
|---|---:|---:|---:|
| Production | 0.101995 | 0.095557 | -0.006437 |
| Regular-(F) shadow | 0.102001 | 0.095528 | -0.006473 |

The larger-domain lapse at its fitted peak is about 6.3% lower. At the common fixed radius (r=9.5), (alpha_{160}-alpha_{80}=+0.004033) for both initializer variants (about +4.2% relative to the smaller-domain value). The sign differs because the fitted peaks occur at different radii and the lapse profile has a strong radial gradient here. These local values establish that lapse profiles differ in the ridge region, but they do **not** establish that the lapse difference causes the residual-amplitude difference.

The largest domain-induced lapse difference remains farther out at (r=39.75), about 0.11614. The first-cell lapse difference alone was not an adequate proxy for the local ridge lapse.

## Duplicate-run status and artifact integrity

At the latest check, workflow [#38114304847](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38114304847) was cancelled. Its uploaded artifact (ID `11693361272`, digest `sha256:03ee8437132cc9b5246f440d545900d26fe3eefa8e567c6f259be2c40495fc91`) is **partial**, containing only the production trajectory through (t=3.0) (13 samples, 800 accepted steps). It contains no completed production case and no shadow case; do not use it as the campaign result.

The intended retained run is [#38114367091](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38114367091). At the latest API check it still reported the numerical-evolution step in progress, with no artifact yet. Both runs use the same script and campaign settings; the cancelled run's partial output is not the authoritative record.

## Interpretation gate

The outer residual ridge is a shared, outward-moving feature in the four completed trajectories, and its late amplitude, fitted radius, and local lapse vary with domain size. The inner residual and global-maximum takeover timing remain initializer-dependent. The causal role of the lapse and the physical interpretation of the ridge remain unresolved. No physics patch is justified by these diagnostics.
