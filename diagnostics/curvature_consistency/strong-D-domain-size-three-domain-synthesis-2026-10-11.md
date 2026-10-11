# Strong-D Domain-Size Audit — Three-Domain Synthesis

**Date:** 2026-10-11  
**Status:** Complete; diagnostic-only; no physics patch justified  
**Branch:** `diag/independent-metric-ricci-audit`

## Purpose

Add the intermediate domain (r_{\max}=120) to the completed (r_{\max}=80) and (160) cases. The spacing was held fixed at \(\Delta r=0.5\), so this is a **domain-size / lapse-boundary sensitivity test**, not a spatial-resolution convergence study.

## Reproducibility and scope

All six trajectories completed to (t=12): production and regular-(F) cubic/Gauss shadow initializers at ((N,r_{\max})=(160,80),(240,120),(320,160)). Each used CFL (=0.0075), (D=10^{-4}), scalar amplitude (0.01), width (7), radiation enabled, and 49 profile samples (initial slice plus saved time samples through (t=12)). Each completed 3,200 accepted steps. Production physics, defaults, gauge, and projection were unchanged. The shadow case changes only the initial-(B) quadrature/initializer.

- Original two-domain workflow: [run #38108719482](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38108719482), artifact ID 11692468908, SHA-256 `98395c9934e950dda46a70e5a011bcc27778641b8ca02b670880bc367e4204b4`.
- Retained intermediate-domain workflow: [run #38114367091](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38114367091), head `8b3bddf08c5a88e5abd02d65b0f67966e147d156`, artifact ID 11693950041, SHA-256 `be2d0f1039ae1753876338ff3fa54eac567d8d8d4d70a5affb9146e2364d7155`.
- The separate duplicate run #38114304847 produced only a partial production trajectory to (t=3) and is excluded from all comparisons.

## Results

### 1. The intermediate domain confirms a smooth domain trend in the late outer residual feature

At (t=12), a quadratic fit through the three neighboring samples of (|H_{\mathrm{vendor}}|) around the maximum in (4\le r\le20) gives:

| (r_{\max}) | Initializer | Grid peak (r) | Fitted peak (r_*) | Fitted (|H_{\mathrm{vendor}}|_*) |
|---:|---|---:|---:|---:|
| 80 | Production | 9.75 | 9.76551 | (2.829204\times10^{-4}) |
| 120 | Production | 9.25 | 9.46998 | (2.668278\times10^{-4}) |
| 160 | Production | 9.25 | 9.34210 | (2.598720\times10^{-4}) |
| 80 | Regular-(F) shadow | 9.75 | 9.76573 | (2.828140\times10^{-4}) |
| 120 | Regular-(F) shadow | 9.25 | 9.46910 | (2.667258\times10^{-4}) |
| 160 | Regular-(F) shadow | 9.25 | 9.34091 | (2.597733\times10^{-4}) |

The fitted radius shifts inward monotonically as the domain grows: about (0.296) from 80 to 120, then (0.128) from 120 to 160; total shift about (0.423). The fitted peak amplitude drops about 5.69% from 80 to 120 and another 2.61% from 120 to 160 (about 8.15% overall). Production and shadow values nearly coincide for this late outer feature at each fixed domain.

The same fitted feature moves outward with time in all domains, but its track is domain-dependent. Production fitted radii at (t=(8,9,10,12)) are approximately:
- (r_{\max}=80): 6.257, 7.103, 8.012, 9.766
- (r_{\max}=120): 6.218, 7.009, 7.858, 9.470
- (r_{\max}=160): 6.200, 6.969, 7.791, 9.342

This supports describing a **shared outward-moving residual feature in the saved profiles**. It does not establish a physical front or prove that (D) itself propagates causally.

### 2. The global off-centre maximum switch is initializer-sensitive

Define the inner competitor as the maximum of (|H_{\mathrm{vendor}}|) over (1.25\le r<4), and the outer competitor over (4\le r\le20). The first *post-initial* sampled crossover to sustained outer dominance (after the inner sector dominates at intermediate earlier samples) is:

| Domain | Production | Regular-(F) shadow |
|---|---:|---:|
| (r_{\max}=80) | (t=8.75) | (t=6.25) |
| (r_{\max}=120) | (t=9.00) | (t=6.25) |
| (r_{\max}=160) | (t=9.00) | (t=6.25) |

These are sampled identity changes, not independently defined front-arrival times. The roughly 2.5–2.75 time-unit initializer dependence means this crossover should not be used as a physical propagation diagnostic.

### 3. The lapse is materially domain-dependent, including in the residual-ridge region

At (t=12), the logged first-cell lapse (cell centre (r=0.25)) is:

| (r_{\max}) | Production | Regular-(F) shadow |
|---:|---:|---:|
| 80 | 0.01655588 | 0.01655614 |
| 120 | 0.01693131 | 0.01693159 |
| 160 | 0.01727752 | 0.01727781 |

This is a smooth rise of about 4.36% from the smallest to the largest domain, with production and shadow nearly identical at each domain. Because the lapse is obtained from an elliptic solve over the domain, this logged first-cell value is not a domain-independent control.

At the fixed common radius (r=9.5), lapse values are (0.0954173, 0.0982247, 0.0994506) for (r_{\max}=80,120,160), respectively—a smooth overall increase of about 4.23%. By contrast, lapse evaluated **at the moving fitted residual peak** is about (0.101995, 0.097483, 0.095557), decreasing about 6.31% overall. These are not contradictory: the fitted peak moves inward across a steep radial lapse gradient. Comparisons at a fixed radius and at a moving peak answer different questions.

The maximum absolute lapse-profile difference over the shared radial grid at (t=12) is:
- (80\) versus (120): (0.07884), near (r=38.25) (shared-grid RMS difference (0.05062));
- (120\) versus (160): (0.03849), near (r=47.25);
- the direct (80\) versus (160) comparison previously measured (0.11614), near (r=39.75).

The 80-versus-120 maximum profile gap grows from about (3.90\times10^{-4}) at (t=0) to (7.88\times10^{-2}) at (t=12). The large difference is not confined to the outer boundary or ridge; it develops through the interior profile as the lapse collapses. The intermediate case lies between the endpoints rather than showing a new discontinuity. This does not by itself establish that the lapse causes the residual amplitude/position change, but it confirms that domain size changes the evolved gauge field enough to remain a serious confounder.

### 4. Initializer-dependent inner residual remains

At (t=12), the central-cell vendor Hamiltonian residual remains roughly (+1.092\times10^{-3}) in production for all three domains, while the regular-(F) shadow is about (-7.5\times10^{-6}). Domain enlargement does not remove this large initializer-associated split. Thus agreement of the late outer feature does **not** imply agreement of the full interior solution or constraint residual.

### 5. Run health and guardrails

All six trajectories reached (t=12) with no evolution failure; cycle-event and handoff counts remain zero. Across all saved samples, (a\gtrsim0.95225), (b\gtrsim0.9973917), and (X\gtrsim0.9353568), remaining finite and positive. Final accumulated (\tau) varies smoothly with domain (approximately 6.91852, 6.91016, 6.90678 for production at (r_{\max}=80,120,160)). This establishes run completion and basic geometric health, not constraint admission or a cycle/bounce/turnaround.

## Conclusion

The intermediate-domain pair confirms a continuous domain-size sensitivity rather than a new anomalous regime. The outer residual feature is robust to the initializer at fixed domain, but its amplitude, fitted radius, and timing of maximum-identity change depend on domain size. Meanwhile the elliptic lapse profile changes substantially throughout the shared interior, and the central/inner constraint residual remains initializer-sensitive.

**No physics patch is justified by this campaign. Do not label the residual feature a physical (D)-front, causal propagation, bounce, or turnaround.**

## Next diagnostic step

Do not launch another domain sweep yet. Use the saved (t=0\ldots12) slices to audit the **actual elliptic lapse equation and boundary contribution** at matched times across the three domains, separating:
1. change already present in the initial slice;
2. the subsequent time growth of the lapse-profile difference;
3. residual changes at fixed radii versus changes caused by the moving residual maximum.

Keep this post-processing/diagnostic-only. Do not alter the lapse equation, boundary condition, gauge, projection, or matter/geometry evolution until that audit identifies a specific defective layer.
