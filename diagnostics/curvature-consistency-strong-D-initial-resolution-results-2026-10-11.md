# Strong-D production initial residual versus spatial resolution

Date: 2026-10-11  
Branch: `diag/independent-metric-ricci-audit`  
Workflow: [run 38099657876](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38099657876)  
Commit: `6b2356e7d50bcb489a86e2dcff5d0682ecf392e4`  
Artifact: [strong-D-production-initial-resolution-trace](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38099657876/artifacts/11687390460)  
Artifact ZIP SHA-256: `9f30e79442237f1d3fbc7aa104af60ada14543bd6fd1555ef374c3c82230a361`

## Configuration

Production `build_initial_data` at fixed (r_{\max}=80), amplitude 0.01, width 7, D amplitude (10^{-4}), radiation on, with the existing centre regularity projection applied unchanged. Resolutions (N=40,80,160,320). Initial slice only; no evolution and no production changes.

## Main measurements

| N | Δr | Raw central (H_{vendor}) | Post-projection central (H_{vendor}) | Post-projection central (C_\Lambda) | Post-projection central curvature gap |
|---:|---:|---:|---:|---:|---:|
| 40 | 2.00 | (3.76922\times10^{-3}) | (6.89618\times10^{-4}) | (-1.82787\times10^{-4}) | (9.09616\times10^{-5}) |
| 80 | 1.00 | (3.98859\times10^{-3}) | (1.00788\times10^{-3}) | (-7.89889\times10^{-5}) | (8.78572\times10^{-5}) |
| 160 | 0.50 | (4.04156\times10^{-3}) | (1.07904\times10^{-3}) | (-3.90041\times10^{-5}) | (8.62473\times10^{-5}) |
| 320 | 0.25 | (4.05429\times10^{-3}) | (1.09606\times10^{-3}) | (-1.94646\times10^{-5}) | (8.59812\times10^{-5}) |

The post-projection central H trends toward a nonzero value near (1.1\times10^{-3}), rather than tending to zero over these resolutions. The central connection constraint decreases by approximately a factor of two per refinement, while the central H residual and central curvature gap level off. This makes the remaining H residual unlikely to be explained solely by the connection constraint or by a vanishing truncation error.

The projection change in H is dominated by its change in the vendor curvature. Its central delta trends from (-3.07960\times10^{-3}) at N=40 to (-2.95823\times10^{-3}) at N=320; other H terms barely change at the centre.

The metric-derived curvature counterfactual also has a nonzero post-projection centre H: (5.98656\times10^{-4}), (9.20027\times10^{-4}), (9.92790\times10^{-4}), (1.01008\times10^{-3}) at N=40,80,160,320. Thus replacing the vendor R by the independent metric-derived R does not remove the plateau.

## Diagnostic interpretation

- Strong evidence of a persistent central initial-data/constraint incompatibility in the production initializer, rather than a residual that vanishes simply by increasing N.
- The centre projection improves H substantially; it does not remove the offset.
- The metric-curvature gap after projection is smaller than the remaining H residual and converges to a small but nonzero value.
- The connection residual decreases with resolution, but H does not track it.
- This does not yet identify whether the mismatch comes from the first-half-cell radial quadrature, the extrapolated centre regularity coefficient, or another detail of the radial reconstruction. No correction has been applied.

## Next targeted experiment

Audit the first half-cell of the radial integral in `engine/v55_initial.py`. The current construction starts its cumulative integral with `0.5 * dr * src[0]`, even though the source has the regular origin form (src(r)=1+r^2F(r)), hence (src(0)=1). On a cell-centred grid, evaluating the whole inner half-cell using (src(r_0)) can misrepresent the (r^2F(r)) contribution by (O(\Delta r^3)) in the integral. Dividing by (r) to form B and then differentiating twice can turn that into a finite first-cell curvature bias. This is a concrete mathematical hypothesis, not yet an admitted cause.

Next: compare the untouched production initial builder with a diagnostic-only shadow reconstruction using a regular-origin expansion for (int_0^{r_0} src(r)dr), while holding all other source equations and settings fixed. Measure H before and after the same centre projection at all four resolutions. Do not patch production code until that controlled comparison discriminates the hypothesis.
