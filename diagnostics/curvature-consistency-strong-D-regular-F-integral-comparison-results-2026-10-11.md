# Strong-D regular-origin radial quadrature comparison — results

Date: 2026-10-11  
Branch: `diag/independent-metric-ricci-audit`  
Workflow: [run 38100069694](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38100069694)  
Commit: `42d9a83e1c910a25e58767338e6d23084441d806`  
Artifact: [strong-D-regular-F-integral-comparison](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38100069694/artifacts/11687406203)  
Artifact ZIP SHA-256: `361f027b0a4c0e55dc7e7bb6235c6e917e489e86f78700dd72ae3ee6d7b81217`

## Design

At (N=40,80,160,320), fixed (r_{\max}=80), (D=10^{-4}), amplitude (0.01), width 7, radiation on, compared:

1. **Production:** untouched `build_initial_data`, cumulative trapezoid of (src(r)).
2. **Origin-only:** analytic integral on ([0,r_0]), trapezoids of (src) over all subsequent centre-to-centre intervals.
3. **Regular-(F):** use (src=1+r^2F(r)). Compute the origin integral using the same regular-origin polynomial extrapolation; on every interval ([r_i,r_{i+1}]), interpolate (F) linearly in (r^2) and analytically integrate (1+r^2F(r)). Apply the *same existing centre projection* to each raw geometry.

Only diagnostic shadow quadrature changes in modes 2/3. The source equations, six B iterations, matter initialization, gauge, projection and production defaults remain unchanged. All term-reconstruction and curvature-replacement diagnostic gates passed.

## Central vendor Hamiltonian residual, before and after existing projection

| N | Production raw | Production projected | Origin-only projected | Regular-(F) raw | Regular-(F) projected |
|---:|---:|---:|---:|---:|---:|
| 40 | (+3.76922\times10^{-3}) | (+6.89618\times10^{-4}) | (+4.21456\times10^{-5}) | (-3.85543\times10^{-5}) | (-7.24234\times10^{-5}) |
| 80 | (+3.98859\times10^{-3}) | (+1.00788\times10^{-3}) | (+3.37230\times10^{-4}) | (-3.93668\times10^{-6}) | (-5.00990\times10^{-6}) |
| 160 | (+4.04156\times10^{-3}) | (+1.07904\times10^{-3}) | (+4.02471\times10^{-4}) | (-2.77626\times10^{-7}) | (-3.05548\times10^{-7}) |
| 320 | (+4.05429\times10^{-3}) | (+1.09606\times10^{-3}) | (+4.18006\times10^{-4}) | (-1.78745\times10^{-8}) | (-1.88776\times10^{-8}) |

At (N=320), the regular-(F) counterfactual reduces (|H_{center}|) from (1.0961\times10^{-3}) (production) to (1.8878\times10^{-8}) after the same projection: roughly five orders of magnitude smaller. Projection changes the regular-(F) central H by only (-1.00\times10^{-9}), unlike the approximately (+3.52\times10^{-4}) shift in the origin-only mode at N=320.

The origin-only correction was not sufficient: it left the post-projection central residual near (4.18\times10^{-4}). The important difference between origin-only and regular-(F) is that regular-(F) also changes each centre-to-centre integration interval from the trapezoid of (src) to the analytic integral of the regular decomposition (1+r^2F(r)). This identifies the interval quadrature—not only the first partial cell—as the key numerical compatibility lead.

## Refinement behaviour in regular-(F) mode

Central post-projection (|H|):
- N=40: (7.24\times10^{-5})
- N=80: (5.01\times10^{-6})
- N=160: (3.06\times10^{-7})
- N=320: (1.89\times10^{-8})

The successive decrease factors are approximately 14.5, 16.4, and 16.2. This is consistent with near-fourth-order convergence of the *central residual* over the finer resolutions, though these four points alone are not a formal convergence proof.

The maximum absolute projected H over the whole grid in regular-(F) mode is (8.05\times10^{-5}), (1.47\times10^{-5}), (3.19\times10^{-6}), (7.65\times10^{-7}) for N=40,80,160,320. It decreases more slowly, with factors around 5.6, 4.6, and 4.2, consistent with a lower-order off-centre residual remaining. At N=320 its maximum is in cell 20 at (r=5.125), not at the centre.

At N=320 after projection in regular-(F) mode:
- central (C_\Lambda=-3.82\times10^{-11});
- central (R_{vendor}-R_{metric}=+6.83\times10^{-9});
- central (H_{vendor}=-1.89\times10^{-8});
- global max (|H_{vendor}|=7.65\times10^{-7}).

## Interpretation

This is strong diagnostic evidence that the production cumulative trapezoid applied directly to (src(r)) is incompatible with the regular-origin structure of the source at the origin and in subsequent cell intervals. Respecting (src=1+r^2F) within the radial integration produces a centre residual that converges toward zero at a much higher rate and leaves the existing centre projection almost non-invasive at high resolution.

It does not yet authorize changing the production initializer: regular-(F) assumes (F) is linear in (r^2) on each interval, and needs an independent quadrature/reference check. No production code or physical model has been changed.

## Next validation

Compare regular-(F) interval integrals against an independent higher-order quadrature applied to the same already-computed (src) samples (or a separately sampled smooth-source reconstruction), retaining the same initial state and source definitions. Check that both quadratures lead to consistent (B(r)), central H convergence, and projection deltas. Only after the quadrature comparison agrees should a separate production patch be considered.
