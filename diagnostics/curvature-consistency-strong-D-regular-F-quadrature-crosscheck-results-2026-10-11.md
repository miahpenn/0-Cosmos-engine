# Strong-D regular-F quadrature cross-check — results

Date: 2026-10-11  
Branch: `diag/independent-metric-ricci-audit`  
Workflow: [run 38100161657](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38100161657)  
Commit: `272b8a398b8b3decfe5b6b33f5c5a6297240ed51`  
Artifact: [strong-D-regular-F-integral-comparison](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38100161657/artifacts/11687920586)  
Artifact ZIP SHA-256: `c73eb1a137df1130ace7b3c73d6ae14ac9ff54df5ae7d9ee101fab441bd5167e`

## Independent formulation tested

This is a follow-up to the linear-(F(r^2)) interval-integral counterfactual. It adds a second shadow that:
- fits local cubic (F) interpolants in (r^2) over four adjacent samples;
- integrates (1+r^2F(r)) over each interval using five-point Gauss-Legendre quadrature;
- evaluates the origin interval using Gauss-5 applied to the first-four-sample cubic extrapolation.

For a cubic in (r^2), the integrand (1+r^2F(r)) is polynomial of degree at most eight in (r), so the five-point Gauss rule is exact for that interpolant in exact arithmetic. The Gauss origin integral agrees with the analytic regular-origin integral to floating-point roundoff (maximum difference in the recorded iteration histories approximately (3\times10^{-17})). This provides an independent implementation check on the integration operation; the interpolant itself remains a sampled-data approximation.

## Central Hamiltonian residual, raw and after the unchanged centre projection

| N | Production projected | Linear-(F) projected | Cubic-(F)/Gauss projected |
|---:|---:|---:|---:|
| 40 | (+6.89618\times10^{-4}) | (-7.24234\times10^{-5}) | (-7.64522\times10^{-5}) |
| 80 | (+1.00788\times10^{-3}) | (-5.00990\times10^{-6}) | (-5.00642\times10^{-6}) |
| 160 | (+1.07904\times10^{-3}) | (-3.05548\times10^{-7}) | (-2.99661\times10^{-7}) |
| 320 | (+1.09606\times10^{-3}) | (-1.88776\times10^{-8}) | (-1.84214\times10^{-8}) |

The independent cubic/Gauss route reproduces the residual convergence of the linear-(F) route at all four resolutions. At N=320, central (|H|) is (1.84\times10^{-8}), about (60{,}000\times) smaller than the production-initialized value (1.096\times10^{-3}). The projection delta at N=320 is only (-2.59\times10^{-10}), rather than (+3.52\times10^{-4}) in production mode or (+3.52\times10^{-4}) in the origin-only variant.

## Whole-grid residual

The cubic/Gauss mode's maximum absolute projected H over the whole grid also decreases with refinement:
- N=40: (7.65\times10^{-5}), at the centre;
- N=80: (5.01\times10^{-6}), at the centre;
- N=160: (3.00\times10^{-7}), at the centre;
- N=320: (1.84\times10^{-8}), at the centre.

This is a particularly useful distinction from the linear-(F) mode, where the whole-grid maximum at N=320 is larger ((7.65\times10^{-7})) away from the centre. The cubic/Gauss interpolant reduces both centre and off-centre initial residuals for this configuration.

At N=320 after centre projection in the cubic/Gauss mode:
- central (C_\Lambda=-3.35\times10^{-11});
- central (R_{vendor}-R_{metric}=+6.82\times10^{-9});
- central (H_{vendor}=-1.84\times10^{-8});
- global max (|H_{vendor}|=1.84\times10^{-8}).

All diagnostic algebraic gates passed. No production source, physics equation, projection, gauge or default was changed.

## What this supports

The stable result under two distinct ways of integrating the regular-source decomposition is strong evidence that the persistent strong-D initial central residual comes from the production radial integration convention: applying a trapezoid directly to (src=1+r^2F) is not compatible with the desired discrete constraint at the centre for this configuration. The origin-only alternative is insufficient; the interval integration of the (r^2F) part also matters.

This remains a diagnostic finding, not a production patch authorization. The two variants share the same source samples and regular-origin decomposition, so an evolution run of the best-behaved shadow and an independent (B(r))/constraint compatibility check are still required.

## Next test

Run the full unchanged strong-D production kernel from (a) the production initializer and (b) the cubic/Gauss shadow initializer at N=160, (r_{\max}=80), D=(10^{-4}), CFL (0.0075), through t=22.5. Record identical time checkpoints, full-grid and centre H, connection constraint, finite-state checks, proper time, cycle/handoff ledger, and term attribution. This will reveal whether the repaired *initial quadrature alone* removes the initial offset without creating a new evolution defect. Keep both trajectories diagnostic-only and do not modify production code.
