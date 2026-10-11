# Strong-D origin-integral counterfactual — results

Date: 2026-10-11  
Branch: `diag/independent-metric-ricci-audit`  
Workflow: [run 38099902239](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38099902239)  
Commit: `cc8cdde781d5ec74df8f0d4efb9df447fdc4410e`  
Artifact: [strong-D-origin-integral-counterfactual](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38099902239/artifacts/11686344335)  
Artifact ZIP SHA-256: `481faf837904b86c23249c405f42ead13f7d4385e323e45dd31f4a39ca7f2b92`

## Design

At (N=40,80,160,320), fixed (r_{\max}=80), amplitude 0.01, width 7, D amplitude (10^{-4}), radiation on, compare:
- untouched `engine.v55_initial.build_initial_data`;
- diagnostic shadow initializer differing only in the first integral from the origin to the first cell centre (r_0=\Delta r/2).

The production first half-cell contribution is (r_0\,src(r_0)=0.5\Delta r\,src(r_0)). The shadow uses the regular-origin form (src(r)=1+r^2F(r)), extrapolates (F) as a polynomial in (r^2) from the first four cell centres, and integrates that polynomial analytically on ([0,r_0]). All other radial integration intervals, six B iterations, matter fields and source terms are unchanged. Both raw and unchanged centre-projected states are measured. This is a diagnostic counterfactual, not a production change.

## Central vendor Hamiltonian residual

| N | Production raw | Production after projection | Shadow raw | Shadow after projection |
|---:|---:|---:|---:|---:|
| 40 | (+3.76922\times10^{-3}) | (+6.89618\times10^{-4}) | (-3.79483\times10^{-5}) | (+4.21456\times10^{-5}) |
| 80 | (+3.98859\times10^{-3}) | (+1.00788\times10^{-3}) | (+4.45421\times10^{-5}) | (+3.37230\times10^{-4}) |
| 160 | (+4.04156\times10^{-3}) | (+1.07904\times10^{-3}) | (+6.25329\times10^{-5}) | (+4.02471\times10^{-4}) |
| 320 | (+4.05429\times10^{-3}) | (+1.09606\times10^{-3}) | (+6.64726\times10^{-5}) | (+4.18006\times10^{-4}) |

At (N=320), the shadow reduces the raw central H residual from (4.0543\times10^{-3}) to (6.6473\times10^{-5}). After the existing centre projection, it reduces H from (1.0961\times10^{-3}) to (4.1801\times10^{-4}). Thus the first-half-cell convention is a strong contributor, but the residual after projection remains nonzero and appears to approach a second plateau.

## What this does and does not establish

- The raw-data residual drops by roughly two orders of magnitude in the shadow. That strongly supports the hypothesis that near-origin integration treatment matters in the production radial reconstruction.
- After projection, the residual drops by about (0.6781\times10^{-3}) at N=320, leaving (0.4180\times10^{-3}). The existing projection / surrounding discrete reconstruction remains implicated as an additional source of mismatch.
- The central connection constraint also decreases under the shadow, to (-3.3659\times10^{-6}) at N=320, but it does not exactly vanish. The residual is not safely attributable to (C_\Lambda) alone.
- Every reported algebraic reconstruction gate passed. No physical source, fitted coefficient, gauge, projection, or production default was changed.

## Next discriminating diagnostic

The current shadow corrects only the interval ([0,r_0]); intervals between cell centres still use a trapezoid on (src(r)). Given (src(r)=1+r^2F(r)), the (r^2F) piece also has a known regular structure. Next compare a second shadow in which:
1. ([0,r_0]) uses the same analytic origin expansion;
2. every interval ([r_i,r_{i+1}]) integrates (1+r^2F(r)) analytically, with (F) interpolated linearly in (r^2) between those sample points.

This isolates the remaining radial-quadrature compatibility question from the existing centre projection. The new method must remain diagnostic-only; even if H improves, a production change would still require separate convergence and compatibility checks.
