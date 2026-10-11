# Strong-D production initialization stage trace — results

Date: 2026-10-11  
Branch: `diag/independent-metric-ricci-audit`  
Workflow: [run 38099589026](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38099589026)  
Commit: `de9107fd1b3f52e8e5d9c98e82f8c619f599a479`  
Artifact: [strong-D-production-initial-stage-trace](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38099589026/artifacts/11687054780)  
Artifact ZIP SHA-256: `80ac0f6c6edd1e21e23c1ad546296cc1c9a9f15d33fbb14d2a9c92af8a057076`

## Exact stage comparison at N=160, r_max=80, D=1e-4

| Stage | Central vendor H | Central vendor R | Central metric-derived R | Central connection constraint |
|---|---:|---:|---:|---:|
| Raw `build_initial_data`, before projection | (4.04156\times10^{-3}) | (6.54869\times10^{-3}) | (4.16096\times10^{-3}) | 0 |
| After existing centre projection | (1.07904\times10^{-3}) | (3.58617\times10^{-3}) | (3.49992\times10^{-3}) | (-3.90041\times10^{-5}) |
| After CMC lapse/gauge setup | (1.07904\times10^{-3}) | (3.58617\times10^{-3}) | (3.49992\times10^{-3}) | (-3.90041\times10^{-5}) |
| Public `kernel.initialize` | (1.07904\times10^{-3}) | (3.58617\times10^{-3}) | (3.49992\times10^{-3}) | (-3.90041\times10^{-5}) |

All term and curvature-replacement identity gates passed at roundoff scale.

## What the stage deltas establish

- The centre projection changes central H by (-2.962523455\times10^{-3}), a reduction of roughly 73.3% from its raw pre-projection value.
- That delta is accounted for by the vendor curvature change: (Delta R=-2.962523447\times10^{-3}). The central extrinsic-A, K and fluid terms are unchanged to the shown precision; scalar contribution changes only by (7.63\times10^{-12}).
- The CMC lapse and gauge reset produce zero measured change to H and the reported constraint terms.
- The public `kernel.initialize` result reproduces the manually staged post-projection/gauge geometry exactly at the recorded precision.
- Before projection, the metric-derived Ricci gap at the centre is (R_{vendor}-R_{metric}=2.38773\times10^{-3}); after projection it is (8.62473\times10^{-5}). Projection therefore changes the curvature-route comparison dramatically as well as reducing the residual.
- The projection moves the centre connection constraint from zero to (-3.90\times10^{-5}). This is an observed consequence of the existing first-cell regularity extrapolation, not evidence by itself that the projection is incorrect.

## Interpretation guardrails

This trace rules out the lapse solve/gauge reset as the cause of the initial H change. It also shows that the centre projection reduces, rather than creates, the central H residual in this configuration. It does not establish that the remaining residual is purely a spatial truncation error or that the projection should be changed: the raw initial geometry is not the regularized production slice, and the Ricci formulations behave differently on it.

Next compare raw and post-projection initial H profiles at N=40,80,160,320 with fixed (r_{\max}=80), fixed physical settings, and the existing projection unchanged. Resolution behavior will help distinguish a shrinking first-cell discretization effect from a persistent mismatch in the radial reconstruction/discrete Hamiltonian constraint.
