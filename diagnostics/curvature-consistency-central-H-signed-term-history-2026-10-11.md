# Central Hamiltonian signed-term evolution comparison

Date: 2026-10-11  
Source: [metric Ricci Hamiltonian residual comparison, run 38096655385](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38096655385)  
Branch: `diag/independent-metric-ricci-audit`

This note reads the central-cell signed terms from the same-state Hamiltonian comparison. These are diagnostic measurements on the candidate initial data and evolved states; they do not indicate a production equation defect by themselves.

## Signed central term changes from t=0 to t=3

Convention:
[
H = R - (A_a^2+2A_b^2) + \frac{2}{3}K^2 - 16\pi\rho
]
The run reports the `extrinsic_A` term as the negative quadratic contribution, `trace_K` as (+2K^2/3), and `matter_source` as (-16\pi\rho).

| Resolution | (\Delta R_{vendor}) | (\Delta(2K^2/3)) | (\Delta(-16\pi\rho)) | Final central H |
|---:|---:|---:|---:|---:|
| N=40 | (-2.17614\times10^{-4}) | (-4.29785\times10^{-6}) | (+1.57785\times10^{-4}) | (-6.41271\times10^{-5}) |
| N=80 | (-1.65325\times10^{-4}) | (-4.27980\times10^{-6}) | (+1.59112\times10^{-4}) | (-1.04920\times10^{-5}) |
| N=160 | (-1.57022\times10^{-4}) | (-4.27868\times10^{-6}) | (+1.59427\times10^{-4}) | (-1.87279\times10^{-6}) |

The change in the (A_a,A_b) quadratic term is negligible at the displayed scale in these runs. The signed terms nearly cancel: as (R) falls, the negative matter-source term becomes less negative. A smaller residual remains from the imperfect cancellation and the decrease in the (K^2) term.

The (K)-term change is remarkably similar across resolutions, approximately (-4.28\times10^{-6}), while the curvature/matter mismatch gets much smaller as resolution increases. The final central constraint residual also falls strongly with resolution.

## Interpretation guardrails

- This signed-term table is accounting on the computed states, not proof of causality or an individual equation being wrong.
- (R_{vendor}) and the total matter density are coupled through the evolution; their changes should not be treated as independent physical mechanisms.
- The independently reconstructed curvature counterfactual shows the curvature-route gap contributes around 30% of the final central vendor residual at (t=3), but it does not remove the remaining residual.
- The current comparison uses fixed (\Delta t=0.0025) across (N=40,80,160); effective CFL therefore varies with resolution.
- No production equations, gauge, sources, projection, or defaults have changed.

## Next discriminating test

Measure the same Hamiltonian term history on (N=80) from a single initial state at CFL (0.03,0.015,0.0075). This tests how much of the signed residual history is timestep-sensitive before considering any production correction.
