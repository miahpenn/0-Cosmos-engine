# Independent metric-Ricci Hamiltonian comparison — results

Date: 2026-10-11  
Branch: `diag/independent-metric-ricci-audit`  
Workflow: [run 38096655385](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38096655385)  
Commit: `0e8b588c41c683904ba80499e9e619468def3db8`  
Artifact: [independent-metric-ricci-hamiltonian-compare](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38096655385/artifacts/11686766093)  
Artifact ZIP SHA-256: `9cc6bd0ef8967e5bfc5b22a434ca5c420ce63095198023614b3a4cbb2ea64dcf`

## Design

On identical evolved states at (N=40,80,160), fixed (r_{\max}=40), and fixed coordinate timestep (\Delta t=0.0025), compute:

- (H_{vendor}=R_{vendor}+T_{A}+T_K-16\pi\rho)
- (H_{metric}=R_{metric}+T_A+T_K-16\pi\rho)

Only the curvature scalar is replaced in the diagnostic counterfactual. Extrinsic-curvature terms and the matter density are identical pointwise. This is not a modified production constraint.

The independent metric-derived curvature uses (R_{areal}=r\sqrt{b}/X), proper radial distance, and does not use evolved `Lambda` or assume areal-coordinate gauge.

## Results at (t=3), centre cell

| N | Signed (H_{vendor}) | Signed (H_{metric}) | (R_{vendor}-R_{metric}) | Gap / (|H_{vendor}|) |
|---:|---:|---:|---:|---:|
| 40 | (-6.41271\times10^{-5}) | (-8.31831\times10^{-5}) | (+1.90560\times10^{-5}) | 29.7% |
| 80 | (-1.04920\times10^{-5}) | (-1.37346\times10^{-5}) | (+3.24264\times10^{-6}) | 30.9% |
| 160 | (-1.87279\times10^{-6}) | (-2.45009\times10^{-6}) | (+5.77301\times10^{-7}) | 30.8% |

At (t=0), the vendor residual is at solver tolerance ((\le 3.4\times10^{-14}) over the grid), as the initializer solves against the vendor discrete constraint. The metric-derived counterfactual has a nonzero residual because its curvature discretization differs from the one used by the initializer. Therefore the (t=0) counterfactual is not a validation of one curvature formula over the other.

At (t=1) and (t=2), the central curvature gap changes sign and its fraction of the vendor residual is non-monotonic; do not extrapolate the (t=3) ~30% fraction to all times. The signed history at (N=160):
- (t=1): (H_{vendor}=-9.8020\times10^{-8}), (H_{metric}=+3.0060\times10^{-7}), curvature gap (-3.9862\times10^{-7}).
- (t=2): (H_{vendor}=-6.3521\times10^{-7}), (H_{metric}=-5.3429\times10^{-7}), curvature gap (-1.0091\times10^{-7}).
- (t=3): (H_{vendor}=-1.8728\times10^{-6}), (H_{metric}=-2.4501\times10^{-6}), curvature gap (+5.7730\times10^{-7}).

## Gates

Workflow succeeded. No diagnostic gate failures. Vendor Hamiltonian-term reconstruction and the curvature-replacement identity closed within floating-point tolerance. The independent curvature mismatch and the Hamiltonian residual both decrease strongly as resolution is increased in this experiment.

## Interpretation

1. The curvature formulation difference is numerically material: at (t=3), replacing only (R) changes the central Hamiltonian residual magnitude by about 30%.
2. It does not explain all the evolved residual; the remaining terms still yield a nonzero (H_{metric}).
3. Since the initial solve targets the vendor discrete constraint, the initialization does not independently adjudicate between curvature discretizations.
4. This is a controlled counterfactual, not authorization to change the vendor formula, production equations, or initial-data acceptance.
5. Fixed timestep means the effective CFL varies with (N), so the resolution factors are indicators, not formal convergence orders.

## Next diagnostic

Hold (N=80) and the initial state fixed while comparing (CFL=0.03,0.015,0.0075), this time measuring both (H_{vendor}) and (H_{metric}) plus their individual curvature, extrinsic-curvature, and matter terms. This isolates timestep sensitivity in the constraint residual itself instead of curvature alone. Keep the entire run diagnostic-only.
