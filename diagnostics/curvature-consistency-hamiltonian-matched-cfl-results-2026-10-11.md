# Matched-CFL Hamiltonian spatial refinement — results

Date: 2026-10-11  
Branch: `diag/independent-metric-ricci-audit`  
Workflow: [run 38097951623](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38097951623)  
Commit: `e8a5a4fcdc2442c7e58c045c454f8fc6c6c57300`  
Artifact: [independent-metric-ricci-hamiltonian-spatial-refine](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38097951623/artifacts/11686677795)  
Artifact ZIP SHA-256: `f699bc6a3a387aebeb8a24025e0bc5474d96ae737b67904a79fd087b43765a85`

## Configuration and gates

- Resolutions (N=40,80,160); fixed (r_{\max}=40)
- Fixed CFL = 0.015; (\Delta t = 0.015\Delta r)
- Amplitude 0.01, width 7, D amplitude (10^{-10}), radiation on
- Samples (t=0,1,2,3)
- All initial data admitted using the diagnostic discrete-consistent initializer.
- Both the vendor-H reconstruction and same-state curvature-replacement identities closed at floating-point tolerance. Zero diagnostic gate failures.
- Diagnostic only: no production physics, gauge, source, projection, or default changes.

## Central residual at t=3

| N | Signed (H_{vendor}) | Signed (H_{metric}) | (R_{vendor}-R_{metric}) |
|---:|---:|---:|---:|
| 40 | (-6.44214\times10^{-5}) | (-8.35424\times10^{-5}) | (+1.91211\times10^{-5}) |
| 80 | (-1.06604\times10^{-5}) | (-1.39482\times10^{-5}) | (+3.28772\times10^{-6}) |
| 160 | (-1.93749\times10^{-6}) | (-2.53183\times10^{-6}) | (+5.94340\times10^{-7}) |

Absolute vendor (H) decreases by factors 6.04 (N40→N80) and 5.50 (N80→N160); the curvature gap decreases by factors 5.82 and 5.53. At (t=3), the center is also the global maximum-(|H|) cell for each resolution. The metric-counterfactual magnitude is larger because the metric-derived curvature difference has the same sign as the negative vendor residual at this sample.

At (t=0), the vendor constraint residual is at the solver tolerance. The metric-derived curvature counterfactual has a nonzero initial residual because it is a different discrete curvature reconstruction; this is not evidence that either expression is physically preferable.

## Cross-check against fixed-coordinate-timestep refinement

The previous matched-configuration run held (\Delta t=0.0025) fixed, while this run held CFL fixed. Both sets show a strong decrease of the vendor Hamiltonian residual as N rises. This is evidence that the observed short-time residual in this weak-D candidate-initial-data experiment is resolution-sensitive and not solely an artifact of one timestep-selection rule.

However, neither experiment is a pure spatial-only or formal convergence-order study:
- fixed (\Delta t) changes CFL as (\Delta r) changes;
- fixed CFL changes (\Delta t) as (\Delta r) changes;
- the initializer solves against the vendor-discrete constraint;
- only three resolutions and short times are sampled.

## Interpretation and next step

The curvature mismatch accounts for a non-negligible portion of the (t=3) vendor H magnitude, but the residual remains after curvature replacement. Therefore the curvature-route gap is contributory in this counterfactual, not a complete explanation or a justified production repair.

The experiment is deliberately weak-D ((10^{-10})) and short-time ((t\le3)). It must not be conflated with the strong-D ((10^{-4})), (r_{\max}=80), (t=22.5) centre-residual investigation. Next, repeat the same-state curvature/Hamiltonian term comparison on that actual production-initialized strong-D run, with the existing repair enabled, so findings address the unresolved operating regime without changing the physical model.
