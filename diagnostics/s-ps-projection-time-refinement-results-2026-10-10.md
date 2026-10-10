# S/PS projection and time-refinement diagnostic — results

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`  
Source commit: \`620e2306c86854cfd6335cdbeb15f8abd2a0124d\`  
Workflow: [run 38082839524](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38082839524)  
Artifact: \`s-ps-projection-refinement\`

## Completion and checks

- Workflow conclusion: success; candidate initializer tests passed.
- Same N=160 candidate trajectory, final coordinate \(t=3\), 400 base steps and 24 additional refinement-branch steps.
- Newton solve: 11 iterations; final max residual \(5.2307550652\times10^{-14}<10^{-13}\); max condition number \(4.4138847648\times10^4<10^{12}\).
- Initial central and full-grid max absolute Hamiltonian residual: \(5.2307550652\times10^{-14}\).
- Final central \(H=-2.2892138974\times10^{-6}\), \(\tau=2.9771256255\), zero cycle events and zero handoffs.
- The symmetric S-field-density split closed exactly at reported precision (maximum closure 0.0).
- No physics equation or production default was changed. Admission remains diagnostic-only.

## Full / half / quarter step comparison

Each row starts from the same saved state and covers the same coordinate-time interval \(\Delta t=0.0075\). The half and quarter paths subdivide only the diagnostic interval; the full-step path continues the main trajectory.

| Probe \(t\) | \(\Delta H\), full | \(\Delta H\), 2 halves | \(\Delta H\), 4 quarters | \(|H_{full}-H_{half}|/|H_{half}-H_{quarter}|\) |
|---:|---:|---:|---:|---:|
| 0.00 | \(-2.15726\times10^{-10}\) | \(-2.16327\times10^{-10}\) | \(-2.19514\times10^{-10}\) | 0.189 (not asymptotic; near-zero change) |
| 0.75 | \(-8.86556\times10^{-10}\) | \(-1.99680\times10^{-9}\) | \(-2.27429\times10^{-9}\) | 4.001 |
| 1.50 | \(-6.65104\times10^{-9}\) | \(-6.56704\times10^{-9}\) | \(-6.54404\times10^{-9}\) | 3.651 |
| 2.25 | \(-1.06085\times10^{-8}\) | \(-9.50814\times10^{-9}\) | \(-9.23343\times10^{-9}\) | 4.005 |

At the three later probes, the full/half versus half/quarter difference ratio is close to 4, consistent with the expected second-order time-step error trend. The \(t=0\) ratio is not informative because the increment is exceptionally small and differences are near the numerical floor. These four local comparisons support time-discretization sensitivity; they do not prove that all of the \(t=3\) residual is temporal error.

## S-sector energy-density change

For the S component, the existing scalar-projection convention is
\[
\rho_S=\tfrac12 P_S^2+\tfrac12(X^2/a)(\partial_r S)^2+\tfrac12 S^2.
\]
Its measured change was split symmetrically into the S/PS field-pair change and the dependence on the metric factor \(X^2/a\):

| Probe \(t\) | Total \(\Delta\rho_S\), full step | Field-pair contribution | Metric contribution |
|---:|---:|---:|---:|
| 0.00 | \(+3.95680\times10^{-10}\) | \(+3.95792\times10^{-10}\) | \(-1.12034\times10^{-13}\) |
| 0.75 | \(+3.86761\times10^{-8}\) | \(+3.86762\times10^{-8}\) | \(-4.90570\times10^{-14}\) |
| 1.50 | \(-1.77334\times10^{-8}\) | \(-1.77334\times10^{-8}\) | \(-2.38416\times10^{-16}\) |
| 2.25 | \(-5.35295\times10^{-8}\) | \(-5.35294\times10^{-8}\) | \(-5.58502\times10^{-14}\) |

The S-sector density changes are extremely stable across full/half/quarter integration. At \(t=0.75\), for example, the full-step versus quarter-step difference in \(\Delta\rho_S\) is only about \(2.33\times10^{-13}\), while the central Hamiltonian increments differ by about \(1.39\times10^{-9}\). The same pattern is present at \(t=2.25\). Thus the time-step sensitivity of the Hamiltonian increment is not explained by a comparably large change in the S-sector energy-density update.

## Interpretation and next focus

This is consistent with the earlier Shapley report: S/PS changes carry the dominant matter-sector contribution, but the residual is the small remainder of cancellations involving the metric, \(\Lambda\), and \(K\). The time-step comparison indicates second-order temporal sensitivity at the later probes. It does not yet identify a faulty equation.

The next narrow diagnostic should symmetrically split the Hamiltonian increment into **all geometry variables versus all matter fields** for the same full/half/quarter endpoints. That will test whether the observed time-step sensitivity lives predominantly in the geometry update or the matter update before inspecting the metric/\(\Lambda\)/\(K\) update operators individually.
