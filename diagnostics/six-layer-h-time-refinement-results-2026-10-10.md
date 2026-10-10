# Six-layer Hamiltonian time-refinement diagnostic — results

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`  
Source commit: \`d63cd01a05f10172c940a333742d53d02750536a\`  
Workflow: [run 38083376102](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38083376102)  
Artifact: \`six-layer-h-refinement\`

## Completion and checks

- Workflow: success; candidate-initializer tests passed.
- Candidate initialization at \(N=160\): 3 Newton iterations; final max residual \(3.3679656286\times10^{-14}<10^{-13}\); max Jacobian condition number \(4.4138843660\times10^4<10^{12}\).
- Initial central Hamiltonian residual: \(-3.2151364904\times10^{-14}\); full-grid maximum absolute residual \(3.3679656286\times10^{-14}\).
- Trajectory: 400 base steps to \(t=3\), plus 24 refinement branch steps; final central \(H=-2.2892141115\times10^{-6}\); \(\tau=2.9771256257\); zero cycle events and handoffs.
- Maximum full-grid closure error across all six-group endpoints: \(1.6940658945\times10^{-21}\).
- Admission remains \`NOT_ADMITTED_DIAGNOSTIC_ONLY\`. No production equation or default changed.

## Full-step versus half-step attribution differences

Values are the component contribution to \(\Delta H_{\mathrm{full}}-\Delta H_{\mathrm{half}}\) at the central cell. Components use the same six-group symmetric Shapley method as the accepted-step audit.

| Probe \(t\) | Net \(\Delta H\) difference | Spatial metric \((a,b,X)\) | \(\Lambda\) | \(K\) | Matter |
|---:|---:|---:|---:|---:|---:|
| 0.00 | \(+5.3114\times10^{-13}\) | \(-2.6781\times10^{-12}\) | \(+3.4325\times10^{-13}\) | \(-1.016\times10^{-15}\) | \(+2.8670\times10^{-12}\) |
| 0.75 | \(+1.1102\times10^{-9}\) | \(+1.1143\times10^{-9}\) | \(+5.4775\times10^{-12}\) | \(-6.1232\times10^{-15}\) | \(-9.5747\times10^{-12}\) |
| 1.50 | \(-8.3967\times10^{-11}\) | \(-9.3549\times10^{-11}\) | \(-1.1381\times10^{-12}\) | \(-1.0752\times10^{-14}\) | \(+1.0730\times10^{-11}\) |
| 2.25 | \(-1.1004\times10^{-9}\) | \(-1.1183\times10^{-9}\) | \(-6.5402\times10^{-12}\) | \(-7.0460\times10^{-15}\) | \(+2.4417\times10^{-11}\) |

The gauge and \(A_a\) contributions are zero or at roundoff at all these central-cell probes. The \(K\) contribution to the *refinement difference* is negligible. \(\Lambda\) is a meaningful part of the absolute per-step balance but contributes only a small part of the difference between step sizes.

At t = 0.75, the full/half versus half/quarter ratio is approximately 4.00 for the net Hamiltonian increment, 4.00 for the metric contribution, and 4.00 for \(\Lambda\), consistent with second-order temporal sensitivity. Similar ratios hold at t = 1.5 and 2.25 (net: 3.66 and 4.01; metric: 3.65 and 4.00). The t = 0 ratio is not informative because the net change is near the numerical floor and results from strong cancellation.

## Interpretation

The time-step sensitivity is localized primarily to the spatial metric group \((a,b,X)\), not to \(K\), \(A_a\), gauge, or \(\Lambda\). The matter contribution largely cancels the metric contribution in the absolute constraint increment, but it changes much less between step subdivisions.

This does **not** establish a metric-equation defect. The three metric variables are algebraically coupled by the regularity projection, so the grouped contribution can mix direct and interaction effects. The next diagnostic separates a, b, and X symmetrically while keeping the other registered groups intact.

## Next step

Run an eight-group Shapley decomposition for \((a)\), \((b)\), \((X)\), lapse/shift/gauge, \(A_a\), \(K\), \(\Lambda\), and matter on the same full/half/quarter endpoints at the four registered probe times. Retain full-grid closure checks and the diagnostic-only admission guardrail.
