# Spatial-metric a/b/X Hamiltonian refinement — results

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`  
Source commit: \`14575a5e0469b2326b93ef508f1b53cfa249eb1d\`  
Workflow: [run 38083630334](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38083630334)  
Artifact: \`spatial-metric-eight-group-refinement\`

## Completion and checks

- Workflow conclusion: success; candidate initializer tests passed.
- Newton solve: 11 iterations; final max residual \(5.2307550652\times10^{-14}<10^{-13}\); max condition number \(4.4138847648\times10^4<10^{12}\).
- Initial central and all-cell maximum \(|H|\): \(5.2307550652\times10^{-14}\).
- Trajectory: 400 base steps and 24 refinement branch steps, final \(t=3\), central \(H=-2.2892138974\times10^{-6}\), \(\tau=2.9771256255\), zero cycle events and handoffs.
- Maximum full-grid attribution closure: \(4.2351647363\times10^{-21}\).
- Admission remains \`NOT_ADMITTED_DIAGNOSTIC_ONLY\`.

## Full-step versus half-step component differences

Central-cell contribution to \(\Delta H_{\rm full}-\Delta H_{\rm half}\):

| Probe \(t\) | Net | \(a\) | \(b\) | \(X\) | \(\Lambda\) | Matter |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | \(+6.012\times10^{-13}\) | \(-9.693\times10^{-13}\) | \(-1.871\times10^{-12}\) | \(+2.326\times10^{-13}\) | \(+3.433\times10^{-13}\) | \(+2.867\times10^{-12}\) |
| 0.75 | \(+1.1102\times10^{-9}\) | \(+3.5179\times10^{-10}\) | \(+6.5847\times10^{-10}\) | \(+1.0409\times10^{-10}\) | \(+5.4775\times10^{-12}\) | \(-9.5747\times10^{-12}\) |
| 1.50 | \(-8.3994\times10^{-11}\) | \(-2.8280\times10^{-11}\) | \(-5.2963\times10^{-11}\) | \(-1.2332\times10^{-11}\) | \(-1.1381\times10^{-12}\) | \(+1.0730\times10^{-11}\) |
| 2.25 | \(-1.1003\times10^{-9}\) | \(-3.5275\times10^{-10}\) | \(-6.6057\times10^{-10}\) | \(-1.0488\times10^{-10}\) | \(-6.5402\times10^{-12}\) | \(+2.4417\times10^{-11}\) |

At the later probes, the magnitude ratio of the full-minus-half to half-minus-quarter differences is near 4 for a, b, and X individually (about 3.6–4.0), consistent with second-order time-step sensitivity. Gauge and \(A_a\) contributions are at zero/roundoff; \(K\) is negligible for the refinement difference; \(\Lambda\) is much smaller than the spatial metric contributions.

## Algebraic-constraint caveat

The production regularity projection enforces \(a\,b^2=1\). Independent a-only or b-only hybrid states violate that relation, so the individual Shapley values are not independent physical forces or a causal ranking. They show where the chosen numerical decomposition allocates the endpoint residual change, and the split closes exactly as an accounting identity. The total \((a,b,X)\) attribution is the robust observation; the individual a/b split requires a constraint-preserving comparison.

## Next diagnostic

Repeat the matched full/half/quarter analysis with \((a,b)\) as one group and \(X\) as a separate group. Every hybrid then takes a and b together from the same endpoint, preserving \(a\,b^2=1\) while retaining the previously registered gauge, \(A_a\), \(K\), \(\Lambda\), and matter groups. Continue to make no physics or production-default changes.
