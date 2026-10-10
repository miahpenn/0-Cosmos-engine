# Constraint-preserving a-b/X Hamiltonian refinement — results

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`  
Source commit: \`714645d2442490dcdce4a39fd6b590025f4e2f81\`  
Workflow: [run 38083927745](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38083927745)  
Artifact: \`constraint-preserving-ab-x-h-refinement\`

## Admission and numerical checks

- Workflow: success; candidate initializer tests passed.
- Newton solve: 3 iterations; final maximum residual \(3.3679656286\times10^{-14}<10^{-13}\); maximum condition number \(4.4138843666\times10^4<10^{12}\).
- Initial central and all-grid maximum \(|H|\): \(3.3679656286\times10^{-14}\).
- Maximum error in \(a\,b^2=1\), across all checked endpoints: \(2.2204460493\times10^{-16}\).
- Trajectory: 400 base steps plus 24 refinement branch steps; \(t=3\), \(\tau=2.9771256257\); zero cycle events and handoffs.
- Maximum full-grid attribution closure: \(2.1175823681\times10^{-21}\).
- Admission remains \`NOT_ADMITTED_DIAGNOSTIC_ONLY\`.

## Refinement result

Central-cell contributions to \(\Delta H_{\rm full}-\Delta H_{\rm half}\):

| Probe \(t\) | Net | Paired \((a,b)\) | \(X\) | \(\Lambda\) | Matter |
|---:|---:|---:|---:|---:|---:|
| 0.00 | \(+5.31\times10^{-13}\) | \(-2.91\times10^{-12}\) | \(+2.33\times10^{-13}\) | \(+3.43\times10^{-13}\) | \(+2.87\times10^{-12}\) |
| 0.75 | \(+1.1102\times10^{-9}\) | \(+1.0103\times10^{-9}\) | \(+1.0406\times10^{-10}\) | \(+5.48\times10^{-12}\) | \(-9.57\times10^{-12}\) |
| 1.50 | \(-8.3967\times10^{-11}\) | \(-8.1244\times10^{-11}\) | \(-1.2304\times10^{-11}\) | \(-1.14\times10^{-12}\) | \(+1.0730\times10^{-11}\) |
| 2.25 | \(-1.1004\times10^{-9}\) | \(-1.0134\times10^{-9}\) | \(-1.0491\times10^{-10}\) | \(-6.54\times10^{-12}\) | \(+2.44\times10^{-11}\) |

At \(t=0.75,1.5,2.25\), the full-minus-half to half-minus-quarter magnitude ratios for the net are approximately 4.00, 3.66, and 4.01. Paired \((a,b)\) and \(X\) show similar second-order refinement ratios. Gauge and \(A_a\) differences are zero/roundoff; \(K\) differences are negligible; \(\Lambda\) differences are small.

## What this does and does not say

The dominant **accounting allocation** in the refinement difference is the paired \((a,b)\) group, with \(X\) smaller. This result preserves \(a\,b^2=1\) in every hybrid, correcting the previous unconstrained split.

It does not prove an \((a,b)\) evolution-equation defect or physical causality. Shapley attribution is an endpoint-residual bookkeeping method; all variables remain coupled in the actual trajectory. The \(t=0\) refinement ratio is uninformative because the net differences are near cancellation/roundoff.

## Next step

Inspect and instrument the actual production RHS/update and regularity projection for \(a,b,X\) at the central cell and first few radial cells. Compare the one-full-step endpoint with two-half and four-quarter endpoints, recording pre-projection RHS, post-projection values, projection deltas, and the \(a\,b^2=1\) residual. This is more direct than further hybrid attribution and can determine whether the step-size sensitivity enters through the RHS or the projection. Keep the diagnostic-only guardrail; do not change the production equations.
