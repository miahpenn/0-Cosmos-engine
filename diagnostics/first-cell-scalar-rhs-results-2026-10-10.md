# First-cell scalar RHS diagnostic — results

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`  
Source commit: \`5e9cd9b3b9b4b7305aed2731463612e1363c1dd2\`  
Workflow: [run 38082419466](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38082419466)  
Artifact: \`first-cell-scalar-rhs-diagnostic\` (report schema \`first_cell_scalar_rhs_and_step_attribution_v1\`)

## Admission and trajectory checks

- Workflow conclusion: success; candidate initializer admission tests passed.
- Solver: 11 iterations; final max residual \(5.2307550652\times10^{-14}\) against tolerance \(10^{-13}\); maximum Jacobian condition number \(4.4138847648\times10^4\), below the \(10^{12}\) limit.
- Initial Hamiltonian residual: central and full-grid maximum absolute value \(5.2307550652\times10^{-14}\), at solver-residual scale.
- Evolution: 400 accepted steps to \(t=3\); final central \(H=-2.2892138974\times10^{-6}\); \(\tau=2.9771256255\); zero cycle events and zero handoffs.
- Four 11-group probes completed. Maximum all-cell component-sum closure across the probed accepted steps: \(3.7965108574\times10^{-20}\).
- Admission remains \`NOT_ADMITTED_DIAGNOSTIC_ONLY\`. The production equations, default initializer, and production branch were not changed.

## Central-H attribution for the next accepted step

All component values below are Shapley attributions of the measured accepted-step increment, not independent physical forces.

| Probe start \(t\) | S/PS | Metric \((a,b,X)\) | \(\Lambda\) | \(K\) | \(\phi/\Pi\) | Net \(\Delta H\) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | \(-1.98947\times10^{-8}\) | \(+2.01593\times10^{-8}\) | \(+6.51807\times10^{-9}\) | \(-1.00701\times10^{-8}\) | \(+3.07501\times10^{-9}\) | \(-2.15726\times10^{-10}\) |
| 0.75 | \(-1.94408\times10^{-6}\) | \(+1.18259\times10^{-6}\) | \(+7.68175\times10^{-7}\) | \(-1.07845\times10^{-8}\) | \(+3.07732\times10^{-9}\) | \(-8.86556\times10^{-10}\) |
| 1.50 | \(+8.91379\times10^{-7}\) | \(-8.16120\times10^{-7}\) | \(-7.23841\times10^{-8}\) | \(-1.14780\times10^{-8}\) | \(+2.75758\times10^{-9}\) | \(-6.65104\times10^{-9}\) |
| 2.25 | \(+2.69068\times10^{-6}\) | \(-1.90705\times10^{-6}\) | \(-7.84342\times10^{-7}\) | \(-1.06402\times10^{-8}\) | \(+2.58790\times10^{-9}\) | \(-1.06085\times10^{-8}\) |

Dark-matter contributions range from \(+1.56\times10^{-12}\) at \(t=0\) to \(-1.60\times10^{-9}\) at \(t=2.25\); baryons and radiation are smaller at these probes. The A_a and lapse/shift/gauge attributions are at or near zero. The D/PD attribution is zero at the report's floating-point precision for all four central-cell probes.

## Direct scalar RHS observations at cell 0

The S/PS RHS changes sign as the field evolves through the initial pulse:

| Start \(t\) | \(S\) | \(\partial_t S\) | \(P_S\) | \(\partial_t P_S\) |
|---:|---:|---:|---:|---:|
| 0.00 | \(+9.99681\times10^{-3}\) | \(0\) | \(0\) | \(-1.13667\times10^{-2}\) |
| 0.75 | \(+6.98707\times10^{-3}\) | \(-7.47593\times10^{-3}\) | \(-7.53102\times10^{-3}\) | \(-7.62252\times10^{-3}\) |
| 1.50 | \(+3.05707\times10^{-5}\) | \(-1.01374\times10^{-2}\) | \(-1.04189\times10^{-2}\) | \(+1.92744\times10^{-4}\) |
| 2.25 | \(-6.82644\times10^{-3}\) | \(-7.24570\times10^{-3}\) | \(-7.30980\times10^{-3}\) | \(+7.74434\times10^{-3}\) |

The D and PD states remain \(O(10^{-10})\) through these probes; their RHS values are nonzero but their contribution to the central Hamiltonian increment is not resolved in this diagnostic. The \(\phi/\Pi\) constraint attributions are only about \(2.6\)–\(3.1\times10^{-9}\) per probe.

## Interpretation

1. The largest terms are the S/PS scalar pair and the geometric metric/\(\Lambda\)/K response. Their substantial cancellation leaves the much smaller net constraint increment.
2. This run does **not** identify D/PD as the source of the central Hamiltonian drift at the sampled steps. It does not rule out D effects elsewhere, later, or in other observables.
3. The result does not prove the S equation or a geometry equation is wrong. Shapley allocation is a symmetric decomposition of the observed discrete constraint increment; closure verifies accounting only.
4. No physics modification, residual repair, damping, floor, clamp, or initializer admission follows from this run.

## Next diagnostic question

Audit the S/PS RHS and its scalar stress-energy projection at cells 0–4 over the same four probe steps. Separate the scalar-field update at fixed pre-step geometry from the full accepted-step change, then compare this with the geometric constraint response under the existing integrator. Use this only to test discrete consistency and expected truncation behavior; do not alter equations unless a reproducible operator-level defect is exposed.
