# 0-Cosmos V5.5 continuation handoff — discrete constraint propagation
**Date:** 2026-10-12  
**Scope:** diagnostics and cross-agent reconciliation; no production physics changes  
**Status:** production execution works for the bounded campaign; scientific admission remains OPEN.

## Executive finding

The current production machine can execute the bounded V5.5 campaign through t=12 at N=40,60,80, but that operational success is not a finding of constraint preservation, mass-balance closure, or physical cycle/bounce behavior. The active defect is constraint growth, led by the momentum constraint in the inner region. Several independent tests now exclude simple explanations, but no single production equation has yet been proved wrong.

The strongest current discriminator is a frozen t=0 state whose discrete Hamiltonian, momentum, and metric-connection constraints close to roundoff. Under the unchanged production semi-discrete RHS, its momentum-constraint tangent is still materially nonzero and agrees closely with the measured one-step quotient. This separates the propagation problem from the initial Hamiltonian residual alone.

**Production source remains frozen at** `6d8b3f2790644591326df57c55908863f994cb0c`. The pinned vendor revision recorded by the tangent report is `d6052d605673ce9d82cdc99b0df79855fa2ea215`. No diagnostic result in this handoff is permission to merge a candidate into production.

## Guardrails for all follow-up work

- Run the full machine when the test calls for evolution; for source attribution, use frozen states and bounded one-step tests first.
- No hand tuning, fitted coefficients, arbitrary damping, floors/clamps, artificial bounce/stop/reset, invented D-to-matter law, or silent physics changes.
- Keep candidate initializers, quadrature, lapse, curvature and RHS variants on diagnostic branches.
- Report signed term budgets and algebraic closure, not norms alone.
- Keep center-adjacent cells (especially i=0,1,2) separate from first admitted inner cells and from fixed physical radii.
- A workflow marked completed or an operational admission flag is not scientific admission.
- Do not infer a formal convergence order from the present short, non-monotone fixed-radius samples.
- No interpretation of bounce, turnaround, or cycle behavior is admitted by these local constraint tests.

## Production campaign and earlier constraint-propagation evidence

Bounded production campaign: [workflow 38208253687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687), artifact `v55-production-t12-n40-n60-n80`, SHA-256 `dbbaa3be94242cdeecd437d0a41c00d4cf3f2d625e1f3aaaad9e1b45833f2047`. It completed through t=12 with strict campaign completion checks, N=40,60,80, r_max=80, CFL=.03, radiation enabled. Constraint/mass-balance admission remains open.

The existing [Lambda/connection audit](https://github.com/miahpenn/0-Cosmos-engine/blob/62d67292054e1bb4f39ff122ccb1923faa0a7255/docs/LAMBDA_CONSTRAINT_PROPAGATION_AUDIT_2026-10-12.md) records:
- Inner H, momentum M and connection residuals grow during evolution, with the largest peaks near the center.
- Halving CFL from .03 to .015 on a fixed grid barely changes early residuals; simple temporal truncation is therefore not the leading explanation (not a formal time-order proof).
- On a fixed physical annulus, spatial refinement reduces volume-weighted residuals, while the peak tracks inward toward the first admitted cells.
- The Lambda connection-rate mismatch is largely predicted by the momentum residual times the known (2\alpha/a) factor. This makes Lambda a downstream witness rather than the leading suspect.
- A one-step trace and semi-discrete tangent agree closely, placing the tendency in the coupled semi-discrete RHS rather than solely in the time integrator.
- A frozen-state regularity projection improves determinant and some first-cell connection/momentum witnesses but worsens H at the first admitted inner cell; it is not a justified repair.

## Diagnostic A — regular-origin half-cell quadrature

Workflow: [38289853406](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289853406).  
Matched one-step propagation: [38290654079](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290654079).  
Diagnostic source: [origin_half_cell_quadrature_regularity_ab.py](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/origin-quadrature-regularity-2026-10-12/tools/origin_half_cell_quadrature_regularity_ab.py).

The baseline initializer integrates the origin-to-first-cell contribution as `(dr/2)*source[0]`. A diagnostic-only counterfactual replaces only that half-cell integral with the integral of the quadratic interpolant in (x=r^2) through the first three cell-centred source values:
`dr*(401/720*source[0] - 31/480*source[1] + 11/1440*source[2])`.
This is exact for source profiles proportional to 1, r² and r⁴. The subsequent inter-centre trapezoids, six nonlinear iterations, matter initialization, CMC lapse solve and other fields were retained.

At t=0 and N=40,60,80,120, the first-cell Hamiltonian residual falls from approximately 3.77e-3, 3.93e-3, 3.99e-3, 4.03e-3 to approximately 3.93e-4, 1.31e-4, 8.17e-5, 6.60e-5, respectively. This is strong evidence that the current half-cell initialization rule contributes materially to the initial central H defect.

Limitations: this is a counterfactual initial-data construction, not an evolution or production patch. Curvature/RHS regularity is not fully restored: the vendor (R_{tf}/r^2) witness does not show a clean finite center limit. The matched one-step changes to H/M/connection are mixed, including nonzero momentum differences. Therefore do not patch the production initializer from this result alone.

## Diagnostic B — smooth manufactured regular-profile operator audit

Concurrent-agent report recorded in [Issue #4](https://github.com/miahpenn/0-Cosmos-engine/issues/4), comment [6130056167](https://github.com/miahpenn/0-Cosmos-engine/issues/4#issuecomment-6130056167). Run: [38288419513](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38288419513). Artifact SHA-256: `83c738421ca014b781d5a07d1218d3de116ddc43bea5bc4d762233bc1dfb2d77`.

Smooth regular spherical manufactured profiles were evaluated without evolution at N=40/60/80/120. On the fixed annulus (3\le r\le6), RMS errors improve from 2.192e-6 to 3.009e-8 for M and 2.651e-6 to 3.658e-8 for C_lambda; descriptive orders are 3.927 and 3.923. Signed-term reconstruction closes at stored precision.

Interpretation: evidence against a generic derivative/parity-sign error or a generic bug in the constraint evaluator. It does not prove that the evolved production RHS satisfies the discrete constraint-propagation identity near the center.

## Diagnostic C — discrete-consistent initial metric and multi-constraint tangent

Workflow: [38290633971](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290633971). Artifact ID `11766405316`; ZIP SHA-256 `2011162ba822bea1b52d2b8c2099284f97f7a0ca5e76e2d2bdc1788c324a7186`. Report: [discrete-initial-constraint-tangent-2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/constraint-tangent-discrete-initial-2026-10-12/diagnostics/discrete-initial-constraint-tangent-2026-10-12.md).

Six frozen t=0 states (baseline/candidate at N=60,80,120) were compared with the actual semi-discrete production RHS. The candidate initial metric is built to satisfy discrete H/M/C_lambda to roundoff. Two symmetric RHS-direction perturbations were used, epsilon=2e-5 and 1e-5, with the CMC lapse re-solved; a separate copy took one ordinary production step.

Candidate initial max |H| is approximately 7.65e-16, 1.29e-15 and 3.16e-15 at N=60,80,120. Candidate H tangent is about 7e-9 in the reported inner/fixed regions. Yet M tangent remains nonzero, and the one-step quotient agrees with it within 0.7% at the first admitted cell:

| N | First admitted r | Candidate Mdot | One-step ΔM/dt |
|---:|---:|---:|---:|
| 60 | 3.333 | -8.7553e-5 | -8.6967e-5 |
| 80 | 2.500 | -1.2354e-4 | -1.2295e-4 |
| 120 | 1.667 | -1.9063e-4 | -1.8993e-4 |

Corrected region maxima from the later reconciliation: admitted-inner max |Mdot| is 8.76e-5 at N=60, 1.24e-4 at N=80, and 1.91e-4 at N=120; first-two-cell maxima are 1.61e-4, 2.05e-4 and 3.06e-4. The first admitted cell moves inward with resolution, so these are not fixed-radius convergence samples.

At N=120, the candidate first-cell signed budget gives geometric momentum tendency about -3.203e-4, scalar-current source about +3.873e-5, net M tendency about -2.816e-4. The geometric sum closes from radial-b, (D_r A_a), and (-(2/3)D_r K) terms to stored precision. The terms do not cancel as the continuum identity would require for this discrete state. Source assembly and geometry reconstruction close algebraically; this is not an arithmetic assembly failure.

At fixed radii, candidate Mdot is sign/resolution-sensitive near r=2, but approaches near-cancellation at r=3–5 for N=120. Do not equate that selected-radius cancellation with center regularity or long-time constraint preservation. The scalar-source derivative is dominated by local S and changes comparatively little between baseline/candidate; most A/B difference comes from geometric response.

Conclusion: correcting the initial H residual does not remove the momentum propagation tendency. The leading next question is discrete compatibility of the metric/A_a/K response with the scalar-current source, particularly near the center. This does not yet identify a unique faulty production term.

## Diagnostic D — Ricci/Aa RHS identity and origin regularity

A separate report, artifact ID `11765776185`, run [38290099699](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290099699), compared vendor Ricci/Aa RHS terms with an independent metric-derived Ricci construction. At r=0.5, the reported production (R_{tf}) is -4.1518e-5 while the independent analytic route gives +2.8828e-5; the production (A_a) RHS is -1.8680e-5 versus +5.7310e-5 from the comparison. The Hessian difference is smaller (~-5.63e-6).

This flags the Ricci/geometry operator and regular-center radial-series identities for further scrutiny. It is not yet proof of an error: formula conventions, geometric definitions, and independent-route assumptions must be checked term-by-term before any conclusion or patch.

## Other counterfactuals that must remain separate

- **Initial-B Newton solve:** closes discrete H at t=0 by construction but does not consistently improve evolved constraints. At N=80, t=2.16 admitted-inner H changes 1.071e-3 to 1.018e-3; at N=120 it worsens 1.539e-3 to 1.674e-3. At N=120 M improves modestly (9.36e-4 to 9.09e-4) and C_lambda (2.56e-4 to 1.84e-4). Runs [38284401548](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548), [38285057887](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38285057887), [38287772729](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38287772729). Not a repair.
- **CMC center-row / production-matched fourth-order lapse counterfactuals:** improve the first-cell K-RHS mismatch algebraically, but both increase admitted-inner M norm after one normal step by roughly 3e-6. Do not patch CMC from the local identity alone.
- **Scalar lapse-gradient inverse-metric audit:** the declared (gamma_{rr}=a/X^2) implies (gamma^{rr}=X^2/a), while the archived RHS uses (X^{-2}). This is a legitimate covariance audit item, but the matched candidate campaign has tiny trajectory deltas and no consistent constraint improvement; it does not explain the present defect alone.
- **Areal-radius odd-parity audit:** effects concentrate in first two cells; diagnostic only, not a repair.
- **Radiation ON/OFF and D-envelope tests:** earlier runs did not change the early D-envelope result. These are separate from the present momentum tangent attribution.
- **Timeline join failure:** run 38285163101 raised `KeyError: (40, 0.72, 'baseline')` in the audit join/report layer. Treat as report coverage/join failure unless independently shown otherwise, not as a production-physics failure.

## Current synthesis

Evidence supports a localized central geometry/constraint-propagation layer:
1. Initial origin quadrature contributes to the initial central H defect.
2. Even a discrete-consistent initial H/M/C_lambda state develops a nonzero momentum tangent under the unchanged production RHS.
3. A generic derivative/parity-sign bug is disfavored by the manufactured-profile audit.
4. Lambda/connection residual mostly follows momentum residual.
5. Signed source and geometry budgets close arithmetically, but the geometric response does not cancel the scalar-current source at central/first-admitted cells.
6. Independent Ricci/Aa comparisons show a possible center-local geometry operator discrepancy, but conventions need validation before calling it a defect.

This remains a diagnosis, not a final causal attribution. No single production equation, coefficient, gauge condition, or boundary condition has been proven wrong.

## Next bounded test — recommended

Use paired frozen states: the baseline initializer and the discrete-consistent candidate at N=60/80/120; add the already saved t=2.16 states after the t=0 accounting is reproduced. No long trajectory initially.

1. Split the actual production (dot A_a) into lapse Hessian, Ricci, anisotropic-stress, nonlinear/extrinsic terms, and split (dot K) with matching signed contributions.
2. Check regular-center radial-series identities at i=0,1,2 and fixed physical radii r=2,3,4,5.
3. For every shared term, record the exact derivative stencil, parity, and geometry convention used by both the constraint and RHS.
4. Keep the scalar-current source budget explicit; test whether its discrete source is compatible with the discrete geometry response, rather than assuming the continuum identity automatically survives discretization.
5. Report signed values, residual after summation, closure error, and separate regions. Verify closure at roundoff before interpreting any mismatch.
6. Re-run the same bounded tangent and one-step comparator after any single demonstrated repair; only then consider a longer campaign.

Acceptance criteria: reproducible signed budget with exact algebraic closure; a minimal failing profile/state distinguishing operator/parity mismatch from physical-source mismatch; no production equation edits until a specific discrete RHS identity fails.

## Repository pointers

- Active issue and live agent reconciliation: [Issue #4](https://github.com/miahpenn/0-Cosmos-engine/issues/4).
- Lambda/connection audit: [docs/LAMBDA_CONSTRAINT_PROPAGATION_AUDIT_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/62d67292054e1bb4f39ff122ccb1923faa0a7255/docs/LAMBDA_CONSTRAINT_PROPAGATION_AUDIT_2026-10-12.md).
- Discrete-consistent tangent report: [diagnostic branch report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/constraint-tangent-discrete-initial-2026-10-12/diagnostics/discrete-initial-constraint-tangent-2026-10-12.md).
- Origin quadrature diagnostic script: [diagnostic branch script](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/origin-quadrature-regularity-2026-10-12/tools/origin_half_cell_quadrature_regularity_ab.py).

## Historical physics locks to preserve

This numerical diagnosis must not silently overwrite the broader GEAR/COSMOS specification: unified total stress-energy and geometry; GEAR-03 kappa=0, g=1; GEAR-63/65 local cubic geometric barrier; G=1; corrected COSMOS potential normalization (legacy V0=1 is control only); no fitted feedback coefficient, imposed clock conversion, D-to-matter conversion law, artificial bounce/stop/reset/ejection/lifetime threshold, or unsupported lambda_cycle. The exact spherical branch does not isotropize archived Bianchi-I shear. Any eventual physical-cycle interpretation must follow admitted invariant witnesses, not be imposed by the diagnostics.
