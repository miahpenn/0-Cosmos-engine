# Production machine continuation — Lambda/connection constraint audit (2026-10-12)

## Objective and admission status

Goal: get the full V5.5 production machine executing reproducibly while locating the earliest source of its evolving constraint failure. No production physics equation, source coefficient, gauge choice, boundary condition, or physical interpretation was changed by the diagnostics below.

The bounded production campaign through t=12 has completed at N=40,60,80, r_max=80, CFL=.03, radiation enabled:
- Production campaign: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687
- Artifact: `v55-production-t12-n40-n60-n80`, SHA-256 `dbbaa3be94242cdeecd437d0a41c00d4cf3f2d625e1f3aaaad9e1b45833f2047`
- Strict campaign completion checks passed. This does not admit constraint preservation or mass-balance closure.

## Findings preserved from the signed-profile / curvature audits

The same-state checkpoint audits covered the production initial slice and the saved t=2,4,6,8,10,12 checkpoints at all three resolutions, without advancing the state.

1. At t=0, the first admitted inner Hamiltonian cell is close to closure; by t=12, the admitted inner H residual is about -1.072e-2 at N=60, r=3.333, and -1.962e-2 at N=80, r=2.5.
2. The same cells have momentum residuals about -9.817e-4 (N=60) and -2.689e-3 (N=80), and metric-connection residuals about -1.1262e-2 and -1.4489e-2. The full-domain largest peaks are in the first two cells and should not be confused with these admitted inner cells.
3. A frozen-state center-regularity projection improves determinant and some first-cell connection/momentum witnesses but makes H worse at the first admitted inner cell. It is not a justified production repair.
4. Independent metric-derived Ricci splits the curvature gap into evolved-connection contribution and remaining metric-route/discretization difference, with roundoff closure. At the t=12 admitted cells, the metric-Ricci counterfactual reduces but does not close H. Therefore neither “Lambda alone” nor the independent metric route alone explains the constraint defect.

Relevant runs:
- Signed source/constraint profiles: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209173946
- Independent metric-Ricci decomposition: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209394006
- Signed Hamiltonian source/increment attribution: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209752406
- One-step full production stage trace: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209727710

## Saved-state Lambda RHS and metric-connection tangent audit

Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210572528
Artifact SHA-256: `8b5f7f675b9af6e82d93e6fce7cf97dbe030d41963cac4b44dd8a4f5b8e5ea75`

The metric-derived connection is
`Lambda_metric = D_r a/(2a^2) - D_r b/(ab) + 2/r (1/b - 1/a)`,
using the pinned fourth-order radial derivative. The analytic directional derivative with respect to the production explicit a,b RHS agrees with a centered finite-difference check at roughly 1e-11 or better. The actual Lambda increment in an observed production step agrees with the archived split `lambda_l2 + lambda_l3` to about 6e-17, so the implementation is applying its declared Lambda RHS consistently.

At t=0, the vector-field difference `Cdot = Lambda_rhs - d_t(Lambda_metric)` is near roundoff. The defect then grows during evolution. At the first admitted inner cell at t=12:
- N=60: `Cdot ≈ -1.9255e-3`
- N=80: `Cdot ≈ -5.2978e-3`

These are not by themselves proof that the Lambda equation is wrong; the connection residual is a propagated Einstein constraint witness and can be sourced by the momentum constraint.

## Critical attribution: the momentum constraint predicts almost all of Cdot

For the pinned production equations, the vendor Lambda L2 equation contains
`lambda_m * alpha * M_geom / a`, while the matter correction in Lambda L3 contributes
`-8*pi*lambda_m*alpha*J/a`. Their sum is exactly
`lambda_m*alpha/a * (M_geom - 8*pi*J)`, i.e. the momentum-constraint residual multiplied by the known factor.

At the same t=12 admitted cells, using the saved checkpoint's actual `2*alpha/a`:
| N | Momentum residual | `(2 alpha/a) × momentum residual` | Measured `Cdot` |
|---:|---:|---:|---:|
| 60 | -9.8167e-4 | -1.9483e-3 | -1.9255e-3 |
| 80 | -2.6894e-3 | -5.3380e-3 | -5.2978e-3 |

The momentum source alone reproduces the measured `Cdot` to within 1.2% and 0.8%. This is the strongest current localization: the connection residual is tracking a growing momentum-constraint defect rather than independently diagnosing a faulty Lambda update.

### One-step stage evidence

The one-step trace at the saved t≈10 states shows the admitted inner H and connection witnesses change during the primary `Aa/K/Lambda` update while a,b,X remain fixed. The same step begins with:
- N=60, i=2: H=-7.9694e-3, connection=-7.5948e-3.
- N=80, i=2: H=-1.2805e-2, connection=-1.2224e-2.

The primary stage pushes both more negative. At t=0, the same one-step control leaves the connection residual only around 1e-7. Thus the defect accumulates rather than appearing as a large first-step jump.

## Fixed-grid temporal refinement

Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38212493285
Artifact SHA-256: `f8ed81f0830ff66322f09b6137a33636d3f8837d953b9e9c380cf1f4d6617a9f`

Matched production runs: N=80, r_max=80, radiation enabled, through t=2.1, CFL=.03 vs .015, with aligned saved times every .15. Both runs passed strict campaign completion.

At t=2.1:
- i=2 connection residual: 3.55561e-4 (CFL .03) vs 3.55797e-4 (CFL .015), difference 2.37e-7.
- i=2 Hamiltonian residual: -9.92681e-4 vs -9.92549e-4, difference 1.33e-7.
- Fixed inner-region peak connection norm (including excluded center cells): 1.25343e-3 vs 1.25095e-3, a difference of about 0.2%.
- Fixed inner-region peak H norm (including excluded center cells): 8.16050e-3 vs 8.15104e-3, about 0.12%.

The same spatial-grid defect is nearly unchanged when halving dt. This strongly argues against simple time-step truncation as the main explanation for early growth, though two CFLs alone are not a formal time-convergence order proof.

## Current next test: same-domain spatial refinement

Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38212665185 (in progress at time of this document update)

Planned configuration: N=40,60,80,120; r_max=80; CFL=.03; radiation enabled; t=2.16; checkpoint interval=.24. Because the nominal steps are .06,.04,.03,.02, all runs reach common sample times without a tiny final step. It reports signed profiles at fixed physical radii and inner/exterior masks using the same r=5 cutoff.

Decision rules:
- If the admitted-inner momentum/H/C residual improves predictably with spatial refinement, quantify that order and retain the production equations while improving spatial resolution/stencil compatibility.
- If it does not improve, inspect the discrete momentum constraint propagation and exact center stencil, not the matter-source coefficients.
- Do not add the center projection as a repair, alter the Lambda equation, add damping, or launch t=160 until this evidence is reviewed.
