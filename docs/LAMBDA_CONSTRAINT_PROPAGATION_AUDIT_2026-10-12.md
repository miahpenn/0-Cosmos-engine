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

## Spatial refinement result details

The fixed physical annulus values clarify that the error is spatially concentrated toward the centre:
- At r=3, signed H is +7.4135e-4 (N40), -1.8587e-4 (N60), -4.1919e-4 (N80), and +9.2175e-6 (N120).
- At r=4, signed H is +1.7024e-4, -3.0519e-4, +1.0214e-4, and -1.1709e-5.
- At r=5, signed H is -4.0086e-4, +1.0104e-4, -1.8623e-5, and -5.4585e-6.
These are fixed-radius interpolated samples; their signs and non-monotonic intermediate values mean a formal spatial order should not be fitted from this short test. Nevertheless, by N=120 the residual at r=3,4,5 is much smaller than at N40/N60/N80.

Volume-weighted RMS H over the fixed physical annulus 2<=r<5 is 7.4135e-4 (N40), 6.5404e-4 (N60), 4.4980e-4 (N80), and 1.2832e-4 (N120). By contrast, the maximum in admitted cells (index>=2 and r<5) shifts inward from r=3.333 (N60, |H|=7.51e-4) to r=2.5 (N80, |H|=1.071e-3) to r=1.667 (N120, |H|=1.539e-3). The full inner-domain peak including the first two excluded cells falls from 9.27e-3 (N60) to 6.97e-3 (N80) and 4.74e-3 (N120).

**Interpretation:** the residual improves at fixed resolved radii and in a volume-weighted inner annulus, while its maximum follows the first admitted cells inward. This is evidence for a shrinking central numerical layer, not a converged central solution. The momentum constraint is also not uniformly converged at the first admitted cell (at N120, |M| reaches about 9.36e-4 inside r<5), so the momentum-source mechanism remains the next target.

## One-step momentum-constraint term decomposition — complete

Latest stage trace:
- Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38213136600
- Artifact: `production-one-step-constraint-stage-trace`
- SHA-256: `f2097b7c1eecafabcfacb81e796ab0aba07d2468785f1c17731953b0c6749686`

The traced momentum constraint is reconstructed from signed terms:
`M_geom = D_r Aa - (2/3) D_r K + 6 Aa chi_r + (Aa-Ab)(2/r + D_r b/b)`,
with `chi_r=-D_r X/(2X)`, `Ab=-Aa/2`. The matter contribution is `8*pi*J`, and `M_residual=M_geom-8*pi*J`. The term closure was exactly zero in every captured stage.

At the first admitted inner cell i=2, the fresh t=0 state begins with M at roundoff. During the Euler predictor, the evolved matter/scalar state changes the current source before the geometric A/K stage responds:
- N=60: M source changes from approximately zero to `-9.6842e-6`, while M_geom remains approximately zero; temporary predictor M residual is `+9.6842e-6`. After the primary geometric update, M_geom is `-4.8340e-6`, source is `-9.6842e-6`, residual `+4.8502e-6`; returned state residual `+4.8491e-6`.
- N=80: predictor source `-6.7068e-6` with M_geom still roundoff, temporary predictor residual `+6.7068e-6`. After primary update M_geom is `-2.2424e-6`, residual `+4.4644e-6`; returned residual `+4.4615e-6`.

These are intermediate split stages, not independent physical slices; do not interpret the predictor-only mismatch as a physical violation. They do show that the stage ordering deliberately advances matter and scalar predictors before geometry responds, with only partial compensation at the primary stage.

At the saved t≈10 state:
- N=60, dt=.04: input M residual `-8.28647e-4`; Euler predictor `-8.33533e-4`; primary stage `-8.41586e-4`; final `-8.41161e-4`. Predictor change is dominated by matter-current source shift (+4.872e-6 in `8πJ`), while M_geom changes only `-1.38e-8`. Across the full step, M_geom worsens by `-7.591e-6` and source moves +`4.924e-6`, giving `ΔM=-1.2514e-5`.
- N=80, dt=.03: input M residual `-6.78261e-4`; Euler predictor `-6.81649e-4`; primary stage `-5.90502e-4`; final `-5.90190e-4`. The predictor change is again almost entirely the matter-current source shift (+3.376e-6 in `8πJ`) with geometry changing only `-1.09e-8`. But the primary update then makes M_geom less negative by +`9.1473e-5`, mostly through the `D_r Aa` term (+`5.5153e-5`) and the radial b term (+`3.9732e-5`), with `D_r K` contributing -`3.3153e-6`; the source changes only +`3.402e-6` over the full step. Net M improves by +`8.8071e-5`.

The primary update makes the momentum residual slightly worse at N=60 but substantially better at N=80. This resolution-dependent sign is important: there is no single universal “primary stage worsens M” story. The geometric terms show strong cancellation between `D_r Aa`, `D_r K`, and the radial b term. The result is consistent with a centre-localized spatial/derivative compatibility problem, but it is not enough to attribute the problem to any one term or matter sector.

## Semi-discrete momentum-residual tangent — complete

Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38213760419
Artifact SHA-256: `d00310158e8d509f34d4ade1c347e2b5fa4542565f51b749e9c781e0b94d862f`

The diagnostic evaluates the actual production semi-discrete RHS on six frozen states: N=60/N=80 at t=0 and saved t≈10/t=12. It computes the directional derivative of every signed geometric momentum term and each source-current term by perturbing the state along the RHS, re-solving the CMC lapse at the perturbed plus/minus states, and comparing two epsilons (2e-5 and 1e-5). It also takes one independent normal production step from each state for a measured ΔM/dt comparator. No long trajectory was launched.

The epsilon check is extremely stable: at i=2, the max absolute difference between the two ε estimates for ẊM is ~1.7e-13 to 6.0e-13 across all six states. The signed momentum term reconstruction and current-sector reconstruction close at 1e-19–1e-21 scale.

| State | Initial M residual at i=2 | RHS-predicted ẊM | Actual one-step ΔM/dt |
|---|---:|---:|---:|
| N60, t=0 | ~0 | +1.21132e-4 | +1.21228e-4 |
| N80, t=0 | ~0 | +1.48735e-4 | +1.48716e-4 |
| N60, t=10 | -8.28647e-4 | -3.21262e-4 | -3.12853e-4 |
| N80, t≈9.99 | -6.78261e-4 | +2.95895e-3 | +2.93569e-3 |
| N60, t=12 | -9.81671e-4 | -8.03541e-6 | -6.34156e-6 |
| N80, t=12 | -2.68944e-3 | -1.47313e-3 | -1.37845e-3 |

At t=0 the semi-discrete RHS itself moves a zero momentum residual away from zero; one step matches that rate to much better than 0.1%. The t≈10 tangent also matches the actual step rate within about 2.6% for N60 and 0.8% for N80. At t=12 the N60 tangent is near cancellation, so relative error is uninformative (absolute mismatch ~1.7e-6); N80 agrees in sign and scale within about 6.4%.

This is the strongest current distinction: **the residual tendency is already present in the semi-discrete coupled RHS, not created solely by choosing a finite time step.** The actual step responds consistently to that tangent. The two-CFL check likewise shows little early change on halving dt.

### Which terms drive the instantaneous residual rate?

At the initially constraint-satisfied i=2 cell:
- N60: geometry-side ẊM contribution `-1.20974e-4`, current-source `d(8πJ)/dt=-2.42106e-4`, net ẊM `+1.21132e-4`.
- N80: geometry-side `-7.48267e-5`, current-source `-2.23561e-4`, net ẊM `+1.48735e-4`.

The local S-field current supplies essentially all of that initial current-source rate: `d(8πJ_S)/dt=-2.42151e-4` at N60 and `-2.23599e-4` at N80. The D and COSMOS-φ current rates, dark matter, baryons, and radiation are orders of magnitude smaller at those cells.

At N60/t=12, the geometry-side rate is `-1.63963e-4` and current-source rate `-1.55928e-4`, nearly canceling; at N80/t=12 the geometry-side rate `-1.61209e-3` dominates the current-source rate `-1.38957e-4`. At N80/t≈10, the geometric rate `+3.06906e-3` is dominated by the `D_r Aa` term `+1.83503e-3` and radial-b term `+1.34463e-3`; the source-current rate is only `+1.10102e-4`.

**Interpretation limit:** a large S-current derivative is not proof the S equation is wrong; it must cancel appropriately against the geometry evolution if the discrete system propagates the momentum constraint. Here the balance differs by a measurable amount, and the geometric contribution changes sign between states/resolutions. No single RHS component is yet identified as the sole defect.

## Saved-state momentum budget across all 21 slices

Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38213685916
Artifact SHA-256: `0db11b923cd5ab5671d5da95739744c5b6fc6d1561de990a2b6e3c95c26324ee`

The exact geometric momentum decomposition, scalar/fluid current assembly, and total residual were rechecked at t=0,2,4,6,8,10,12 for N=40,60,80. Geometry and total-source closures all pass at roundoff. At t=12/i=2:
- N60: `M_geom=-9.71025e-4`, residual `M=-9.81671e-4`, scalar source term `-8πj_S=-1.082e-5`; DM source term about `+1.55e-7`, with baryon/radiation smaller.
- N80: `M_geom=-2.68327e-3`, residual `M=-2.68944e-3`, scalar source term `-8πj_S=-6.317e-6`; DM source term about `+1.31e-7`, with baryon/radiation smaller.

So the *value* of the late momentum residual is overwhelmingly geometric, even though the time derivative of the S current is important in its instantaneous propagation balance. This is why “the scalar source causes it” would be too strong: the residual depends on the difference between geometric and total-current evolution.

## Immediate next diagnostic: spatial behavior of the semi-discrete tendency

Evaluate the same no-trajectory momentum-residual tangent at N=40,60,80,120 on (a) fresh t=0 states and (b) saved t≈2.16 states from the successful same-domain refinement. Compare ẊM, each signed geometric rate term, and each scalar/fluid-current rate at common physical radii (r=2,3,4,5,6) as well as in fixed inner/exterior regions. Keep CMC, matter, boundary, and coefficients unchanged.

This will distinguish a central spatial layer in the instantaneous defect from a resolution-independent source/geometry inconsistency. No production equation change is authorized until the fixed-radius tangent itself has been checked for spatial refinement.

## Matched-CFL spatial refinement — complete

Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38212879394
Artifact SHA-256: `181e807240a54126db39a35fdd22f865509b96027db0a724853512225dd8b6e8`

Configuration: N=40,60,80,120; r_max=80; CFL=.03; radiation enabled; t=2.16; checkpoint interval=.24. All four runs reached exactly t=2.16, with 36,54,72,108 accepted steps respectively; all state and diagnostic values remained finite. The common time targets are exact multiples of each nominal dt (.06,.04,.03,.02).

At t=2.16, max `|H|` on the fixed exterior region r>=5 decreased strongly:
- N=40: 4.0086e-4
- N=60: 1.6886e-5
- N=80: 8.9180e-6
- N=120: 5.4585e-6

For the fixed physical annulus 2<=r<5, the volume-weighted RMS Hamiltonian residual fell from 7.4135e-4 (N40) to 6.5404e-4 (N60), 4.4980e-4 (N80), and 1.2832e-4 (N120). The volume weight uses `sqrt(a)*b/X^3 * r^2 dr` from the saved production metric; this is a norm diagnostic, not an acceptance threshold.

The maximum in the first **admitted** inner cells is not uniformly convergent in the sup norm: H peaks are about 7.51e-4 at r=3.333 (N60), 1.071e-3 at r=2.5 (N80), and 1.539e-3 at r=1.667 (N120). The peak moves inward with grid spacing, whereas the fixed-radius samples and volume-weighted annular norm improve. The full inner peaks, including the first two excluded cells, decline from about 9.27e-3 (N60) to 6.97e-3 (N80) and 4.74e-3 (N120).

This supports a **central error layer that narrows with refinement**, but does not yet establish pointwise regular convergence at the centre. N=40 has no cell with index>=2 and r<5, so it cannot supply the admitted-inner max for that mask.

The exterior convergence is encouraging; central regularity and the momentum-constraint source remain open. Do not add the center projection as a repair, change the Lambda equation, add damping, or launch t=160 based on these results alone. The next step is a one-step signed decomposition of the momentum constraint into its exact geometric derivatives and matter-current source across the existing initial and t≈10 stage captures.
