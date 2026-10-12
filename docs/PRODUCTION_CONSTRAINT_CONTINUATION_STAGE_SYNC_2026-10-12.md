# Production constraint continuation — stage-sync and regularized-derivative probes (2026-10-12)

## Scope and admission

This report records two bounded, diagnostic-only tests against the frozen production source commit:

- Production source: \`6d8b3f2790644591326df57c55908863f994cb0c\`
- Pinned vendor: \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`
- Base campaign through \(t=12\): [run 38208253687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687), resolutions N=40/60/80, \(r_{\max}=80\), CFL=.03, radiation enabled.
- No production equation, coefficient, gauge, boundary, or campaign setting was changed. These probes live on \`diag/lambda-metric-connection-rhs-2026-10-12\`; \`main\` is untouched.

“Success” below means the diagnostic workflow completed, not that the full production equations are constraint-preserving.

## New test A — Lambda predictor stage synchronization

Run: [38277888835](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38277888835)  
Artifact: \`lambda-pirk-stage-sync-counterfactual\`  
Tool: \`tools/test_lambda_predictor_stage_sync.py\`

Four frozen one-step inputs were compared: N=60/80 at the initial slice and saved \(t\approx10\). The counterfactual recomputed the predictor \(\lambda_{L2}\) term after advancing the predictor's \(A_a,K\) primary block. No trajectory was advanced beyond one step per path.

For N=80 saved \(t\approx10\), the predictor \(\lambda_{L2}\) changed locally by at most \(3.95\times10^{-5}\), but returned-state constraint changes were negligible:
- admitted-stencil H max-norm change: approximately \(7\times10^{-12}\);
- admitted-stencil M-direct max-norm change: approximately \(7\times10^{-10}\);
- admitted-stencil connection C max-norm change: approximately \(5\times10^{-12}\).

N=60 gives the same qualitative result. The production H residual near \(r=2.5\) remains approximately \(-1.288\times10^{-2}\); the connection residual remains approximately \(-1.226\times10^{-2}\). The stage-sync counterfactual does not move the state appreciably toward constraint closure.

**Decision:** demote stale predicted-\(A_a/K\) evaluation inside the predictor \(\lambda_{L2}\) call as the leading cause. Do not edit production stage ordering from this test.

## New test B — Tensor-regularized \(D_r A_a\) counterfactual

Run: [38279405856](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38279405856)  
Artifact: \`regularized-aa-derivative-counterfactual\`  
Tool: \`tools/test_regularized_aa_derivative_counterfactual.py\`

The pinned vendor includes \`regularized_aa_prime\`, which differentiates a regular trace-free component using \(A_a=(2/3)r^2 A_\lambda\). The diagnostic substituted this operator only inside the momentum-constraint expression and Lambda L2 RHS, then compared one normal production step on N=60/80 initial and saved \(t\approx10\) states. It did not project any state field.

### Result 1 — no evolution change

All compared evolution state fields and H/C diagnostics were bit-identical between baseline and alternative paths for the tested single steps. H, \(C_\Lambda\), and M evaluated with the original direct-derivative definition had exactly zero reported delta.

This is explained by the pinned Lambda L2 equation at \(\lambda_m=2\). Its explicit term \(-2\alpha D_rA_a/a\) cancels the \(D_rA_a\) part of \(+2\alpha M_{\rm geom}/a\) after the momentum constraint is expanded. Replacing that derivative coherently in both places therefore cannot change the Lambda RHS. This operator substitution cannot repair the evolved constraint drift.

### Result 2 — the reported M constraint is derivative-sensitive, but not uniformly improved

At the saved \(t\approx10\), cell \(i=2\):
- N=60, \(r=3.333\): direct M changes from \(-8.2865\times10^{-4}\) on input to \(-8.4116\times10^{-4}\) after one baseline step; using the regularized derivative for the *diagnostic definition* gives \(-5.7323\times10^{-5}\) on input and \(-8.4007\times10^{-5}\) on output.
- N=80, \(r=2.5\): direct M is \(-6.7826\times10^{-4}\) on input and \(-5.9019\times10^{-4}\) after the step; the regularized-definition values are \(-5.9378\times10^{-4}\) and \(-5.1374\times10^{-4}\).

But the result is not a global improvement: at N=60, neighboring first cells move the other way, and for an initially constraint-satisfied state the regularized derivative makes the post-step \(i=2\) M witness *larger* (about \(6.7\times10^{-5}\) versus \(4.85\times10^{-6}\) direct). H and C are unchanged.

**Decision:** do not switch the production derivative or quote its smaller N=60 \(i=2\) residual as a fix. The alternative changes the discrete diagnostic meaning, does not change the state, and is not uniformly better across cells/resolutions.

## Reconciled production findings

The previous stage / fixed-radius audits remain the current authority:

- [Signed one-step production stage trace](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38213136600): the constraint budget reconstructs to roundoff and localizes the growing residual to geometry/matter RHS balance, not a failure in the diagnostic algebra.
- [Saved-state Lambda RHS tangent](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210572528): the semi-discrete \(C_\Lambda\) tangent is near roundoff at \(t=0\), then grows; at the saved late-time admitted cell, the momentum-constraint addition accounts for nearly all of \(\dot C_\Lambda\).
- [Semi-discrete M tangent](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38213760419) and [per-step M observer](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38214554266): the M residual tendency is already present in the semi-discrete coupled RHS; one-step rates agree closely with it, and the residual is strongly time-dependent/oscillatory.
- [Fixed-radius geometric M decomposition](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38215858158): at \(t=2.16\), resolved fixed-radius H norms improve toward N=120, while the innermost cells and the geometric M terms remain non-monotonic. The dominant geometric combination shifts between \(D_rA_a\), \(-(2/3)D_rK\), and \(3A_a/r\); no single term has a consistent signed error across radii/resolutions.
- [Centre-series RHS rate audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38215595339): the earliest evolved central feature is a large and resolution-sensitive \(D_rK/r\) ratio in the first cell. At \(t=2.16\), for example, the N=80 values across the first three cells are approximately \(-4.22\times10^{-3}, -1.196\times10^{-3}, +1.264\times10^{-4}\), and N=120 gives \(-9.558\times10^{-3}, -2.717\times10^{-3}, +2.776\times10^{-4}\). This localizes a central K-gradient layer, but does not establish whether lapse Hessian, Ricci, or A/K evolution causes it.
- [Centre parity/operator audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38214547118): manufactured parity-polynomial tests pass at about \(2\times10^{-15}\); no generic sign/parity bug was found in the shared derivative.

The aggregate evidence still favors a central spatial/derivative compatibility issue. The full trajectory is not yet admitted as constraint-preserving and should not be interpreted as a physical result.

## Literature check

The relevant numerical-relativity literature supports the *methodological caution*, not a conclusion about this particular code. Montero & Cordero-Carrión demonstrate stable spherical BSSN evolutions using PIRK without an origin-regularization algorithm ([Phys. Rev. D 85, 124037 (2012)](https://doi.org/10.1103/PhysRevD.85.124037)). Baumgarte et al. discuss covariant/reference-metric BSSN and regularity handling in spherical polar coordinates ([Phys. Rev. D 87, 044026 (2013)](https://doi.org/10.1103/PhysRevD.87.044026)). Neither result validates the present CMC-gauged, matter-coupled split implementation; these are comparison targets for operator ordering and regularity formulation, not a basis for inserting an untested projection.

## Next step

Do not run a longer campaign and do not alter the Lambda equation. The next bounded investigation should split the initial/evolved \(K\) RHS into the *separately signed* lapse-Laplacian and nonlinear/extrinsic-curvature contributions while recording the CMC lapse's own elliptic residual and the same-cell \(D_rK/r\) series. Compare initial and \(t\approx2.16\) saved states at N=40/60/80/120 and fixed physical radii. The purpose is to identify whether the central K-gradient layer originates in the CMC lapse solve, in the \(K\) RHS decomposition, or in their discrete interface. Do not project or edit fields; only a source-localized residual that survives this operator-level test can justify a production change.


## New CMC/K-RHS frozen-state finding — center row is inconsistent (2026-10-12)

Workflow: [38279973263](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38279973263), artifact \`cmc-K-rhs-discrete-compatibility\`.  
Frozen-state center-row comparison: [38280503385](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38280503385), artifact \`cmc-lapse-center-row-counterfactual\`.  
Counterfactual code is on this branch: \`tools/test_cmc_lapse_center_row_counterfactual.py\`. No production state was advanced in either audit.

### The specific mismatch

The production grid is cell-centered, \(r_i=(i+1/2)\Delta r\). However, \`engine/cmc_gauge.py\` replaces the first CMC equation with \(\alpha_0=\alpha_1\). That equality puts a zero slope across the first *two cell centers*; it is not the regular even-parity ghost rule at the origin.

On the same frozen state, the present CMC solve’s second-order rows close their own discrete equation to roundoff for \(i=1,\ldots,N-2\). But the actual production \(K\) RHS uses the pinned fourth-order BSSN lapse Laplacian. Comparing the production \(K\) RHS to the target \( \dot K_{\rm CMC}\) exposes a much larger first-cell inconsistency:

- Across the initial, \(t=2.16\), and \(t\approx10\) samples, the current first-cell production-\(K\)-RHS residual is approximately \(1.2\)–\(1.6\times10^{-3}\) in magnitude at N=60/80/120.
- At N=80, it is \(-1.5501\times10^{-3}\) initially, \(+1.3262\times10^{-3}\) at \(t=2.16\), and \(+1.4985\times10^{-3}\) at \(t\approx10\).
- At \(i=1\), the fourth-vs-second-order operator gap is much smaller but still around \(4\)–\(6\times10^{-5}\); at \(i=2\) it is around \(0.8\)–\(2.0\times10^{-5}\) over these samples. The operator mismatch thus does not end at the special center row.

### Center-row counterfactual

The frozen-state counterfactual held the CMC target \(\dot K\) fixed and changed only the first linear-system row to a cell-centered even-parity ghost, \(\alpha_{-1}=\alpha_0\), while solving the lapse PDE at \(i=0\). It compared that against (a) the existing second-order solve and (b) a full production-matched fourth-order solve.

- For N=80 initially, the first-cell production-\(K\)-RHS residual fell from \(-1.5501\times10^{-3}\) to \(+2.4547\times10^{-5}\) under the even-parity second-order row.
- At N=80, \(t=2.16\), it fell from \(+1.3262\times10^{-3}\) to \(-1.0868\times10^{-5}\).
- At N=80, \(t\approx10\), it fell from \(+1.4985\times10^{-3}\) to \(-2.7669\times10^{-6}\).
- The lapse profile changes by up to roughly \(5\times10^{-4}\)–\(7\times10^{-4}\) for N=80 and about \(1.2\times10^{-3}\) at N=60/t≈10. Yet the admitted \(i\ge2,\ r<5\) K-RHS mismatch remains approximately \(10^{-5}\); the center-row correction alone does not eliminate the interior fourth-/second-order gap.
- The fully fourth-order lapse solve makes the production K-RHS residual roundoff-sized by construction. Treat this as an operator-consistency control, **not** independent evidence that its gauge evolution is the correct repair.

### Decision

This is the first concrete defect candidate in the active CMC/K-RHS interface: the solver’s first row is not the even-parity boundary closure for the cell-centered grid. It can plausibly seed the narrow center layer, but it is not yet the established cause of the later Hamiltonian/momentum constraint drift. Two cautions remain:

1. Changing the first row does not remove the bulk mismatch between the second-order elliptic operator and the fourth-order evolution RHS.
2. The first two cells are excluded from the main admitted-inner constraint gate; we need evidence that the revised lapse changes the evolved \(H\), \(M\), and \(C_\Lambda\) at admitted cells/fixed radii, not just \(K\)'s first-cell algebraic residual.

The bounded next test is a single normal production step on identical inputs, with only the CMC center row changed and all output constraints evaluated using the same original production diagnostic. Workflow: [CMC center-even lapse one-step counterfactual](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38280666579). Until that output is reviewed, do not patch \`engine/cmc_gauge.py\`, do not replace the production derivative operator, and do not launch a long campaign.

## Research comparator

Montero & Cordero-Carrión’s spherical BSSN/PIRK implementation uses cell-centered radial grids and careful origin treatment while demonstrating regular evolutions without a separate origin-regularization algorithm ([Phys. Rev. D 85, 124037 (2012)](https://doi.org/10.1103/PhysRevD.85.124037)). The paper supports treating parity and center-row closure as structural parts of the discretization; it does not directly validate the CMC row or matter-coupled split used here. The current counterfactual is therefore a useful interface test, not an imported repair.


### One-step result for the even-parity center-row alternative

Run: [38280666579](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38280666579), artifact \`cmc-center-even-lapse-one-step\`.  
Tool: \`tools/test_cmc_center_even_lapse_one_step.py\`. One ordinary production step was compared from identical states, using the same original H/M/\(C_\Lambda\) diagnostics.

The alternative does **not** improve the complete constraint set. Across all six states (N=60/80 initial; N=80/120 at \(t=2.16\); N=60/80 at \(t\approx10\)), the admitted-inner momentum residual norm rises by approximately \(3.0\times10^{-6}\) to \(3.8\times10^{-6}\) in max norm after one step. The RMS increase is roughly \(1.4\times10^{-6}\) to \(2.7\times10^{-6}\). Connection changes are small but mixed; the Hamiltonian norm improves slightly, by roughly \(2\times10^{-7}\) to \(3\times10^{-7}\) RMS. At \(t=0\), the alternative moves an initially near-constraint-satisfied M residual farther from zero at admitted cells.

Representative admitted-inner changes (alternative minus baseline, max norm):
| Frozen input | \(\Delta\max |H|\) | \(\Delta\max |M|\) | \(\Delta\max |C_\Lambda|\) |
|---|---:|---:|---:|
| N80, \(t=0\) | \(-5.16\times10^{-7}\) | \(+3.79\times10^{-6}\) | \(+1.16\times10^{-7}\) |
| N80, \(t=2.16\) | \(-3.98\times10^{-7}\) | \(+3.03\times10^{-6}\) | \(-8.88\times10^{-8}\) |
| N120, \(t=2.16\) | \(-4.08\times10^{-7}\) | \(+3.12\times10^{-6}\) | \(+5.98\times10^{-8}\) |
| N80, \(t\approx10\) | \(-4.27\times10^{-7}\) | \(+3.53\times10^{-6}\) | \(+1.02\times10^{-7}\) |

The center-even row correctly reduces the algebraic first-cell lapse/K-RHS mismatch, but it is not sufficient as a production change and does not improve momentum propagation. Preserve it as a diagnostic finding only.

### Next discriminator: the full CMC operator/order compatibility

Workflow started: [CMC production-matched one-step counterfactual](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38280957254), tool \`tools/test_cmc_production_matched_one_step.py\`.

This next experiment uses the exact production fourth-order even-parity derivative matrices in the CMC elliptic solve, retaining the same stage-wise target functional and outer alpha normalization. It compares one ordinary production step against the current solve on the same frozen inputs. Its elliptic residual is expected to close by construction, so the decision is based on whether the *evolved, original-definition* H, M, and connection witnesses improve—not on the lapse-equation residual alone. No long trajectory is planned unless this short comparison demonstrates a clear multi-constraint benefit.

The original local GEAR archive also documents a central finite-volume divergence mistake in a different radial scalar operator (GEAR-70C.1). That is a useful reminder to test radial origin terms cell by cell, but it must not be transplanted as a patch to this BSSN/CMC machine: its operator and variable representation differ. The current one-step evidence specifically rejects a standalone center-row replacement.



### Production-matched CMC one-step result — no constraint improvement

Run: [38280957254](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38280957254), artifact \`cmc-production-matched-one-step\`.  
Tool: \`tools/test_cmc_production_matched_one_step.py\`.

The CMC solver was rebuilt diagnostically with the exact fourth-order, even-parity derivative matrices used by the production BSSN K RHS, then substituted only for the lapse solves on the same six frozen inputs. This version solves its own K-RHS operator identity by construction; the relevant comparison is the evolved state measured by the original production H/M/\(C_\Lambda\) constraints.

The result closely repeats the even-parity second-order counterfactual:
- At N=80, \(t=0\), admitted-inner max-norm changes after one step are \(\Delta\max|H|=-5.16\times10^{-7}\), \(\Delta\max|M|=+3.79\times10^{-6}\), and \(\Delta\max|C_\Lambda|=+1.15\times10^{-7}\).
- At N=80, \(t=2.16\), they are \(-3.99\times10^{-7}\), \(+3.06\times10^{-6}\), and \(-9.12\times10^{-8}\).
- At N=120, \(t=2.16\), they are \(-4.10\times10^{-7}\), \(+3.15\times10^{-6}\), and \(+6.18\times10^{-8}\).
- At N=80, \(t\approx10\), they are approximately \(-4.25\times10^{-7}\), \(+3.61\times10^{-6}\), and \(+1.07\times10^{-7}\).

The full fourth-order lapse operator and the center-row-only fix therefore both make the Hamiltonian norm slightly smaller but worsen the admitted-inner momentum residual at the same scale. They do not change the established large late-time H/M/C defects meaningfully.

**Decision:** downgrade the CMC center row and the CMC/K-RHS stencil mismatch as the primary cause of the evolved constraint failure. Keep this as a genuine local operator-consistency finding, but do not patch \`engine/cmc_gauge.py\`: the one-step constraint evidence does not support it as a repair.

## New equation-level consistency lead: scalar lapse-gradient contraction

The next diagnostic will audit the scalar momentum equation’s \(\nabla^i\alpha\,\nabla_i f\) term against the physical radial inverse metric used elsewhere in the same production code.

The production metric convention is \(\gamma_{rr}=a/X^2\), so \(\gamma^{rr}=X^2/a\). The scalar Laplacian and stress projection both use \(X^2/a\), but the extra lapse-gradient term in \`engine/scalar_system.py\` is currently coded as \(X^{-2}\,\partial_r\alpha\,\partial_r f\). The same coefficient appears in the retained \`engine/reference_pirk_unified.py\`, so provenance alone does not establish correctness. The covariant contraction under the current declared metric variables indicates \(X^2/a\); the next step is an isolated one-step comparison to quantify whether this term can account for the S-current/momentum-constraint tangent mismatch. No production change is authorized yet.
