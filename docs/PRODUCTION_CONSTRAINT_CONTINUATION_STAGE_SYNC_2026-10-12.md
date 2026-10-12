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
