# Primary-Predictor Boundary \(A_a\) Counterfactual — CMC/RHS Response

**Date:** 12 October 2026 (UTC)  
**Status:** Diagnostic-only. No production law change or scientific admission is authorized by this result.  
**Production main:** \`2035acc1210db18da0cededea18750b7ed813260\` (unchanged).  
**Diagnostic branch:** \`diag/radiation-primary-boundary-cf-2026-10-12\`  
**Matched counterfactual run:** [38321467866](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321467866) — clean rerun successful.  
**Artifact:** [11779103151 — radiation-primary-boundary-counterfactual](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321467866/artifacts/11779103151)  
**Artifact ZIP SHA-256:** \`3edf5a3c8500a904d51fd6b5af7deea459b1ab1a205571f3032034fe5d43fd09\`  
**Embedded report SHA-256:** \`7077442ccef26aac4a1f16a6f0b9b1bcc3079adb0116b1e7a1146c6b40c90164\`  
**Tested instrumentation head:** \`667cf399df480b901c72423ada01fa98129506bf\`  
**Rate-corrected source baseline:** \`1c63090d4f42396f899421a2379f70320689d5ad\`  
**Pinned vendor:** \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`

## Question and bounded method

The preceding completed-state boundary increment was already tested in a frozen-\(\alpha\) direct-channel calculation. It modestly worsened the next predictor margin but did not restore admissibility when removed. This experiment tested the *distinct earlier call*: the primary-predictor outer-light-boundary operation.

The source-matched replay advanced the N=80, \(R_{\max}=80\), radiation-on run with fixed \(dt=0.03\) to the accepted state at \(t=46.92\), \(\tau=21.998447614469036\). Two branches then followed the normal coupled kernel through \(t=46.98\):
- **Actual:** normal primary-predictor boundary formula and normal completed-state boundary.
- **Counterfactual:** execute the primary-predictor boundary formula for measurement, then restore only the pre-call outer \(A_a\) at cell 79 for that one primary-predictor call. All other data/stages remain on the normal path; in particular, the completed-state boundary is still enabled.

Both arms then form the next explicit geometry/matter predictor at \(t=47.01\) from their own accepted \(t=46.98\) states. The diagnostic stops before CMC primitive recovery rejects the predictor; it does not attempt retry, clipping, projection, or physical evolution past the failure. The process-local wrapper is restored after the measurement.

CI passed the full repository test suite (**94 passed in 10.54 s**), Python compilation, exact source/vendor guard and this two-step comparison. The CI success is a successful diagnostic capture, not a successful \(t=50\) physical trajectory.

## Primary-boundary correction and CMC response

At the target primary-predictor boundary call, outer cell \(i=79\), \(r=79.5\):

| Quantity | Value |
|---|---:|
| Outer \(A_a\) before boundary formula | \(0.1539084555761373\) |
| \(A_a\) returned by boundary formula | \(0.15395152911917892\) |
| Boundary correction \(\Delta A_a\) | \(+4.307354304161626\times10^{-5}\) |
| \(A_a\) used by actual arm | \(0.15395152911917892\) |
| \(A_a\) used by counterfactual at this call | \(0.1539084555761373\) |

The CMC solver responds measurably to this outer \(A_a\) difference:

| CMC witness | Actual | No-primary-correction arm | Actual minus control |
|---|---:|---:|---:|
| Minimum lapse | 0.29441369052911975 | 0.294412690174016 | \(+1.00035510376\times10^{-6}\) |
| Maximum lapse | 1.737514814723544 | 1.737510093899299 | \(+4.720824245\times10^{-6}\) |
| Outer lapse (normalization) | 1.0 | 1.0 | 0 |
| CMC target \(K\)-rate | 0.22737526399668426 | 0.2273744914238959 | \(+7.72572788355\times10^{-7}\) |
| Maximum full lapse-profile difference | — | — | \(4.86488604090\times10^{-6}\) |

So the boundary operation is not mathematically invisible to the gauge response. The key question, however, is the scale of its effect on the imminent cone crossing.

## Stage-1 RHS and geometry response

The matched stage-1 recomputation produced nonzero but small changes across the coupled RHS and geometry-stage blocks. Selected outer-cell differences are:

| Stage-1 witness at \(i=79\) | Actual | No-primary-correction arm | Absolute difference |
|---|---:|---:|---:|
| Explicit \(\dot a\) | -0.18928658690385228 | -0.18923362709244376 | \(5.29598114085\times10^{-5}\) |
| Explicit \(\dot b\) | 0.19634418310840085 | 0.19628924867760644 | \(5.49344307944\times10^{-5}\) |
| Explicit \(\dot X\) | 0.08650935394552563 | 0.08650935394552563 | 0 |
| Primary \(K\) RHS block | 0.06114681342678672 | 0.06112692249633438 | \(1.98909304523\times10^{-5}\) |
| Primary \(A_a\) RHS block | 0.042310251489555355 | 0.04229837447568554 | \(1.18770138698\times10^{-5}\) |
| \(\Lambda\) L2 block | 0.2216145302332118 | 0.221560219476208 | \(5.43107570038\times10^{-5}\) |
| Radiation energy RHS | 0.00018758892668780292 | 0.0001875758380127368 | \(1.30886739501\times10^{-8}\) |
| Radiation momentum RHS | 0.00010941971496171925 | 0.00010941898771073017 | \(7.27250989085\times10^{-10}\) |

The difference affects the CMC lapse profile, the stage-1 RHS and geometry, and some evolved conservative variables. The maximum RHS changes should not be confused with the final predictor margin: they are rates/blocks, and the fully coupled updates and subsequent boundary operation still run normally in both arms.

## Next explicit predictor: decisive scale comparison

Both arms reproduce the same failure cell, stage, and step. The actual predictor agrees with the previously archived rate-corrected witness to the frozen criterion \(5\times10^{-11}\) in \(C\).

| Witness at the next predictor, \(t=47.01\), \(i=79\) | Actual | No-primary-correction arm |
|---|---:|---:|
| \(a\) | 0.6032317132323858 | 0.6032323947016255 |
| \(X\) | 0.9465358623871615 | 0.9465358555930193 |
| \(\gamma^{rr}=X^2/a\) | 1.485217237641921 | 1.4852155384764436 |
| Densitized cone margin \(C=U_E-\sqrt{\gamma^{rr}}|U_r|\) | \(-1.6997425987812665\times10^{-8}\) | \(-1.7003753622897404\times10^{-8}\) |
| \(|S|/E\) | 1.0000573283857723 | 1.0000573497601783 |

The actual-minus-control margin difference is only
\[
\Delta C = +6.3276350847386864\times10^{-12}.
\]
This is about **0.037% of the magnitude of the approximately \(1.70\times10^{-8}\) cone overshoot**. Both values of \(C\) remain negative. Removing the primary-boundary correction does not admit the predictor; at this bounded horizon, the actual branch is microscopically *less* outside the cone than the no-primary-correction branch.

That sign is a reminder not to infer causal direction from the sign of a boundary correction alone: CMC, the staged geometry RHS, matter RHS, and the next boundary all interact. The important quantitative result is the tiny final margin response despite measurable stage-level differences.

## Joint accepted-state witnesses at \(t=46.98\)

The one-step accepted endpoint comparison shows the same limited effect:

| Witness | Actual | No-primary-correction arm |
|---|---:|---:|
| Full-domain Hamiltonian maximum | 0.012112228777012048 | 0.012112225911846165 |
| Full-domain momentum maximum | 0.001565767277511094 | 0.0015690094428324403 |
| Full-domain connection maximum | 0.002861315886176663 | 0.0028613149458070475 |
| Determinant-constraint maximum | \(2.30854016205\times10^{-5}\) | \(2.30853929806\times10^{-5}\) |
| Last-interior CMC residual | \(-4.06667745233\times10^{-4}\) | \(-4.06673812167\times10^{-4}\) |
| Misner–Sharp mass at diagnostic surface | 5.063415687415684 | 5.063414360331468 |
| \(H_{\rm eff}\) | -0.0963988519355134 | -0.0963988352080843 |
| Accepted outer radiation ratio | 0.998375557016913 | 0.998375564869101 |

These are finite, non-identical responses. They do not close the joint constraints or conservation gates: the Hamiltonian, momentum, connection, determinant and CMC residuals remain nonzero, and the next explicit predictor violates radiation admissibility in both arms.

## Interpretation and disposition

**Established:** the primary-predictor boundary correction changes the stage data and CMC solution, but suppressing it for one call changes the next predictor cone margin by only \(6.33\times10^{-12}\) and does not remove the radiation failure.

**Not established:** this bounded two-step counterfactual is not a long alternate trajectory, and does not prove the boundary formula is universally correct or incorrect. It says this specific primary-predictor correction is not a sufficient explanation or repair for the captured late-time cone crossing. It does not rule out other cumulative effects of the boundary, interior source coupling, centre regularity, or metric/matter stage synchronization.

**Decision:** do not alter or remove the boundary law; do not add a clamp, floor, damping, cone relaxation, or fitted correction. Keep the primary and completed-state boundary hypotheses separate from the independently open centre-series, H/M/\(C_\Lambda\), source reciprocity, and Misner–Sharp tests. The most productive next move is to reconcile this report with the parallel stacked-source boundary-history artifact before asking for another test. No production code, \`main\`, pinned vendor, or canonical GEAR archive was changed.
