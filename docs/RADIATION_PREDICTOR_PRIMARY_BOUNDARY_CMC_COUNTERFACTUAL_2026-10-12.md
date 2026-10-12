# Primary-Predictor Boundary \(A_a\) Counterfactual — CMC/RHS Response

**Date:** 12 October 2026 (UTC)  
**Status:** Diagnostic-only. No production law change or scientific admission is authorized by this result.  
**Production main:** \`2035acc1210db18da0cededea18750b7ed813260\` (unchanged).  
**Diagnostic branch:** \`diag/radiation-primary-boundary-cf-2026-10-12\`  
**Matched counterfactual run:** [38321444722](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321444722) — successful capture.  
**Artifact:** [11778259794 — radiation-primary-boundary-counterfactual](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321444722/artifacts/11778259794)  
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

CI passed the full repository test suite (**94 passed in 7.61 s**), Python compilation, exact source/vendor guard and this two-step comparison. The CI success is a successful diagnostic capture, not a successful \(t=50\) physical trajectory.

## Primary-boundary correction and CMC response

At the target primary-predictor boundary call, outer cell \(i=79\), \(r=79.5\):

| Quantity | Value |
|---|---:|
| Outer \(A_a\) before boundary formula | \(0.15390845494664215\) |
| \(A_a\) returned by boundary formula | \(0.15395152848912053\) |
| Boundary correction \(\Delta A_a\) | \(+4.307354247837236\times10^{-5}\) |
| \(A_a\) used by actual arm | \(0.15395152848912053\) |
| \(A_a\) used by counterfactual at this call | \(0.15390845494664215\) |

The CMC solver responds measurably to this outer \(A_a\) difference:

| CMC witness | Actual | No-primary-correction arm | Actual minus control |
|---|---:|---:|---:|
| Minimum lapse | 0.29441369046699944 | 0.29441269011190896 | \(+1.00035509048\times10^{-6}\) |
| Maximum lapse | 1.7375148191265442 | 1.7375100983023495 | \(+4.7208241947\times10^{-6}\) |
| Outer lapse (normalization) | 1.0 | 1.0 | 0 |
| CMC target \(K\)-rate | 0.22737526310717676 | 0.22737449053440156 | \(+7.72572775198\times10^{-7}\) |
| Maximum full lapse-profile difference | — | — | \(4.86488599116\times10^{-6}\) |

So the boundary operation is not mathematically invisible to the gauge response. The key question, however, is the scale of its effect on the imminent cone crossing.

## Stage-1 RHS and geometry response

The matched stage-1 recomputation produced nonzero but small changes across the coupled RHS and geometry-stage blocks. Selected outer-cell differences are:

| Stage-1 witness at \(i=79\) | Actual | No-primary-correction arm | Absolute difference |
|---|---:|---:|---:|
| Explicit \(\dot a\) | -0.1892865863947715 | -0.1892336265839812 | \(5.29598107903\times10^{-5}\) |
| Explicit \(\dot b\) | 0.19634418216714966 | 0.19628924773711212 | \(5.49344300375\times10^{-5}\) |
| Explicit \(\dot X\) | 0.08650935394552563 | 0.08650935394552563 | 0 |
| Primary \(K\) RHS block | 0.06114681293728212 | 0.061126922007171264 | \(1.98909301109\times10^{-5}\) |
| Primary \(A_a\) RHS block | 0.04231025115132725 | 0.042298374137658946 | \(1.18770136683\times10^{-5}\) |
| \(\Lambda\) L2 block | 0.2216145292574766 | 0.22156021850111876 | \(5.43107563578\times10^{-5}\) |
| Radiation energy RHS | 0.00018758892957665933 | 0.00018757584090202343 | \(1.30886746359\times10^{-8}\) |
| Radiation momentum RHS | 0.00010941971803392654 | 0.00010941899078295786 | \(7.27250968675\times10^{-10}\) |

The difference affects the CMC lapse profile, the stage-1 RHS and geometry, and some evolved conservative variables. The maximum RHS changes should not be confused with the final predictor margin: they are rates/blocks, and the fully coupled updates and subsequent boundary operation still run normally in both arms.

## Next explicit predictor: decisive scale comparison

Both arms reproduce the same failure cell, stage, and step. The actual predictor agrees with the previously archived rate-corrected witness to the frozen criterion \(5\times10^{-11}\) in \(C\).

| Witness at the next predictor, \(t=47.01\), \(i=79\) | Actual | No-primary-correction arm |
|---|---:|---:|
| \(a\) | 0.6032317141253851 | 0.603232395594617 |
| \(X\) | 0.9465358619432119 | 0.9465358551490698 |
| \(\gamma^{rr}=X^2/a\) | 1.4852172340500571 | 1.485215534884605 |
| Densitized cone margin \(C=U_E-\sqrt{\gamma^{rr}}|U_r|\) | \(-1.699741681811858\times10^{-8}\) | \(-1.700374445276964\times10^{-8}\) |
| \(|S|/E\) | 1.0000573283549097 | 1.0000573497293141 |

The actual-minus-control margin difference is only
\[
\Delta C = +6.3276346510578174\times10^{-12}.
\]
This is about **0.037% of the magnitude of the approximately \(1.70\times10^{-8}\) cone overshoot**. Both values of \(C\) remain negative. Removing the primary-boundary correction does not admit the predictor; at this bounded horizon, the actual branch is microscopically *less* outside the cone than the no-primary-correction branch.

That sign is a reminder not to infer causal direction from the sign of a boundary correction alone: CMC, the staged geometry RHS, matter RHS, and the next boundary all interact. The important quantitative result is the tiny final margin response despite measurable stage-level differences.

## Joint accepted-state witnesses at \(t=46.98\)

The one-step accepted endpoint comparison shows the same limited effect:

| Witness | Actual | No-primary-correction arm |
|---|---:|---:|
| Full-domain Hamiltonian maximum | 0.012112228676626889 | 0.012112225811460924 |
| Full-domain momentum maximum | 0.0015657672595360958 | 0.0015690094248118204 |
| Full-domain connection maximum | 0.002861315858133817 | 0.0028613149177650343 |
| Determinant-constraint maximum | \(2.30854014656\times10^{-5}\) | \(2.30853928258\times10^{-5}\) |
| Last-interior CMC residual | \(-4.06667742098\times10^{-4}\) | \(-4.06673809033\times10^{-4}\) |
| Misner–Sharp mass at diagnostic surface | 5.0634156458115545 | 5.063414318727431 |
| \(H_{\rm eff}\) | -0.09639885154075516 | -0.09639883481332635 |
| Accepted outer radiation ratio | 0.9983755569474224 | 0.9983755647996045 |

These are finite, non-identical responses. They do not close the joint constraints or conservation gates: the Hamiltonian, momentum, connection, determinant and CMC residuals remain nonzero, and the next explicit predictor violates radiation admissibility in both arms.

## Interpretation and disposition

**Established:** the primary-predictor boundary correction changes the stage data and CMC solution, but suppressing it for one call changes the next predictor cone margin by only \(6.33\times10^{-12}\) and does not remove the radiation failure.

**Not established:** this bounded two-step counterfactual is not a long alternate trajectory, and does not prove the boundary formula is universally correct or incorrect. It says this specific primary-predictor correction is not a sufficient explanation or repair for the captured late-time cone crossing. It does not rule out other cumulative effects of the boundary, interior source coupling, centre regularity, or metric/matter stage synchronization.

**Decision:** do not alter or remove the boundary law; do not add a clamp, floor, damping, cone relaxation, or fitted correction. Keep the primary and completed-state boundary hypotheses separate from the independently open centre-series, H/M/\(C_\Lambda\), source reciprocity, and Misner–Sharp tests. The most productive next move is to reconcile this report with the parallel stacked-source boundary-history artifact before asking for another test. No production code, \`main\`, pinned vendor, or canonical GEAR archive was changed.
