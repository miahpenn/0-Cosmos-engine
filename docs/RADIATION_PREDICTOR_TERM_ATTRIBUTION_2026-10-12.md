# Radiation Predictor Metric: Frozen Term Attribution and Control Synthesis

**Date:** 2026-10-12 (UTC)  
**Status:** Diagnostic-only. No production equation change or repair is justified by this report.  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Current production main:** `2035acc1210db18da0cededea18750b7ed813260` (verified unchanged during this review).  
**Method:** Offline inspection and calculations from completed GitHub Actions artifacts; no new trajectory was launched for this report.

## Executive conclusion

The prior high-resolution first-cell failure and the late-time outer-cell failure are different events.

- The N=240/320 first-cell failure was a false alarm in the diagnostic wrapper: the old RHS query attempted a hypothetical unit-sized update and then applied the finite-update radiation cone guard. The rate-only RHS correction removes that mismatch while keeping the guard on actual finite-time updates.
- With that rate-only correction active, a full bounded N=80 trace still captures the late-time radiation-cone violation at the outermost cell during the first predictor CMC lapse solve. A separate stage ledger reproduces the same transition and attributes it to the explicit Euler geometry predictor, before algebraic regularity projection.
- At the previous accepted step, the raw metric predictor narrows the margin substantially but remains admissible. One step later, the same decomposition crosses the cone. The projection marginally improves the margin; it does not initiate the crossing.
- In the terminal metric change, both the a update and the X update raise the inverse radial metric \(\gamma^{rr}=X^2/a\), which narrows the radiation cone. The b update has no direct term in \(\gamma^{rr}\), although b participates in the determinant/regularity projection that also adjusts a.
- Previously tried accepted-metric recovery, adaptive retries, and domain-extension controls do not establish a repair. The failure still occurs at an outermost cell near central proper time \(\tau\sim22\), with the larger-domain run failing at a later coordinate time. That correlation is a diagnostic clue, not proof of a continuum or boundary-independent physical effect.

The evidence points to the **explicit metric-predictor / matter-stage synchronization** as the next engineering layer to investigate. It does not prove the underlying physical geometry equations are wrong and does not justify changing the radiation cone inequality, clipping the state, or adding a floor, damping term, or fitted correction.

## 1. Rate-corrected trace: failure is real and reproduced

### Latest trace

- Workflow: [38310237842](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842)
- Result: successful diagnostic capture, **not** a successful physical horizon.
- Artifact: [ZIP artifact 11775230036](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842/artifacts/11775230036)
- Verified artifact ZIP SHA-256: `f64298e2819e0dd9c4c2df712a4d1250e4855a5b0e467430fe7aa419a533e49f`
- Embedded trace SHA-256: `53a1f9f188e0eea4575da4458d09891f50144b9f89f945c995cd9980221fb619`
- Frozen source used before trace instrumentation: `1c63090d4f42396f899421a2379f70320689d5ad`
- Instrumentation commit: `51335a56873d65b8ed28863b7658ce02793ae396`
- Pinned numerical vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`

Configuration: \(N=80\), \(R_{\max}=80\), \(dr=1\), \(CFL=0.03\), nominal \(dt=0.03\), radiation ON, amplitude \(0.01\), width \(7\), D amplitude \(10^{-10}\); accepted-metric recovery OFF, radiation retries OFF.

Captured failure:

| Quantity | Value |
|---|---:|
| Last accepted coordinate time | \(t=46.98000000000102\) |
| Central proper time | \(\tau=22.01604689460664\) |
| Last accepted step index | 1566 |
| Attempted coordinate time | \(t=47.01000000000102\) |
| Stage | explicit predictor, first CMC lapse |
| Cell | \(i=79,\ r=79.5\) |
| Radiation \(E\) | \(3.97836632884\times10^{-8}\) |
| \(|S|\) | \(3.97859440216\times10^{-8}\) |
| \(|S|/E\) | \(1.00005732839\) |

The workflow passed its tests and captured the expected failure after 1,566 steps; it had 524 stage snapshots, 66 predictor budgets, 66 completed-state budgets, 67 geometry-predictor budgets, zero rejected attempts, and zero instrumentation errors. Another rate-corrected trace, [run 38310216296](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310216296), reproduced the same cell/stage/time pattern with a slightly different terminal ratio. The small numeric variation is retained rather than claiming bitwise identity across instrumentation/source revisions.

### Ordered matter/metric budget

The relevant conserved-variable cone margin and scale-invariant ratio are

\[
C=U_E-\sqrt{\gamma^{rr}}\,|U_r|,\qquad
\frac{|S|}{E}=\frac{\sqrt{\gamma^{rr}}\,|U_r|}{U_E}.
\]

For the latest rate-corrected trace at cell 79, step 1566, \(dt=0.03\):

| Ordered budget stage | \(C\) | Change from previous stage |
|---|---:|---:|
| Accepted radiation and accepted metric | \(+4.72908995\times10^{-7}\) | — |
| Flux update, evaluated on accepted metric | \(+1.48403208\times10^{-6}\) | \(+1.01112309\times10^{-6}\) |
| Flux plus geometric-source update, still on accepted metric | \(+2.21975975\times10^{-6}\) | \(+7.35727666\times10^{-7}\) |
| Same updated conservative state, predictor metric | \(-1.69974260\times10^{-8}\) | \(-2.23675717\times10^{-6}\) |

Budget closure error is exactly zero; the predictor conservative energy and momentum reconstruct from the recorded RHS with zero recorded error. Thus the RHS update alone, assessed on the accepted metric, increases the margin. The metric used at predictor recovery removes slightly more margin than is available and flips the sign.

The captured \(\gamma^{rr}\) rises from \(1.4628938913\) to \(1.4852172376\), a change of \(+0.0223233464\). In that trace's decomposition, \(+0.0139806007\) comes from updating a first and \(+0.0083427455\) from updating X; the algebraic regularity projection contributes zero to that trace's recorded \(\gamma^{rr}\) change. A separate, earlier-source stage ledger resolves the pre-projection state explicitly below.

## 2. Frozen stage-ledger decomposition at the prior and terminal step

### Provenance and validation

- Stage-ledger workflow: [38309949383](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309949383)
- Artifact: [ZIP artifact 11773874630](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309949383/artifacts/11773874630)
- Artifact ZIP SHA-256: `356f067a583ec3596c3b15e8b5d17a18058cc7bd7a2bbc7fd379b6b1e2acbe93`
- Trace JSON SHA-256: `d2e9d507f8300598c4dd5ac45b4bcd3f77d205adc388e460e6a75ae614a214c9`
- Stage-ledger JSON SHA-256: `35af5f350fe2888323cc851ab05190c5ddb981ef5042d656f4f45d4a803398bf`
- Configuration SHA-256: `1315fb36f2fd8e0ee3d064bc1d91cb9449829dea75f6f7863ae99e223e4503e7`
- Frozen source baseline: `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`
- Instrumentation commit: `527b44931b059f4d118a86c898b91d79ca4d96f2`
- Vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`

The full regression suite reported 141 passed. The ledger captured 100 accepted steps × cells 75–79; its stage-order validation passed with no errors. The captured source snapshot differs from the rate-corrected trace's source revision, so its detailed state values are reported as a separate frozen audit, not silently mixed with the rate-corrected trajectory.

### Hold the step's actual conservative radiation state fixed while changing only the metric

The saved ledger's \(\gamma^{rr}=X^2/a\) identity and explicit RHS reconstruction close at zero error for a and at \(2.8\times10^{-17}\) for X. At zero shift and with KO disabled in this diagnostic, the explicit identities match \(\dot a=-2\alpha a A_a\) and \(\dot X=\alpha X K/3\). For component attribution below, the computation order is: accepted metric → apply the a update only → apply the X update with predictor a held → apply the actual algebraic regularity projection. Because \(X^2/a\) is nonlinear, the split between the a and X increments depends slightly on this chosen ordering; the total metric change and cone-margin change do not.

| Ledger step | Stage | \(\gamma^{rr}\) | \(C\) | \(|S|/E\) |
|---:|---|---:|---:|---:|
| 1565, \(t=46.95\) | Accepted metric | 1.4413770169 | \(+2.58748840\times10^{-6}\) | 0.9911161620 |
| 1565 | After a update only | 1.4548222763 | \(+1.24424510\times10^{-6}\) | 0.9957280303 |
| 1565 | After X update, raw Euler metric | 1.4628580325 | \(+4.44398220\times10^{-7}\) | 0.9984742108 |
| 1565 | After regularity projection / first CMC inputs | 1.4628266993 | \(+4.47512736\times10^{-7}\) | 0.9984635175 |
| 1566, \(t=46.98\) | Accepted metric | 1.4630613182 | \(+2.21713013\times10^{-6}\) | 0.9925229671 |
| 1566 | After a update only | 1.4770455270 | \(+8.13948905\times10^{-7}\) | 0.9972550448 |
| 1566 | After X update, raw Euler metric | 1.4853905759 | \(-2.02334730\times10^{-8}\) | 1.0000682352 |
| 1566 | After regularity projection / first CMC inputs | 1.4853571835 | \(-1.69002052\times10^{-8}\) | 1.0000569941 |

At step 1565, the explicit metric update reduces C by \(2.14309018\times10^{-6}\), but the margin remains positive; projection restores \(3.11451627\times10^{-9}\). At step 1566, the a update reduces C by \(1.40318122\times10^{-6}\), and the subsequent X update reduces it by another \(8.34182252\times10^{-7}\). The total raw-predictor margin change is \(-2.23736360\times10^{-6}\), just enough to put the radiation state outside the cone.

The terminal projection changes a by \(+1.35603248\times10^{-5}\), b by \(+2.89450505\times10^{-5}\), and X by zero. Its net effect is to reduce \(\gamma^{rr}\) by about \(3.33923\times10^{-5}\), increase C by \(3.33326779\times10^{-9}\), and lower \(|S|/E\) by about \(1.1241\times10^{-5}\). The corrected state remains inadmissible. Across all 500 captured cell-step pairs, the ledger reports that projection increased the cone margin and reduced \(|S|/E\) in every pair; exactly one raw-predictor pair and one post-projection pair are outside the cone, both the terminal cell of the terminal step.

**Admitted attribution:** the raw explicit Euler geometry predictor is already outside the cone; projection does not create the crossing and slightly improves it. The a and X updates both contribute to increasing \(\gamma^{rr}\) under the stated decomposition order. b has no direct appearance in \(\gamma^{rr}\), but its value enters the algebraic determinant constraint \(ab^2=1\) and can affect the projection-induced adjustment of a. This is a numerical stage attribution, not proof of a defective physical law.

## 3. Existing controls: what they do and do not establish

| Control | Existing evidence | Result and interpretation |
|---|---|---|
| Accepted-metric radiation recovery, no retries | [Run 38018878024](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38018878024), \(N=80,\ R_{\max}=80\) | Does not reach target; the failure moves to the completed-state outer light-boundary stage at \(t=46.92,\ \tau\approx21.9985\), cell 79, with \(|S|/E\approx1.00115\). Not a repair. |
| Accepted-metric recovery plus adaptive admissibility retries | [Run 38023932307](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38023932307) | Eventually underflows at \(dt\approx6.82\times10^{-15}\) at the completed-state outer boundary around \(t=46.93255,\ \tau=22.00210\). Smaller rejected trials do not establish a completed physical step. |
| Larger domain, \(R_{\max}=120\), target \(t=50\) | [Run 38024259654](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38024259654) | Reaches \(t=50\), \(\tau=21.16458\), without failure; outermost-cell \(|S|/E\approx0.34766\). This horizon is earlier in central proper time than the failing runs. |
| Larger domain matched to late proper time | [Runs 38025556389](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38025556389) and [38028306184](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184) | At \(R_{\max}=120\), the outer cell \(i=119,\ r=119.5\) reaches the cone boundary near \(t=61.0465,\ \tau=21.91307\) during the explicit predictor first CMC lapse solve; retries underflow. The boundary audit's post-boundary P-minus mismatch is \(\sim1.0\times10^{-17}\), with its A-a correction formula error \(\sim1.2\times10^{-17}\). This weakens a simple allegation that the explicit boundary correction formula itself is algebraically wrong; the failing predictor stage occurs before that boundary update. |
| Rate-only matter RHS, frozen-dt trace | [Run 38310237842](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842) | Correctly gets past the false unit-step check at this configuration, then independently captures the real late-time outer-cell predictor crossing described above. |

These are not one perfectly matched convergence suite: source revisions and retry configurations differ. They do support two carefully bounded conclusions: accepted-metric recovery and shrinking the time step have not yielded an admitted continuation, and extending the outer boundary from 80 to 120 alone does not eliminate the late-time outer-cell crossing at a comparable late proper time. They do **not** prove the failure is a physical singularity, nor do they prove it is solely a boundary artifact.

## 4. Decision and next gate

1. **Do not change** the cone inequality \(E\geq |S|\), add a floor/clip/damping term, or insert a recovery metric override into production. Existing controls show those are not an admitted repair.
2. Keep the matter-RHS rate extraction correction separate from this event. Its purpose is to query the semi-discrete tendency without running a fictitious \(dt=1\) finite update; it retains validation on genuine finite-time updates.
3. Treat the current terminal cause as a stage-interface candidate: explicit metric predictor and updated radiation conservative state are not simultaneously primitive-admissible, even though the matter budget is admissible on the accepted metric.
4. Before another trajectory, perform one bounded frozen-input audit of the exact production source used by the prospective candidate, checking the signed a/X geometry RHS, metric/conserved-state stage synchronization, outer-cell source/flux closure, and relation to the existing H/M/\(C_\Lambda\), centre-series, and Misner–Sharp gates. Use the two existing terminal/preceding-step witnesses above as regression targets. Any counterfactual must be labelled diagnostic-only and must not change the equations used by the production run.
5. Preserve the no-bounce/no-turnaround guardrail. This radiation-cone failure is a numerical admission failure; it is not evidence for a bounce, turnaround, cycle completion, or a D-to-matter identification law.

No long campaign was started for this report. No production code, canonical GEAR file, or main branch was changed.


## 5. Source-path and upstream-boundary refinement (cross-agent reconciliation)

A complementary source-path audit is recorded in [RADIATION_PREDICTOR_TRACE_AFTER_RATE_FIX_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/RADIATION_PREDICTOR_TRACE_AFTER_RATE_FIX_2026-10-12.md), alongside the fixed-dt artifact. It refines what “predictor CMC failure” means and adds same-source field counterfactuals.

### Failure site is primitive recovery, not the CMC linear solve

The exact call path is: V55ProductionKernel.step forms the Euler geometry predictor and radiation conservative predictor, then calls _solve_lapse(grid, g_explicit1, s1, m1). At the start of solve_cmc_lapse, total_matter_projection(...) evaluates the supplied scalar/fluid sources. Radiation primitive recovery checks the invariant \(E\geq |S|\) and raises before the CMC matrix is assembled/solved to completion. Thus the exception label explicit_predictor_first_CMC_lapse identifies the calling stage; it is **not** evidence that the CMC linear system itself diverged.

### Same-source one-field counterfactuals

These values are all from the corrected-rate trace [run 38310237842](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842), with each step's predictor conservative radiation state held fixed while varying the metric. Unlike the historical ledger above, these are a single-source comparison.

| Accepted time | Metric used | \(\gamma^{rr}=X^2/a\) | Cone margin \(C\) | \(|S|/E\) |
|---:|---|---:|---:|---:|
| 46.92 | Accepted metric | 1.4203544029 | \(+2.92753094\times10^{-6}\) | 0.9897538451 |
| 46.92 | a update only | 1.4332809124 | \(+1.64361212\times10^{-6}\) | 0.9942474718 |
| 46.92 | X update only | 1.4280219498 | \(+2.16525438\times10^{-6}\) | 0.9924217601 |
| 46.92 | Both explicit updates | 1.4410182410 | \(+8.77874710\times10^{-7}\) | 0.9969274995 |
| 46.95 | Accepted metric | 1.4412159444 | \(+2.58989744\times10^{-6}\) | 0.9911068669 |
| 46.95 | a update only | 1.4546577813 | \(+1.24701279\times10^{-6}\) | 0.9957180348 |
| 46.95 | X update only | 1.4491752752 | \(+1.79398118\times10^{-6}\) | 0.9938398667 |
| 46.95 | Both explicit updates | 1.4626913466 | \(+4.47393488\times10^{-7}\) | 0.9984637500 |
| 46.98 | Accepted metric | 1.4628938913 | \(+2.21975975\times10^{-6}\) | 0.9925132639 |
| 46.98 | a update only | 1.4768744921 | \(+8.16949732\times10^{-7}\) | 0.9972446175 |
| 46.98 | X update only | 1.4711576616 | \(+1.38976919\times10^{-6}\) | 0.9953126300 |
| 46.98 | Both explicit updates | 1.4852172376 | \(-1.69974260\times10^{-8}\) | 1.0000573284 |

The same-source table adds a useful discriminant: **neither the a nor the X change alone crosses the cone on the terminal step; their joint explicit update does**. The state is increasingly close to the cone over the three steps, rather than showing a one-step discontinuity in the reconstruction identity.

For the terminal updated conservative state, the exact threshold is \(\gamma^{rr}_{\mathrm{crit}}=(U_E/|U_r|)^2=1.4850469621\). The predictor reaches \(1.4852172376\), exceeding that threshold by \(1.70276\times10^{-4}\), or \(0.01147\%\) of the threshold. The overshoot is \(0.7628\%\) of the total joint \(\gamma^{rr}\) increment. That quantifies how close the update was to preserving admissibility; it does not prescribe a timestep tweak or establish the proper repair.

### What the signed geometry equations establish

The pinned numerical operator's explicit BSSN equations are in [vacuum_bssn.py at vendor commit d6052d6](https://github.com/RonBibb/bb-palatini-unified-r0/blob/d6052d605673ce9d82cdc99b0df79855fa2ea215/vacuum_bssn.py#L268-L295). At the captured cell, shift is zero and KO strength is zero, so

\[
\dot a=-2\alpha a A_a,\qquad
\dot X=\frac{\alpha X K}{3}.
\]

On the corrected-rate trace at \(t=46.98\), \(\alpha=1\), \(A_a=0.157772386214\), \(K=0.282048234945\), \(a=0.608996684842\), and \(X=0.943873683323\). Substitution reproduces the saved rates \(\dot a=-0.192165720328\) and \(\dot X=+0.088739302131\); the recorded formula residuals are both zero. This establishes faithful execution of the pinned explicit equations at this cell. It does not establish that the incoming \(A_a,K\) satisfy all continuum equations, constraints, or boundary/CMC coupling requirements.

### Upstream boundary influence remains open

The boundary routine for the *attempted* \(t=47.01\) predictor has not run when radiation primitive recovery fails. But the accepted state at \(t=46.98\) already contains all prior completed-state updates, including any previous outer light-constraint boundary correction. At zero shift the next a-RHS uses the incoming \(A_a\) directly, so “the current boundary call occurs later” does **not** exonerate earlier boundary influence.

The existing expanded-domain trace [run 38028306184](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184), configuration \(N=120,R_{\max}=120\) and a different source revision, records a completed-state outer \(A_a\) correction of \(+1.185133219\times10^{-4}\) at \(t\approx60.99\) (\(\tau\approx21.89305\)). Its characteristic target/formula residual closes near \(10^{-17}\), but the *size of the correction* is not zero. Holding that correction alone in the zero-shift equation gives the local first-order estimate \(\delta\dot a=-2\alpha a\,\delta A_a\approx-1.62214\times10^{-4}\); over \(dt=0.03\), \(\delta a/a\approx-7.1108\times10^{-6}\), implying \(\delta\gamma^{rr}/\gamma^{rr}\approx+7.1108\times10^{-6}\) and \(\delta(|S|/E)/(|S|/E)\approx+3.5554\times10^{-6}\) at fixed conservative state. This is a scale estimate from a different domain/source, not a counterfactual replay of the N=80 trajectory.

Therefore the current-step boundary condition is not the direct operation at the failing substage, and the boundary formula's algebraic closure is strong; nevertheless, prior boundary-conditioned \(A_a\) remains a possible contributor to the later explicit metric RHS. The corrected-rate N=80 trace does not record matched pre/post-boundary \(A_a\) snapshots, so no causal boundary attribution is admitted. Any next study should compare those quantities on a matched source and state, while also inspecting \(A_a/K\) constraint/source evolution. The archived boundary trace is not sufficient to justify changing the boundary law.



## 6. Cumulative near-null approach and the boundary admission distinction

The full corrected-rate trace contains ten completed steps leading into the failed predictor, providing an existing (no-replay) trend:

- At outer cell \(i=79\), \(\gamma^{rr}\) increases monotonically from \(1.278974\) to \(1.462894\), a 14.38% rise.
- Its accepted densitized cone margin falls from \(3.7051\times10^{-6}\) to \(4.7291\times10^{-7}\), an 87.2% reduction.
- At every one of those completed steps, flux plus geometric source raise the margin when measured against the accepted metric, but the metric-stage contribution is more negative; the accepted margin therefore decreases step by step.
- At the attempted terminal predictor, cells \(i=75\ldots78\) remain inside the cone with margins approximately \(4.74, 5.10, 5.78, 8.60\times10^{-6}\), while cell \(i=79\) reaches \(-1.70\times10^{-8}\). Cells 75–78 actually have larger fractional \(\gamma^{rr}\) changes (about 2.06–2.81%) than cell 79 (about 1.53%). The distinguishing condition is the much smaller pre-existing outer-cell radiation margin, not simply the largest metric change.

Source: the existing ten-step table in [RADIATION_PREDICTOR_TRACE_AFTER_RATE_FIX_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/RADIATION_PREDICTOR_TRACE_AFTER_RATE_FIX_2026-10-12.md), corrected-rate run [38310216296](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310216296). All recorded ordered budgets close algebraically. This is evidence for a persistent near-null-state / metric-stage coupling in the captured discretization; it does not establish a singularity, physical null limit, defective continuum geometry equation, or bounce.

### Boundary residual is not boundary admission

The pinned numerical module bssn_characteristic_boundary.py explicitly states that it does **not** implement a complete constraint-preserving boundary. The V5.5 outer routine sets only Aa[-1] using a frozen second-order reconstruction of the incoming \(P_-\) light/constraint characteristic. In addition, the production CMC lapse solve normalizes \(\alpha[-1]=1\) instead of imposing its lapse PDE at that outer row.

The archived \(R_{\max}=120\) boundary artifact reports a post-correction \(P_-\) residual near \(10^{-17}\), but the routine solves for Aa[-1] from its own target, so that residual is an algebraic consistency check, not independent evidence that the target is the correct continuum incoming datum or compatible with the CMC/zero-shift gauge. The archived run also shows a material correction \(\Delta Aa\approx+1.1851\times10^{-4}\) at a late step. Because it is a different source revision and domain, that magnitude cannot be used as a causal replay of the current \(N=80\) failure.

Therefore:
1. the boundary correction at the *attempted failing step* has not yet run when primitive recovery throws;
2. previous accepted-step boundary-conditioned Aa still enters the next zero-shift \(\dot a=-2\alpha a Aa\), so upstream boundary influence remains possible;
3. the current \(N=80\) corrected trace does not record the matched pre-/post-boundary Aa values needed to quantify that pathway;
4. no boundary-law edit is justified from either a tiny post-residual or the nonzero correction magnitude alone.

The next existing diagnostic capture is on frozen historical source 9eed1cf95c52f001c2b87afbbe6da37faa56d50e, not the corrected-rate source 1c63090d4f42396f899421a2379f70320689d5ad. It may clarify how the historical Aa/K sources feed its predictor, but it must not be presented as a source-matched explanation of the corrected-rate trajectory. Compare formulas and structure, not raw trajectory values.

