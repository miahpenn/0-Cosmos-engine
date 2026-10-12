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
