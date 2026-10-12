# Matched Outer-Boundary \(A_a\) Pair Audit — Radiation Predictor Failure

**Date:** 12 October 2026 (UTC)  
**Status:** Diagnostic-only; production admission remains open.  
**Production main:** 2035acc1210db18da0cededea18750b7ed813260 (unchanged).  
**Diagnostic branch:** diag/radiation-boundary-aa-pair-2026-10-12  
**Trace workflow:** [38320777962](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38320777962) — workflow success means regression and diagnostic capture succeeded, not that the physical trajectory reached its requested horizon.  
**Artifact:** [radiation-boundary-aa-pair, ID 11778836819](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38320777962/artifacts/11778836819)  
**Artifact ZIP SHA-256:** 30d8a69179caba7640bad5d0cb1244aaf7ad49b9954394ff8f5d059d3b96f0c1  
**Embedded trace SHA-256:** 86e63096dec6ccc35f9f071c044222bec583ab2a5375aef2b806db49e731fba7  
**Embedded summary SHA-256:** c63c1be61892c4f026edf617168d4cb2ac106d6d72f1c85b823c9c6c64b04ca0  
**Instrumented source head:** 4b07ee7ed373f8414a5f20ca3640e846d106bc7a  
**Rate-extraction-corrected source baseline:** 1c63090d4f42396f899421a2379f70320689d5ad  
**Pinned numerical vendor:** d6052d605673ce9d82cdc99b0df79855fa2ea215

## Executive result

The paired recorder captures outer-five-cell fields immediately before and after **both** production outer-light-boundary calls: the primary predictor call and the completed-state call. It matches each completed-state post-boundary value to the accepted state at the next step; this is the relevant pair when asking whether that boundary-written outer Aa influences the next explicit a-predictor.

The corrected-rate trace again captures the known late-time radiation-cone violation. At the last accepted state \(t=46.98000000000102\), central proper time \(\tau=22.016046893100416\), the next step's explicit predictor fails at \(t=47.01000000000102\), cell \(i=79\), \(r=79.5\). The actual predictor has \(C=-1.6997416818\times10^{-8}\) and \(|S|/E=1.00005732835\).

**The preceding completed-boundary Aa correction makes this cell's next predictor margin slightly worse in a frozen, direct-channel counterfactual, but removing that correction does not restore admissibility.** The counterfactual changes only the incoming accepted Aa used by the local zero-shift identity for \(\dot a\); it holds accepted \(\alpha,K,X\), the predicted radiation conserved state, and the measured regularity-projection \(\Delta a\) fixed. It is not a replay of the coupled step and excludes indirect changes through the CMC solve or other stage feedback.

**Decision:** the present artifact does not justify changing or removing the outer boundary law. The direct \(Aa\to\dot a\to\gamma^{rr}\to C\) influence is measurable but insufficient to explain or repair the late crossing by itself. Continue to treat the explicit metric/matter-stage compatibility issue as open.

## Run and recorder admission

Configuration: \(N=80\), \(R_{\max}=80\), \(\Delta r=1\), fixed \(dt=0.03\), radiation ON, scalar amplitude \(0.01\), scalar width \(7\), \(D\) amplitude \(10^{-10}\), no accepted-metric recovery, and no radiation admissibility retries. The trace source is based on the rate-extraction-corrected engine and the pinned vendor listed above.

The workflow:
- checked that engine/vendor/dependency changes relative to the frozen rate-fix source were limited to diagnostic instrumentation;
- verified the vendor gitlink;
- compiled the recorder;
- passed **94 repository-wide tests**;
- recorded the expected radiation predictor failure and uploaded the trace;
- reported **132 boundary-pair records**, **335 matched completed-boundary outer-cell comparisons**, **0 Aa mismatches** against the associated accepted state, and **0 instrumentation errors**.

The first trial of this new test had two exact-float assertions and failed at test level (93 passed, 1 failed); these were replaced with numeric tolerances. That first run did not start the trajectory and was superseded by the successful full gate above. It was not a physics failure.

The trace has status expected_radiation_failure_captured, not a completed \(t=50\) evolution. There are 1,566 accepted steps before the captured rejection and no retry attempts.

## Paired boundary measurements at the outermost cell

In the table, the displayed target time is the time of the accepted state that receives the completed-state boundary result. \(\Delta A_a\) is post-boundary minus pre-boundary.

| Accepted target \(t\) | \(A_a\) pre-boundary | \(A_a\) post-boundary | \(\Delta A_a\) |
|---:|---:|---:|---:|
| 46.92000000000102 | 0.1502468015826905 | 0.1503137469129523 | \(+6.6945330\times10^{-5}\) |
| 46.95000000000102 | 0.1539393306852895 | 0.1540091532150225 | \(+6.9822530\times10^{-5}\) |
| 46.98000000000102 | 0.1576999717970437 | 0.1577723855707334 | \(+7.2413774\times10^{-5}\) |

At these three points the boundary correction increases monotonically in magnitude, while the accepted radiation cone margin decreases. Correlation alone is not causal attribution. The paired instrumentation confirms that the recorded post-boundary value equals the Aa used in the next step at the same cell exactly at stored precision; all matched comparisons passed the explicit Aa consistency check.

The boundary routine itself directly changes \(A_a[-1]\), not \(a,b,X\) or the conservative radiation variables. Consequently, the immediate before/after snapshots have zero direct \(\Delta C\) and zero direct \(\Delta\gamma^{rr}\). Its possible next-step effect is assessed through the explicit geometry equation \(\dot a=-2\alpha a A_a\) at zero shift and zero KO strength.

## Frozen counterfactual at the failing predictor

The actual step starts at \(t=46.98000000000102\), with cell \(i=79\). The immediately preceding completed-boundary pair is from the step targeting that accepted state.

| Witness | Actual accepted \(A_a\) | Replace \(A_a\) by pre-boundary value (frozen counterfactual) |
|---|---:|---:|
| \(A_a\) entering failing step | 0.1577723855707334 | 0.1576999717970437 |
| \(\dot a\) from \(-2\alpha aA_a\) | \(-0.1921657198214381\) | \(-0.1920775203250831\) |
| Predicted \(a\) after retaining measured projection delta | 0.6032317141253851 | 0.6032343601102758 |
| \(X\) held fixed | 0.9465358619432119 | 0.9465358619432119 |
| Predictor \(\gamma^{rr}=X^2/a\) | 1.4852172340500571 | 1.485210719397344 |
| Predictor cone margin \(C=U_E-\sqrt{\gamma^{rr}}|U_r|\) | \(-1.6997416818\times10^{-8}\) | \(-1.6347122273\times10^{-8}\) |
| \(|S|/E\) | 1.00005732835491 | 1.00005513506184 |

The boundary increment changes the local a rate by \(-8.8199496\times10^{-5}\), and the one-step predicted a by \(-2.6459849\times10^{-6}\) relative to the counterfactual. Holding the other listed fields fixed, that raises the predicted inverse radial metric by about \(6.514653\times10^{-6}\) and lowers the cone margin by \(6.5029455\times10^{-10}\).

This is approximately **0.029%** of the full accepted-metric-to-predictor-metric cone-margin loss (\(2.2367572\times10^{-6}\)) and about **3.8%** of the magnitude of the final negative margin. Most importantly, the no-boundary-increment margin remains negative (\(-1.6347\times10^{-8}\)); the same predictor is still inadmissible. The direct boundary increment is therefore a small contributor in this bounded calculation, not a sufficient standalone cause or repair target.

## What this establishes—and what it does not

**Established for the captured source and state:**
1. The completed-state boundary call adds \(7.2414\times10^{-5}\) to the outermost \(A_a\) at the accepted \(t=46.98\) state.
2. That post-boundary value is the accepted value consumed by the following explicit a-predictor; no matching discrepancies or instrumentation errors occurred.
3. Replacing only that increment with its measured pre-boundary \(A_a\) slightly improves the cone margin, but does not restore \(E\ge |S|\).
4. The genuine predictor failure remains at the same time, cell, and stage pattern seen in the prior corrected-rate trace.

**Not established:**
- This does not measure a fully coupled alternate trajectory; \(\alpha\), the CMC solution, the next \(K/A_a\) source response, or other stage dependencies are intentionally held fixed.
- It does not prove the boundary condition is correct or incorrect, and it does not rule out wider upstream influence from boundary-conditioned evolution over earlier steps.
- It does not validate centre regularity, H/M/\(C_\Lambda\) propagation, Misner-Sharp balance, or full production readiness.
- It is not evidence for a physical bounce, turnaround, or completed cosmological cycle.

## Disposition and next diagnostic

Keep the \(P_-\) boundary formula and production physics unchanged. Do not add a clamp, floor, damping, retry rule, fitted correction, or radiation-cone relaxation.

The next useful bounded question is no longer “does the previous completed-boundary \(\Delta A_a\) by itself cause the cone crossing?”—this test says no. Reuse this exact trace to examine the **source/stage coupling** that produces the accumulating \(A_a,K\) and \(a/X\) trajectory, particularly whether the primary-predictor boundary increment changes the next explicit RHS or the CMC/lapse response in a way missed by the frozen-alpha test. That should be a separately declared counterfactual with its own state/source hashes and joint H/M/\(C_\Lambda\), centre-series and cone witnesses. Do not launch another long trajectory or alter the boundary law unless a specific signed identity fails and a bounded counterfactual supports the change.

No production code, main, pinned vendor, or canonical GEAR archive was modified by this audit.
