# Completed-State Boundary History — Follow-up Stage Budget

**Date:** 12 October 2026 (UTC)  
**Classification:** Evidence addendum only; no production source or boundary-law change.  
**Source snapshot:** `cdcec8282bff0443f2fe4674b15929822ff38736` (stacked radiation RHS/rest-current candidate).  
**Vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`.  
**Workflow:** [38321456272](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321456272).  
**Artifact:** [11779635426 — boundary-cmc-radiation-audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321456272/artifacts/11779635426).  
**Artifact ZIP SHA-256:** `448f8184571a7d0d0e4603e53531ef3148f07c0b565acfa250e31a230a8377d0`.  
**Embedded trace SHA-256:** `06138c2b31db6d6bd2f7807a4cdb655d535eb48ec392c22aaae897ff855391ce`.  
**Run outcome:** expected radiation predictor rejection reproduced and captured; 90 repository-wide tests passed before trace execution.

## Why this follow-up exists

The earlier [boundary-history report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/boundary-cmc-coupling-audit-2026-10-12/docs/BOUNDARY_CMC_HISTORY_COUNTERFACTUAL_2026-10-12.md) lacked the completed-state trapezoidal matter-margin budget. This follow-up retained it. It complements, but must not be numerically pooled with, the primary-predictor boundary counterfactual on the rate-corrected source baseline described in [the adjacent report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/radiation-primary-boundary-cf-2026-10-12/docs/RADIATION_PREDICTOR_PRIMARY_BOUNDARY_CMC_COUNTERFACTUAL_2026-10-12.md).

The run captures 522 stage snapshots, 65 predictor budgets, 66 geometry-predictor budgets and 66 completed-state budgets. Both geometry-identity and completed-state-budget error counters are zero; there are zero instrumentation errors. Conservative radiation variables reconstructed from RHS0 and from the trapezoidal RHS0/RHS1 update are identical at stored precision.

## Failure and completed-boundary history counterfactual

The actual accepted state reaches \(t=47.130000000001026\), \(\tau=22.020654806515967\). The next attempted step (\(dt=0.03\), target \(t=47.16000000000103\)) fails at the outer cell \(i=79,\ r=79.5\), during the **explicit-predictor** radiation guard, before that attempted step's P-minus boundary call.

| Witness | Actual | Restore pre-boundary \(A_a[-1]\) only |
|---|---:|---:|
| Accepted outer \(A_a[-1]\) at \(t=47.13\) | 0.15026677150036508 | 0.15025766140928215 |
| Preceding completed-state boundary increment | \(+9.11009108293\times10^{-6}\) | removed for this one stored component |
| Predictor \(|S|/E\) at \(i=79\) | 1.0004012501699084 | 1.0004012440766743 |
| Radiation admissibility | **FAIL** | **FAIL** |

The ratio difference is \(6.0932341\times10^{-9}\), only about **0.00152% of the ratio excess above one**. The counterfactual exception text identifies the explicit-predictor guard; the trace's legacy `failure_stage` field says `accepted_start_CMC_lapse` because the exception is raised outside the wrapped lapse solver. The exception text and `actual_failure_stage` field are the authoritative witnesses; this label mismatch is instrumentation metadata, not a physical result.

The P-minus residual changes from \(-3.6806777384\times10^{-6}\) before the boundary operation to \(6.94\times10^{-18}\) after it, by design of the boundary solve. This checks the formula's internal discrete target closure; it is not independent proof of a complete constraint-preserving boundary.

## Full source/metric budget on the preceding positive-margin step

At accepted \(t=47.10\), the predictor targeting \(t=47.13\), cell \(i=79\), the predictor-only budget is:

| Stage | Cone margin \(C=U_E-\sqrt{\gamma^{rr}}|U_r|\) |
|---|---:|
| Accepted state and metric | \(+6.33426202561\times10^{-7}\) |
| After flux increment on accepted metric | \(+1.37273304763\times10^{-6}\) |
| After geometric/source increment on accepted metric | \(+1.99264074702\times10^{-6}\) |
| Same predicted matter on explicit predictor metric | \(+2.66313538935\times10^{-7}\) |

The signed changes are \(\Delta C_{\rm flux}=+7.39306845066\times10^{-7}\), \(\Delta C_{\rm source}=+6.19907699398\times10^{-7}\), and \(\Delta C_{\rm metric}=-1.72632720809\times10^{-6}\). Recorded budget closure is zero. Flux and sources restore positive margin when evaluated against the accepted metric; the explicit metric update consumes more than that restoration.

The complete RHS0/RHS1 trapezoidal budget then checks the *completed* state on the same step:

| Completed-state witness | Value |
|---|---:|
| Accepted-state \(C\) | \(+6.33426202561\times10^{-7}\) |
| Completed conservative state on accepted metric | \(+2.01945659426\times10^{-6}\) |
| Completed conservative state on RHS1 spatial metric | \(+2.93680795754\times10^{-7}\) |
| Completed conservative state on final completed metric | \(+2.80644372659\times10^{-7}\) |
| Total accepted-to-completed margin change | \(-3.52781829903\times10^{-7}\) |

The metric conversion from accepted to RHS1 stage contributes \(-1.72577579850\times10^{-6}\); the RHS1-to-completed spatial-metric change contributes another \(-1.30364230953\times10^{-8}\). The signed flux/source and metric pieces close to the observed final margin with zero stored reconstruction error. This shows that the completed slice remains admissible but leaves a small reserve for the following explicit predictor, where the actual guard rejects the state.

## Exact geometry identity and attribution limit

At the preceding predictor cell \(i=79\), the pinned zero-shift equations are executed exactly:
\[
\dot a=-2\alpha a A_a,\qquad \dot X=\frac{\alpha X K}{3}.
\]
The signed residuals against the captured explicit geometry RHS are both zero. The inverse-metric decomposition is:

| Quantity | Value |
|---|---:|
| Accepted \(\gamma^{rr}\) | 1.50991851671 |
| After \(a\)-update only | 1.52337881154 |
| After \(X\)-update, before any regularity projection | 1.53138262359 |
| Contribution from \(a\) | \(+0.01346029483\) |
| Contribution from \(X\) | \(+0.00800381205\) |
| Contribution from regularity projection | 0 |
| Total \(\Delta\gamma^{rr}\) | \(+0.02146410688\) |

This identifies the numerical mechanism: the explicit \(a/X\) predictor increases \(\gamma^{rr}\) and consumes most of the radiation margin on this near-saturated outer state. The current evidence still does **not** distinguish an integrator/stage compatibility problem from an incorrect upstream continuum or boundary-conditioned trajectory; the formula identities demonstrate faithful evaluation of the pinned equations, not their complete constraint/conservation compatibility.

## Decision and next step

- Do not change or disable the P-minus boundary based on this counterfactual; restoring the preceding pre-boundary \(A_a\) does not rescue the predictor.
- Preserve the exact \(E\geq |S|\) condition. No clamp, floor, radiation-cone relaxation, damping, fitted coefficient or forced retry rule is justified.
- Keep the boundary law, primary-predictor CMC response, centre-series/\(A_a,K\) compatibility and metric–matter stage synchronization as separate questions.
- Reuse this budget in the next source/geometry review. A new long run is not justified until a specific signed stage identity fails and a bounded counterfactual distinguishes the implicated equation or split.

No production main, physics equation, pinned vendor, or canonical GEAR archive was changed by this follow-up diagnostic.
