# Boundary–CMC history counterfactual and radiation predictor budget

**Date:** 12 October 2026 (UTC)  
**Classification:** diagnostic-only; no production equations, matter source laws, gauge conditions, boundary law, or model parameters changed.  
**Source baseline:** `cdcec8282bff0443f2fe4674b15929822ff38736`  
**Pinned geometry vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Run:** [38321399848](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321399848)  
**Artifact:** [11778478114 — boundary-cmc-radiation-audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321399848/artifacts/11778478114)  
**Artifact ZIP:** `boundary-cmc-radiation-audit-v2.zip`  
**Trace JSON SHA-256:** `f0bcfd73b6911b4693e4f43483974b19b036deeed7dd4c2d0d29f3464c7bf5a7`  
**Trace workflow source SHA:** `c5c225d6e01c265522d4f27291f6f1492ee12b73`

## Executive finding

The one-step counterfactual restoring the *preceding completed-boundary call's pre-update* outer (A_a[-1]) does **not** rescue the next radiation predictor. Both actual and counterfactual one-step updates reject radiation at the same outer cell and at the same explicit-predictor admissibility gate.

| Witness | Actual accepted state | Counterfactual accepted state |
|---|---:|---:|
| Outer (A_a[-1]) at (t=47.13) | (0.15026677099199126) | (0.15025766090063600) |
| Explicit-predictor (|S|/E) at (r=79.5) | (1.000401250165977) | (1.000401244072738) |
| Radiation admission | FAIL | FAIL |

The tested perturbation is (Delta A_a[-1]=-9.1100913553\times10^{-6}). It changes the failing ratio by only (6.093239\times10^{-9}), around **0.00152% of the ratio's excess above one**. Therefore the immediately preceding boundary update, by itself, is not a sufficient cause of this observed failure. This single-step test does **not** exclude cumulative effects from earlier boundary conditioning or prove the boundary's full dynamic compatibility.

### Actual rejecting event

- Last accepted coordinate time: (t=47.130000000001026).
- Attempted step: (t=47.16000000000103), (\Delta t=0.03).
- First failing gate: explicit-predictor radiation admissibility, before the predictor's CMC source/lapse solve and before this attempted step's P-minus boundary call.
- Cell: (i=79, r=79.5).
- Actual-stage values: (E=3.4249001439887531\times10^{-8}), (|S|=3.4262743857399831\times10^{-8}), (|S|/E=1.000401250165977).

The counterfactual exception message also identifies the explicit-predictor guard at (i=79), (r=79.5), with (E=3.424899145209403\times10^{-8}), (|S|=3.426273365691143\times10^{-8}), and ratio (1.000401244072738). The artifact's counterfactual `failure_stage` field incorrectly retains its first CMC-call label because this guard is thrown directly by the production step, not inside the wrapped lapse solver. The exception text is the authoritative stage witness. The tracer's stage-label parser was corrected subsequently in commit `e7e9d4d713633c5de8fc5f6e0e5490b2749ab8cd`; this metadata fix does not change the captured numerical values.

## Boundary call witness

At the preceding completed-state P-minus call (outer cell (i=79, r=79.5)):

| Witness | Pre-boundary | Post-boundary |
|---|---:|---:|
| (A_a[-1]) | (0.150257660900636) | (0.15026677099199126) |
| Zero-shift (dot a[-1]=-2\alpha a A_a) | (-0.17805883526202831) | (-0.17806963093291875) |
| P-minus residual | (-3.6806778524\times10^{-6}) | (6.94\times10^{-18}) |
| Radiation (|S|/E) on same completed slice | (0.9988576966) | (0.9988576966) |

The P-minus residual closing after the update is expected: the boundary value is solved from this same discrete relation. It demonstrates algebraic target closure, not independent proof of a complete constraint-preserving boundary or compatibility with the CMC-normalized time evolution. The current attempted step fails before its own boundary call, so prior-step influence is the only relevant boundary pathway for the immediate event.

## Source-separated cone-margin budget

The corrected-rate/source-capture trace contains a complete, exactly closing cone-margin budget at the preceding positive-margin predictor (accepted (t=47.10), predictor attempt (t=47.13), outer cell (i=79)):

| Stage / contribution | (C=U_E-\sqrt{\gamma^{rr}}|U_r|) |
|---|---:|
| Accepted state and accepted metric | (+6.33426220995\times10^{-7}) |
| After flux increment on accepted metric | (+1.37273306437\times10^{-6}) |
| After geometric/source increment on accepted metric | (+1.99264075334\times10^{-6}) |
| Same predicted matter evaluated on predictor metric | (+2.66313550817\times10^{-7}) |

The signed increments are (Delta C_{flux}=+7.3930684337\times10^{-7}), (Delta C_{source}=+6.1990768897\times10^{-7}), and (Delta C_{metric}=-1.7263272025\times10^{-6}), with exact recorded closure error (0.0). This reinforces the stage-level result: flux and source improve the available margin on the accepted metric; the explicit predictor metric consumes more than those increments restore. The next attempted predictor crosses the cone.

This trace version did not invoke the completed-state trapezoidal budget writer (`completed_state_budget_count=0)); its predictor budgets do report `source_capture_missing=false`, and the recorded predictor margin budget closes exactly. A follow-up instrumentation run was queued to retain the completed-state budget as well. No new long physical trajectory or production change follows from that instrumentation enhancement.

## Disposition and next step

1. Do not modify or disable the P-minus boundary on the basis of this result. The immediately preceding single update is too small to explain/rescue the measured crossing in the matched one-step counterfactual; cumulative boundary effects remain untested.
2. Keep the exact (E\ge|S|) criterion and current no-clamp/no-floor rule.
3. Continue the equation-level audit of the explicit predictor metric update, CMC slice/outer normalization, and conserved radiation state as one coupled stage. Keep any additional test bounded and stage-matched; record raw (a,X,A_a,K), lapse, (gamma^{rr}), (C), and source/flux contributions.
4. Do not merge diagnostic branches or claim scientific production admission from this run. `main` and the GEAR archive remain unchanged.

The run completed as an **expected failure capture** with the full repository regression suite passing before the trace. That status means the diagnostic reproduced and recorded the numerical failure, not that the production trajectory passed.


---

## Follow-up accounting-complete run — source-separated completed-state budget

**Run:** [38321456272](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321456272) — completed successfully as an expected radiation-failure capture after the full repository suite passed.  
**Artifact:** [11779635426 — boundary-cmc-radiation-audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321456272/artifacts/11779635426)  
**Artifact ZIP SHA-256:** `448f8184571a7d0d0e4603e53531ef3148f07c0b565acfa250e31a230a8377d0`  
**Trace JSON SHA-256:** `06138c2b31db6d6bd2f7807a4cdb655d535eb48ec392c22aaae897ff855391ce`  
**Summary JSON SHA-256:** `64638067e3e30d6d379e6e0a4b38cbe19b58905981b85a0cf8d53da6a5cf67ec`  
**Instrumentation commit:** `bc7cf1ab052644457d7c3a1b5fd404d736605ec0`; same source/vendor baseline as above.

This run restores the completed-state budget call that was absent from the earlier boundary-history run. It captured:
- 65 source-separated explicit-predictor cone-margin budgets;
- 66 accepted-to-completed trapezoidal budgets;
- 66 geometry predictor budgets;
- zero instrumentation errors;
- zero missing source-capture cells in the completed-state budgets;
- maximum completed-state margin-budget closure error (5.421010862427522\times10^{-20});
- exact (zero recorded) predictor RHS reconstruction errors for energy and radial momentum, exact (zero recorded) completed trapezoidal reconstruction errors, zero recorded explicit-metric identity closure error, and zero recorded predictor cone-budget closure error.

At the final completed candidate (t=47.130000000001026), cell band (i=75\ldots79), the saved accepted-state (C) at (i=79) is (2.806443726585477\times10^{-7}). The accepted-to-completed margin budget remains consistent: matter flux and source increments increase the margin on the accepted metric, while the accepted-to-RHS1 predictor metric change is negative. The next attempted explicit predictor fails at (t=47.16000000000103), (i=79, r=79.5), with (|S|/E=1.0004012501699084).

The same matched preceding-boundary history counterfactual reproduces the result: resetting only the accepted outer (A_a[-1]) to its immediate pre-boundary value still fails at (i=79, r=79.5), with (|S|/E=1.0004012440766743). This is an insignificant change to the threshold excess and does not rescue admissibility. The exception text establishes the explicit-predictor guard as the failed stage; the counterfactual object's `failure_stage` field remains a stale internal label (`accepted_start_CMC_lapse`) on this tested workflow commit. No inference is drawn from that metadata field.

This follow-up closes the previous diagnostic's missing completed-state accounting. It does not promote the boundary, establish covariant conservation, or admit the trajectory. A smaller-(dt) one-step comparison from the saved accepted state is now the next bounded discriminator for the predictor-stage failure; any successful smaller step remains a diagnostic-only observation, not permission to add an adaptive production policy automatically.
