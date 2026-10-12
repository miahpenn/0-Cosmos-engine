# PIRK Lambda Predictor Stage-Order Audit — 2026-10-12

## Disposition

**A source-level stage-order discrepancy is confirmed and a diagnostic-only correction passes its regression and repository test gates. Its one-step global effect is extremely small; it does not explain or fix the saved late radiation-cone crossing by itself. No promotion to `main`.**

- Candidate branch: `diag/pirk-lambda-predictor-stage-order-2026-10-12`.
- Tested source/workflow head: `ce9a49560b04ccbe33d5b51ac86219a5c4a3c2f9`.
- Code correction commit: `a725730b8f536cda8f086b286114f6003ec187ca`.
- Test: `tests/test_pirk_lambda_predictor_stage_order.py`.
- Bounded A/B driver: `tools/diag_pirk_lambda_predictor_stage_order_ab.py`.
- Workflow: [38318309196](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196) — success.
- Targeted stage-order test: **1 passed**.
- Full repository test suite: **36 passed**.
- Pinned vendor checked in CI: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- A/B artifact: [11777443806](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196/artifacts/11777443806).
- Artifact ZIP SHA-256: `1ec986b7cbbe97ba23184c2d448aba26199257f6ae9d6c4182490c81af3867fe`.
- Embedded `report.json` SHA-256: `283c5b337a90670c9a397c9584d19447497017cc44a640ca96bc7303a03186ef`.

## Exact ordering discrepancy

The pinned source is `RonBibb/bb-palatini-unified-r0` at the hash above. In `moving_puncture_pirk.py`, the sequential PIRK2 predictor performs:

1. Construct the explicit `U1) stage, holding the old primary variables (V_0=(A_a,K)).
2. Evaluate the primary L2 blocks on the initial and `U1,V0` states.
3. Form (V_1=(A_{a,1},K_1)).
4. Evaluate `lambda_l2_rhs(grid, stage1)` on ((U_1,V_1)), then form predictor (Lambda_1) with the trapezoidal L2 pair plus the initial L3 term.

The previous production code computed `gterms_pred = geometry_stage_terms(..., g_explicit1, ...)` before assigning the new `g1.Aa` and `g1.K`, then reused `gterms_pred["lambda_l2"]` for (Lambda_1). That evaluated the (Lambda) L2 block on ((U_1,V_0)), not ((U_1,V_1)).

The discrepancy matters algebraically because the pinned `lambda_l2_rhs` includes the terms involving (A_a\partial_r\alpha), (alpha\partial_r A_a), the spherical (A_a-A_b) contribution, and the momentum constraint. It is not a metric-only function of (U_1).

The isolated correction now computes `vacuum.lambda_l2_rhs(grid, g1, lambda_m=LAMBDA_M)` **after** the predictor (A_a,K) assignments and before the existing outer-boundary conditioning and subsequent CMC lapse solve. It leaves the source terms, CMC gauge, boundary algorithm, timestep, cone criterion, and all coefficients unchanged. The new unit test spies on the actual call ordering and verifies the predictor (Lambda) equation uses the value from that updated stage.

## Matched one-step A/B

The audit ran one ordinary production step for each arm at N=80 and N=160, `r_max=80`, (mathrm{CFL}=0.015), using bitwise-identical initial states. The legacy arm uses the cached second-call (Lambda)-L2 vector (explicit (U_1) with old (A_a,K)) at the third call to reproduce the prior ordering; the candidate arm uses the actual updated-stage vector. This isolates only the predictor (Lambda)-L2 input. No long trajectory was launched.

| Quantity | N=80 | N=160 |
|---|---:|---:|
| (dt) | 0.015 | 0.0075 |
| max (|lambda_{L2}(U_1,V_1)-lambda_{L2}(U_1,V_0)|) | (5.3217372642\times10^{-8}) | (1.3313374252\times10^{-8}) |
| Expected predictor (max|\Delta\Lambda_1|) from this block | (3.9913029482\times10^{-10}) | (4.9925153444\times10^{-11}) |
| Observed predictor (max|\Delta\Lambda_1|) | (3.9913029479\times10^{-10}) | (4.9925153421\times10^{-11}) |
| Final-state (max|\Delta A_a|) | (1.6330\times10^{-12}) | (1.9526\times10^{-13}) |
| Final-state (max|\Delta\Lambda|) | (6.94\times10^{-17}) | (2.17\times10^{-18}) |
| Final-state (max|\Delta a|) | (2.22\times10^{-16}) | 0 |
| Final-state (max|\Delta X|) | 0 | 0 |
| Change in full-domain (max|H|) | (-2.43\times10^{-17}) | 0 |
| Change in full-domain (max|M|) | (-4.2540\times10^{-12}) | (-1.0631\times10^{-12}) |
| Change in full-domain (max|C_\Lambda|) | (-6.94\times10^{-17}) | (-2.17\times10^{-18}) |

All four steps (legacy/candidate at each resolution) completed and remained radiation-cone admissible. Scalar and conservative matter fields were bitwise identical at the final stage in this short comparison. The initial-state final cone ratios were small (max full-domain (|S|/E\approx0.00501)); this is not a stress test of the nearly saturated late-time state.

The stage value itself changes by the expected (dt/2) times the L2-vector difference, which confirms the correction is active and isolated. Most of that difference does not survive as a final (Lambda) difference because the later split and completed-stage update follow. The one-step constraint and metric deltas are correspondingly tiny. There is no evidence from this A/B that the ordering discrepancy is a dominant driver of the radiation-cone event.

## Relation to the late radiation failure

The saved corrected-rate trace rejects the attempt from accepted (t=46.98000000000102) toward (t=47.01000000000102), at the **first explicit predictor before the first CMC solve**, cell (i=79). That rejection happens while the explicit metric predictor (a,b,X) and matter predictor are being formed, before the predictor (A_a,K,Lambda) stage corrected here. Therefore this change cannot directly rescue that same step's explicit-predictor cone crossing. It could only change later states indirectly through accumulated evolution and the existing boundary/gauge coupling; that has not been tested.

The separately archived geometry/matter cone budget remains the stronger local explanation of the actual rejected stage: the radiation update improves the cone margin on the accepted metric, and the explicit predictor metric then consumes more than that margin. The stage-order correction is a standards/conformance repair, **not** a substitute for resolving that geometric predictor crossing.

## Gate and next decision

- Retain this isolated branch and its test as the reference-sequence correction candidate.
- Do not merge it to `main` solely on this one-step result. A source-matched replay would be required to establish any effect on the (t\approx46.98) event; that is intentionally not run here.
- The outer-boundary/CMC compatibility audit already exists in the cross-agent research records; do not duplicate it as a generic follow-up. See [the corrected-rate trace report, including its boundary/CMC sections](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/RADIATION_PREDICTOR_TRACE_AFTER_RATE_FIX_2026-10-12.md) and [the P-minus constraint-combination derivation](https://github.com/miahpenn/0-Cosmos-engine/blob/research/radiation-predictor-term-attribution-2026-10-12/docs/RADIATION_PREDICTOR_TERM_ATTRIBUTION_2026-10-12.md). The remaining gap is the exact matched-state boundary-to-CMC/history attribution at the corrected-rate failure; the recorded failing explicit predictor precedes the boundary operation in that step. Keep this as an open coupling question, not a proven boundary defect.
- Keep centre-series, H/M/C\(_\Lambda\), radiation admissibility, and Misner–Sharp criteria independent. No conservation tolerance is inferred by this one-step A/B.

No physical law, damping, floor, clamp, projection, source coefficient, or boundary condition was changed by the candidate.


## Final confirmed A/B values and checkpoint — 12 October 2026

The bounded A/B and its checksums are now fully recorded in [the experiment report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/pirk-lambda-predictor-stage-order-2026-10-12/docs/PIRK_LAMBDA_PREDICTOR_STAGE_ORDER_AB_2026-10-12.md). Exact retained provenance:
- A/B workflow [38318309196](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196), success at `ce9a49560b04ccbe33d5b51ac86219a5c4a3c2f9`.
- Stage-order CI [38318119471](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318119471), success: 1 focused test and 36 full-suite tests passed.
- Artifact [11777443806](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196/artifacts/11777443806), ZIP SHA-256 `1ec986b7cbbe97ba23184c2d448aba26199257f6ae9d6c4182490c81af3867fe`; report JSON SHA-256 `283c5b337a90670c9a397c9584d19447497017cc44a640ca96bc7303a03186ef`.
- Exact pinned vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`; production `main` remains `2035acc1210db18da0cededea18750b7ed813260`.

| Measured difference, corrected minus stale | N=80 | N=160 |
|---|---:|---:|
| Predictor Lambda delta | 3.9913029479e-10 | 4.9925153421e-11 |
| Completed Lambda max-array delta | 6.9388939039e-17 | 2.1684043450e-18 |
| Completed Aa max-array delta | 1.6329877488e-12 | 1.9526410357e-13 |
| Completed H max-absolute delta | -2.4286128664e-17 | 0.0 |
| Completed M max-absolute delta | -4.2540321716e-12 | -1.0631220974e-12 |
| Completed C_Lambda max-absolute delta | -6.9388939039e-17 | -2.1684043450e-18 |

The raw Lambda-L2 stage arrays differ by 5.3217372642e-8 (N=80) and 1.3313374252e-8 (N=160), confirming the stage-order discrepancy. After the full step, the max `|S|/E` is identical between arms (0.00501354074617448 at N=80; 0.005015096011912016 at N=160), as are the radiation cone margins (2.2067153165e-9 and 5.5163762656e-10 respectively). All four one-step arms completed.

**Decision:** preserve the sequential stage-order correction as a discrete conformance candidate, but do not promote it or infer it fixes the accepted-t≈46.98 outer cone crossing. The full-run causal question, centre regularity, constraint propagation, Misner–Sharp balance, convergence and production admission remain open.
