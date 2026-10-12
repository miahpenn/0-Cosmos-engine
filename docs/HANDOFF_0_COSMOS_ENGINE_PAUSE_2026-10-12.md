# 0-Cosmos-engine — detailed handoff (pause point)

**Date:** 2026-10-12  
**Repository:** https://github.com/miahpenn/0-Cosmos-engine  
**Handoff branch:** `research/production-readiness-review-2026-10-12`  
**Status:** Investigation paused at an evidence-gathering boundary. No production physics change is authorized by this handoff.

## 1. Executive summary

The current investigation is testing whether the late radiation/metric failure is rooted in a local numerical defect, a source-term/stage-timing mismatch, a boundary-to-global-CMC coupling, or a genuine dynamical interface. These remain distinct hypotheses. No conclusion that the engine has a physical bounce, turnaround, cycle, or other specific global interpretation is admitted by the evidence summarized here.

Recent work improved diagnostic capture and separated two calls that had risked being conflated: the integrator's final explicit-state primary-L2 RHS and a post-completion diagnostic primary-L2 RHS. The diagnostic-only ledger retains phase identity and is intended to fail closed if call count or phase labels are unexpected. The latest observed ledger workflow runs were still in progress when inspected. This handoff therefore does **not** claim the latest workflow passed or its output was admitted.

A separate offline analysis was saved as `docs/RADIATION_PMINUS_GLOBAL_CMC_WEIGHT_BOUND_2026-10-12.md`. It derives conservative upper bounds on the global CMC target perturbation from partial outer-grid captures. Those bounds are informative but do not replace a matched corrected-rate failure-state full-grid calculation.

## 2. Active branches and key evidence

### 2.1 Radiation predictor metric-stage ledger

- Branch: `diag/radiation-predictor-metric-stage-ledger-2026-10-12`
- Source/capture work examined at commit `1777b5bda7c6cfd46a4fc7580d69ffbfa14879da`; a workflow-only trigger revision followed on this branch.
- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/workflows/radiation-predictor-metric-stage-ledger.yml
- Recent runs inspected: [38313377947](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38313377947), [38313377854](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38313377854), [38313351496](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38313351496).

The run/job inspection showed these runs were in progress at the time checked. Fetch final outcomes, job logs, artifact list, and hashes before using their evidence. The recorder labels primary-L2 captures `integrator_final_l2_rhs` and `post_completion_diagnostic_l2_rhs`, storing the first two non-nested calls under separate stable keys. A third call is retained under an explicit unexpected-call key. The validator expects two calls on non-failure steps and validates phase/index, cell indices, and required finite fields.

**Instrumentation limitation:** the L2 snapshot contains only the last five cells (`CELLS = 5`). That is appropriate for the outer-cone/stage-timing question but cannot reconstruct a full-grid global CMC target or full lapse response.

### 2.2 Corrected-rate radiation work and trace

- Branch: `research/production-readiness-review-2026-10-12`
- Corrected-rate source commit previously inspected: `1c63090d4f42396f899421a2379f70320689d5ad`.
- Report: `docs/RADIATION_PREDICTOR_TRACE_AFTER_RATE_FIX_2026-10-12.md`.
- Related investigation: `docs/RADIATION_PREDICTOR_TERM_ATTRIBUTION_2026-10-12.md`.

The corrected-rate trajectory has an unresolved radiation/metric failure; this is not cleared by a passing capture harness or unrelated regression test. The production-readiness question remains whether the recorded metric cone/radiation RHS defect is produced by the actual production source path and at which time stage the causal divergence begins. Keep conservative radiation state, metric reconstruction, geometric source terms, primitive recovery, and outer-boundary treatment separate.

### 2.3 Historical partial artifacts and global CMC boundary bound

Saved offline analysis: `docs/RADIATION_PMINUS_GLOBAL_CMC_WEIGHT_BOUND_2026-10-12.md`.

Historical evidence inspected:
- [Run 38311614754](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38311614754), artifact `11774209967`, N=80 outer five-cell ledger near t=46.95.
- [Run 38028306184](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184), artifact `11661332148`, historical N=120 trace near t=60.99.
- [Run 38210163221](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221), artifact `11729296911`, initial-state all-grid CMC profile; **not** the late radiation-failure state.

The partial-band method uses positive proper-volume weights to upper-bound the outer cell's full-grid share. For the historical N=80 capture, the bound on the full-grid target shift was approximately 6.90e-6. For the historical N=120 capture, it was approximately 2.02e-6. These are different historical trajectories/source contexts. Neither is the exact corrected-rate failure-state target perturbation, and neither quantifies the solved-lapse response.

## 3. Established vs. open

### Established within inspected artifacts/code

1. The diagnostic recorder distinguishes two non-nested primary-L2 calls by phase.
2. The current ledger is outer-band instrumentation, not a complete-grid checkpoint.
3. Under the pinned spherical BSSN relation Ab = -Aa/2, the local boundary-induced K-RHS increment is calculable from saved pre/post Aa values: ΔLK[-1] = (3/2) alpha[-1] (Aa_post² - Aa_pre²).
4. A positive partial sum of proper-volume weights gives a conservative upper bound on the boundary cell's share of the complete-grid weight, provided omitted weights are non-negative.
5. Existing partial artifacts do not supply a complete late-time corrected-rate failure state, so they cannot yield the exact global ΔKdot_CMC or paired late-time solved-lapse response.
6. Passing diagnostic infrastructure does not admit the physics result, full trajectory, or production readiness.

### Open questions

1. On the corrected-rate source and exact failure trajectory, which captured A/B/C stage first shows a meaningful radiation-cone/metric inconsistency?
2. Does divergence originate in centre-face flux, geometric source, primitive recovery, or a coupling/order mismatch? Planned observational comparison: N=80/160 × beta=0/-0.04.
3. Is the outer P-minus boundary operation materially coupled to the global CMC target and solved lapse at the failure time?
4. Does a matched full-grid capture show that the CMC/lapse effect is decision-relevant compared with the local radiation/metric residual?
5. Are the ledger's captured stages complete, phase-correct, and deterministically admitted by CI? Latest workflow outcomes were not yet verified.
6. What does the full bidirectional machine do when matter, radiation, scalar, and geometry lanes are coupled without artificial interface rules? Do not infer the answer from an isolated lane.

## 4. Next steps — in order

### Step 1 — close current workflow admission

Fetch final status, job logs, artifact list, and hashes for runs `38313377947`, `38313377854`, and `38313351496`. For each candidate result, record exact commit SHA/branch, run URL/conclusion, artifact name/ID/size/SHA-256, admission/validation lines and stage labels, corrected-rate source provenance, and whether the run is an instrumentation control or physical trajectory. Do not use an in-progress run as evidence. If it fails, diagnose the harness/test/physics layer specifically and change only that layer.

### Step 2 — reconcile A/B/C capture specification with code

Keep capture strictly observational: no changes to fluxes, sources, primitive recovery, guards, state updates, or boundary rules. Four cases remain N=80/160 × beta=0/-0.04, with three named moments A/B/C. Ensure the report says exactly what each moment means in actual integrator call order. Preserve same-run source and state hashes so comparisons cannot silently cross revisions.

### Step 3 — identify first divergence from admitted captures

Compare at the same named stage and cell range. Separate centre-face numerical flux, geometric source terms, conservative radiation state, primitive recovery, reconstructed metric/cone invariant, boundary operation, CMC target, and lapse solve. Use existing commit history and archived data first; do not rerun covered cases merely to regenerate a familiar symptom. If a new full run is needed, run only one expensive campaign at a time.

### Step 4 — decide whether full-grid CMC sensitivity capture is justified

If the CMC boundary pathway remains decision-relevant, add a diagnostic-only capture at the existing production boundary call on the same snapshot, containing the minimum full-grid data needed to reproduce the CMC target and lapse solve. Capture pre/post Aa, all-grid a,b,X,alpha,K, proper-volume weights, pointwise K-RHS, and matter/scalar terms consumed by the pinned solver. Run paired frozen-state solves with only boundary Aa changed, recomputing the target consistently. Record full-grid residuals and verify all other fields are identical.

Do not approximate missing inner-grid weights by extrapolation. Do not substitute initial-state lapse profiles for late-time failure states. Do not modify the boundary rule or production CMC definition to make the response small.

### Step 5 — reconcile with blueprint and full engine

Use repository GEAR files and authoritative REV16 criterion as blueprints, and inspect commits/results before proposing a repeat test. Preserve the standards: full bidirectional machine rather than isolated toy success; no hand tuning, arbitrary damping, artificial floors/clamps, fitted feedback coefficients, invented bounce law, or silent physics; diagnostics/admission before interpretation; diagnostic/campaign changes separate from physics changes; exact-commit CI success; machine evidence over preferred interpretation.

## 5. Guardrails

- Do not merge the diagnostic branch into `main` merely because its tests pass.
- Do not claim the latest stage-ledger workflow passed until final run and artifact are checked.
- Do not call partial-band CMC estimates exact full-grid target changes.
- Do not claim solved-lapse sensitivity from a target-shift bound alone.
- Do not mix historical N=120/source data with corrected-rate N=80 results without explicit provenance labels.
- Do not change physics to make a numerical result pass.
- Do not declare a bounce, cycle, turnaround, or physical interface from these diagnostic traces alone.
- If evidence conflicts, preserve both provenance trails and reconcile source/commit before interpretation.

## 6. Session changes

1. Inspected the current radiation stage-ledger implementation and validator at the pinned diagnostic commit.
2. Inspected recent ledger workflow runs; they were still in progress at the time of inspection.
3. Recovered historical radiation artifacts and computed conservative boundary-to-global-CMC weight bounds from captured outer bands.
4. Added `docs/RADIATION_PMINUS_GLOBAL_CMC_WEIGHT_BOUND_2026-10-12.md` with the offline result.
5. No production physics equation, source term, flux, recovery, guard, or boundary rule was changed.

## 7. Resume checklist

- [ ] Verify final results of runs 38313377947 / 38313377854 / 38313351496.
- [ ] Record artifact hashes and admission lines.
- [ ] Confirm A/B/C capture semantics against actual call order.
- [ ] Review the N=80/160 × beta=0/-0.04 observational matrix without changing physics.
- [ ] Identify the first divergence stage from admitted same-source captures.
- [ ] Decide, based on evidence, whether full-grid late-time CMC/lapse sensitivity is necessary.
- [ ] Update this handoff with exact commits, outputs, and falsifiable conclusions before any physics modification.

---
**Handoff disposition:** paused, diagnostic-only. Next action is evidence admission and stage reconciliation, not another speculative physics change.