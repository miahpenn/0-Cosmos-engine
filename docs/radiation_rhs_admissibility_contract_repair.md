# Radiation RHS admissibility contract repair

> **Historical candidate note — 13 October 2026.** This document records the PR #10-era guard-placement candidate and its narrow checks. Its embedded statements that PRs #8/#9/#10/#11 remain draft or unmerged are historical and are not the live repository status. Use [docs/STATUS.md](STATUS.md) for current production/admission state and [Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135) for active work. Preserve the failed-run and candidate evidence below; do not use this document as a current implementation or merge-status checklist.

Status: implementation regression gates pass; independent implementation review and admission remain pending.

## Exact repair validation

### Original guard-placement candidate

- Tested implementation head: `5a5a9657e5d0e45c899addd5851ab2bd8db1cfe0`.
- CI run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317104269 (success).
- Focused contract / real-step regressions: **4 passed**.
- The command `python -m pytest -q tests` reported **26 passed**, covering only the `tests/` subtree at that commit; it was not the repository-wide suite.

### Latest repository-wide gate on the repair branch

- CI/workflow-only head: `d1338a8473aa8bc2474019f99b13001edde221a6`.
- CI run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318333422 (success).
- Focused repair tests: **4 passed**.
- Repository-wide `python -m pytest -q`: **83 passed**.
- This commit changes the workflow's suite command and does not change the engine implementation. The run checks out the pinned vendor submodule `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

Commit `13f76f9d1bc711db4c784629e649c21791c38b0f` / run `38316799110`, and the intermediate documentation-tested head `5a5a9657e5d0e45c899addd5851ab2bd8db1cfe0` / run `38317104269`, remain useful historical results. They are not the newest full-repository-suite pair.

## Defect boundary

The matter RHS wrapper calls the finite-volume operator with `dt=1.0` and subtracts the original conserved state to extract the rate. Applying the physical update guard to that intermediate mistakenly validates `U + L(U)`, rather than the timestep-scaled state the production integrator admits.

## Repair

- `species_rhs` calls `evolve_species` with `validate_admissibility=False`, because the unit-step result is only an algebraic RHS extraction.
- Direct `evolve_species` callers retain validation by default.
- A shared, non-mutating radiation-state validator enforces the existing `E >= 0` and `|S| <= E` bound and rejects non-finite values.
- The production kernel checks the actual initial radiation state, timestep-scaled explicit predictor, primary predictor, and completed final radiation state against the corresponding stage metric.
- No projection, clipping, retries, damping, flux/source/recovery changes, or new physics are introduced.

## Related implementation candidate requiring comparison

PR #8 (https://github.com/miahpenn/0-Cosmos-engine/pull/8) is an earlier main-based alternative. Instead of bypassing only the guard during RHS extraction, it separates the semi-discrete flux/source calculation into `_species_rhs_components` and returns rates directly. Its main-based regression run `38307704487` passed 2 focused tests and the then-current repository-wide suite (**81 passed**); its matched seven-resolution diagnostic is `38307235876`. Neither PR #8 nor PR #10 has a formal review submission at the time of this handoff.

These approaches must be compared before selecting the implementation to promote. PR #8 is a more direct rate-extraction API, but it reorganizes arithmetic into a rate followed by multiplication by `dt`, so floating-point operation order and regression equivalence deserve explicit scrutiny. PR #10 keeps the existing update implementation and moves guard use to actual state stages, but must prove its skip is confined to algebraic extraction. This note does not choose between them. Review both contracts, tests, and production-kernel call paths.

## Regression gates

1. At N=160, radiation RHS rates must be extractable and finite.
2. A constructed real candidate with `|S|/E = 2` must be rejected by the same validator.
3. A full baseline N=160 production step must complete with stage checks and the final state must pass the validator.
4. The latest repair branch must pass focused tests and the repository-wide suite.

## Post-repair four-force counterfactual (separate verification composition)

The follow-up diagnostic combines PR #10's guard-placement repair with the diagnostic-only `beta_rad` counterfactual from draft PR #9 on a temporary verification branch. It is not a production branch and must not be merged.

- Verification head: `bd48e1307e4ffb5872c88edbea77595c35dca803`.
- CI run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317257752 (success).
- Artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317257752/artifacts/11777401783.
- JSON SHA-256: `58ee12cfdd548b363b3308db0515a68c83686dacde665f85406ecf3a318c471c`.
- All four one-step cases passed: N=80/160, each with `beta_rad=0.0` and `-0.04`.

The centre geometric source and first interior HLL flux were unchanged at captured precision between null and candidate. Global endpoint differences were extremely small in this one-step probe. The Misner–Sharp current remains unavailable (NaN), so this does not provide complete term-by-term conservation attribution or establish physical validity of `beta_rad=-0.04`. See [the full post-repair result](https://github.com/miahpenn/0-Cosmos-engine/blob/verify/post-repair-radiation-counterfactual-2026-10-12/docs/radiation_four_force_post_repair_counterfactual_result.md).

## Separate representation issue: radiation rest-current

Draft PR #11 (https://github.com/miahpenn/0-Cosmos-engine/pull/11) proposes routing radiation HLL/outer-face fluxes through existing radiation-specific helpers to retain the exact-zero radiation `rest` slot. Its one-step A/B saw no change to E/S, cone observables, geometry, or H/M/C_Lambda at captured precision. This is a separate representation-invariant repair, not evidence for or against the RHS contract fix; it is stacked on this PR and remains unreviewed/unmerged. The late-time cone crossing and Misner–Sharp balance remain open.

## Admission boundary

The current PR #10 branch has passed the focused gate and full repository-wide suite, and the post-repair four-force probe completes its one-step matrix. These are narrow software/execution results, not full production readiness or candidate-physics validation. The independent comparison of PR #8 and PR #10, the separate PR #11 representation review, longer-horizon tests, and Misner–Sharp balance remain open.

- PR #10 remains draft/unmerged pending independent review.
- PR #8 remains draft/unmerged pending comparison.
- PR #11 remains draft/unmerged and stacked on PR #10.
- PR #9 remains draft/unmerged and diagnostic-only.
- No coupling promotion, no physics claim, no guard relaxation, no clipping or damping.
