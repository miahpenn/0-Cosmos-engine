# Radiation RHS admissibility contract repair

Status: implementation regression gate passes; independent review and admission pending.

## Exact repair validation

- Tested repair head: `5a5a9657e5d0e45c899addd5851ab2bd8db1cfe0`.
- CI run on that head: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317104269 (completed successfully).
- Focused contract/real-step tests: **4 passed**.
- Existing test suite: **26 passed**.
- The passing workflow includes the N=160 real production-step regression and direct-update guard-retention test.

Commit `13f76f9d1bc711db4c784629e649c21791c38b0f` and run `38316799110` are earlier validation history; they are not the latest head/run pair. The subsequent documentation-only fix is included in tested head `5a5a9657e5d0e45c899addd5851ab2bd8db1cfe0`.

## Defect boundary

The matter RHS wrapper calls the finite-volume operator with `dt=1.0` and subtracts the original conserved state to extract the rate. Applying the physical update guard to that intermediate mistakenly validates `U + L(U)`, rather than the timestep-scaled state the production integrator admits.

## Repair

- `species_rhs` calls `evolve_species` with `validate_admissibility=False`, because the unit-step result is only an algebraic RHS extraction.
- Direct `evolve_species` callers retain validation by default.
- A shared, non-mutating radiation-state validator enforces the existing `E >= 0` and `|S| <= E` bound and rejects non-finite values.
- The production kernel checks the actual initial radiation state, timestep-scaled explicit predictor, primary predictor, and completed final radiation state against the corresponding stage metric.
- No projection, clipping, retries, damping, flux/source/recovery changes, or new physics are introduced.

## Regression gates

1. At N=160, radiation RHS rates must be extractable and finite.
2. A constructed real candidate with `|S|/E = 2` must be rejected by the same validator.
3. A full baseline N=160 production step must complete with stage checks and the final state must pass the validator.
4. The existing full test suite must pass.

## Post-repair four-force counterfactual (separate verification composition)

The follow-up diagnostic combines this repair with the diagnostic-only `beta_rad` counterfactual from draft PR #9 on a temporary verification branch. It is not a production branch and must not be merged.

- Verification head: `bd48e1307e4ffb5872c88edbea77595c35dca803`.
- CI run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317257752 (completed successfully).
- Artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317257752/artifacts/11777401783.
- JSON SHA-256: `58ee12cfdd548b363b3308db0515a68c83686dacde665f85406ecf3a318c471c`.
- All four one-step cases passed: N=80/160, each with `beta_rad=0.0` and `-0.04`.

The centre geometric source and first interior HLL flux were unchanged at captured precision between null and candidate. Global endpoint differences were extremely small in this one-step probe. The Misner–Sharp current remains unavailable (NaN), so this does not provide complete term-by-term conservation attribution or establish physical validity of `beta_rad=-0.04`. See [the full post-repair result](radiation_four_force_post_repair_counterfactual_result.md).

## Admission boundary

This establishes the narrow API-contract regression and one real N=160 production step. It does not establish full production readiness, validate the candidate coupling, or close the Misner–Sharp balance gap.

- PR #10 remains draft/unmerged pending independent review.
- PR #9 remains draft/unmerged and diagnostic-only.
- No coupling promotion, no physics claim, no guard relaxation, no clipping or damping.
