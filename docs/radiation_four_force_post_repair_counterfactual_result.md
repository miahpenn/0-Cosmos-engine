# Post-repair radiation four-force counterfactual result

Status: diagnostic result only; no physics admission or promotion.

## Exact execution

- Verification branch: `verify/post-repair-radiation-counterfactual-2026-10-12`
- Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317257752
- Tested head: `bd48e1307e4ffb5872c88edbea77595c35dca803`
- JSON artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317257752/artifacts/11777401783
- JSON SHA-256: `58ee12cfdd548b363b3308db0515a68c83686dacde665f85406ecf3a318c471c`
- Focused radiation null-identity tests: passed.
- Existing full test suite: passed.
- Four one-step cases: **4 passed, 0 failed**.

This branch is a temporary verification composition: it combines the guard-placement repair from draft PR #10 with the diagnostic-only beta_rad counterfactual from draft PR #9. It is not a production branch and should not be merged.

## Fixed case matrix and outcome

| Resolution | beta_rad | dt | Outcome |
|---:|---:|---:|---|
| 80 | 0.0 | 0.015 | passed |
| 80 | -0.04 | 0.015 | passed |
| 160 | 0.0 | 0.0075 | passed |
| 160 | -0.04 | 0.0075 | passed |

The N=160 baseline now completes a real production step without the former unit-step RHS extraction rejection. The real-state validator remains active for initial, predictor, and final states.

## Observed null-versus-candidate deltas

The centre recovered primitive remains at rest with `|S|/E = 0` in cells 0–2 before the step. The centre geometric source and first interior HLL momentum flux are identical between null and candidate at the captured precision for each resolution. This supports the narrow statement that the prior centre failure was the guard/API contract, not a candidate-specific centre flux/source change.

The global endpoint witnesses differ only at very small scales for this one-step probe. Examples:

- N=80: `H_eff` changes from `0.0056686720884344755` to `0.005668672088434430`; maximum admitted Hamiltonian residual changes by about `-4.02e-14`.
- N=160: `H_eff` changes from `0.005668820081124814` to `0.005668820081124803`; maximum admitted Hamiltonian residual changes by about `-1.95e-14`.
- The admitted momentum residual is about `4.2e-6` after the step in both arms.

These endpoint changes are observations, not evidence of a meaningful or correct physical coupling. The artifact does not supply a complete term-by-term conservation/balance attribution, and the Misner–Sharp current remains explicitly unavailable (NaN). No conclusion about the physical validity of `beta_rad=-0.04` follows.

## Disposition

- PR #10 remains draft/unmerged pending independent review.
- PR #9 remains draft/unmerged and diagnostic-only.
- No coupling promotion, no physics claim, no guard relaxation, no clipping or damping.
- The next step is review of the narrow repair and the captured artifact; a broader production-readiness claim still requires the full established admission gates.
