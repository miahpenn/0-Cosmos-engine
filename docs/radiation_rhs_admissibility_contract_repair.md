# Radiation RHS admissibility contract repair

Status: implementation and CI regression gates pass; independent review and admission pending.\nValidated head: `13f76f9d1bc711db4c784629e649c21791c38b0f`.\nCI run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38316799110 This repair is intentionally separate from PR #9 and does not add or promote a radiation four-force.

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
4. The existing full test suite must pass.\n\n## Current result\n\n- Focused repair tests: **4 passed** (including a direct-update guard-retention check and one full N=160 production step).\n- Existing test suite: **26 passed**.\n- No physics validation is inferred from these tests; the four-force counterfactual remains in draft PR #9.

A passing single step is a narrow regression result, not full production readiness. The diagnostic PR #9 remains draft and independent.
