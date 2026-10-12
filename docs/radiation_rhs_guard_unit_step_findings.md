# Radiation RHS unit-step guard finding

Status: diagnostic finding; repair not yet applied  
Branch: `diag/radiation-four-force-counterfactual-2026-10-12`  
PR: https://github.com/miahpenn/0-Cosmos-engine/pull/9 (must remain draft)  
Capture run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38316101510  
Structured rejection run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38316279301\nLatest JSON artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38316279301/artifacts/11776975591\nLatest JSON SHA-256: `f6ab04a363bfbde8195f278f44445e93638822b1e5139904b6315b1630230a8f`  
Parent artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38316101510/artifacts/11776544112

## Finding

The N=160 centre rejection is produced while extracting the matter RHS, not after applying the caller's physical timestep.

The call chain is:

1. `V55ProductionKernel.step(state, dt)` calls `self._rhs(...)`.
2. `_rhs` calls `species_rhs(...)` for radiation.
3. `species_rhs` calls `evolve_species(..., dt=1.0)` and then returns `advanced - state` as the derivative.
4. The radiation admissibility guard currently runs inside `evolve_species`, before `species_rhs` subtracts the original state. Therefore it tests the unit-step state `U + L(U)`, even though the integrator will use `U + dt L(U)`.

This explains why the guard can reject a derivative-extraction call before the production integrator has applied its actual timestep. This is a guard-placement / API-contract defect candidate, not evidence that the radiation four-force candidate is physically valid or invalid.

## Captured N=160 centre values

Both `beta_rad=0` and `beta_rad=-0.04` begin with a recovered centre primitive at rest and `|S|/E = 0`. The captured centre-face boundary flux is the existing zero vector. At the first interior face, the radial-momentum HLL flux is approximately `2.3399337052591616e-10`.

For the null case, the unit-step operator produces:
- flux momentum tendency: approximately `-9.35973482103665e-10`;
- geometric momentum source: approximately `+7.48761506593284e-10`;
- unit-step rejected ratio: `|S|/E ~= 1.3585922116`.

For `beta_rad=-0.04`, the corresponding rejected ratio is approximately `1.3584029299`. The centre source and first-face flux are the same to the displayed precision in the two cases. The small ratio difference does not establish a coupling effect on the centre transport mechanism.

## Actual-timestep scale check

Using the recorded unit-step operator difference as `L(U)`, and the configured caller timestep `dt=0.0075`, the diagnostic first-stage candidate

`U_candidate = U + dt * (U_unit_step - U)`

has `|S|/E = 0.0101313642296932` for `beta_rad=0` and `0.0101313537021683` for `beta_rad=-0.04`, below the existing `|S| <= E` bound. These structured values are recorded under `capture_C.dt_scaled_rhs_candidate` in the latest artifact. This is an observational, first-stage timestep-scaled estimate using the pre-step metric. It is **not** the completed PIRK predictor or final state, and it does not prove that the complete stage remains admissible.

## Gate status

- Focused null-identity tests: passed (3 tests).
- Existing test suite: passed (25 tests).
- Four-case runner: expected to exit nonzero while the N=160 unit-step guard remains in place; artifacts are still uploaded.
- No guard has been relaxed, no flux/source/recovery equation has been changed, and no coupling or physics interpretation has been promoted.

## Next repair boundary

Keep the diagnostic probe and physics counterfactual separate from the repair. A separate, reviewed repair should align the admissibility check with the actual state being admitted: either keep RHS extraction free of a unit-step candidate guard and validate the real timestep-scaled predictor/final states, or define an explicit validation mode so the RHS wrapper cannot accidentally invoke the update guard. Preserve a regression that proves the initial RHS can be evaluated and that actual-stage candidates are checked without projection, clipping, damping, or altered flux/source terms.

Do not merge PR #9 as a fix. The current result diagnoses the failure mechanism; it does not validate the production machine.
