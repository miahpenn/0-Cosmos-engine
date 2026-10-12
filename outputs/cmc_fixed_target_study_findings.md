# CMC fixed-target and self-consistency iteration study

Date: 2026-10-12  
Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38204727992  
Artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38204727992/artifacts/11726897293  
Diagnostic script: https://github.com/miahpenn/0-Cosmos-engine/blob/84f0480d2d36c6a3eee300384d60092c66ead3cf/diagnostics/cmc_fixed_target_study.py

## Scope

N=64, r_max=32, localized S amplitude 1e-3, D amplitude 0, radiation off. Geometry other than alpha, scalar fields and matter fields were frozen. No production equations were changed.

## Results

- The fixed-target control behaves exactly as expected: first solve shifts the input lapse by RMS 4.3787449e-9, and the next three solves reproduce the first output to reported RMS difference 0.0. The tridiagonal solver is not accumulating the shift when its target is held fixed.
- With the production target recalculated from each returned lapse, the iteration did **not** converge within 200 iterations.
- The per-iteration lapse correction stayed near 4.375e-9 RMS; its last contraction ratio was 1.00000495, not below 1 in a meaningful way.
- The output lapse drifted from the initial lapse by 8.7539e-7 RMS after 200 updates.
- The target drifted from its initial value by about -7.1755e-9 after 200 updates. The test remained finite and the fixed-target control passed.

## Interpretation

This is strong evidence that ordinary repeated substitution (Picard iteration) of the existing target/solve pair is not converging at this weak-field initial slice. Holding the target fixed removes the repeated shift. The measured feedback is therefore associated with recomputing the target from the evolving lapse input, not with a changing linear system or a failed linear solve.

This does **not** prove that no self-consistent solution exists. A non-convergent fixed-point iteration can fail even when another solution method would find a root, and the result applies only to this weak-field frozen slice. In particular, do not conclude that the gauge has no solution or extrapolate to strong-D proper time without additional tests.

## Next test

Before any production change, characterize the residual function directly around the candidate lapse: evaluate the scalar fixed-point defect F(alpha)=solve_cmc_lapse(alpha-dependent target)-alpha along a small controlled family of lapse perturbations, using a diagnostic-only target/solve wrapper. This will estimate the local response/slope and distinguish an affine near-constant correction from a possible root with a poorly conditioned Picard map. Repeat on a representative strong-D snapshot only after identifying a robust way to capture that state.
