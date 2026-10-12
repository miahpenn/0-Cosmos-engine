# CMC lapse same-slice residual audit

Date: 2026-10-12  
Diagnostic branch: `diag/bond-sweep-self-consistent`  
Diagnostic commit: `2b3c2a8c7b924729d2f1a47500a9743b99712a4d`  
Workflow run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38204378856  
JSON artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38204378856/artifacts/11726318438  
Script: https://github.com/miahpenn/0-Cosmos-engine/blob/2b3c2a8c7b924729d2f1a47500a9743b99712a4d/diagnostics/cmc_lapse_residual_audit.py

## Scope

The diagnostic uses the actual `V55ProductionKernel` initialized at N=64, r_max=32, localized S amplitude 1e-3, D amplitude 0, radiation off. It freezes non-lapse geometry, scalar and matter fields, and performs eight repeated CMC lapse solves, feeding each returned lapse into the next call. Production equations were not changed. The script reconstructs the exact tridiagonal operator assembled in `engine/cmc_gauge.py`, including the center and pinned outer boundary rows.

## Measured result

- The lapse shift per solve is effectively constant: RMS range **4.3786015e-9 to 4.3787613e-9** across eight calls.
- Maximum RMS change between consecutive shift profiles is **6.8952e-14**. Consecutive profile cosine similarity is approximately **1.0**.
- The radial profile is pinned to zero at the outer cell and reaches about **-6.0124e-9** near the center.
- Reconstructing the linear solve with the target actually used by that solve gives a maximum relative residual of **2.997e-15**. The tridiagonal solve itself is accurate to floating-point precision for its supplied right-hand side.
- Recomputing `target_kdot` using the returned lapse changes the target by approximately **-3.607e-11 per call**. With that recomputed target, the same output lapse has a maximum discrete equation residual of about **3.608e-11**; the residual does not shrink across the repeated calls.

## Interpretation

This isolates the effect from linear-solver inaccuracy: the solve accurately returns the lapse for the target computed from its *input* lapse, but the CMC target changes when evaluated on the *output* lapse. The coupled solve is therefore not a fixed point after one pass, and repeated calls at this state do not converge toward one; they add almost the same radial correction each time.

That establishes a measurable target-dependence/fixed-point residual in this small inhomogeneous case. It does **not** yet establish whether this is an implementation defect or a consequence of the intended target definition, nor does it measure how the effect scales in the strong-D regime or explain the strong-D proper-time difference. Those remain separate tests.

## Next discriminating test

Keep the production equations unchanged. Compare (A) the residual and `target_kdot` response with the lapse included in the target evaluation, against (B) a diagnostic-only evaluation where the target is calculated from a controlled, unchanged lapse input while the linear solve is reconstructed separately. Record the per-call target, equation residual, lapse correction, and boundary rows. Only after that comparison should we decide whether the CMC target/solve coupling needs a production-layer repair.
