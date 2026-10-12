# CMC target-closure rank audit — 2026-10-12

**Run:** [38294955148](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294955148) — successful  
**Artifact:** [11768056222 — cmc-target-closure-rank-probe](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294955148/artifacts/11768056222)  
**Artifact ZIP SHA-256:** 64cc1666a2010d73b462cac48c78706ac48be1a1d02fea7de343bd09172c4e24  
**Diagnostic source:** diagnostics/cmc_target_closure_rank_probe.py  
**Diagnostic commit:** be7dac470eca22a1267dbbeac7776719afbd88d6  
**Inputs:** archived t=0 R80 and R160 profiles; both ZIPs hash-verified by workflow.  
**Production source/physics:** unchanged; no evolution run.

## Main result

The CMC target-closure residual is effectively flat across nine deliberately varied fixed target values. On the archived t=0 R80 slice, with the R160 lapse value imposed at the R80 outer cell:

- Trial Kdot spans 6.1294749e-5 to 6.3028915e-5.
- Recomputed target-output slope dKdot_out/dKdot_in = 0.9999999997573141.
- Closure residual slope dH/dKdot = -2.4268635370e-10, where H = Kdot_out - Kdot_trial.
- H ranges only from -8.60513e-14 to -8.56135e-14, a total range 4.37882e-16.
- The returned lapse is affine in the trial target to numerical precision: max RMS deviation from the endpoint-interpolated profile 4.48e-15; max pointwise deviation 6.66e-15.
- The elliptic PDE rows close with max residual 3.71e-15. The maximum discrepancy between H and the weighted mean pointwise K-RHS residual is 1.14e-20.
- The outer normalization cell's K-RHS residual is approximately -4.58e-12, and that cell carries 1.86291% of the projection weight.

## Interpretation

This is strong numerical evidence that the current global target closure is **nearly redundant with the interior elliptic equation** on this frozen slice. The CMC solve enforces the pointwise K RHS equal to the imposed target on every interior PDE row; the target functional then volume-averages that same RHS. The remaining H about -8.58e-14 is quantitatively accounted for by the outer normalization cell, whose lapse is prescribed rather than solved by an elliptic PDE row but whose K RHS still participates in the average.

This does **not** prove an exact continuum family of solutions, and it does not establish that the production target definition is physically wrong. It establishes that the present discrete target equation supplies almost no independent rank to select Kdot over the tested target interval. Therefore the target/lapse closure cannot be treated as a robust independent root condition without first deriving the intended global CMC condition and its discrete boundary contribution.

## Relationship to previous CMC findings

The earlier pointwise closure audit found, for three targets, H about -8.58e-14, with H matching the weighted mean pointwise residual to 1.6e-21–4.0e-21. Interior PDE rows closed to 3.6e-15–5.4e-15 max; the largest residual was at the outer normalization cell, about -4.5882e-12. This nine-target rank probe now shows the same seam persists across a target interval, rather than being an isolated near-root.

The initial-slice target/lapse iteration audit also found tiny steps (about 6.4e-11 RMS on boundary-matched R80 and 2.2e-12–2.6e-12 on native R160), but could not reach its strict 1e-13 alpha-step threshold. The current result gives a structural reason not to interpret that iteration as a conventional well-conditioned root search: the target map is almost an identity.

## What this result does not decide

- It does not identify a faulty local centre stencil.
- It does not prove the continuum CMC equation or the chosen target average is incorrect.
- It does not justify dropping the outer cell from the target average; that would change the discrete definition without authority.
- It does not justify a fitted correction, damping, centre projection, or any physics/gauge patch.
- It does not replace the separate central regularity and radiation-admissibility gates.

## Next discriminating work

1. Complete the read-only CMC lapse centre-row audit on frozen t=0 R80/R160 profiles: report the solved alpha_i, regular-even-series fit, exact discrete elliptic row residual, independent K-RHS residual, and lapse-Laplacian bracket at the first cells and fixed radii.
2. Compare the actual centre row and boundary normalization row with the authoritative CMC derivation, explicitly deriving whether the global target must include the normalization cell and how Kdot is determined.
3. If the derivation reveals a rank-deficient discrete system, design a diagnostic-only alternative formulation and compare all constraint/admissibility witnesses before any production change.
4. Do not resume long production trajectories while centre-series and high-resolution radiation admissibility remain unresolved.

## Reproduction

Workflow .github/workflows/cmc-target-closure-rank-probe.yml runs diagnostics/cmc_target_closure_rank_probe.py on the pinned archived profiles. It samples nine fixed targets and reports alpha-profile affinity, the closure slope, interior PDE residuals, normalization-cell residual, and weighted residual identity. The artifact ZIP is preserved in the local library under ARTIFACTS_2026-10-12/cmc-target-closure-rank-probe.zip.
