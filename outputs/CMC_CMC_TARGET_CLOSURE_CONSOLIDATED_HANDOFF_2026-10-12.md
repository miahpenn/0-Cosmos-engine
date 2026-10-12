# CMC target–lapse closure investigation — consolidated handoff

**Prepared:** 2026-10-12 (UTC)  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Active branch:** `diag/independent-metric-ricci-audit`  
**Scope:** diagnostic-only analysis of initial-slice CMC target/lapse coupling. No production equations, kernel, or physics were changed by these audits.

## Executive result

The diagnostics provide strong numerical evidence that, for the tested frozen initial slice and discretization, the global CMC target condition is nearly redundant with the elliptic lapse PDE rows plus the outer lapse-normalization row. This is a structural hypothesis supported by several independent numerical probes; it is **not yet a proof of exact algebraic dependence**, not proof of multiple continuum roots, and not by itself evidence of a production physics defect.

The next required step is to derive and evaluate the coupled discrete residual/Jacobian (or Schur complement) for the complete system, then compare it to the authoritative CMC derivation. Do not patch production or spend another expensive evolution run until that derivation is resolved.

## 1. Nine-point target-closure rank probe

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294955148
- Artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294955148/artifacts/11768056222
- Report: `outputs/cmc_target_closure_rank_probe_2026-10-12.md`
- Source commit: `1985daaee34e89fe1da68968f0518249f79cd3bc`
- Workflow commit: `be7dac470eca22a1267dbbeac7776719afbd88d6`
- Artifact ZIP SHA-256: `64cc1666a2010d73b462cac48c78706ac48be1a1d02fea7de343bd09172c4e24`
- Extracted report JSON SHA-256: `8ce2d07e90f0a760fd0683efb72e0af6530893b13f03776dc1ed67efb240b1ab`

Setup/results:
- Frozen R80 initial state; non-lapse fields fixed.
- R80 outer alpha matched to R160 at `r=79.75`, with `alpha_out = 1.0002206712250323`.
- Nine trial Kdot values ranged from `6.1294749186e-5` to `6.3028914972e-5`, a span wider than the native R80/R160 target gap.
- Solved target output versus input slope: `0.999999999757`.
- Closure `H = target_kdot(alpha_solved) - Kdot_trial` spanned only `4.378821524125831e-16`; H values were approximately `-8.6051e-14` to `-8.5613e-14`.
- Linear-fit slope `dH/dKdot = -2.426863537036872e-10`.
- Lapse profiles were affine in Kdot to RMS error `4.4823e-15`.
- PDE-row K RHS residual maximum `3.7093e-15`; outer normalization-cell residual maximum `4.5993e-12`.
- H matched the weighted mean pointwise residual to `1.14e-20`.

Interpretation: over the sampled range, the target closure barely changes while trial Kdot varies substantially, consistent with the global projection reproducing the imposed Kdot through the PDE and boundary closure.

## 2. Fixed-target boundary/target cross-check

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38293992203
- Artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38293992203/artifacts/11767094230
- Report: `outputs/cmc_fixed_target_boundary_crosscheck_initial_2026-10-12.md`
- Diagnostic script commit: `cc664bdfcc408f79f8605f1ecf35f58c237061d8`
- Workflow commit: `a5ffbe468024cc76b20b47b768dd2e73c54dbcc8`
- Artifact ZIP SHA-256: `dbf718c489a6813e16cf4c208c1f055f2b6bedd4810203b3b7ac9565dfbf871e`
- Report JSON SHA-256: `013c349a5cb6296cccf18e135f39449716e11d8929390b0bd7f7a76289075f6a`

Results:
- R80 with R160 fixed target and matched boundary reproduces the R160 lapse on the shared core: RMS `2.4499692141e-12`, maximum `2.6325608360e-12`.
- Keeping that same boundary but switching only to the native R80 target reproduces the archived lapse mismatch: RMS `2.5835384621e-4`, maximum `3.5712617679e-4` (archived mismatch RMS `2.5835378535e-4`, maximum `3.5712609202e-4`).
- Changing boundary normalization separately has an effect of RMS `2.11866561e-4`, maximum `2.20671225e-4`.
- Archived R80 was already boundary-matched to R160, so boundary normalization alone does not explain the archived shared-core mismatch.
- Native R80 versus R160 target gap: `-3.4683315716587076e-7`.

Interpretation: the initial shared-core lapse mismatch is directly explained by the target difference in this controlled comparison. This does not yet explain why the targets differ.

## 3. Initial-slice target–lapse fixed-point audit

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687
- Artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687/artifacts/11728232378
- Report: `outputs/cmc_target_lapse_fixed_point_initial_2026-10-12.md`
- Diagnostic commit: `04172fadbea50c4315054fec16bf41419e8e7f3b`
- Artifact ZIP SHA-256: `21018a09c7f30cf57915e4f1aebd3be62135efcd04941f2ece1e4efcc38f1f87`
- Report JSON SHA-256: `b951aae61ba16b3c1850d26534a28ce40a6ee09b69178c45f847b74ebc99c7fb`

Results:
- Undamped frozen-state map `alpha -> solve_CMC_lapse(target_kdot(alpha))` ran 20 iterations on R80 and R160.
- Both reported `max_iterations_reached` because the requested alpha RMS-step threshold was `1e-13`.
- Final alpha step: R80 `6.39264e-11`; R160 `2.23689e-12`.
- Target drift: about `-8.58e-14` on R80 and `-9e-16` on R160.
- Ratios fluctuated around 1. Do not describe this as divergent; do not describe it as formally converged under the stated tolerance.
- Fixed-target control with three repeated solves returned identical arrays (RMS difference exactly zero).
- Initial target R80 `6.198841550010287e-5`; R160 `6.233524865732698e-5`; gap `-3.4683315722410597e-7`.

## 4. Independent pointwise residual audit

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221
- Artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221/artifacts/11729296911
- Artifact ZIP SHA-256: `cfb1b7a71d4267c15f2b07dd123ecac2f83195b9be012a4e6a84cf1f86966c33`

Results:
- PDE-row residuals are around `1e-15`.
- Outer normalization cell residual is around `-4.5882e-12`.
- H matches weighted mean pointwise residual to about `1.6e-21`–`4e-21`.
- Three trial Kdot values (native R80, midpoint, archived R160) yield H around `-8.58e-14`; range only `2.29e-17` over a target gap of `3.468e-7`; endpoint slope about `1.4e-11`.
- `engine/cmc_gauge.py`: `target_kdot` computes a proper-volume-weighted average of raw K RHS `l2["K"] + l3["K"]`; `solve_cmc_lapse` solves PDE rows on cells `0..N-2`, with the last cell enforcing lapse normalization.

Interpretation: this implementation structure makes the near-redundancy plausible. The authoritative derivation still must decide whether this is expected by construction or signals a discrete closure issue.

## 5. Related trajectory/counterfactual findings

- At `t=12`, `tau80 - tau160 = +0.01062109`.
- Integrated central lapse difference is `0.01063961`; the ~0.18% difference is a bookkeeping identity associated with `tau = integral(alpha dt)`, not a causal explanation.
- Frozen-slice lapse re-solve correction integrates to only about `2.5e-6` of the original proper-time discrepancy.
- Frozen-state target swap substantially improves radial profile error but does not reproduce the center-lapse history.
- Forcing the full R160 target history into R80 closed 99.9941% of the final proper-time gap; every-other-snapshot forcing closed 99.9813%, with interpolation residual still relevant.
- R160 outer shell contributes about 87.5% of proper-volume weight, but full versus inner-only target differs by at most `7.1518e-12`; simple outer-shell averaging does not explain the target gap.
- A separate evolving target-history counterfactual workflow `38206154450` was previously noted as in progress; its final status/result must be checked before relying on it.

## 6. Interpretation and constraints

- The target gap is the direct cause of the controlled initial shared-core lapse mismatch.
- The target/lapse closure looks nearly rank-deficient or algebraically dependent under the tested discrete construction, but exact dependence remains unproven.
- Do not infer a bounce, turnaround, or new physics from these diagnostics.
- Do not introduce arbitrary damping, clamps, floors, fitted feedback coefficients, hand-tuned target values, or other manufactured closure rules.
- No production patch is justified by the present evidence.

## 7. Required next test

Build the full coupled discrete residual system:
1. Elliptic PDE rows `F_i(alpha, Kdot) = 0`, `i = 0..N-2`.
2. Outer lapse-normalization row `B(alpha) = alpha[N-1] - alpha_out = 0`.
3. Global target condition `H(alpha, Kdot) = target_kdot(alpha) - Kdot = 0`.

Compute the Jacobian/rank or Schur complement, keeping row/unknown ordering explicit. Compare with the authoritative CMC derivation and document whether H contributes an independent equation after the PDE and boundary rows. Include numerical conditioning and grid dependence. Keep it diagnostic-only; no expensive evolution or production change until the equation count and intended gauge closure are settled.

## 8. Repository and Library status

The four focused diagnostic reports are present in the active diagnostic branch:
- `outputs/cmc_target_closure_rank_probe_2026-10-12.md`
- `outputs/cmc_fixed_target_boundary_crosscheck_initial_2026-10-12.md`
- `outputs/cmc_target_lapse_fixed_point_initial_2026-10-12.md`
- `outputs/cmc_strongD_cmc_lapse_handoff_2026-10-12.md`

The Library already contains earlier fixed-point and pointwise-residual artifacts and the rank-probe ZIP. The newest boundary-cross-check archive/report and standalone rank/boundary Markdown reports had not been verified as uploaded to Library at the time of this handoff. Uploads must be verified individually; earlier quota failures mean do not assume they were saved just because they exist in the repository or working directory.

## 9. Workflow/branch housekeeping note

A workflow file was briefly created on `main` because a branch argument was omitted, then deleted and recreated on the diagnostic branch. The `main` tree hash was verified identical before and after, so there was no net file-tree change and no production code change, though the `main` commit history includes an add/delete pair.
