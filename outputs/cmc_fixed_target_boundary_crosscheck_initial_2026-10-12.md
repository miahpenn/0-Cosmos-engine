# Fixed-target CMC boundary/target cross-check — initial slice

**Date:** 2026-10-12  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Branch:** `diag/independent-metric-ricci-audit`  
**Workflow:** [38293992203 — success](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38293992203)  
**Artifact:** [11767094230 — cmc-fixed-target-boundary-crosscheck-initial](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38293992203/artifacts/11767094230)  
**Diagnostic script commit:** `cc664bdfcc408f79f8605f1ecf35f58c237061d8`  
**Workflow commit:** `a5ffbe468024cc76b20b47b768dd2e73c54dbcc8`  
**Artifact ZIP SHA-256:** `dbf718c489a6813e16cf4c208c1f055f2b6bedd4810203b3b7ac9565dfbf871e`  
**Extracted report SHA-256:** `013c349a5cb6296cccf18e135f39449716e11d8929390b0bd7f7a76289075f6a`

## Executive result

This frozen `t=0` test independently varies (a) the fixed global target and (b) the R80 outer lapse boundary. It advances no time and changes no production code.

**The archived shared-core lapse difference is reproduced by changing the target while keeping the R80 boundary matched to R160.** Solving R80 with the archived R160 target and boundary alpha `1.0002206712250323` matches the R160 fixed-target lapse on the shared 160 cells to RMS `2.4499692141e-12`, maximum absolute difference `2.6325608360e-12`. With the same boundary but R80’s native target, the lapse difference is RMS `2.5835384621e-4`, maximum `3.5712617679e-4`—essentially the archived mismatch (RMS `2.5835378535e-4`, max `3.5712609202e-4`).

This is strong evidence that the **different target values** account for the initial shared-core lapse discrepancy in this matched-boundary comparison. It is not a claim that the target difference has a physical cause or that the coupled nonlinear system has been proven to possess two distinct roots.

## Method and admission

The workflow restored exact archived `t=0` full snapshots, froze all non-lapse fields, and ran the production CMC lapse solve under fixed target values. R80 was solved using both outer boundary choices: (1) R160 alpha at the shared-core edge `r=79.75`, `1.0002206712250323`; and (2) native alpha `1.0`. Two target values were tested: R80 native target and archived R160 target. R160 itself was solved with the archived target and native outer alpha 1.0. All alpha outputs were finite and positive; archived R160 target reconstruction differed by only `1.675e-17`.

Targets:
- R80 native: `6.198841550014436e-5`
- R160 archived: `6.233524865732698e-5`
- R80 − R160: `-3.4683315716587076e-7`

## Numerical comparisons

| Comparison | RMS alpha difference | Maximum absolute difference |
|---|---:|---:|
| Archived R80 vs archived R160 shared core | `2.5835378535e-4` | `3.5712609202e-4` |
| R80, R160 target + matched boundary vs R160 fixed-target solution | `2.4499692141e-12` | `2.6325608360e-12` |
| R80, R160 target + alpha=1 boundary vs R160 shared-core solution | `2.1186655855e-4` | `2.206712223998e-4` |
| R80 boundary effect at fixed R160 target | `2.1186656100e-4` | `2.2067122503e-4` |
| R80 target effect with matched boundary | `2.5835384621e-4` | `3.5712617679e-4` |
| R80 target effect with alpha=1 boundary | `2.5835384622e-4` | `3.5712617680e-4` |

The first row is the archived discrepancy. The second row shows that matching the target as well as the boundary removes it down to roundoff-level profile differences. The fifth row shows that switching only the target back to R80’s native value recreates it. Boundary normalization matters when changed, but it does **not** explain the archived R80/R160 discrepancy by itself, because the archived R80 case was already boundary-matched to R160 at the shared edge.

## Relation to the undamped map audit

The prior 20-iteration map had target drift only around `8.6e-14` for R80 and `9e-16` for R160, with alpha steps around `6.4e-11` and `2.2e-12`, respectively. Those are tiny relative to the `3.47e-7` target gap. Taken together, these observations are compatible with each archived target/lapse pair being nearly stationary under its own target, despite different target values and lapse profiles. However, the map did not meet its strict `1e-13` alpha-step threshold, and no direct coupled Jacobian/residual analysis has yet proved multiplicity or uniqueness. Use wording such as “two distinct near-stationary candidates” until the coupled system is analyzed.

## Relation to pointwise residual audit

The independent pointwise K-RHS/lapse residual audit (run [38210163221](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221), artifact [11729296911](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221/artifacts/11729296911)) found PDE-row residuals around `1e-15`, while the outermost lapse-normalization cell has a residual around `-4.5882e-12`. The weighted pointwise residual mean matches the scalar closure `H(Kdot)` to about `1.6–4.0e-21`. This points to a small boundary-row/global-average seam worth auditing against the authoritative CMC derivation; it does not justify changing the production target definition without that derivation.

## Interpretation limits

- Snapshot-only; no evolution, no changed production source, no physics claim.
- The test separates boundary normalization from target choice for the archived initial states, but does not establish why the target functional has distinct values on the two lapse profiles.
- The fixed-target solve is repeatable; that is not a proof of uniqueness of the coupled equation.
- Do not introduce damping, fitted feedback, clamps, floors, or a manufactured law to force one branch.

## Next step

Construct an independently evaluated residual/Jacobian for the coupled discrete system: the CMC elliptic equation on PDE rows plus `Kdot - target_kdot(alpha) = 0`, with the outer normalization cell treated explicitly and not silently included/excluded. Compare the R80-native and R160-target near-stationary candidates under the **same** shared-core boundary condition. Separately audit the pointwise outer-cell contribution against the authoritative CMC derivation. Keep all work diagnostic-only and do not launch another expensive trajectory yet.


## Target-closure flatness across three fixed-target trials

The pointwise residual artifact [11729296911](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221/artifacts/11729296911) provides three fixed-target R80 solves with the same frozen non-lapse state and matched outer boundary. The trial Kdot values span the entire native-R80 to archived-R160 target gap:

| Trial | Kdot trial | H = target_kdot(alpha_solved) − Kdot |
|---|---:|---:|
| Native R80 | (6.198841550010287\times10^{-5}) | (-8.581045487891792\times10^{-14}) |
| Midpoint | (6.216183207871493\times10^{-5}) | (-8.582846618750833\times10^{-14}) |
| Archived R160 | (6.233524865732698\times10^{-5}) | (-8.580558952166889\times10^{-14}) |

Across a Kdot span of (3.46833157\times10^{-7}), H varies by only (2.28767\times10^{-17}); the endpoint finite-difference slope is about (1.40\times10^{-11}). For every trial, H matches the proper-volume weighted mean of the pointwise K-RHS residual to (1.6\text{–}4.0\times10^{-21}). PDE-row residuals are at ~(10^{-15}); the outermost normalization cell carries a residual of about (-4.5882\times10^{-12}).

**Interpretation:** the coupled target-closure residual is nearly flat with respect to trial Kdot over this interval after solving the lapse equation. This is strong evidence that the global target equation supplies little independent selection of Kdot beyond the elliptic PDE rows and the boundary closure, and it is consistent with distinct near-stationary target/lapse candidates. It is not yet a proof of a continuum of exact roots or of a production defect. The next diagnostic must derive the discrete residual/Jacobian and treat the normalization cell explicitly, then check the authoritative CMC derivation before any production change.


## Structural reading of the target/lapse closure

In `engine/cmc_gauge.py`, `target_kdot` computes a proper-volume weighted average of the production raw K RHS `l2["K"] + l3["K"]`. `solve_cmc_lapse` solves the CMC elliptic equation with the chosen Kdot on the first (N-1) cells and uses the last cell only to impose lapse normalization. The pointwise audit confirms PDE-row raw K RHS agrees with trial Kdot to about (10^{-15}); only the normalization cell remains at about (4.59\times10^{-12}) residual.

This makes the near-flat closure understandable: after solving the lapse for a trial Kdot, re-averaging the same K RHS almost necessarily returns that Kdot. The global projection may therefore be algebraically dependent on the PDE rows, apart from the normalization-cell contribution, rather than an independent equation selecting Kdot. Across the tested (3.4683\times10^{-7}) target interval, closure H varies only (2.29\times10^{-17}), slope approximately (1.40\times10^{-11}).

**This is a structural hypothesis supported by the discrete residuals, not a declaration that the production gauge is wrong.** The authoritative CMC derivation must decide whether Kdot is meant to be fixed by a separate condition or whether the current projection is intentionally a gauge convention. Next derive the discrete coupled system explicitly and determine whether the target equation adds an independent constraint. Do not patch the target average or remove the normalization cell by intuition.
