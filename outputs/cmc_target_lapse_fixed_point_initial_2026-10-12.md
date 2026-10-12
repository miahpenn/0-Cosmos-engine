# Initial-slice CMC target–lapse fixed-point audit

**Date:** 2026-10-12  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Branch:** `diag/independent-metric-ricci-audit`  
**Run:** [38209333687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687) — successful  
**Artifact:** [11728232378 — cmc-target-lapse-fixed-point-initial](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687/artifacts/11728232378)  
**Diagnostic commit:** `04172fadbea50c4315054fec16bf41419e8e7f3b`  
**Artifact ZIP SHA-256:** `21018a09c7f30cf57915e4f1aebd3be62135efcd04941f2ece1e4efcc38f1f87`  
**Report schema:** `cmc_target_lapse_fixed_point_initial_v1`

## Executive result

The initial-slice audit passed its workflow and input-integrity gates. It does **not** report formal fixed-point convergence within 20 iterations because the stopping tolerance was set to an alpha RMS change of (10^{-13}), below the observed iteration-to-iteration numerical floor. However, the map steps are already extremely small: roughly (6.4\times10^{-11}) RMS for boundary-matched R80 and (2.2\text{–}2.6\times10^{-12}) for native R160. The target drift between the target evaluated on input alpha and the target evaluated on the returned alpha is only about (8.6\times10^{-14}) on R80 and (9\times10^{-16}) on R160. The increments do not show a reliable contracting sequence; their ratios fluctuate around one, consistent with the changes being at or near the numerical floor.

The fixed-target control reproduces each CMC lapse output **bit-for-bit** across three repeats on both domains (RMS differences exactly zero). Thus the elliptic solve is repeatable for a held target. The result is not evidence of a large unresolved target/lapse oscillation on the initial slice. It is evidence that this simple undamped Picard iteration cannot certify a (10^{-13}) fixed point because its residual step is dominated by very small numerical changes. Do not label the map “divergent”; do not label it “converged” under the declared stopping rule.

## 1. What was tested

The script `diagnostics/cmc_target_lapse_fixed_point_initial.py` restores the archived (t=0) states from the boundary-matched R80 profile and R160 reference profile. All non-lapse fields are frozen. It applies the map

[
\alpha_{n+1}=\mathrm{solve\_CMC\_lapse}\left(\dot K[\alpha_n]\right),
\qquad \dot K[\alpha_n]=\mathrm{target\_kdot}(\alpha_n).
]

There is no damping, clipping, parameter fit, altered production source, or time evolution. The iteration cap is 20; the declared stop threshold is alpha RMS change (le 10^{-13}).

- R80: (N=160, r_{\max}=80, \Delta r=0.5); outer lapse boundary fixed to the R160 archived value at (r=79.75), (\alpha_{\mathrm{out}}=1.0002206712250323).
- R160: (N=320, r_{\max}=160, \Delta r=0.5); native outer lapse boundary (\alpha_{\mathrm{out}}=1).
- The fixed-target control uses the archived R160 full-domain target, (\dot K=6.233524865732698\times10^{-5}), for three repeated lapse solves.

## 2. Input and admission checks

- Workflow status: **success**; all setup, exact-profile download/hash verification, diagnostic execution, and artifact upload steps passed.
- Both restored states are finite and have the expected cell counts and radial grids.
- R160 archived target reconstruction error: exactly (0).
- Initial native targets:
  - R80: (6.198841550010287\times10^{-5})
  - R160: (6.233524865732698\times10^{-5})
  - R80 − R160: (-3.4683315722410597\times10^{-7})
- Alpha remains positive on every recorded map step.
- The report is diagnostic-only; `production_physics_changed=false` and `production_kernel_source_changed=false`.

## 3. Iteration summary

| Quantity | R80 boundary-matched | R160 native |
|---|---:|---:|
| Map iterations | 20 | 20 |
| Declared status | max_iterations_reached | max_iterations_reached |
| First alpha-step RMS | (6.16437\times10^{-11}) | (2.44491\times10^{-12}) |
| Last alpha-step RMS | (6.39264\times10^{-11}) | (2.23689\times10^{-12}) |
| Last alpha-step max abs | See artifact JSON | See artifact JSON |
| Input-to-output target drift, typical | (-8.58\times10^{-14}) | about (-9\times10^{-16}) |
| Alpha minimum | (1.0002202830\) (last step) | (0.9999989453\) (last step) |
| Alpha maximum | (1.0144660371\) (last step) | (1.0148231649\) (last step) |
| Successive-step RMS ratio | around 1.00 | fluctuates around 1.00 |

R80’s step RMS is about 26 times R160’s, but both are extremely small in absolute terms. R80’s successive-step ratios span approximately 0.9969–1.0030 after the first few iterations; R160’s are noisier, approximately 0.91–1.10. Neither sequence provides a stable contraction factor that can be interpreted as a meaningful convergence rate. The small changes also drift slowly in one direction in the stored target values; see every iteration in the JSON artifact rather than inferring convergence from the last row alone.

### Correct reading of “max_iterations_reached”

The code uses a strict (10^{-13}) alpha-RMS stop threshold. It did not hit that threshold, so the status is correctly `max_iterations_reached`. This is **not** equivalent to saying the target–lapse map has a physically meaningful instability. The observed changes are already at the (10^{-11}) / (10^{-12}) level, and target drift is smaller still. A future solver should report a residual-based stopping criterion tied to the actual coupled equations and floating-point/discretization scale, not loosen this run’s tolerance post hoc and call that a proof.

## 4. Fixed-target control

For both R80 and R160, three repeated solves with the same fixed target return exactly the same alpha profile:

- R80 repeat 1 vs 2 RMS: (0); repeat 2 vs 3 RMS: (0).
- R160 repeat 1 vs 2 RMS: (0); repeat 2 vs 3 RMS: (0).
- All fixed-target outputs are positive.

This control confirms deterministic repeatability of the elliptic solve when its target is held fixed. It does not independently prove the full coupled target/lapse equation has a unique root.

## 5. Interpretation alongside the earlier diagnostics

The target/lapse domain-size investigation already found:

1. At (t=0), R80 and R160 non-lapse state fields on the shared core agree to roundoff, while alpha differs (max (3.57126\times10^{-4}), RMS (2.58354\times10^{-4})).
2. Replacing only R80’s alpha with R160 alpha in the frozen target functional nearly eliminates the initial target gap: (-3.46833\times10^{-7}) becomes (-5.88335\times10^{-14}). The reverse transplant barely changes R160’s full-domain target.
3. Full R160-target-history forcing into R80 closed 99.9941% of the final proper-time gap. Retaining every other target sample still closed 99.9813%, so the broad effect is robust to that tested temporal coarsening; the tiny remaining difference is interpolation-sensitive.
4. The R160 outer shell has about 87.5% of the proper-volume weight but shifts its own full-domain target relative to inner-only averaging by at most (7.1518\times10^{-12}) over saved snapshots. The observed R80/R160 target difference is therefore not explained by simple outer-shell averaging.

The present audit adds that the *initial* undamped target-to-lapse map produces only tiny per-iteration changes and a repeatable fixed-target solve. It does not yet establish whether the initial R80/R160 alpha difference comes from boundary normalization, a different branch/normalization of the coupled CMC solution, a discretization/closure convention, or a defect elsewhere in the production target/lapse construction. Do not jump from alpha–target sensitivity to a physical mechanism.

## 6. Limitations / unresolved checks

- This is a frozen (t=0) slice. It does not test evolution-stage target interpolation or feedback over the full trajectory.
- The two domains use deliberately different outer boundary normalizations for the native comparison (R80 is boundary-matched to R160 at (r=79.75); R160 uses alpha=1 at its own outer edge). The comparison is intentionally diagnostic, but boundary normalization remains a possible contributor to the different initial alpha profiles.
- The target map’s small iteration steps are not a substitute for evaluating a coupled nonlinear residual. The script records target drift and alpha step size, but no independently derived Newton/Jacobian solve of the full coupled discrete equations has yet been run.
- A tolerance of (10^{-13}) is below the empirical iteration noise for this map. Keep the declared result as “threshold not reached”; do not retcon the status.
- No production physics or kernel was edited, and no trajectory was launched by this diagnostic.

## 7. Recommended next steps

1. **Audit the exact initial lapse boundary and normalization convention** in the two archived runs and the production CMC solver. Reconstruct the same shared-core boundary value with a fixed target, then compare the resulting profiles against archived alpha at full precision.
2. **Compute a direct coupled residual** for the discrete system consisting of the CMC elliptic equation and (\dot K-\mathrm{target\_kdot}(\alpha)=0). Keep this diagnostic-only; first derive the exact discrete residual and Jacobian rather than introducing damping or fitted feedback.
3. **Run a fixed-target boundary cross-check**: solve R80 with both its archived boundary-matched alpha and the native alpha=1 boundary, using the same fixed target; solve R160 with the same shared-core target and compare only the common radial interval. This separates boundary normalization from target functional feedback.
4. Only after these checks, decide whether a coupled solver is needed. No more full trajectories are warranted for this question until the initial-slice algebra is understood.

## 8. Reproducibility and provenance

- Diagnostic source: `diagnostics/cmc_target_lapse_fixed_point_initial.py`
- Workflow: `.github/workflows/cmc-target-lapse-fixed-point-initial.yml`
- Commit: `04172fadbea50c4315054fec16bf41419e8e7f3b`
- Run: [38209333687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687)
- Artifact: [11728232378](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687/artifacts/11728232378)
- Artifact SHA-256: `21018a09c7f30cf57915e4f1aebd3be62135efcd04941f2ece1e4efcc38f1f87`
- Extracted JSON report SHA-256: `<report-sha256-to-fill>`
- Exact input ZIP checksums are recorded in the workflow and input-file hashes in the JSON report.
- `main` is not the target branch for this diagnostic. Keep all further work on `diag/independent-metric-ricci-audit`; production source remains untouched.
