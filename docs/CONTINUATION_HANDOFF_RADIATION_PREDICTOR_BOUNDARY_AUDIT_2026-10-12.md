# Continuation handoff — radiation predictor / P-minus boundary audit

**Date:** 12 October 2026  
**Status:** Diagnostic research only; production admission remains open.  
**Canonical rule:** preserve the exact radiation conservative cone \(E\ge |S|\). No clamp, relaxed cone, artificial floor, fitted coefficient, or manufactured physics law.

## Main result

The old high-resolution early radiation RHS failure was a false-positive guard check: `species_rhs()` used `evolve_species(dt=1.0)` to extract a semi-discrete RHS, and the cone guard incorrectly treated that artificial unit-time increment as a production update. The isolated rate-only correction fixed that software defect; the full repository suite passed.

The separate, late radiation predictor failure remains. In the corrected-rate fixed-`dt=0.03` trace, 1,566 accepted steps complete, ending at \(t=46.98000000000102\), central proper time \(\tau=22.016046893100416\). The attempted predictor at \(t=47.01000000000102\) crosses the exact radiation cone at the outermost cell \(i=79, r=79.5\): \(|S|/E=1.0000573283549097\), densitized margin \(C=-1.6997416818\times10^{-8}\).

At that failing same-cell stage, flux plus geometric source *increase* the cone margin on the accepted metric from \(4.7291\times10^{-7}\) to \(2.2198\times10^{-6}\). The predictor metric then contributes \(-2.2368\times10^{-6}\), turning the final margin negative. The explicit BSSN updates at zero shift are \(\dot a=-2\alpha a A_a\) and \(\dot X=\alpha X K/3\); the combined Euler metric update increases \(\gamma^{rr}=X^2/a\) by about 0.0223233. Updating either \(a\) or \(X\) alone leaves the same conservative state admissible; their combined update crosses the cone. The centre-regularity projection changes nothing at that outermost failing cell.

The last ten completed steps in the corrected-rate artifact show the outer inverse radial metric rising 1.278974 → 1.462894 (+14.38%) while the accepted cone margin declines 3.7051e-6 → 4.7291e-7 (−87.2%). Every step's matter flux/source increases the margin on the accepted metric, while the geometry stage removes more than the matter update adds. Neighboring cells tolerate even larger fractional metric changes because their radiation states are farther from null.

## Primary P-minus boundary correction: new admitted cross-agent ledger

Workflow [38311614754](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38311614754) completed successfully with the full test suite, source/vendor pin checks, and stage-ledger admission. Artifact [11774209967](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38311614754/artifacts/11774209967) SHA-256:

`989b51ffec4c9153d228ce3e8f80ef38873337985b3a4271ef904d2c1ff6b30b`

Internal `stage_ledger.json` SHA-256: `2c56e08e0331b7b06dcf7017f53ffb47413eb647658c6a17bd075a8f902a7bb6`. The stage summary SHA-256 is `1f2a749a6b8cb131a0d6d626d0f9f40380cad15501813911c3d744be306284fa`. It admits 100 captured stages for cells 75–79, has no validation errors, exact `primary_l3_total - primary_l3_vacuum - primary_l3_matter_delta` closure, and completed-to-next-accepted geometry handoff error at most 2.22e−16.

**Source caveat:** the ledger is frozen to historical source `9eed1cf95c52f001c2b87afbbe6da37faa56d50e` and vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. It is *not* a matched A/B against corrected source `1c63090d4f42396c899421a2379f70320689d5ad`; its near-identical late event supports stage-order interpretation but its exact values must remain source-qualified.

Reconstruction of the actual primary-predictor PIRK formula at outer cell 79 gives a direct P-minus boundary correction to \(A_a[-1]\), growing over the last ten completed steps from +2.53865e−5 at t=46.68 to +5.75540e−5 at t=46.95. Across the last five captured cells, the correction is exactly zero at cells 75–78 and applies at cell 79 only, consistent with the boundary routine assigning only the outermost `Aa[-1]`. At t=46.95, this is about 0.0365% of post-boundary \(A_a\), and about 1.6% of the unconditioned primary-predictor \(A_a\) increment.

The boundary assignment also shifts the geometric momentum-constraint term at that same primary predictor. For the pinned fourth-order outer derivative, the exact same-slice jump is

\[
\Delta M = \Delta A_a\left[\frac{25}{12\,\Delta r}+6\chi+\frac32\left(\frac2r+\frac{D_r b}{b}\right)\right],\quad \chi=-\frac12\frac{D_r X}{X}.
\]

At t=46.95, \(\Delta M=+1.22523\times10^{-4}\); across the last ten completed steps this grows monotonically from +5.39832e−5 to +1.22523e−4. The formula and directly recomputed geometric-M difference agree within 4.3e−17 for all 99 completed primary-predictor stages. Since this boundary operation changes only \(A_a[-1]\), not the metric or matter state, this is also the signed jump in the full production M residual at fixed same-slice matter. The *absolute* full M residual change cannot be assessed from this ledger alone because it does not contain a complete total-\(j\) profile in the needed stage form. Therefore, this proves a growing, non-negligible boundary-conditioned contribution to the outer momentum-constraint equation—not that the target is wrong and not that it caused the radiation crossing.

## Boundary interpretation guardrail

The exact frozen P-minus target is inherited from the pinned numerical BSSN operator, not prescribed independently by the archived GEAR physical law. Its module says it is not a complete constraint-preserving boundary. Near-roundoff post-boundary residual demonstrates algebraic closure of the chosen spatial target, not independent validation of the continuum incoming datum, dynamic characteristic evolution, all constraint modes, or compatibility with the production CMC/zero-shift gauge (`alpha[-1]=1`).

An older expanded-domain run (N=120, Rmax=120, historical source) achieved post-boundary target/formula residuals near 1e−17 yet corrected outer \(A_a\) by +1.1851e−4 at t≈60.99; the eventual radiation cone failure moved to r=119.5 but was not eliminated. This neither convicts nor exonerates the boundary. The corrected-rate failing predictor throws during total stress-energy construction for the first predictor CMC lapse, before the current predictor's P-minus boundary call. Prior boundary-conditioned accepted \(A_a\) remains a plausible pathway into the next \(\dot a\), and the new historical ledger confirms a growing contribution, but attribution to the corrected-rate trajectory remains unmeasured.

## Relevant repo documents and coordination

- Main radiation trace report: [RADIATION_PREDICTOR_TRACE_AFTER_RATE_FIX_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/RADIATION_PREDICTOR_TRACE_AFTER_RATE_FIX_2026-10-12.md), current commit `96661619cec681c2e4dd35261414401ec60eb32e`.
- Combined CMC/radiation/Misner–Sharp one-step A/B: [COMBINED_CMC_RADIATION_MS_AB_ARTIFACT_FORENSIC_READ_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/COMBINED_CMC_RADIATION_MS_AB_ARTIFACT_FORENSIC_READ_2026-10-12.md), commit `048a31147b8cd17d4fd655fd445a02e18bd7279a`.
- Cross-agent coordination thread: [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5), latest comment `6133175935`.
- Corrected-rate trace [run 38310216296](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310216296), artifact `11774592255`, ZIP SHA-256 `d2f0cdc037701f7b03cb76b26b5ccefa935a77fadfbd62b23969328ea7bb7c53`.
- Historical-source metric ledger [run 38311614754](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38311614754), artifact `11774209967`.

No `main`, production physics, or canonical GEAR archive changes were made. No physical interpretation is admitted: no bounce/turnaround/cycle inference; no full-trajectory or production pass.

## Next steps — offline first, no new trajectory yet

1. Derive the frozen spatial P-minus target's compatibility (or incompatibility) with the actual production CMC and zero-shift evolution equations. Keep algebraic target closure, characteristic transport, momentum/Hamiltonian/connection constraints, and radiation cone margins as separate checks.
2. Reconcile the saved Aa/K source-stage values and boundary-induced M jump with the full total-stress-energy constraint budget. The L3 matter-minus-vacuum correction at the outer terminal source in the historical ledger is small relative to the L3 totals (Aa −1.404e−4 vs +4.234e−2; K +2.519e−4 vs +6.119e−2), but this does not measure all matter influence because the geometry and CMC solve have already responded to total Tμν.
3. If a matched corrected-source boundary attribution is still essential after this derivation, make the instrumentation record the pre/post `Aa[-1]` and total M terms directly at both boundary call sites, fail-closed pin all runtime/physics source blobs, and preregister the decision metric before any long run. Do not launch a second expensive trajectory simply to repeat the existing failure.
4. Preserve exact `E>=|S|`; no clamping, weakened cone, accepted-metric recovery, or projection tweak. Keep all production changes separate from the diagnostic branches until a specific defective layer is identified and a matched regression shows improved constraint/regularity budgets.
