# Build Status — 2026-10-12

## Objective

Get the unified V5.5 production machine running through a controlled, synchronized multi-resolution campaign, while leaving the locked physics unchanged until numerical admission is satisfied.

## Current production state graph

The canonical repository contains the production-facing state graph:

1. corrected archive COSMOS background equations;
2. local GEAR-03 S/D sector with GEAR-63/65 barrier;
3. scalar S/D/COSMOS evolution and stress projections;
4. conservative spherical dark-matter, baryon, and radiation matter sectors;
5. locked beta DM/scalar four-force exchange;
6. one total Einstein stress-energy assembly;
7. pinned reference-metric BSSN/PIRK numerical geometry;
8. synchronized geometry/scalar/matter state evolution;
9. derived proper time, H_eff, e-folds, trapping, and Misner-Sharp observables;
10. cycle/handoff ledgers without a synthetic lambda_cycle;
11. checkpoint and long-horizon multi-resolution campaign infrastructure.

## Physics locks and scope exclusions

- Bianchi-I shear remains outside the exact spherical local source.
- No D-to-matter identification law.
- No shell EOS/action or phenomenological interface source.
- No fitted coefficient, imposed clock conversion, bounce/reset/branch flip, lapse floor, ejection threshold, or physical stop.
- A finite campaign horizon is a computation boundary, not physics.
- Diagnostic-only forced histories are never promoted to the production equations.

## Validation state

Independent review reported the full test suite passed 79/79 on 2026-10-11, and a short production smoke campaign (`resolution=40`, `r_max=40`, `final_time=1.0`) completed without NaNs or assertion failures.

A separate strict-admission branch has now rerun the full repository suite and the same bounded full-machine smoke successfully:
- Branch: `fix/strict-campaign-admission-2026-10-12`
- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207394186
- Branch changes are not yet merged into `main`.

The branch adds a fail-closed campaign gate: all requested resolutions must be present, every resolution must report `completed`, and every trajectory must reach the requested final coordinate time within floating-point endpoint tolerance. The campaign summary records admission status and errors before a failure is raised; the campaign workflow is configured to upload available output even if admission fails. This prevents an early numerical stop from being mistaken for a successful campaign.

Current classification:

**IMPLEMENTED / UNIT-TESTED / SHORT-SMOKE-TESTED / NOT YET MULTI-RESOLUTION-CAMPAIGN-VALIDATED**

## Recent diagnostic evidence

### CMC target counterfactual

The diagnostic-only R80 evolution with the archived R160 CMC `Kdot` history completed through `t=12`:
- Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450
- Commit: `51045b8e0dfa8a00d167f7df78da9f07cfa772c5`
- 3,200 accepted steps, 49 saved slices, 12,801 lapse-target calls.
- Target delivery error was exactly zero; final central lapse was approximately 0.0172773 and final proper time approximately 6.90677991.

In the autonomous/baseline R80 comparison, final proper time differed from the archived R160 trajectory by about `+1.1743e-2`, with maximum saved-slice absolute difference about `1.5702e-2`. Jointly prescribing the archived R160 outer-lapse boundary history and full-domain CMC `Kdot` target reduced those differences to about `+6.30e-7` at the final time and `1.2192e-5` maximum over saved slices. Most core geometric and scalar profile discrepancies were also strongly reduced.

**Interpretation limit:** this is a full evolving, jointly forced gauge/boundary counterfactual, not a self-consistent R160 solution on an R80 domain. It demonstrates strong sensitivity to the combined outer-lapse boundary and global CMC target history; it does **not** isolate the `Kdot` contribution by itself, prove that the archived histories are the unique correct gauge, or establish a physical cause. The target history is interpolated. The separate Kdot interpolation audit passed 24 held-out times with maximum absolute error `5.29e-8`; the archived outer-lapse history has a larger reported held-out absolute interpolation discrepancy (about `1.5e-4` over a 0.5 time spacing), so stage-level target accuracy remains a limitation.

### N=640 signed-pole audit

The workflow at https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38183314088 was **cancelled at the 240-minute execution limit** after reaching approximately `t=8.375`, short of the requested `t=12`. Its partial artifact was uploaded and reproducible; this is not a completed through-`t=12` result. The early sign crossings disappeared by approximately `t=5.125`, while a later Hamiltonian residual ridge grew afterward. Term reconstruction closes to roundoff, but that is not independent constraint satisfaction.

### Other current audits

- Successful CMC Kdot interpolation audit: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868
- Strong-D CMC lapse Jacobian audit: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205542397
- Strong-D CMC lapse time series: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205709871
- R80/R160 CMC target counterfactual comparison: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205934746
- Bond sweep is diagnostic evidence only; it does not by itself explain the full evolving constraint discrepancy.

## Open validation gates

1. Run a controlled production campaign through `t=12` for `N=40,60,80`, `r_max=80`, `CFL=0.03`, radiation enabled, and locked physics. Preserve all partial ledgers if a resolution fails.
2. Require strict admission: all requested resolutions present, all runs completed, exact campaign horizon reached within endpoint tolerance, and no non-finite state.
3. Review per-resolution Hamiltonian, momentum, connection, determinant, trapped-root, Misner-Sharp current, proper-time, H_eff, matter-exchange, and boundary-flux histories. Completion alone is not physics admission.
4. Compare the same-domain resolution sequence before drawing conclusions from cross-domain comparisons.
5. Continue the N=320 reconstruction audit using the frozen `1e-8` per-field relative criterion. Do not widen tolerances to force a pass. PR #2 remains draft because earlier host/runtime-dependent outer-Hamiltonian discrepancies and preregistration-source provenance remain unresolved.
6. Only after the intermediate campaign is admitted should the longer multi-resolution horizon be run.

A numerical/gauge failure is not a physical interpretation. No production physics change is authorized by unit tests, a short smoke run, or a target-forced counterfactual alone.

## Reproducibility note

Record commit SHA, pinned vendor revision, settings, runtime, checkpoint hashes, and all failing diagnostics for every admitted campaign. Preserve preregistered tolerances and report host/runtime variance explicitly.
