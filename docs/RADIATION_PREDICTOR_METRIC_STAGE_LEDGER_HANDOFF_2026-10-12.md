# Radiation Predictor Metric / Regularity-Projection Attribution — SI Handoff

**Date:** 12 October 2026  
**Purpose:** Preserve exact evidence and next admission gate for the late-time outer-cell radiation cone crossing.  
**Status at report creation:** Existing late-time failure is reproduced and algebraically attributed to the predictor-metric change as a whole; raw Euler metric versus algebraic projection remains the remaining split. The dedicated stage-ledger workflow is running; its final artifact has not yet been available at report time.  
**Diagnostic branch:** `diag/radiation-predictor-metric-stage-ledger-2026-10-12`  
**Frozen physical-source baseline:** `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Diagnostic workflow commit:** `e9eb4b2b6625008f5837532b0862888ba4751ded`  
**Stage-ledger run:** [38308235431](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38308235431)

## Executive summary

Keep this late-time outer-cell failure separate from the high-resolution first-cell false failure fixed by the matter RHS rate-extraction correction.

The saved N=80 trace reached 1,566 accepted steps at coordinate time (t=46.98) and central proper time (	au=22.016115393070738). The next explicit predictor, at attempted time (t=47.01), failed during the first predictor CMC lapse input evaluation at outer cell (i=79, r=79.5). The updated conservative radiation state is admissible when evaluated on the accepted metric, but inadmissible on the predictor metric.

The saved ordered cone budget proves the flux and geometric-source increments are not the proximate cause in that one-step budget: evaluated on the accepted metric they improve the cone margin. The change from accepted metric to predictor metric then subtracts (2.2340303316\times10^{-6}) from the densitized cone margin and makes it slightly negative. This attributes the terminal crossing to the **combined predictor metric update**, not yet to the raw geometry RHS versus the algebraic regularity projection.

The pinned kernel explicitly calls the reference `enforce_algebraic_regularity` projection after the explicit Euler geometry predictor. The original trace does not record the pre-projection geometry, so its last metric budget cannot isolate that projection's share. The dedicated stage-ledger workflow is intended to record: accepted geometry/RHS → unprojected Euler predictor → post-projection metric → first predictor CMC inputs. Until its artifact validates, no claim should be made about the relative sizes of those terms.

## 1. Existing reproduced event and provenance

- Original trace workflow: [38068938182](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938182)
- Original artifact ID: `11676486591`
- Original artifact ZIP SHA-256: `bef03745ad4ea1ab62d8a42de9b957c57210571d933cf2728c96a6196b970c43`
- Internal trace SHA-256: `642a82b8b97ddc849b08f70b4379672b08e81f842d3ae8bdeaaa1d5374afd8a3`
- Trace schema: `0star_radiation_predictor_trace_v1`
- Configuration: N=80; (r_{max}=80); (dr=1); CFL=0.03; (dt=0.03); amplitude=0.01; width=7; (D)-amplitude (10^{-10}); radiation enabled.
- Source baseline: `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`; trace instrumentation `527b44931b059f4d118a86c898b91d79ca4d96f2); vendor submodule `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

Failure facts:
- Last accepted (t=46.98000000000102), (	au=22.016115393070738), step index 1566.
- Attempted (t=47.01000000000102), (dt=0.03), stage `explicit_predictor_first_CMC_lapse`.
- Cell (i=79), (r=79.5), (alpha=1), (eta=0).
- (U_E=2.965253923984548\times10^{-4}), (U_r=2.433162925154048\times10^{-4}).
- Accepted metric: (a=0.6089553401215667), (X=0.9438956524107626), (gamma^{rr}=1.4630613181946635).
- Predictor metric recorded by the old trace: (a=0.6032035002114748), (X=0.946558319470313), (gamma^{rr}=1.4853571835116133).
- Predictor cone margin (C=-1.6900205171508503\times10^{-8}); (|S|/E=1.0000569941246342). This is a small but real failure under the unchanged cone criterion, not a diagnostic sign convention error.

## 2. Ordered cone budget on the terminal step

Use the actual invariant:
[
C=U_E-\sqrt{X^2/a}\,|U_r|,qquad
\frac{|S|}{E}=\frac{\sqrt{X^2/a}\,|U_r|}{U_E}.
]

| Terminal cell budget stage | Cone margin (C) | Interpretation |
|---|---:|---|
| Accepted radiation + accepted metric | (+4.6971156703\times10^{-7}) | Inside cone, close to boundary |
| Flux-only update, evaluated on accepted metric | (+1.4812240207\times10^{-6}) | Margin increases by (+1.0115124536\times10^{-6}) |
| Geometric-source update, still on accepted metric | (+2.2171301265\times10^{-6}) | Further increase (+7.3590610578\times10^{-7}) |
| Same predictor conservative state, predictor metric | (-1.6900205172\times10^{-8}) | Crossing caused by metric change in this split |

With the predictor conservative variables held fixed:
- On the accepted metric, (|S|/E=0.9925229670601998).
- On the predictor metric, (|S|/E=1.0000569941246342).
- The change in (gamma^{rr}) is (+1.5239187203\%); its square-root factor rises by (+0.7590783603\%).
- Over the same step, (U_E) rises by (1.8453523\%), (U_r) rises by (1.2471916\%), and (|U_r|/U_E) falls by (0.5873225\%).

So the conservative matter update, viewed on the accepted metric, improves the margin; the predictor metric change is large enough to overwhelm that improvement. This is not proof that the physical geometry law is wrong: a predictor stage can leave the primitive domain even when the accepted state and the finite-step conservative source update are individually admissible.

## 3. Why the projection split matters

In the frozen production kernel `9eed1cf...`, `V55ProductionKernel.step` first regularizes/copies the accepted geometry, calculates the explicit geometry RHS, then constructs
[
g_{\mathrm{Euler}}=g_0+\Delta t\,R_{\mathrm{explicit}}(g_0),
]
and applies `_enforce_center_regularity(grid, g_explicit1)` before the first predictor CMC lapse solve.

The pinned vendor's `vacuum_bssn.enforce_algebraic_regularity` enforces (a b^2=1) over all cells by replacing
[
a\to(a/b)^{2/3},\qquad b\to(a/b)^{-1/3},
]
and also enforces the spherical-centre series identities at the first cell. It therefore alters the outer-cell (a,b) pair if the raw Euler predictor departs from the algebraic determinant constraint; since (gamma^{rr}=X^2/a), any change in (a) alters the cone diagnostic even though (X) and the radiation conserved state are unchanged by that projection.

At the terminal outer cell, the accepted and recorded projected predictor states both have (ab^2\approx1). The old trace did not save raw pre-projection (a,b), so this fact alone cannot say whether projection was negligible or contributed materially to the final (gamma^{rr}) shift.

## 4. New diagnostic-only stage-ledger run

The recorder has been wired on this isolated branch through [workflow run 38308235431](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38308235431). It checks the frozen source blobs and pinned vendor before running; then runs the complete repository tests and the bounded trace. The workflow records, over the outer five cells:
1. accepted geometry and the explicit RHS;
2. unprojected Euler-predictor values reconstructed as accepted state plus (dt\times RHS);
3. the actual post-`enforce_center_regularity` geometry;
4. the actual geometry/radiation values handed to the first predictor CMC solve.

It verifies stage ordering, reconstructed (\gamma^{rr}=X^2/a), the Euler update, projection determinant, cone values and exact correspondence between projected values and first CMC inputs. The frozen physical source files (`production_kernel.py`, CMC gauge, scalar/matter source, Valencia and adapter modules) are compared blob-for-blob to the original source baseline, and the vendor submodule is pinned.

At the time this report was written, the source check, compilation and full repository test step had passed; the stage-capture step was still in progress and no final artifact had yet been listed. Therefore **the raw-RHS-versus-projection numerical result is not yet admitted by this report**. Update this section and the handoff once the run artifact is available. Do not infer raw projection impact from the old trace.

## 5. Separate correction already made: artificial unit-step RHS

This late-time event must not be conflated with the separate N=240/N=320 one-step *wrapper* failures caused by `species_rhs` calling `evolve_species(dt=1.0)` to extract a rate. A clean `main`-derived numerical correction is preserved as draft [PR #8](https://github.com/miahpenn/0-Cosmos-engine/pull/8). Its push and PR validation runs passed; an independent SI audit reproduced the artificial unit-step failure on frozen states and verified the real small-(dt) update remains admissible:
- Main-based PR validation: [38307704487](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307704487) and [38307788533](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307788533).
- Matched one-step A/B after correction: [38307235876](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307235876).
- Independent frozen-state confirmation: [38307129998](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307129998).

PR #8 remains draft for independent code review. It changes neither the radiation cone criterion nor the physical laws, and it does not fix the late-time predictor-metric crossing.

## 6. Still-open production gates / guardrails

- Do not merge the CMC centre-row candidate based on the radiation metric audit. Its centre (\dot A_a/r_0^2) does not yet demonstrate a regular (O(r^2)) continuum limit.
- Do not call the one-step Misner–Sharp balance residual a pass; the residuals are nonzero and no tolerance has been preregistered.
- Keep constraint propagation (H/M/(C_\Lambda)), spherical-centre regularity, metric/matter stage synchronization, radiation admissibility, conservation, resolution/CFL/domain convergence, and long-horizon behavior as independent gates.
- Do not relax (E\geq|S|), clip the conservative state, add arbitrary damping/floors, fit a feedback coefficient, or alter the light-boundary/CMC equations without a failed algebraic identity and a bounded counterfactual test.
- Do not infer a bounce, turnaround or completed cycle from a numerical failure or a short diagnostic.
- Preserve the GEAR archive and REV16 scorecard; the currently recorded `0*_eff` criterion remains false until its independent gates close.

## 7. Next decision rule

Once the stage-ledger artifact is available, compare the terminal-cell cone ratio/margin at all four stages. Then classify the crossing:
- If the **unprojected** Euler metric already makes the actual updated radiation state inadmissible, investigate the explicit geometry RHS / stage coupling and its finite-volume metric synchronization.
- If the unprojected predictor remains admissible but the **projection** crosses the cone, quantify how the exact algebraic projection changes (a) and (gamma^{rr}) across the outer band before considering any repair.
- If both are admissible but first CMC inputs differ, investigate instrumentation/order consistency first.
- In every branch, reproduce the signs and witnesses on saved data before any new physics change or long campaign.

No physics patch is justified by the current terminal budget alone.


## 8. Recorder implementation attempts — provenance and correction (12 Oct 2026)

The first stage-ledger attempt did not produce a scientific ledger. [Run 38308235431](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38308235431) passed compilation and **137 repository tests** before the recorder accessed a nonexistent `SphericalCellGrid.resolution` attribute. Artifact `11774360345` (921 bytes) contains provenance/hash files only; it contains no stage trace or ledger and must not be interpreted as an admitted capture.

A second attempt, [run 38308988809](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38308988809), passed compilation and **138 repository tests** in 30.75 seconds, but an overly broad textual correction had generated the invalid expression `self.latest_state.len(grid.centers)`. At the first captured step it failed with `AttributeError: 'ProductionState' object has no attribute 'len'`. This is also an instrumentation-only failure; artifact `11774360345` belongs to the prior run, and run 38308988809 has no admitted stage ledger.

The recorder has since been corrected to use one tested helper:
```python
def _grid_size(grid):
    return len(grid.centers)
```
The row builder and output serializer call `_grid_size(...grid)`, the regression test imports and exercises this helper on the actual pinned `SphericalCellGrid(n=80, r_max=80)`, and the source fingerprint list now names the active workflow file. The physics-source SHA/blob pins and vendor pin are unchanged.

- Grid-size helper source commit: `b76fbdfe7ff43e2878d532526a15123fcc30c02a`
- Test commit: `8a69e5dc16e3d3e4cf39fa0a91fc8ce312f55ba1`
- Active-workflow fingerprint commit: `9c3df9f85bc6dfd6c5ea01bb7cde871f5ea046b6`
- Corrected tagged workflow commit: `4d9fcb2bf1d7faccbe2013e1c71c70ea7c32bbd9`
- Corrected run: [38309326002](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309326002) (queued when this addendum was written)

Do not accept the corrected run until all of the following hold: (1) frozen source/vendor guard passes; (2) full repository test suite passes; (3) capture reaches the target failure and writes `stage_ledger.json` plus its summary; (4) ledger `validation.ok=true`, with all four stages for cells 75–79 and zero stage/Euler/projection mismatches; and (5) artifact hashes verify. A green test suite without a complete ledger is not a successful attribution.


### Additional capture-order correction discovered before admission

A focused reread of the frozen `V55ProductionKernel.step` found that after `_rhs` returns, production evaluates `geometry_stage_terms` first on the accepted state (`gterms0`), then again on the explicit predictor (`gterms_pred`), before the next call to `_rhs` increments the RHS counter. Therefore a recorder that keys only on `rhs_count == 1` can overwrite the accepted-state RHS with predictor-state terms. This is instrumentation contamination; it must not be interpreted as a geometry-law discrepancy.

A third, focused recorder guard is now in the branch:
- `_capture_first_geometry_terms_call(step)` tracks the first geometry-term evaluation during the capture window and accepts only its first call with `rhs_count == 1`.
- New tests verify the first call records, subsequent predictor/completed-state calls do not, and calls before (t=44) do not consume the counter.
- Source commit: `4f2e22681bfd086cbd35a22c58ac8e1452adace6`.
- Regression commit: `1f1f50b5f4d56cbe6297722f6149195e065cb36d`.

The currently running [run 38309326002](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309326002) was triggered before this last guard. Any ledger it produces must be checked for `validation.ok`; if the RHS identity fails or accepted geometry does not match the real Euler starting state, it is not admitted. The next tagged run must include this guard and its two new tests before any raw-versus-projection conclusion is recorded.


## 9. Final validated stage-ledger result — 12 October 2026

**This section supersedes the earlier “still running / not yet admitted” notes above.** The corrected run completed successfully, including the full trace capture and ledger validation.

- Run: [38309949383](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309949383)
- Workflow commit: `0b7b012e3d79880dc4b8bd74ffadb2891f6c513b`
- Frozen physics source baseline: `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`
- Pinned vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- Complete repository suite: **141 passed** in 30.70 seconds.
- Artifact ID: `11773874630`; archive: `radiation-predictor-metric-stage-ledger`.
- Downloaded artifact ZIP SHA-256: `356f067a583ec3596c3b15e8b5d17a18058cc7bd7a2bbc7fd379b6b1e2acbe93`.
- `stage_ledger.json` SHA-256: `35af5f350fe2888323cc851ab05190c5ddb981ef5042d656f4f45d4a803398bf`.
- `stage_ledger_summary.json` SHA-256: `6f2085abe8462d9ac33c71ef6c9c07bdb5fc63a360e73bc74b0776860b70526c`.
- `trace.json` SHA-256: `d2e9d507f8300598c4dd5ac45b4bcd3f77d205adc388e460e6a75ae614a214c9`.
- `summary.json` SHA-256: `79970f812f4fc45791b6d1b63951e4da1065c686554a47e3884a1f72b5778d82`.

### Admission checks

The ledger reports `status=radiation_admissibility_failure_captured`, `validation.ok=true`, 100 steps (indices 1467–1566), exactly five outer cells (75–79), and **zero validation errors**. Each of the 100 steps contains all four required stages. The frozen-source/vendor guard passed. The Euler predictor reconstructed as accepted geometry plus (dt\times\) the recorded explicit RHS with **zero absolute difference** across the captured fields. The post-projection determinant condition (ab^2=1) has maximum absolute error (3.33\times10^{-16}). Projected geometry and first predictor CMC inputs match exactly for (a,b,X).

This means the diagnostic ledger is now admitted as a faithful stage capture. The physical evolution remains failed at the recorded predictor stage; successful capture is not successful evolution.

### Validated terminal-cell witness

At the failing outer cell (i=79), (r=79.5), (dt=0.03), on the same updated conservative radiation state:

| Stage | (a) | (X) | (gamma^{rr}=X^2/a) | Cone margin (C) | (|S|/E) |
|---|---:|---:|---:|---:|---:|
| Accepted metric | 0.6089553401 | 0.9438956524 | 1.4630613182 | (+2.2171301\times10^{-6}) | 0.9925229671 |
| Unprojected Euler predictor | 0.6031899399 | 0.9465583195 | 1.4853905759 | (-2.0233473\times10^{-8}) | 1.0000682352 |
| After algebraic regularity projection | 0.6032035002 | 0.9465583195 | 1.4853571835 | (-1.6900205\times10^{-8}) | 1.0000569941 |
| First predictor CMC inputs | 0.6032035002 | 0.9465583195 | 1.4853571835 | (-1.6900205\times10^{-8}) | 1.0000569941 |

The original captured failure is unchanged: last accepted (t=46.98000000000102), central proper time (	au=22.016115393070738); attempted (t=47.01000000000102); the first predictor CMC lapse fails at (i=79).

### Cross-window result: projection is not the initiating crossing

Across all 500 cell-step pairs (100 steps × five outer cells):
- Algebraic projection **increased** the cone margin in all 500 pairs; there were zero decreases and zero ties.
- Projection reduced (|S|/E) in all 500 pairs.
- Exactly one raw Euler-predictor pair was outside the cone, and exactly one post-projection pair was outside it: the terminal cell of the terminal captured step.
- At that pair, projection improves the margin by (+3.3332677935\times10^{-9}) and reduces the ratio by (1.1241\times10^{-5}), but does not restore admissibility.

**Admitted numerical conclusion:** the cone crossing is already present on the unprojected explicit Euler metric predictor. The following algebraic regularity projection marginally improves the cone margin; it does not create the crossing. The proximate numerical failure is therefore the explicit geometry predictor (or synchronization of that predictor with the updated radiation state), not the projection itself. This is a stage attribution only, not proof that the underlying physical geometry equations are wrong.

### Next bounded discriminator

Do not change the cone condition and do not launch a long campaign yet. The next task is a frozen-state term attribution at the same final step: decompose the recorded explicit (dot a,dot b,dot X) by its existing signed RHS components and quantify their contributions to (Deltagamma^{rr}), (C) and (|S|/E), holding the actual updated conservative radiation state fixed. Reconstruct the same values at the prior accepted step as a control. Any proposed correction must then preserve the existing equations and improve the entire gate set (centre-series regularity, H/M/(C_\Lambda), radiation cone, outer boundary, and Misner–Sharp balance), not merely this terminal ratio.

The separate main-based matter-RHS rate-extraction fix remains in draft [PR #8](https://github.com/miahpenn/0-Cosmos-engine/pull/8); it addresses the artificial unit-step wrapper failure, not this validated late-time geometry-predictor crossing. `main` remains at `2035acc1210db18da0cededea18750b7ed813260`.
