# Radiation Predictor Metric Split — Preliminary Reconstruction

**Date:** 2026-10-12  
**Classification:** PRELIMINARY / UNADMITTED. This note preserves an independently reconstructed signal from a stage ledger that failed its admission validator. It must not be treated as a successful physics run.  
**Frozen production source baseline:** `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Active diagnostic branch:** `diag/radiation-predictor-metric-stage-ledger-2026-10-12`  
**Current intended validation run:** [38309949383](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309949383)

## Question

At the known late-time outer-cell radiation-cone failure, is the algebraic metric regularity projection itself pushing the radiation state out of the cone, or is the unprojected explicit Euler metric already inadmissible?

This is distinct from the high-resolution false failure caused by rate extraction through an artificial `dt=1` call; that independent correction is preserved in draft [PR #8](https://github.com/miahpenn/0-Cosmos-engine/pull/8).

## Frozen reference event

The original saved N=80 trace records the last accepted state at:
- (t=46.98000000000102),
- (	au=22.016115393070738),
- accepted step index 1566, with (dt=0.03).

The next explicit predictor fails at (t=47.01000000000102), in `explicit_predictor_first_CMC_lapse`, outer cell (i=79, r=79.5). The accepted metric cone ratio is below 1; the post-projection predictor metric ratio is (1.0000569941246342), so the unchanged radiation admissibility requirement is violated slightly.

Existing ordered budget: the radiation conservative flux and geometric-source increments, evaluated using the accepted metric, improve the cone margin. Changing to the predictor metric then drives that state outside the cone. That alone identifies a metric-stage contribution but does not split raw geometry evolution from projection.

## Attempt 3 artifact — rejected admission, valuable raw rows

Workflow [38309326002](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309326002) passed source/vendor checks, compilation and 138 tests, then reproduced the terminal failure and wrote 100 steps × five outer cells for all four stages. The run ended nonzero because `stage_ledger_summary.json` correctly reported `validation.ok=false`.

- Artifact ID: `11774216769`
- Archive: downloadable from the Artifacts section of [workflow 38309326002](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309326002); this unadmitted archive is also preserved in the current session as a local copy.
- ZIP SHA-256: `7d5c59bf585d3a7733f1cddccf09e876b4705acafc36eb88e0942b77151e4758`
- Internal `stage_ledger.json` SHA-256: `bb7944cb21c5f9f8f0c2b85b066a49c6c82243489cebf7e083d6064b3bbcb3f8`
- Trace SHA-256 recorded in the package: `52d65eec16aa39ad0d7e58277ee4a8715be0309012145b9286ec0881fb3574f5`
- Captured failure: `ValueError`, (t=47.01000000000102), cell 79, radiation ratio (1.0000569941246342); capture status `radiation_admissibility_failure_captured`.

Why admission failed:
1. The accepted-state geometry RHS record was overwritten by a subsequent predictor geometry-term call. The frozen step evaluates `gterms0` on accepted geometry and `gterms_pred` on predictor geometry before the next matter RHS call increments its counter. The original capture guard did not distinguish these calls. This caused 1,485 Euler-reconstruction mismatches.
2. The intermediate accepted/raw/projected metric rows did not carry their own `U_E` and `U_r` fields/cone reconstructions. This caused missing-field errors across the four stages (5,970 validation errors total).
3. Thus the report is not admitted even though it contains valid recorded raw/post-projection metric arrays and exact first-CMC radiation values.

## Provisional cone split, reconstructed from actual recorded values

For every captured step, we used the exact (U_E,U_r) from that step's first predictor CMC input record and independently reconstructed the cone on the corresponding raw and projected metrics. This bypasses only the missing metadata fields for the calculation; it does **not** cure the failed Euler identity or upgrade the artifact's admission status.

At terminal cell (i=79):
| Quantity | Unprojected Euler predictor | Post-projection predictor |
|---|---:|---:|
| (a) | captured raw value | captured projected value |
| (gamma^{rr}=X^2/a) | (1.4853905758553785) | (1.4853571835116133) |
| (ab^2) | (0.9999325599738486) | approximately (1) |
| Cone margin (C=U_E-|U_r|sqrt{X^2/a}) | (-2.0233472965\times10^{-8}) | (-1.6900205172\times10^{-8}) |
| (|S|/E) | (1.0000682352118357) | (1.0000569941246342) |

The projection changes the margin by (+3.3332677935\times10^{-9}): it **improves** the margin but does not restore admissibility. In the provisional sample of 100 steps × five outer cells, the projected-minus-raw cone-margin delta was positive in all 500 cell-step pairs. Only the terminal outer cell at the terminal attempted step crossed the cone, and it was already outside before projection.

**Inference, still provisional:** the metric projection does not appear to initiate the crossing. The leading target becomes the raw explicit geometry update or its synchronization with the matter predictor. Do not change physical equations until the next ledger has passed its algebraic admission contract.

## Instrumentation fixes on the current branch

The current recorder revision adds:
- `_grid_size(grid) -> len(grid.centers)`, tested against the actual pinned `SphericalCellGrid`; no `.resolution` attribute is assumed.
- `_capture_first_geometry_terms_call(step)`, so only the first accepted-state geometry-term call is stored as the explicit RHS; predictor geometry evaluations cannot overwrite it.
- `_attach_radiation_state(values, radiation)`, putting exact `U_E`, `U_r`, (gamma^{rr}), cone margin and ratio into accepted, raw and projected records for independent reconstruction.
- Focused regression tests for grid access, first-call capture order, and radiation conservative metadata.
- A workflow guard that keeps physical source files byte-identical to the frozen baseline and verifies the vendor commit.

Run [38309617374](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309617374) was launched after the RHS-capture guard but before the radiation metadata helper landed; its validation rejection is expected and documented. Run [38309949383](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309949383) is the intended run with both fixes. At this note's creation, its test suite is green and the bounded capture is running. No additional long run should be launched unless this run fails a specific admission check that is first repaired.

## Admission contract for the next result

A ledger is accepted only if:
1. source/vendor pins pass;
2. compile and the full repository suite pass;
3. the known terminal failure is reproduced on the frozen production source;
4. `stage_ledger.json` includes all four stages for five cells on every step in the capture window;
5. `validation.ok=true`, with zero raw Euler identity mismatches, zero reconstructed-cone mismatches, the post-projection determinant satisfied, and exact match between projected (a,b,X) and the first CMC inputs;
6. trace, ledger, summary and source hashes verify.

A green suite or a captured failure alone is insufficient.

## Production-readiness boundary

No physical law, cone condition, radiation EOS, CMC gauge, outer-light boundary, damping or limiter was changed. The full machine is still not admitted as production ready: centre-series regularity, constraint propagation, conservation tolerances, resolution/CFL/domain convergence and long-horizon behavior remain separate gates. The canonical GEAR archive and REV16 scorecard remain authoritative; no bounce, turnaround, reset or closure should be inferred from this numerical failure.
