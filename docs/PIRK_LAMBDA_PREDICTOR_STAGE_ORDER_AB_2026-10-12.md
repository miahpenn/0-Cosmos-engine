# PIRK Lambda predictor stage-order A/B — 12 October 2026

**Status:** Diagnostic-only one-step A/B completed at N=80 and N=160. The test confirms a stale-vs-sequential stage-value mismatch in Lambda-L2, but the measured one-step state and radiation-cone impact is negligible at these settings. It does not explain the separate late-time cone crossing.

## Provenance

- Workflow: [38318309196](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196) — success.
- Exact audit commit: `ce9a49560b04ccbe33d5b51ac86219a5c4a3c2f9`.
- Preceding stage-order regression: [38318119471](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318119471) — success.
- Artifact: [11777443806 — pirk-lambda-predictor-stage-order-ab](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196/artifacts/11777443806).
- Artifact ZIP SHA-256: `1ec986b7cbbe97ba23184c2d448aba26199257f6ae9d6c4182490c81af3867fe`.
- Embedded report SHA-256 (also in report.sha256): `283c5b337a90670c9a397c9584d19447497017cc44a640ca96bc7303a03186ef`.
- Pinned geometry/PIRK vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Tests: 1 focused stage-order test and 36 full-suite tests passed before the A/B capture.
- Main was unchanged at `2035acc1210db18da0cededea18750b7ed813260`.

## Experiment

The stale arm replays the former behavior: the third Lambda-L2 call returns the Lambda-L2 array captured on the explicit U1 geometry before Aa/K are updated. The sequential candidate recomputes Lambda-L2 on stage1 after the predictor Aa/K updates and before the predictor boundary and the next CMC lapse solve. Both arms use bit-identical initial states, the same pinned vendor, the same radiation setting, and exactly one ordinary kernel step. All other production equations and update order remain shared.

| Quantity | N=80 | N=160 |
|---|---:|---:|
| dt | 0.015 | 0.0075 |
| max difference between corrected and stale Lambda-L2 arrays | 5.3217372642e-8 | 1.3313374252e-8 |
| expected/observed predictor Lambda difference | 3.9913029479e-10 | 4.9925153421e-11 |
| max completed-state Lambda-array difference | 6.9388939039e-17 | 2.1684043450e-18 |
| max completed-state Aa-array difference | 1.6329877488e-12 | 1.9526410357e-13 |
| completed-state H full-domain max-absolute delta | -2.4286128664e-17 | 0.0 |
| completed-state M full-domain max-absolute delta | -4.2540321716e-12 | -1.0631220974e-12 |
| completed-state C_Lambda full-domain max-absolute delta | -6.9388939039e-17 | -2.1684043450e-18 |
| max radiation |S|/E, legacy vs candidate | 0.00501354074617448 / identical | 0.005015096011912016 / identical |
| min radiation cone margin, legacy vs candidate | 2.2067153165e-9 / identical | 5.5163762656e-10 / identical |

All four arm outcomes (legacy/candidate x N=80/160) completed. The radiation full-domain and outer-band witnesses are admissible in every arm. The report captures full diagnostic summaries, candidate-minus-legacy array deltas, the stage L2 delta, and source/vendor provenance.

## Interpretation

1. The pinned sequential stage order does require the Lambda-L2 evaluation to follow the predictor Aa/K update. The regression test verifies that call sequence and confirms the stale and recomputed L2 arrays differ above the test's tolerance.
2. The bounded A/B confirms that this order difference propagates into the predictor Lambda value by about 4.0e-10 at N=80 and 5.0e-11 at N=160 in this initial state.
3. After a complete step, the metric-array difference is at roundoff for Lambda and no measurable change appears in the radiation E/S cone witnesses. H and C_Lambda deltas are roundoff-sized; the largest reported full-domain M max-absolute delta is about 4.25e-12 at N=80 and 1.06e-12 at N=160.
4. Thus the ordering mismatch is a legitimate discrete stage-order defect candidate, but this one-step A/B gives no evidence that it materially drives the late-time outer-cell radiation crossing. Do not claim it cures the crossing or use it to justify physical changes.

## Admission limits and next step

This remains a diagnostic branch, not a merged production repair. The one-step result does not establish long-time convergence, constraint propagation, Misner-Sharp balance, or production readiness. If proposing a source correction, keep it in a separate reviewable draft branch, preserve the pinned-vendor stage-order test, and demand a matched longer bounded trace only after source-level analysis identifies how a per-step roundoff-scale final-state difference could accumulate. Do not start a long campaign or change physics, boundary conditions, gauge, or radiation admissibility from this result alone.

## Required follow-up to the late cone investigation

The late corrected-rate N=80 crossing remains the higher-priority unresolved blocker: accepted t is about 46.98, attempted t about 47.01, cell i=79; the simultaneous explicit a/X update pushes the predictor beyond the radiation cone threshold. That source-matched finding remains independent of this Lambda one-step A/B. Preserve both records without conflating their causes.