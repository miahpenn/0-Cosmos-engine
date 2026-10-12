# Radiation Predictor Failure: Geometry–Source Margin Budget

**Date:** 2026-10-12  
**Status:** offline attribution from a previously captured corrected-rate trace; no evolution rerun and no physics modification.

## Provenance

- Canonical trace workflow: [38310237842](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842).
- Artifact: [11775230036 — radiation-trace-after-rate-fix](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842/artifacts/11775230036).
- Artifact ZIP SHA-256: `f64298e2819e0dd9c4c2df712a4d1250e4855a5b0e467430fe7aa419a533e49f`.
- Embedded `trace.json` SHA-256: `53a1f9f188e0eea4575da4458d09891f50144b9f89f945c995cd9980221fb619`.
- Embedded `summary.json` SHA-256: `fd3f5a81992277974d675b0908fec5c68aa950e25c4c0e0a06746034bb939627`.
- Rate-fix source: `1c63090d4f42396f899421a2379f70320689d5ad`.
- Instrumentation commit: `51335a56873d65b8ed28863b7658ce02793ae396`.
- Pinned vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Configuration: N=80, (r_{\max}=80), (Delta r=1), (dt=0.03), radiation on, (D) amplitude (10^{-10}), current production outer light-constraint boundary and radiation closure. Predictor radiation recovery uses the actual predictor metric; no retries, clipping, floors, damping, or candidate source changes.

The run completed 1,566 accepted steps and captured the expected rejection when attempting (t=47.01000000000102) from accepted (t=46.98000000000102), at the first explicit predictor CMC-lapse solve, cell (i=79), (r=79.5). The source capture and all saved reconstruction identities are present for this failing step.

## Exact failing-step budget (cell 79)

The conservative cone margin is (C=U_E-\sqrt{\gamma^{rr}}\,|U_r|), with the source's densitized conservative variables and metric for the indicated stage.

| Stage / decomposition | (C) or increment |
|---|---:|
| Accepted state on accepted metric | (+4.7290899469088366\times10^{-7}) |
| After predictor matter flux, still on accepted metric | (+1.4840320824389746\times10^{-6}) |
| After predictor matter sources, still on accepted metric | (+2.2197597480853805\times10^{-6}) |
| Change from accepted to predictor metric at fixed predicted matter state | (-2.236757174073193\times10^{-6}) |
| Final explicit-predictor margin on predictor metric | (-1.6997425987812665\times10^{-8}) |

The matter update improves the cone margin on the accepted metric. The budget's separate integrated contributions are:

- Flux contribution to (C): (+1.011123087748091\times10^{-6}).
- Geometric/source contribution to (C): (+7.357276656464059\times10^{-7}).
- Combined matter contribution: (+1.7468507533944968\times10^{-6}).
- Geometry-metric contribution: (-2.236757174073193\times10^{-6}).
- Net accepted-to-predictor change: (-4.899064206786963\times10^{-7}), closing exactly to the saved (C) values in the artifact (closure error (0.0)).

Thus the failed margin is not explained by the captured matter flux/source update alone. It is the same predicted (U_E,U_r) evaluated with the explicit predictor metric that crosses the existing cone.

## Threshold and geometry decomposition

At the failed predictor:

| Quantity | Value |
|---|---:|
| Predictor (U_E) | (2.964923180513671\times10^{-4}) |
| Predictor (U_r) | (2.4330069621024254\times10^{-4}) |
| (|S|/E) | (1.0000573283857723) |
| (a) accepted (	o) predictor | (0.6089966848422357 	o 0.6032317132323858) |
| (X) accepted (	o) predictor | (0.943873683323241 	o 0.9465358623871615) |
| (\gamma^{rr}=X^2/a), accepted (	o) predictor | (1.4628938912877238 	o 1.485217237641921) |
| (Delta\gamma^{rr}) | (+0.022323346354197193) ((+1.5259716708\%)) |
| Critical (\gamma^{rr}_{\mathrm{crit}}=(U_E/|U_r|)^2) for this predicted matter state | (1.4850469620709852) |
| Predictor metric excess above that threshold | (+0.00017027557093585166) ((+0.0114660059\%) relative to threshold) |

The predictor is only slightly outside the cone: the recorded (|S|/E) excess is about (5.73\times10^{-5}) relative to unity. But the matter state is near the cone boundary, so the actual metric move exhausts the available margin.

The recorded (\gamma^{rr}=X^2/a) decomposition attributes the inverse-metric increment to:

| Contribution | (Delta\gamma^{rr}) | Share of increment |
|---|---:|---:|
| (a) update | (+0.01398060076501273) | 62.6277% |
| (X) update | (+0.008342745589184464) | 37.3723% |
| Centre-regularity projection | (0.0) | 0% |
| Total | (+0.022323346354197193) | 100% |

The actual explicit RHS values recorded at this cell are (dot a=-0.19216572032833137) and (dot X=+0.08873930213068544) at the accepted metric/stage. The raw Euler candidate matches the recorded predictor fields. The centre-regularity projection changes neither (a), (b), nor (X) here.

## What the artifact verifies — and what it does not

The budget reconstruction is internally consistent over the saved window:

- 67 geometry-predictor budgets: maximum (\gamma^{rr})-decomposition closure error (0); maximum (a)-RHS formula residual (0); maximum (X)-RHS formula residual (2.7755575615628914\times10^{-17}); maximum (X) projection change (0).
- 66 radiation predictor budgets: maximum conservative reconstruction/closure errors (0).
- 66 completed-state budgets: maximum budget/source-split closure error (5.421010862427522\times10^{-20}); conservative-variable predictor/completed reconstruction errors (0).
- No instrumentation errors were recorded.

This rules out the obvious bookkeeping explanation that the traced cone budget was assembled inconsistently, and shows precisely how the existing predictor calculation reaches the rejected state. It does **not** prove the (a) or (X) evolution law is wrong, determine whether the chosen stage/timestep/gauge coupling is appropriate near a saturated cone, or establish which layer should change. It is a same-step sensitivity attribution, not an alternative-physics comparison.

The immediate unresolved question is now narrower: what production-admissible time-integration or coupled geometry/matter treatment is justified when the correct explicit predictor metric consumes the remaining cone margin? The answer must preserve the exact (E\ge |S|) admission rule and be derived/tested without arbitrary limiter, clamp, floor, damping, or fitted coefficient.

## Next gate

1. Independently audit the pinned-vendor (a/X) equations and the production adapter's split/order against the archived PIRK blueprint. The recorded formula identities pass; this is a mathematical/partition compatibility audit, not permission to change the RHS.
2. Use the critical (\gamma^{rr}_{\mathrm{crit}}) value above as a preregistered scalar witness for any offline candidate comparison on this frozen failing snapshot. Report the complete metric and matter state together, with the same source revision.
3. Do not change the outer boundary, radiation transport, CMC projection, or physical source law based on this attribution alone.
4. Do not launch another long trajectory until a candidate-specific, bounded question and fail-closed criterion are approved. The one-step `beta_rad` counterfactual report is a separate experiment and does not explain this late event.

No trajectory was rerun for this analysis.
