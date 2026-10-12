# 0-Cosmos Engine — Current-State Reconciliation and Next Diagnostic Gate

**Date:** 12 October 2026 (UTC)  
**Type:** research / handoff only; no production source changes  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Production branch observed:** `main` at `2035acc1210db18da0cededea18750b7ed813260`  
**Radiation trace source:** `1c63090d4f42396f899421a2379f70320689d5ad`  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Disposition:** continue bounded, frozen-state/source-level diagnosis. Do not start another long trajectory or merge a diagnostic candidate on this evidence.

## Executive reconciliation

Several recent results must be kept separate:

1. **The N=240/320 first-cell radiation “failure” was not a real production-step failure.** The earlier diagnostic called `evolve_species(dt=1.0)` to extract a RHS and applied the finite-time radiation-cone check to that artificial `U + L(U)` state. The corrected contract bypasses that *post-update* guard only for derivative extraction; real direct finite-time updates and the actual production stages remain validated. The one-normal-step capture passes at N=240 (`dt=0.01)) and N=320 (`dt=0.0075)) at initial, explicit-predictor, primary-predictor, and final stages. This is one-step admission only, not multi-step stability.
   - Run: [38318958640](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318958640)
   - Report: [high-resolution actual-stage check](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/high-resolution-radiation-real-step-after-rhs-contract-2026-10-12/docs/HIGH_RESOLUTION_RADIATION_REAL_STEP_AFTER_RHS_CONTRACT_2026-10-12.md)
   - Exact tested head: `a5ca89e57599048503cf259d0fb10f6bfd518700`; report JSON SHA-256 `21742fcefebb3155b61aab4963086a0e31875756fe17c1af03cce09146f07453`.

2. **The separate N=80 late-time radiation cone crossing is real in the captured corrected-rate trajectory.** At accepted (t=46.98000000000102), the next explicit predictor (attempted (t=47.01000000000102), (dt=0.03)) reaches cell (i=79, r=79.5) with (|S|/E=1.0000573283549097) and densitized cone margin (C=-1.699741681811858	imes10^{-8}). The saved budget closes exactly. On the accepted metric, matter flux and geometric-source increments together improve (C) by (+1.7468507316	imes10^{-6}); evaluating the predicted matter on the explicit predictor metric then contributes (-2.2367571625	imes10^{-6}), putting the state just outside the cone.
   - Trace: [38310216296](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310216296), artifact [11774592255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310216296/artifacts/11774592255)
   - Geometry/source attribution: [report](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/RADIATION_PREDICTOR_FAILURE_GEOMETRY_SOURCE_ATTRIBUTION_2026-10-12.md), artifact ZIP SHA-256 `f64298e2819e0dd9c4c2df712a4d1250e4855a5b0e467430fe7aa419a533e49f`.

3. **The cell-centred CMC row is a strong numerical-discretization defect candidate, but not a complete repair.** With (r_i=(i+1/2)Delta r), the old `alpha[0] = alpha[1]` row equates two samples at distinct positive radii and omits the first-cell CMC PDE. The new three-arm one-step comparison separates row choice from stencil order:
   - A: old equality row;
   - B: second-order first-cell CMC PDE row;
   - C: parity-aware fourth-order PDE row matched to the production K RHS.
   
   On N=80/120/160, A→B produces the large momentum/connection improvement; B→C constraint differences are around (10^{-11}) or smaller on fixed-annulus witnesses. Thus the dominant benefit is restoring the missing first-cell PDE row, not fourth-order accuracy by itself. Full-domain H responds mixedly and the required (\dot A_a/r_0^2) centre-series criterion remains open. Keep PR #7/CMC candidate isolated.
   - Run: [38304709682](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304709682)
   - Artifact: [11771896263](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304709682/artifacts/11771896263), ZIP SHA-256 `59918d53bc883fbc796bf79057666692a01c49a0a355bd26f89cbac22fbb08d8`
   - Analysis in [production-readiness review](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/PRODUCTION_READINESS_DERIVATION_REVIEW_2026-10-12.md).

4. **The broader constraint/balance gate remains open.** A discrete-consistent initial H/M/(C_\Lambda) slice still develops a material nonzero momentum-constraint tangent matching the one-step quotient. Manufactured regular-profile operators give near-fourth-order errors on a fixed annulus, so a generic parity-sign or constraint-evaluator bug is disfavoured; no unique bad production RHS term has been proven. The t=12 N=40/60/80 campaign is execution evidence only, not scientific admission. The Misner–Sharp residual remains nonzero and needs a preregistered tolerance.

## Source-level audit: what the current trace can and cannot attribute

I inspected the corrected-rate trace source `1c63090d4f42396f899421a2379f70320689d5ad`, the production stage order in `engine/production_kernel.py`, the adapter `engine/v55_pirk_adapter.py`, and the exact pinned vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215` (`moving_puncture_pirk.py`, `vacuum_bssn.py`, and `bssn_characteristic_boundary.py`).

### 1. The failing attempted step does not reach its own boundary call

The zero-shift explicit predictor first advances (a,b,X) from the accepted slice and updates the conservative matter/scalar state. The kernel then calls the CMC lapse solve for that explicit predictor. The saved exception occurs during its total-matter/radiation primitive recovery, before that solve returns. The P-minus outer-boundary operation for the **current attempted step** is later: after the primary (A_a,K,\Lambda) update, and before its CMC re-solve; it is applied again after the completed final slice.

Therefore the current-step boundary call is **not the direct operation that crosses the cone at this substage**. But this does not exclude a *history pathway*: the accepted (A_a[-1]) entering the current explicit RHS may already have been set by the previous step's boundary operation. The corrected trace does not save the pre-boundary (A_a[-1]) for that accepted state, so the existing trace alone cannot measure that counterfactual. The older 100-stage P-minus ledger is source-pinned to historical commit `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`, not the corrected-rate trace source, and must not be treated as a matched causal A/B.

### 2. The geometry-to-cone sensitivity is exact for the recorded state

At zero shift and with the recorded `KO_EPSILON=0`, the pinned explicit equations reduce to
[
\dot a=-2\alpha a A_a,\qquad \dot X=\alpha X K/3.
]
Since (\gamma^{rr}=X^2/a),
[
\frac{\dot\gamma^{rr}}{\gamma^{rr}}
=2\frac{\dot X}{X}-\frac{\dot a}{a}
=2\alpha\left(A_a+K/3\right).
]
Thus both (A_a) and (K) affect the predictor inverse radial metric. This identity is also directly verified cell-by-cell in the trace artifact (RHS formula residual zero to stored precision); it is not evidence of an incorrect evolution law.

At the failing cell, the critical metric for the predicted conservative state is
[
\gamma^{rr}_{\rm crit}=(U_E/|U_r|)^2
=1.4850469620709852,
]
while the explicit predictor has (\gamma^{rr}=1.4852172340500571), only (1.7027557\times10^{-4}) above critical (about (0.01147\%)). The inverse-metric increment is (+0.0223233464); the saved decomposition attributes 62.63% to (a)'s update and 37.37% to (X)'s update. The update is reproducible and the budget closes, but the trace does not prove whether the timestep/stage coupling, boundary history, or another continuum-compatible effect is ultimately responsible.

### 3. P-minus algebraic closure is not yet CMC/zero-shift evolution compatibility

The pinned boundary module forms (\Omega_-), chooses
[
\Omega_-[-1] = \Omega_-[-3]+2\Delta r\,S_-[-2],
]
and solves this relation for (A_a[-1]). Its residual is computed from the same ((\Omega_-[-1]-\Omega_-[-3])/(2\Delta r)-S_-[-2]) construction. A near-roundoff residual therefore establishes the algebraic closure of this frozen discrete target, not independent proof that the target is a complete constraint-preserving boundary, dynamically compatible with the CMC elliptic lapse normalization (\alpha[-1]=1), or compatible with every constraint mode. The vendor module explicitly disclaims a complete constraint-preserving boundary. It corrects (A_a[-1]) only, and the production kernel then re-solves CMC and recomputes the secondary stage.

**No boundary defect is established here.** The next evidence should decide compatibility from the actual coupled equations and matched state data rather than assuming either that the P-minus closure proves the boundary good or that the late cone crossing proves it bad.

## Bounded next gate — no new long run

1. **Derive first:** reconstruct the continuum and frozen-discrete evolution identity for (\Omega_-) under the actual production gauge ((\beta^r=0), CMC lapse with outer normalization), including the matter terms in the total (T_{\mu\nu}). State precisely which field equations the P-minus source assumes and which terms survive or cancel under this gauge. Keep target closure, characteristic direction/speed, H/M/(C_\Lambda) propagation and radiation-cone margin as separate tests.
2. **Then do one frozen-state A/B, if the derivation leaves a live ambiguity:** on an identical saved accepted slice, record pre/post-boundary (A_a[-1]), (Omega_-), (a_t), (X_t), (\gamma^{rr}_{\rm pred}), (C_{\rm pred}), CMC target/lapse, H/M/(C_\Lambda) and Misner–Sharp terms. The boundary-on/off arm is strictly diagnostic and must not become a production alternative. Pin all source and vendor blobs. If no full-grid accepted snapshot exists at the corrected trace's failing time, say so and instrument the *next targeted capture* to preserve one; don't pretend the historical artifact is matched.
3. Join the CMC-row candidate only after it passes the same-state centre-series, radiation stage-cone and Misner–Sharp checks. Do not merge #7, #8, #10 or #11 from unit tests alone. For RHS extraction, PR #8's direct rate terms avoid constructing `U+L(U)` with a fake dt; PR #10's narrower flag suppresses only the derivative-operation output guard while retaining real-update checks. These implementations need independent operation-order/regression review before choosing one canonical path. PR #11's zero radiation rest-current channel is orthogonal and its stated bit-identical A/B is not a cure for the late metric-driven crossing.
4. Only after the bounded identities pass: agree preregistered tolerances, freeze an exact candidate commit, run exact-commit CI, then one controlled multi-resolution full trajectory. No t=160 campaign yet.

## Guardrails

- Preserve exact radiation (E\ge |S|); no clipping, floors, relaxed cone, accepted-metric recovery, fitted limiter, damping or projection.
- No production physics, EOS, gauge, boundary or coefficient change without a demonstrated equation-level failure and matched regression.
- Keep the canonical GEAR archive and REV16 criterion unchanged. Numerical progress does not make the outstanding D-to-matter, local-to-global (Q), clock-ownership or scientific cycle admission gates true.
- Keep branch-specific test counts branch-specific. A successful diagnostic artifact or unit suite is not full production admission.

**Bottom line:** the software false-positive in high-resolution RHS extraction is isolated and corrected in tested candidates; the real late radiation crossing is a predictor-metric margin loss at a near-saturated state. The cell-centred CMC row issue is independently supported but incomplete. The immediate task is a source-derived boundary/gauge compatibility decision and, only if needed, a matched frozen-state test—not another long evolution.
