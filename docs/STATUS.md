# Build Status — 2026-10-12

## Objective and current classification

Objective: get the unified V5.5 production machine running through a controlled, synchronized multi-resolution campaign without changing locked physics before numerical admission.

The canonical state graph is implemented. The first synchronized production campaign at N=40,60,80, r_max=80, CFL=0.03, radiation enabled, through t=12 has completed, and the fail-closed completion gate passed. This is **horizon completion, not physical/numerical admission**.

Current classification:

**IMPLEMENTED / 83-TESTED / SHORT-SMOKE-TESTED / t=12 MULTI-RESOLUTION RUNS COMPLETE / CONSTRAINT AND BALANCE ADMISSION STILL OPEN**

## Locked physics and scope

- Corrected archive COSMOS background equations.
- Local GEAR-03 S/D sector with GEAR-63/65 barrier and corrected COSMOS potential normalization V0=8.242522415500654e-5.
- Conservative spherical dark matter, baryons and radiation, with the locked beta DM/scalar exchange.
- One total Einstein stress-energy assembly and pinned reference-metric BSSN/PIRK numerical geometry.
- No fitted coefficient, imposed clock conversion, D-to-matter law, shell EOS, phenomenological interface source, artificial bounce/reset/branch flip, lapse floor, ejection threshold or physical stop.
- Bianchi-I shear remains outside the exact spherical source. Diagnostic-only forced CMC histories are not part of production physics.

## Engineering admission fix and checks

Draft PR #3: https://github.com/miahpenn/0-Cosmos-engine/pull/3

Branch: fix/strict-campaign-admission-2026-10-12

- The campaign runner now writes an admission result and raises if requested resolutions are missing, a resolution did not complete, or a trajectory fails to reach the requested horizon.
- Available campaign outputs are uploaded even when the run fails.
- The runner checks numeric geometry, scalar and matter state fields and stored diagnostic values for non-finite entries.
- The initial t=0 slice is now in the history, and an endpoint discrepancy within floating-point roundoff is recorded and canonicalized without a tiny extra evolution step.
- The strict smoke workflow passed. Latest strict smoke: full repository suite passed and N=40, r_max=40, CFL=0.03, t=1 finished and produced its ledger. Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253695

## Synchronized production campaign through t=12

Latest fixed-radius campaign run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687

- N=40: 200 steps, tau_final=11.9559993814.
- N=60: 300 effective steps after endpoint correction, tau_final=11.9284259877.
- N=80: 400 effective steps after endpoint correction, tau_final=11.9184114053.
- All requested resolutions reached t=12 and the state/diagnostic values were finite.
- Endpoint canonicalization in the rerun corrected only +7.9936057773e-14 for N=60 and +6.3948846218e-14 for N=80; no additional physical evolution step was inserted.
- Post-correction Misner-Sharp current residual max/RMS: N=40 2.3899e-3 / 1.2621e-3; N=60 5.5462e-4 / 2.3152e-4; N=80 6.0062e-4 / 2.1438e-4. The enormous old N=60/N=80 maxima were endpoint-sliver artifacts. The remaining current balance is smaller and improves strongly in RMS from N=40, but its resolution behavior is not yet a completed balance proof.

### Constraint-location finding

The earlier norm called hamiltonian_max excluded the first two *cells*, not the same physical region at each resolution. After correcting the recorded peak location, that maximum sat exactly on the first included cell: r=5.0 for N=40, r=3.3333 for N=60, r=2.5 for N=80. This means the former cross-resolution comparison used a moving near-center cutoff and is not a valid convergence norm.

At t=12, the full-domain Hamiltonian residual maxima were:
- N=40: 4.0417e-2 at r=1.0.
- N=60: 7.9698e-2 at r=0.6667.
- N=80: 1.3494e-1 at r=0.5.

Those full-domain maxima include the first two center cells, where the pinned vendor's formal constraint stencil is not admitted. The admitted inner-region comparison excluding those two cells still shows a substantial N=60/N=80 Hamiltonian residual at t=12 (about 1.072e-2 at r=3.333 and 1.962e-2 at r=2.5), so the concern is not solely the excluded center cells. The fixed physical exterior domain r >= 5.0 is much better behaved at N=60/N=80.

### Signed checkpoint profile audit and initial-slice comparison

Post-processing-only audit succeeded on 21 states: the production initial slice and six scheduled checkpoints for N=40,60,80. It reloads saved states, computes signed residual profiles without evolution, and preserves the actual accepted state time separately from the scheduled checkpoint label.

- Audit workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209173946
- Artifact: `checkpoint-constraint-profiles-t12`, SHA-256 `af215b6a59105100777c1f709039157d9f2d94091a95d45cad313df4ebd6b`.
- At t=0, the H residual is about 3.8e-3 to 4.0e-3 in the first center cell, but at the first admitted inner cell i=2 it is only about 1e-5. Initial momentum is at roundoff and the connection constraint is zero by construction.
- By t=12, the admitted inner H residual has grown to about 1.072e-2 (N=60, r=3.333) and 1.962e-2 (N=80, r=2.5). At those same points the raw geometric Hamiltonian term is about -9.41e-3 and -1.81e-2 respectively, while the assembled total matter source (16\pi\rho) is about +1.31e-3 and +1.52e-3. The mismatch is therefore primarily on the geometric side of the production constraint expression, not just a matter-source magnitude error. This does not yet identify which geometric derivative or stage is responsible.
- The growing Hamiltonian and connection errors co-locate toward the center. At N=80, t=12, r=2.5, the connection residual is about -1.449e-2 alongside H about -1.962e-2. This suggests testing the BSSN connection/Ricci pathway on the same state, but is not by itself causal proof.

### Frozen center-regularity projection sensitivity

The pinned vendor exposes `enforce_algebraic_regularity`. We applied it to a copy of each saved geometry at t=0 and t=12, reassembled the matter source, and recomputed constraints. No projected state was evolved or inserted into production.

- It repairs the determinant identity to roundoff and reduces the t=12 inner connection maxima from 2.07e-2 to 7.28e-3 (N=40), 2.52e-2 to 1.14e-2 (N=60), and 3.09e-2 to 1.46e-2 (N=80).
- It also reduces the first-cell momentum maxima, but it does **not** remove the admitted inner Hamiltonian defect: at i=2 it changes H from -1.072e-2 to -1.174e-2 for N=60, and from -1.962e-2 to -2.142e-2 for N=80. The exterior H residual is essentially unchanged at N=60/N=80.
- Decision: do **not** add this projection to the production step as a proposed fix. Its frozen-state effect shows that center algebraic regularity matters for some witnesses, but projection alone does not cure the admitted Hamiltonian residual and worsens it at the first admitted cell.

### Independent metric-Ricci decomposition on the same 21 states

The saved-state curvature audit is complete:
- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209394006
- Artifact: `checkpoint-constraint-profiles-t12`, SHA-256 `71306184a4980bea7fdeb0376307efe3ee09d90d34717484e80fef5686b269e6`.
- All algebraic split, connection-formula and Hamiltonian-gap closure checks pass at roundoff scale (worst reported closure below 3e-17).
- It independently reconstructs spatial scalar curvature from the physical metric using `R_areal = r*sqrt(b)/X` and proper radial distance, without using evolved Lambda; production BSSN Ricci remains the comparator.
- At the first admitted inner cell at t=0, H is about -6.25e-6 (N=60) and +1.77e-5 (N=80); the curvature-route gap is only +3.42e-5 and +2.65e-5.
- At t=12 the same cells have H = -1.072e-2 (N=60, r=3.333) and -1.962e-2 (N=80, r=2.5). The metric-derived Ricci counterfactual reduces these to -7.09e-3 and -1.277e-2, respectively, but does not close them.
- At those locations the Ricci gap is -3.63e-3 (N=60) and -6.86e-3 (N=80). Its exact split is a positive evolved-connection contribution (+5.20e-3, +9.57e-3) plus a larger negative remaining metric-route/discretization difference (-8.83e-3, -1.643e-2). The two components partly cancel in the total curvature gap. Neither term alone explains the whole residual.
- At the first center cell, the metric-derived curvature route reduces H substantially (N=80: -0.135 to -0.038), but that cell is excluded from the admitted constraint stencil. The admitted i=2 comparison is more relevant.
- The same-state H gap equals the Ricci gap to roundoff because all other source and extrinsic-curvature terms are unchanged in the counterfactual. This validates the accounting split; it does not show that metric-derived Ricci should replace BSSN Ricci in production.

**Next diagnosis:** inspect per-stage geometry/connection evolution and derivative/operator consistency around the first admitted inner cells. Do not replace the production Ricci equation or add a center projection just because either counterfactual reduces a residual.

The determinant constraint remains small in the unprojected trajectory (roughly 2.8e-9, 4.4e-10, and 7.7e-10 for N=40,60,80), and CMC residuals are small, but that does not compensate for the central Hamiltonian growth or establish full constraint satisfaction.

## CMC target counterfactual and causal limits

R80 with jointly prescribed archived R160 outer-lapse boundary history and full-domain CMC Kdot target completed through t=12:
https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450

It reduced the final proper-time discrepancy relative to archived R160 from about +1.1743e-2 to +6.30e-7 and maximum saved-slice discrepancy from about 1.5702e-2 to 1.2192e-5. This is a **joint gauge/boundary counterfactual**, not a self-consistent R160 solution on R80. It does not isolate Kdot.

The separate frozen-state Kdot target swap (https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205934746) did not recover the R160 central lapse at t=12; its residual difference was about -3.793e-3. Therefore the near-match in the evolved run cannot be assigned to Kdot alone. If needed, separately evolved counterfactuals must vary outer-lapse boundary and Kdot independently.

The interpolation audit passed 24 held-out times with maximum absolute Kdot interpolation error 5.29e-8; the archived outer-lapse history has a larger reported held-out absolute interpolation discrepancy around 1.5e-4 over half-time-unit spacing.

## N=640 signed-pole audit

Run https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38183314088 was cancelled at the 240-minute limit near t=8.375, short of t=12. Its partial artifact uploaded. Early sign crossings disappeared by about t=5.125; a later residual ridge grew afterward. Term reconstruction closes to roundoff, but this is not independent constraint satisfaction.

## Next gates

1. Use the completed same-state Ricci decomposition to trace how the geometry/connection residual first grows across a production step, starting from t=0 and then an early evolved checkpoint. Keep the matter/scalar equations unchanged and preserve signed H/connection witnesses at every stage.
2. Use the geometric audit to identify whether the admitted inner H drift is tied to the connection/Ricci pathway. Only then design a single bounded numerical-layer intervention with preregistered pass/fail metrics.
3. Check matter and radiation boundary flux direction and local conservation terms from raw signed terms, not ratios whose denominators may be near zero.
4. Preserve frozen N=320 reconstruction tolerance 1e-8 and report host/runtime variance; do not relax thresholds after seeing results.
5. Do not launch the t=160 campaign until the t=12 fixed-domain constraint behavior and balance residual are reviewed.
6. **Completed:** the production t=12 workflow now uses per-branch concurrency cancellation and no longer triggers on edits to its own workflow YAML. It triggers only on engine, tests or dependency changes, preventing documentation/workflow-maintenance pushes from launching a fresh expensive campaign.

A passing unit suite or a completed run does not itself prove physical turnaround, ejection, or multi-cycle behavior. Numerical/gauge failures remain numerical/gauge results unless invariant, independent and resolution-consistent evidence supports a physical interpretation.
