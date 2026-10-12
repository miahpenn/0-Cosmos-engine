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

## New source-budget and single-step stage localization

### Signed source-budget / Shapley attribution

The post-processing source-budget report ran on the original admitted campaign checkpoints only; it executed zero production steps and made no source or physics edits.
- Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209752406
- Artifact: `v55-production-hamiltonian-source-budget`
- SHA-256: `358b40fd3c9337ee9010752ed2386a1496b84fb3b6ccfc0ad70f68f86dae9565`
- It contains 21 source-budget snapshots and 18 interval-by-interval signed Hamiltonian increments. Shapley component attribution reconstructs each ΔH to a maximum error of about `5.5e-17` (roundoff).

At the first admitted inner cell i=2:
- N=60, t=10→12: ΔH = `-2.7518e-3`. Component contributions: Lambda `-3.6971e-3`, a `+3.0591e-3`, b `-2.1217e-3`, X `-1.63898e-4`, matter `+1.7401e-4`; Aa and K contributions are each below `3e-6`.
- N=80, t≈9.99→12: ΔH = `-6.8203e-3`. Component contributions: Lambda `-6.6889e-3`, a `+2.2994e-3`, b `-2.3291e-3`, X `-3.2892e-4`, matter `+2.3379e-4`; Aa and K contributions together are only a few `1e-6`.
- These same-cell contributions show that the matter-source increment is much smaller than the total drift and that the Lambda/connection pathway is strongly implicated in the evolving *constraint residual*. The components are a nonlinear same-state attribution, not independent causes, and do not yet prove why the Lambda equation departs from the metric connection.

### One-step production-stage trace

The trace ran exactly one production step from each of four states: N=60/N=80 production initializer at t=0 and saved N=60/N=80 states near t=10. No campaign-length evolution was performed.
- Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209727710
- Artifact: `production-one-step-constraint-stage-trace`
- SHA-256: `0c35e93f0279837bb4fb95c2a036c97388b88c59c6072cac485af30c83dbedba`

At N=60, the saved t=10 state (dt=0.04) has i=2 H `-7.96937e-3`. The explicit predictor first moves it slightly toward zero (`-7.95412e-3`); the primary update then moves it to `-8.01458e-3`. During that primary update, a,b,X are unchanged, the metric-derived H counterfactual remains approximately `-5.2008e-3`, while Lambda moves from `-1.16352e-3` to `-1.18921e-3` and the connection residual changes from `-7.63333e-3` to `-7.65902e-3`.

At N=80, saved near t=9.99 (dt=0.03), i=2 H moves from `-1.280465e-2` at input to `-1.280097e-2` at the explicit predictor, then `-1.288246e-2` at the primary stage. During that primary update, a,b,X remain fixed, the metric-derived H counterfactual changes by less than `8e-8`, Lambda changes by `-2.57535e-5`, and the connection residual changes by the same amount (`-1.2237166e-2` to `-1.2262919e-2`). The returned final state is `H=-1.288321e-2`.

At t=0, the one-step control keeps the connection residual near zero (final i=2 values `1.96e-7` at N=60 and `1.36e-7` at N=80). Thus the defect is not injected as a large jump by the first step; it accumulates through evolution. At the later saved state, one sampled step localizes a significant added H defect to the primary `Aa/K/Lambda` update, specifically the Lambda/connection-vs-metric-Ricci channel because the metric-derived H remains essentially fixed there. This is a stage-level numerical association, not yet proof that the continuum Lambda equation is wrong: intermediate PIRK stages are not independent physical slices.

### Next targeted audit

Before touching production evolution, compare the full-step change in `Lambda - Lambda_metric(a,b)` against the change predicted by the metric connection from the explicit a/b rates and the pinned PIRK Lambda RHS, cell-by-cell at the initial and t≈10 states. Report the separate `lambda_l2`, `lambda_l3`, and metric-connection derivative contributions. This should test whether the discrete connection constraint is propagated consistently by the existing split scheme. Keep all physics terms and coefficients unchanged.

## Lambda/connection propagation and time-step audit

The follow-on diagnostic series now ties the connection-constraint growth to the evolving momentum-constraint residual.

- Saved-state Lambda RHS audit (21 states, no evolution): https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210572528
- One-step Lambda connection budget: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210399733
- One-step stage trace: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209727710
- Matched CFL refinement: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38212493285

The analytic time derivative of the metric-derived connection agrees with an independent centered directional-difference check to roughly (10^{-11}) or better. The production step's actual Lambda increment agrees with the declared split `lambda_l2 + lambda_l3` to about (6	imes10^{-17}). At t=0 the vector-field difference between Lambda RHS and the metric-connection time derivative is near roundoff; a single production step creates only an (O(10^{-7})) connection residual at the first admitted inner cell.

At t=12, the momentum-residual term in the exact Lambda equation predicts nearly all of the connection-residual rate at the first admitted inner cell:
- N=60, r=3.333: momentum residual (-9.8167	imes10^{-4}); ((2alpha/a)M approx -1.9483	imes10^{-3}); measured (dot C_Lambda approx -1.9255	imes10^{-3}).
- N=80, r=2.5: momentum residual (-2.6894	imes10^{-3}); ((2alpha/a)M approx -5.3380	imes10^{-3}); measured (dot C_Lambda approx -5.2978	imes10^{-3}).

The match is within 1.2% and 0.8%. Current working diagnosis: the Lambda residual is predominantly **tracking the growing momentum constraint**, not independently proving the Lambda equation defective.

A matched fixed-grid CFL refinement at N=80, r_max=80, through t=2.1 (CFL .03 and .015) passed strict campaign completion. At t=2.1 the first admitted inner cell has (C_Lambda=3.5556	imes10^{-4}) vs (3.5580	imes10^{-4}), and H=-(9.9268	imes10^{-4}) vs -(9.9255	imes10^{-4}). The fixed inner peak norms differ by only about 0.2% for (C_Lambda) and 0.12% for H. This argues against simple time-step truncation as the dominant cause of the early defect.

### Completed same-domain spatial refinement

Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38212879394
Artifact SHA-256: `181e807240a54126db39a35fdd22f865509b96027db0a724853512225dd8b6e8`

The bounded production run at N=40,60,80,120, r_max=80, CFL=.03, radiation on, t=2.16 passed strict completion at every resolution (36,54,72,108 accepted steps; no tiny terminal steps).

Spatial findings at t=2.16:
- Fixed exterior `max|H|` for r>=5: 4.0086e-4 (N40), 1.6886e-5 (N60), 8.9180e-6 (N80), 5.4585e-6 (N120).
- Volume-weighted RMS H over the fixed annulus 2<=r<5: 7.4135e-4, 6.5404e-4, 4.4980e-4, 1.2832e-4.
- At fixed r=3, signed H is +7.4135e-4, -1.8587e-4, -4.1919e-4, +9.2175e-6 respectively; at r=5 it falls from -4.0086e-4 to -5.4585e-6 by N120.
- The maximum in admitted cells (index>=2 and r<5) moves inward and grows from |H|=7.51e-4 at r=3.333 (N60), to 1.071e-3 at r=2.5 (N80), to 1.539e-3 at r=1.667 (N120). The full inner peak, including the two cells excluded by the constraint stencil, decreases from 9.27e-3 (N60) to 6.97e-3 (N80) to 4.74e-3 (N120).

The fixed-radius and volume-weighted metrics improve, but the admitted-inner maximum follows the moving centre and does not converge in the sup norm. Working characterization: a narrowing central numerical layer; pointwise central regularity remains unadmitted. N40 has no cells with index>=2 and r<5, so it cannot provide an admitted-inner max for this mask.

Combined with the time-refinement result, the evidence does not support a simple CFL error or a broad-domain failure. It points to a centre-localized spatial/momentum-constraint propagation problem. Preserve physics; the next diagnostic decomposes the momentum-constraint residual term-by-term at each production stage. Do not add damping, center projection, a custom Lambda equation, or a physical source as a shortcut.

## Momentum-constraint term budget and semi-discrete tangent (2026-10-12)

The central momentum residual has now been audited as an exact sum of signed geometric derivative terms and separately projected scalar/DM/baryon/radiation currents at all 21 initial/checkpoint states, followed by an RHS tangent test and one-step comparison.

Saved-state term budget:
- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38213685916
- Artifact SHA-256: `0db11b923cd5ab5671d5da95739744c5b6fc6d1561de990a2b6e3c95c26324ee`
- Geometry, current assembly, and total momentum residual closures all pass at roundoff.
- At t=12/i=2, N60 has M_geom=-9.71025e-4 vs M=-9.81671e-4; scalar current term -8pi*j_S=-1.082e-5 and DM about +1.55e-7. N80 has M_geom=-2.68327e-3 vs M=-2.68944e-3; scalar term -6.317e-6 and DM about +1.31e-7. Late M values are dominated by geometry, not a large fluid-current source.

Semi-discrete tangent:
- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38213760419
- Artifact SHA-256: `d00310158e8d509f34d4ade1c347e2b5fa4542565f51b749e9c781e0b94d862f`
- On six frozen states (N60/N80 at t=0, t≈10, t=12), the tangent uses the actual production RHS, perturbs by ±epsilon, re-solves the CMC lapse, and finite-differences every signed geometry/source term. Two epsilons (2e-5 and 1e-5) agree for dM/dt within absolute 1e-12; geometric/source assembly closes at about 1e-19 or smaller.
- At t=0, M starts at roundoff but the RHS predicts dM/dt=+1.21132e-4 (N60) and +1.48735e-4 (N80). One actual production step yields +1.21228e-4 and +1.48716e-4 respectively, agreeing within 0.1%. This means the residual tendency exists in the semi-discrete coupled RHS; it is not generated solely by finite-dt error.
- At N60/t≈10, dM/dt=-3.21262e-4 vs actual step -3.12853e-4; at N80/t≈10, +2.95895e-3 vs +2.93569e-3. At t=12, N60 is near cancellation (-8.04e-6 predicted vs -6.34e-6 actual), while N80 is -1.47313e-3 predicted vs -1.37845e-3 actual.
- On the initial slice, the source-current tendency is overwhelmingly the local S-field current: d(8pi*j_S)/dt=-2.42151e-4 (N60) and -2.23599e-4 (N80); D, COSMOS-phi and all fluid-current rates are orders smaller. Geometry contributes -1.20974e-4 / -7.48267e-5 respectively, leaving the measured nonzero dM/dt. At N80/t≈10 the geometric rate (+3.06906e-3, mainly D_r Aa and radial-b terms) dominates the current-source rate (+1.10102e-4).
- Interpretation: no single scalar or geometry term is isolated as the sole error. The defect is a mismatch in the *balance* between the semi-discrete geometric momentum constraint and evolving source current. The fixed-grid CFL refinement hardly changes the early defect, while the N=40/60/80/120 t=2.16 refinement shows a narrowing central spatial layer.

### Spatial refinement of the semi-discrete momentum tangent — completed

Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38214065011
Artifact SHA-256: `6c19e84cb197dbd4ed248791e8a06bcc42bd982339d196ea380cf9c58a43884b9`

The RHS tangent audit sampled the initial slice and saved t=2.16 states for N=40,60,80,120. It used the actual production semi-discrete RHS, plus/minus directional perturbations with CMC lapse re-solved, two epsilons, and an actual one-step comparison. It did not advance a trajectory beyond one local step.

At fixed radii r=3,4,5,6,8, the initial-slice ẊM becomes small as resolution increases: for example at r=3 it goes -4.80e-4 (N40), -9.59e-5 (N60), +7.59e-5 (N80), +1.91e-7 (N120); at r=5 it decreases from +8.77e-5 (N40) to +4.84e-7 (N120). At t=2.16, the N120 tangent at r=3,4,5,6,8 is -7.81e-5, -2.01e-5, +1.65e-6, -8.25e-8, +9.25e-8. Coarse resolutions remain non-monotonic and differ in sign at some fixed radii.

The inner r=2 tangent remains the exception: it is non-monotonic and strongly resolution-sensitive at both t=0 and t=2.16. Thus the defect narrows spatially, but pointwise central regularity has not been admitted.

**Operational decision:** keep all production physics, scalar/matter sources, CMC gauge, and boundaries locked. Do not merge PR #3 yet and do not launch t=160. Next diagnostic is an operator-level check of regular-centre parity and discrete derivative compatibility in the first cells, including the directional derivatives of Aa, K, X, a and b and the independently decomposed matter-current tangent. No centre projection, damping, fitted coefficient, or custom Lambda equation should be introduced as a shortcut.

## Centre parity / derivative-operator audit (2026-10-12)

Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38214547118
Artifact SHA-256: `dc4959fc117cb428ffcfc180db143b0c74adfd7d60f5e9db0d356361c7191981`

The pinned cell-centred fourth-order derivative was verified against manufactured even/odd regular polynomials at N=40,60,80,120. First-five-cell first- and second-derivative errors are at floating-point scale; the worst is about 2.0e-15. No basic parity sign or polynomial-exactness bug was found in the generic operator.

A frozen-state local even-in-r polynomial probe changes the t=2.16 momentum residual substantially in r<=5, but it is not evidence that the production stencil is wrong: the fit is a diagnostic alternate derivative, its evolved-field fit residuals reach 1e-3 to 2e-3 (especially for K), and its support contains different counts of cells at different N. It must not be promoted into production.

This audit therefore narrows the issue but does not resolve it: production profiles develop a centre-localized spatial layer, and the connection/momentum constraints are sensitive to derivative reconstruction there. The next test is direct evaluation of regular-centre series ratios `Aa/r^2`, `Lambda/r`, `(a/b-1)/r^2`, and `D_r f/r` for the even fields `K,X,a,b,S,D,phi` at t=0 and t=2.16. This will compare regularity at the first cells across resolutions without imposing a replacement stencil.

## Latest centre-series and exact per-step momentum audits (2026-10-12)

### Centre-series ratio audit

Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38214837536
Artifact SHA-256: `d5af2d743068000a354506eba3dc34ff1be5f9a9fba6482727045e5178be43d4`

At t=0 the centre ratios are bounded and `Aa/r²`, `D_r K/r` are roundoff-close to zero. The first-cell `Lambda/r` is about 2.85e-3 (N40), 2.97e-3 (N60), 3.02e-3 (N80), 3.05e-3 (N120), but it drops by roughly a factor two from first to second cell across every N. The production `C_lambda` is zero by construction on the initializer.

At t=2.16 the first cell becomes strongly resolution-sensitive:
- `Lambda/r` first-to-second cell: N40 3.52e-3 to 1.28e-3; N60 4.49e-3 to 1.78e-3; N80 5.70e-3 to 2.23e-3; N120 9.06e-3 to 3.26e-3.
- `C_lambda/r` first cell: +1.64e-3 (N40), -2.99e-3 (N60), -1.95e-3 (N80), +3.59e-3 (N120), with signs changing across resolutions.
- `Aa/r²` first cell changes sign with resolution: -7.72e-4 (N40), -9.96e-4 (N60), +1.029e-3 (N80), +4.429e-4 (N120).

The direct ratio audit supports a first-cell regularity layer that grows relative to the adjacent-cell behavior as the mesh is refined. It does not find a generic parity-stencil implementation error: the manufactured even/odd polynomial checks pass to ~2e-15. Do not replace the production stencil or apply a fitted centre projection.

### Full replay with per-step signed momentum budget

Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38214554266
Artifact SHA-256: `223901a0fbd2ce4ddbfb39cce2cb2b3c88eb24286f7d249c62cc3ecad2f9b569`

The N80 observer replayed all 400 accepted production steps through t=12 and reproduced the archived endpoint exactly field-by-field, including proper time. Constraint-term closures remain roundoff-level.

The first admitted inner momentum residual at r=2.5 is **oscillatory, not monotonic**: zero crossings occur near t=1.86, 3.93, 4.92, 7.08, 8.01, 10.26, 11.07, with prominent negative excursions near t≈3,6,9.12,12 and partial recovery between them. Over the complete record, RMS is about 1.083e-3 for M, 9.119e-3 for H, and 6.816e-3 for C_lambda. M tracks the geometric momentum term strongly (correlation 0.9975); scalar-current RMS is about 6.49e-5.

At t=12, i=2, the geometric M terms are:
- D_r Aa: +2.6496e-4
- -(2/3)D_r K: -6.3120e-4
- 6 Aa chi_r: +1.0166e-5
- radial 2(Aa-Ab)/r: -2.30978e-3 (dominant)
- b-gradient term: -1.74225e-5
- M_geom: -2.68327e-3; total M residual: -2.68944e-3.

The sparse two-time-unit tangent comparison must not be read as instantaneous causality because the trajectory reverses several times between those samples. The per-step observer is the more faithful reference for temporal variation. Current diagnosis remains: a centre-localized spatial/derivative-compatibility issue, with a strongly oscillating constraint mode, rather than a simple CFL error or an oversized matter-current source. This remains a numerical diagnosis; no bounce/ejection claim is admitted.

## Next gates

1. Compare first-cell regularity-series evolution rates implied by the production RHS against the ratios now recorded at t=0/t=2.16: Aa/r², Lambda/r, (a/b-1)/r² and D_r f/r for K,X,a,b,S,D,phi. Use the exact per-step N80 replay to select signed intervals around the large M excursions. Keep this as an operator-level diagnostic; no modified trajectory.
2. Use the geometric audit to identify whether the admitted inner H drift is tied to the connection/Ricci pathway. Only then design a single bounded numerical-layer intervention with preregistered pass/fail metrics.
3. Check matter and radiation boundary flux direction and local conservation terms from raw signed terms, not ratios whose denominators may be near zero.
4. Preserve frozen N=320 reconstruction tolerance 1e-8 and report host/runtime variance; do not relax thresholds after seeing results.
5. Do not launch the t=160 campaign until the t=12 fixed-domain constraint behavior and balance residual are reviewed.
6. **Completed:** the production t=12 workflow now uses per-branch concurrency cancellation and no longer triggers on edits to its own workflow YAML. It triggers only on engine, tests or dependency changes, preventing documentation/workflow-maintenance pushes from launching a fresh expensive campaign.

A passing unit suite or a completed run does not itself prove physical turnaround, ejection, or multi-cycle behavior. Numerical/gauge failures remain numerical/gauge results unless invariant, independent and resolution-consistent evidence supports a physical interpretation.
