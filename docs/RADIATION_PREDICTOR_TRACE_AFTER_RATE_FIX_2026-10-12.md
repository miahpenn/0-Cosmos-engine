# Radiation predictor failure after matter-RHS rate-extraction correction

**Date:** 12 October 2026  
**Classification:** Diagnostic capture; production admission remains open.  
**Decision:** The former high-resolution early-step radiation RHS false failure was a unit-increment guard artifact. A separate, real late-time predictor-metric cone crossing remains after rate extraction was corrected. Do not clamp or relax the radiation invariant. No production/main/GEAR merge is authorized by this evidence.

## Provenance

### Corrected-rate, fixed-production-dt trace

- Workflow: [38310216296](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310216296) — completed successfully as a trace-capture workflow; all repository tests passed.
- Artifact: [radiation-trace-after-rate-fix, artifact 11774592255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310216296/artifacts/11774592255)
- Artifact ZIP SHA-256 (verified locally): `d2f0cdc037701f7b03cb76b26b5ccefa935a77fadfbd62b23969328ea7bb7c53`
- Embedded `trace.json` SHA-256: `1a34e0a46b5be64e7133dae5f835ac13fa478bf68ad121699f7c6f11cfd1c9aa`
- Embedded `summary.json` SHA-256: `886a05e478375074766f6193b6e7d403c987d9ec310b709dd9f2d153fa9930cf`
- Branch: `diag/radiation-trace-after-rate-fix-2026-10-12`
- Trace head: `d7a0c975800d62965df59bb2dba7a807bc39a6af`
- Corrected rate-only source base: `1c63090d4f42396f899421a2379f70320689d5ad`
- Pinned vendor submodule: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- The fail-closed workflow proved the adapted diagnostic file was the only allowed `engine/` difference from the corrected rate-fix source; the vendor pin matched. Python 3.12.15, NumPy 2.5.4. Regression suite, compile, trace capture, hash collection, and artifact upload all succeeded.
- Configuration: N=80, Rmax=80, dr=1, CFL=0.03, fixed dt=0.03, D amplitude 10^-10, radiation ON, scalar profile amplitude 0.01 and width 7. Target t=50; no admissibility retries, accepted-metric recovery disabled.

### Independent raw-geometry stage ledger

- Workflow: [38309949383](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309949383)
- Artifact: [radiation-predictor-metric-stage-ledger, artifact 11773874630](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309949383/artifacts/11773874630)
- ZIP SHA-256 (verified locally): `356f067a583ec3596c3b15e8b5d17a18058cc7bd7a2bbc7fd379b6b1e2acbe93`
- Stage-ledger JSON SHA-256: `d2e9d507f8300598c4dd5ac45b4bcd3f77d205adc388e460e6a75ae614a214c9`
- This ledger uses a different frozen historical source baseline, `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`, rather than the corrected rate-fix source. Its exact values are corroborating stage-order evidence, not a matched numerical A/B. It validates 100 ordered ledger records and reports the same failing time/cell.

## Corrected-rate trace: observed event

The workflow captured 1,566 accepted steps and 66 same-cell radiation predictor budgets, with no rejected-step retries or instrumentation errors. The trace did not reach its target time.

- Last accepted coordinate time: t = 46.98000000000102.
- Central proper time at that state: tau = 22.016046893100416.
- Next explicit predictor attempt: t = 47.01000000000102, dt = 0.03.
- First failure stage: `explicit_predictor_first_CMC_lapse`, entering `V55ProductionKernel._solve_lapse`.
- Cell: i = 79, r = 79.5, the outermost cell.
- Failure: exact radiation cone test E ≥ |S| fails. Logged physical values are E = 3.9783663187660695 × 10^-8, |S| = 3.9785943919623527 × 10^-8, and |S|/E = 1.0000573283549097. The densitized cone margin C = U_E − sqrt(gamma^rr) × |U_r| is −1.699741681811858 × 10^-8.

This is a real finite-dt predictor failure on the corrected-rate source, not the former dt=1 surrogate check. The exception is raised when the first predictor CMC lapse attempts to invert the predictor radiation state, before that predictor’s first CMC solve returns. A successful workflow means expected failure capture and artifact validation succeeded; it is not an evolution pass.

## Exact call path at failure

The stage name can be misleading if read as “the CMC linear solve failed.” On the pinned candidate’s source path, `V55ProductionKernel.step` first forms the explicit predictor geometry (a, b, X) and the first matter predictor from their accepted-state RHS values, then calls `_solve_lapse(grid, g_explicit1, s1, m1)`. Inside `solve_cmc_lapse`, the first operation after building the metric is `total_matter_projection(...)`; that sums the scalar and fluid stress projections. The fluid projection calls radiation primitive recovery, which enforces the exact conservative cone. The failure is thrown there, while constructing the matter source supplied to CMC, before the lapse PDE is assembled/solved to completion. See the pinned [production predictor stage](https://github.com/miahpenn/0-Cosmos-engine/blob/1c63090d4f42396f899421a2379f70320689d5ad/engine/production_kernel.py#L243-L275), [CMC source assembly](https://github.com/miahpenn/0-Cosmos-engine/blob/1c63090d4f42396f899421a2379f70320689d5ad/engine/cmc_gauge.py#L52-L85), and [radiation primitive admissibility check](https://github.com/miahpenn/0-Cosmos-engine/blob/1c63090d4f42396f899421a2379f70320689d5ad/engine/matter_system.py#L145-L183).

The explicit predictor stage does not call the outer light-constraint boundary before this CMC source assembly; the boundary operations occur later in the primary predictor and completed-state sequence. So the boundary is **not a direct operation at the failing substage**. However, the current accepted geometry is the output of prior completed steps, where the boundary has acted on the outer A_a state. Since A_a directly sets the next step's a-dot at zero shift, earlier boundary influence cannot be excluded from this stage-level finding.

## Same-cell cone-margin budget

The recorded budget separates (1) accepted state plus accepted metric, (2) flux-only matter increment evaluated on the accepted metric, (3) geometric-source increment evaluated on the accepted metric, and (4) the updated conservative state assessed on the predictor metric.

For cell 79 in the failing attempt:

| Quantity | Recorded value |
|---|---:|
| Accepted cone margin C₀ | +4.729090141987707 × 10^-7 |
| After flux-only increment, accepted metric | +1.484032098292450 × 10^-6 |
| After flux plus geometric source, accepted metric | +2.219759745727783 × 10^-6 |
| Metric-change contribution to C | −2.236757162545902 × 10^-6 |
| Final predictor cone margin | −1.699741681811858 × 10^-8 |
| Budget closure error | 0 |

At this cell the conservative energy flux changes U_E by −4.1564015995 × 10^-7, and radial momentum U_r by −1.1796291609 × 10^-6. The geometric source then adds +5.7872482608 × 10^-6 to U_E and +4.1765310727 × 10^-6 to U_r. Taken together, flux plus source leave a positive margin on the accepted metric. The same predictor conservative values assessed on the updated metric cross the cone. The algebraic stage budget closes exactly under its defined order, and predictor reconstruction residuals for U_E and U_r are zero.

This localizes the immediate trigger: the predictor spatial metric, not the unit-time radiation RHS wrapper, produces the final crossing in this step. It does not identify a unique underlying term in the geometry RHS as the ultimate defect and is not proof of a missing physical law.

## Geometry update decomposition

The rate-fixed trace records the following at cell 79:

| Cell-79 metric quantity | Value |
|---|---:|
| Accepted a | 0.6089966857200283 |
| Explicit RHS for a | −0.1921657198214381 |
| a after Euler update | 0.6032317141253851 |
| Accepted X | 0.9438736828908237 |
| Explicit RHS for X | +0.08873930174627159 |
| X after Euler update | 0.9465358619432119 |
| gamma^rr = X²/a, accepted | 1.462893887838753 |
| gamma^rr, after a update | 1.476874488513230 |
| gamma^rr, after X update | 1.485217234050057 |
| Change from a update | +0.013980600674477 |
| Change from X update | +0.008342745536827 |
| Change from centre-regularity projection at this cell | 0 |
| Decomposition closure error | 0 |

The explicit decrease in a, followed by increase in X, both increase gamma^rr = X²/a. Together they account for +0.022323346211304 in gamma^rr. The centre-regularity projection leaves a, b and X unchanged at this outer cell, so it is not the direct cause of this crossing. This is limited to the recorded failing cell and step; it does not exclude earlier boundary influence or all regularity effects elsewhere.

### Exact origin of the explicit metric rates

The pinned numerical operator states the explicit geometry equations in [`vacuum_bssn.py` at vendor commit `d6052d6`](https://github.com/RonBibb/bb-palatini-unified-r0/blob/d6052d605673ce9d82cdc99b0df79855fa2ea215/vacuum_bssn.py#L268-L295):

- `a_dot = beta*a_r + 2*a*beta_r - (2/3)*a*div_beta - 2*alpha*a*Aa + KO[a]`
- `X_dot = beta*X_r + X*(alpha*K - div_beta)/3 + KO[X]`

The captured outer cell has beta = 0; the explicit operator uses its default KO strength of zero, so the actual rates reduce to `a_dot = -2*alpha*a*Aa` and `X_dot = alpha*X*K/3`. Substitution of the recorded accepted values alpha = 1, Aa = 0.1577723855707334 and K = 0.2820482338520797 reproduces the saved rates `a_dot = -0.1921657198214381` and `X_dot = +0.08873930174627159`. This is the expected extrinsic-curvature-driven explicit BSSN metric update at zero shift, not evidence that a missing direct stress-energy term belongs in these two explicit equations. The local geometry-stage adapter separately places total-stress-energy terms into the Aa/K primary RHS blocks; that upstream coupling and the constraints still need their own admission.

To isolate the metric factors without modifying any state, I evaluated the *same predictor matter values* on counterfactual metrics formed by changing only one of a or X:

| Metric used for the same predictor radiation state | gamma^rr used | Cone margin C | Admissibility |
|---|---:|---:|---|
| Accepted metric, no explicit metric update | 1.462893887838753 | +2.2197597457 × 10^-6 | Admissible |
| a updated; X held at accepted value | 1.476874488513230 | +8.1694973651 × 10^-7 | Admissible |
| X updated; a held at accepted value | 1.471157658083646 | +1.3897691933 × 10^-6 | Admissible |
| Both explicit updates applied | 1.485217234050057 | −1.6997416818 × 10^-8 | Inadmissible |

Thus neither metric-field change alone crosses the cone in this one step; their combined increase in gamma^rr does. Using the recorded a-then-X decomposition, the a update consumes 1.4028100092 × 10^-6 of margin and the X update consumes a further 8.3394715333 × 10^-7. This ordering attribution is path-dependent because gamma^rr = X²/a is nonlinear, but the separate-field counterfactuals show the joint update is necessary for the observed crossing at this stage.

## Approach to the crossing

The last three fixed-dt attempts show the predictor margin shrinking toward zero:

| Accepted time | C on accepted metric | C on predictor metric |
|---:|---:|---:|
| 46.92000000000102 | +1.3010870 × 10^-6 | +8.7787474 × 10^-7 |
| 46.95000000000102 | +9.0393351 × 10^-7 | +4.4739351 × 10^-7 |
| 46.98000000000102 | +4.7290901 × 10^-7 | −1.6997417 × 10^-8 |

On the failing step, flux and source lift C to about +2.22 × 10^-6 on the accepted metric, but the predictor metric contributes about −2.24 × 10^-6. This confirms rapidly growing metric sensitivity at the outer edge for this frozen configuration. It does not yet distinguish geometry RHS discretization, constraint propagation, boundary-conditioned geometry, gauge/CMC coupling, or a deeper physical model requirement.

## Ten-step late-window trend from the same corrected-rate trajectory

The artifact contains ten *completed* production steps leading into the failing predictor, not just three final snapshots. On each completed step, the trapezoidal radiation flux and source increase the cone margin when evaluated on the accepted metric, yet the evolved spatial metric takes away more margin than the matter update adds. Over this ten-step window, the outer-cell inverse radial metric rises monotonically from 1.278974 to 1.462894 (+14.38%), while the accepted cone margin falls from 3.7051e−6 to 4.7291e−7 (−87.2%). This is a persistent accumulation, not a one-step unit-Δt artifact.

| Accepted interval | C at interval start | C after flux + source on accepted metric | C at completed state | Metric change contribution through RHS1 | Net completed-step ΔC |
|---|---:|---:|---:|---:|---:|
| 46.680 → 46.710 | +3.7051e−6 | +4.8577e−6 | +3.4502e−6 | −1.3931e−6 | −2.5488e−7 |
| 46.710 → 46.740 | +3.4502e−6 | +4.6642e−6 | +3.1852e−6 | −1.4639e−6 | −2.6495e−7 |
| 46.740 → 46.770 | +3.1852e−6 | +4.4629e−6 | +2.9091e−6 | −1.5381e−6 | −2.7609e−7 |
| 46.770 → 46.800 | +2.9091e−6 | +4.2526e−6 | +2.6207e−6 | −1.6156e−6 | −2.8849e−7 |
| 46.800 → 46.830 | +2.6207e−6 | +4.0317e−6 | +2.3184e−6 | −1.6964e−6 | −3.0230e−7 |
| 46.830 → 46.860 | +2.3184e−6 | +3.7988e−6 | +2.0007e−6 | −1.7805e−6 | −3.1767e−7 |
| 46.860 → 46.890 | +2.0007e−6 | +3.5506e−6 | +1.6646e−6 | −1.8678e−6 | −3.3610e−7 |
| 46.890 → 46.920 | +1.6646e−6 | +3.2774e−6 | +1.3011e−6 | −1.9574e−6 | −3.6350e−7 |
| 46.920 → 46.950 | +1.3011e−6 | +2.9721e−6 | +9.0393e−7 | −2.0486e−6 | −3.9715e−7 |
| 46.950 → 46.980 | +9.0393e−7 | +2.6344e−6 | +4.7291e−7 | −2.1414e−6 | −4.3102e−7 |

The C values are densitized cone margins at cell i=79. The metric contribution shown is the accepted-to-RHS1 stage change; the completed-stage correction is smaller but also negative (about −1.4e−8 to −2.0e−8 in these rows). All saved completed-state budgets report zero algebraic closure error. Thus the matter update is not what drives the late-window *loss of margin* in this discrete trace; the repeatedly changing spatial metric is. That still does not make the metric RHS or the physical geometry “wrong” by itself: the radiation state at i=79 is already much closer to null than its neighbors, and a valid continuum solution could in principle approach that boundary. The failed predictor means the current discretized stage cannot represent that attempted next state while preserving the exact radiation primitive condition.

The same last-step trace gives another useful nuance: the *fractional* rise in γ^rr on cells i=75–78 (about 2.06–2.81%) is larger than on i=79 (about 1.53%), but those neighboring radiation states keep C comfortably positive. The outermost cell crosses because the margin there is orders smaller than adjacent-cell margins; this is a coupled near-null-state / metric-stage event, not evidence that the largest fractional metric update alone predicts failure. The measured final predictor margins are i=75: +4.7428e−6, i=76: +5.0984e−6, i=77: +5.7824e−6, i=78: +8.6027e−6, i=79: −1.6997e−8.

## Cross-check against the older stage ledger

The historical-source ledger also captured t = 47.01, cell 79. In that run gamma^rr rose from 1.4628939 to 1.4853572 across the raw Euler geometry update. The margin remained negative after the centre-regularity stage (C = −1.6900205 × 10^-8, ratio = 1.0000569941); the raw metric already crossed the cone and projection only improved it slightly. Because the old ledger used source base `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`, its exact numbers must not be treated as a paired comparison with the corrected-rate trace. It corroborates the stage ordering: the predictor metric is inadmissible before the lapse rejection.

## Duplicate trigger cross-check

A workflow-only label/comment commit queued a second run of the same frozen-dt trace; this was redundant work, not a new candidate or parameter change. It completed with the same 1,566 accepted steps and the same failure at t = 47.01000000000102, cell 79, r = 79.5.

- Repeat run: [38310237842](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842), same corrected source base `1c63090d4f42396f899421a2379f70320689d5ad` and pinned vendor.
- Repeat artifact: [artifact 11775230036](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842/artifacts/11775230036)
- Repeat artifact ZIP SHA-256: `f64298e2819e0dd9c4c2df712a4d1250e4855a5b0e467430fe7aa419a533e49f`
- Repeat `trace.json` SHA-256: `53a1f9f188e0eea4575da4458d09891f50144b9f89f945c995cd9980221fb619`
- Repeat failure ratio: 1.00005732838577255 versus 1.00005732835490968 in the first capture; predictor margin −1.6997425987812665 × 10^-8 versus −1.699741681811858 × 10^-8.
- Repeat central proper time differs by only about 1.5 × 10^-9, and the computed predictor inverse metric by about 3.6 × 10^-9.

This is an independent trigger of essentially the same result, not a scientifically distinct ablation. Do not commission another repetition; next work should be offline analysis of the captured same-step geometry/source and constraint budgets.

## Existing expanded-domain boundary audit (historical source; not a matched A/B)

Before drawing any conclusion that the outer metric behavior is independent of boundary treatment, I inspected the archived [proper-time/domain control run 38028306184](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184), [artifact 11661332148](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184/artifacts/11661332148). Its ZIP SHA-256 is `f4f06a54bf66f0e5eed597efe0b6b0e605e31375edd7303bcb4edd4cdb97fac5`; its trace SHA-256 is `f418a023c44dcd51802696e6468e3e09691d01ea44ecbb5a9463881be6b14e76`. It is the expanded N=120, Rmax=120 run at an older frozen source baseline (`a7170c12ed08b4fdb3c52ca98f02e106a3f39b82)), so all numerical statements here are historical cross-checks, not a quantitatively paired comparison with the current N=80/rate-fixed trace.

That older run captured 174 pre/post boundary updates, with no boundary-instrumentation errors. The post-boundary P-minus target mismatch maxed at 1.0408341 × 10^-17 and the A_a correction-formula error at 1.2397586 × 10^-17, so the implemented correction follows its declared algebraic target to roundoff. The maximum pre-boundary target mismatch was 1.1540766 × 10^-4 at t=60.99000000000155; the completed-state boundary changed outer A_a from 0.1187022562 to 0.1188207695 (delta +1.1851332 × 10^-4). Corrections remained at +6.3026 × 10^-5 at t=61.02000000000155 and +3.2333 × 10^-5 at t=61.03500000000155, then collapsed as the adaptive retries drove dt toward floating-point underflow. At the eventual failure, the accepted state was t=61.046503364247606, central tau=21.9130694389674, and the cone failure was at the new outermost cell r=119.5. Thus expanding the domain relocates the observed failure outward but does not eliminate it.

This does **not** prove the P-minus boundary law is incorrect: the residual audit shows it achieves its stated target. It also does not prove the boundary is innocent: the imposed outer value changes A_a, and that value feeds the next explicit metric RHS. The current N=80 rate-fixed artifact does not record same-step pre/post boundary states, so exact boundary attribution at t≈46.98 is unmeasured. Keep this as an open coupling question rather than a production fix or an exonerated layer.

## Boundary/gauge compatibility: what the source itself guarantees

The boundary reconstruction used by V5.5 is not part of the archived GEAR physical law. It comes from the pinned numerical operator [`bssn_characteristic_boundary.py` at vendor commit `d6052d6`](https://github.com/RonBibb/bb-palatini-unified-r0/blob/d6052d605673ce9d82cdc99b0df79855fa2ea215/bssn_characteristic_boundary.py). The module's own top-level note says it does **not** implement a complete constraint-preserving boundary; it isolates characteristic/gauge behavior for a bounded diagnostic. In the production path, `apply_light_constraint_boundary` sets only `Aa[-1]` from a frozen second-order radial reconstruction of the incoming light/constraint characteristic (P_-). Separately, the production CMC solve makes `alpha[-1] = 1` its outer normalization row rather than imposing the lapse PDE at that cell.

This matters for interpreting the old proper-time-matched boundary artifact. Its post-boundary residual reaches ~1e−17 because the routine has just solved for `Aa[-1]` to make the prescribed target residual zero by construction:
`target = omega[-3] + 2*dr*source[-2]`, followed by the assignment of `Aa[-1]` that sets `omega[-1] = target`. The test also confirms the correction matches that same algebraic formula. Those are good implementation/accounting checks, but they are **not independent confirmation** that the chosen interior-extrapolated target approximates the correct continuum incoming datum, that it is compatible with this CMC/zero-shift gauge, or that all incoming constraint/gauge modes are controlled.

The historical R=120 audit's pre-boundary mismatch of up to 1.1541e−4 and `Aa[-1]` correction of +1.1851e−4 are therefore evidence of a material correction being applied, not a proof that the correction is either wrong or correct. The current rate-fixed R=80 capture does not record the exact same-state pre-/post-boundary `Aa` and characteristic values around its last accepted update. Because the current failing predictor's direct call path has not yet applied the boundary, we cannot assign this event to the boundary routine. Conversely, the earlier `Aa[-1]` values enter the next explicit `a_dot = -2*alpha*a*Aa` at zero shift, so historical boundary influence remains a plausible *numerical pathway* that has not been quantitatively isolated on a matched state.

There is an additional equation-level limit to this diagnostic: the implemented target is spatial,
`omega[-1] = omega[-3] + 2*dr*source[-2]`, and the reported residual is
`(omega[-1]-omega[-3])/(2*dr) - source[-2]`. It tests the finite-difference radial reconstruction. It does not by itself test the *time-dependent* incoming-characteristic evolution, a no-incoming-wave condition, or all constraint-preserving boundary modes. The pinned module explicitly calls the operation “frozen.” Because production fixes `beta=0` and solves `alpha` from a CMC elliptic equation (with `alpha[-1]=1)), while the vendor's literal principal matrix is formulated for a broader gauge subsystem, the exact characteristic completeness and compatibility of this frozen radial target in the production CMC system must be derived or independently admitted. Do not assume the frozen target is wrong; equally, do not treat its near-zero spatial residual as proof that the dynamic boundary problem is solved.

This is a gauge/constraint-compatibility gate, not a license to replace the boundary law. The archived GEAR blueprint and unified V5.5 architecture require one evolved spatial metric and total stress tensor, but they do not make the vendor's exact (P_-) target a physical law. Keep the following checks distinct: (1) residual closure of the coded reconstruction, (2) compatibility with the CMC outer normalization and current zero-shift evolution, (3) full incoming characteristic/constraint control, and (4) impact on radiation-cone and Hamiltonian/momentum/connection budgets. Only (1) has been demonstrated to roundoff here.

## Prior canonical R120 ledger: radiation face and CMC outer-layer clues

An older, canonical production-path R120 proper-time ledger (N=120, Rmax=120, CFL=0.03; source baseline predates the rate-extraction correction) adds a separate boundary clue. It must not be treated as a matched numerical A/B against the N=80 corrected-rate trace, but it helps set the next audit priorities:

- The terminal outer radiation characteristic speeds are both outward, with `radiation_outer_vhat≈0.9999999999999996`, `s_minus≈0.9999999999999984`, `s_plus≈0.9999999999999997`, and `radiation_outer_all_outgoing=1`. This weakens the simple claim that an inward radiation characteristic directly drives the endpoint.
- In the final accepted row, the outer-face HLL energy flux is 0.0019801631 and the inner face of the same outermost cell is 0.0024670764; the analogous radial-momentum fluxes are 0.0019186033 and 0.0021722137. The reported relative differences (19.74% in energy, 11.68% in momentum) are the **difference between two distinct faces of one finite-volume cell**. The production source explicitly computes the outer-cell divergence as `-(F_outer-F_inner)/dr), exactly these flux differences with the prescribed sign. They are not two estimates at the same face and therefore are **not a boundary-closure residual** or evidence of a boundary defect. A true closure check must compare admissible reconstructions at the same interface or an independently defined integrated conservation budget.
- The CMC outer normalization is `alpha[-1]=1`, while `alpha[-2]≈1.456444`, `∂r alpha[-1]≈−0.537565`, and the last-interior CMC residual is `−3.72689e−4`. Geometry/lapse values are strongly nonuniform.
- The largest explicit `a` RHS and vacuum `Aa/K` RHS magnitudes peak at (r=106.5), 13 cells inside the R=120 edge (r=119.5); the largest absolute Hamiltonian residual peaks at (r=89.5). So the episode contains a broad outer geometry/lapse structure, not merely a bad value in the last cell.
- The same legacy ledger's P-minus pre-boundary mismatch/correction was recorded in the separate boundary audit: the chosen formula closes algebraically to roundoff, but prior outer (A_a) is measurably changed by that condition.

These legacy observations are consistent with coupled near-null radiation, outer-cell flux divergence, and CMC/geometry sensitivity; they do not independently demonstrate a radiation boundary defect. They do not transfer to the corrected-rate source as exact values and do not identify one unique cause. Keep three pieces of boundary treatment distinct: (1) radiation outer-face flux/conservation at the same interface; (2) P-minus (A_a) reconstruction; and (3) CMC lapse outer normalization. A proper comparison must audit their own terms and the total constraint budget, not combine them into one generic “boundary bug.”

The local GEAR V5.5 exact-spherical prototype also uses a different outer treatment: its `fill_boundaries` copies the fourth-from-edge value into the last three positions under a comment describing smooth/asymptotically constant outer values for the short domain. That is not the same algorithm as the pinned frozen P-minus reconstruction and is not, by itself, proof that either choice is physically admitted. Treat the blueprint/operator difference as a question to reconcile, not as authority to copy either numerical boundary unexamined.

## Historical-source boundary correction measured from the admitted Aa/K stage ledger

The separate [Aa/K source-stage ledger run 38311614754](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38311614754) has completed and passed its fail-closed checks. Its [artifact 11774209967](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38311614754/artifacts/11774209967) has ZIP SHA-256 `989b51ffec4c9153d228ce3e8f80ef38873337985b3a4271ef904d2c1ff6b30b`; internal `stage_ledger.json` SHA-256 `2c56e08e0331b7b06dcf7017f53ffb47413eb647658c6a17bd075a8f902a7bb6`; `stage_ledger_summary.json` SHA-256 `1f2a749a6b8cb131a0d6d626d0f9f40380cad15501813911c3d744be306284fa`. The source/vendor pins passed, the complete test suite passed, all 100 captured steps had the required raw RHS/Euler/projection/first-CMC stages for cells 75–79, the signed L3 total−vacuum split closed exactly, and the completed-to-next-accepted handoff closed to (2.22\times10^{-16}). There were no stage or source-admission errors.

**Source qualification:** this recorder is frozen to historical production source `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`, with vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. It is not an A/B against the corrected rate-extraction source `1c63090d4f42396f899421a2379f70320689d5ad`. The terminal event is close in time and state but the artifact's terminal τ is `22.016115394562938` and its radiation conservatives differ slightly. Retain all exact figures under the historical-source provenance.

### Isolated measurement of the primary-predictor outer boundary correction

For each completed step, the saved accepted-state stage and explicit-predictor stage contain the exact (A_a) RHS pieces used by the production PIRK update. The pre-boundary primary predictor at cell (i) can therefore be reconstructed as

[
A_{a,\mathrm{pre}}^{(1)}=A_{a,0}+\Delta t\left(\tfrac12 L_{2,0}^{A_a}+\tfrac12 L_{2,\mathrm{pred}}^{A_a}+L_{3,0}^{A_a}\right).
]

At the outer cell this reconstruction differs from the captured post-boundary primary-predictor (A_a) only by the incoming-(P_-) boundary assignment; the centre-regularity projection does not alter (A_a) there. I evaluated all 99 completed primary-predictor stages in the admitted 100-step window.

| Accepted coordinate time | Reconstructed (A_a) before primary boundary | Captured (A_a) after primary boundary | Boundary-stage Δ(A_a) |
|---:|---:|---:|---:|
| 46.680 | 0.126391027654 | 0.126416414193 | +2.53865387e−5 |
| 46.710 | 0.129578770133 | 0.129607274080 | +2.85039467e−5 |
| 46.740 | 0.132840238989 | 0.132872035203 | +3.17962139e−5 |
| 46.770 | 0.136175182709 | 0.136210430679 | +3.52479700e−5 |
| 46.800 | 0.139583169121 | 0.139622003056 | +3.88339350e−5 |
| 46.830 | 0.143063575059 | 0.143106101470 | +4.25264105e−5 |
| 46.860 | 0.146615593583 | 0.146661885930 | +4.62923477e−5 |
| 46.890 | 0.150238229133 | 0.150288318479 | +5.00893457e−5 |
| 46.920 | 0.153930296253 | 0.153984160929 | +5.38646760e−5 |
| 46.950 | 0.157690420189 | 0.157747974203 | +5.75540139e−5 |

Across the five captured cells, the algebraically reconstructed correction is exactly zero at (i=75,76,77,78), and all of the boundary correction sits at the imposed outer point (i=79), as expected for the existing boundary routine that assigns only (A_a[-1]). Over the listed ten steps the correction grows by a factor of about 2.27. At (t=46.95), Δ(A_a=5.7554\times10^{-5}), about (0.0365\%) of post-boundary (A_a), and about (1.6\%) of the step's unconditioned primary-predictor (A_a) increment. This is not numerically zero, and it demonstrates a growing boundary-conditioned contribution in the same late window where the radiation cone is eroding.

### The same primary-boundary assignment shifts the discrete momentum-constraint term

This can be quantified without interpreting the boundary law. The pinned vendor's geometric momentum-constraint expression is

[
M_{\mathrm{geom}}=D_r A_a-\tfrac23 D_r K+6 A_a\chi
+\tfrac32 A_a\left(\tfrac2r+\frac{D_r b}{b}\right),
\qquad \chi=-\tfrac12\frac{D_r X}{X}.
]

The production diagnostic uses (M=M_{\mathrm{geom}}-8\pi j_{\mathrm{total}}). The boundary operation assigns only (A_a[-1]); it does not alter the metric or the matter state, so at the same primary-predictor slice the change in the *full* diagnosed (M) is exactly the change in (M_{\mathrm{geom}}). On the pinned fourth-order outer stencil, (Delta(D_r A_a)=25\Delta A_a/(12\Delta r)). Substitution of the recorded same-slice (b,X,K) values yields the following terminal-band change at cell 79:

| Accepted time | Boundary Δ(A_a[-1]) | (M_{\mathrm{geom}}) before boundary | (M_{\mathrm{geom}}) after boundary | Same-slice Δ(M) |
|---:|---:|---:|---:|---:|
| 46.680 | +2.53865e−5 | +8.78449e−4 | +9.32432e−4 | +5.39832e−5 |
| 46.710 | +2.85039e−5 | +9.42838e−4 | +1.00346e−3 | +6.06203e−5 |
| 46.740 | +3.17962e−5 | +1.01159e−3 | +1.07923e−3 | +6.76311e−5 |
| 46.770 | +3.52480e−5 | +1.08488e−3 | +1.15987e−3 | +7.49829e−5 |
| 46.800 | +3.88339e−5 | +1.16285e−3 | +1.24547e−3 | +8.26220e−5 |
| 46.830 | +4.25264e−5 | +1.24562e−3 | +1.33611e−3 | +9.04894e−5 |
| 46.860 | +4.62923e−5 | +1.33330e−3 | +1.43182e−3 | +9.85149e−5 |
| 46.890 | +5.00893e−5 | +1.42598e−3 | +1.53259e−3 | +1.06608e−4 |
| 46.920 | +5.38647e−5 | +1.52371e−3 | +1.63837e−3 | +1.14657e−4 |
| 46.950 | +5.75540e−5 | +1.62654e−3 | +1.74907e−3 | +1.22523e−4 |

For all 99 captured completed primary-predictor stages, the formula-based Δ(M) agrees with the directly recomputed geometric-constraint difference to within (4.3\times10^{-17}). At (t=46.95), the boundary raises the geometric momentum term by (1.2252\times10^{-4}), about 7.5% of its pre-boundary value at that cell. This is a *measured same-slice shift in the constraint residual*, not just an arbitrary comparison of (A_a).

The signed effect on the full (M) residual is positive because the matter term is unchanged by this (A_a)-only assignment. Whether that worsens the *absolute* full constraint residual depends on the pre-boundary sign and the exact (-8\pi j_{\mathrm{total}}) value; the ledger does not capture the complete total-(j) term in a way that admits an absolute-residual statement here. Therefore the allowed conclusion is: the frozen primary boundary imposes an increasing, quantitatively non-negligible change on the outer momentum-constraint equation. That warrants a boundary/CMC/constraint compatibility audit, but is **not** proof the P-minus target is invalid and does not establish it caused the later radiation crossing.

### The P-minus correction also perturbs the *global* CMC time-rate projection

The pinned production CMC code supplies a third, nonlocal coupling. In `engine/cmc_gauge.py`, `target_kdot` forms the actual (K)-RHS array as `primary_l2_rhs["K"] + primary_l3["K"]), then computes its proper-volume weighted mean using
[
w_i=r_i^2\sqrt{a_i}\,b_i/X_i^3\,\Delta r,\qquad
\dot K_{\rm CMC}=\frac{\sum_i w_i L_K[i]}{\sum_i w_i}.
]
In the pinned BSSN (K) RHS, the (A_a)-dependent term is (\alpha A_2=\tfrac32\alpha A_a^2) for (A_b=-A_a/2). Therefore, on the same primary-predictor slice, the (A_a[-1])-only boundary assignment changes the outer (K)-RHS by
[
\Delta L_K[-1]=\tfrac32\alpha[-1]
\left(A_{a,post}^2-A_{a,pre}^2\right),
\qquad
\Delta\dot K_{\rm CMC}=
\frac{w_{-1}}{\sum_i w_i}\,\Delta L_K[-1].
]
At (t=46.95), (\alpha[-1]=1) and this local RHS change is (+2.72321186\times10^{-5}). The exact global Δ\(\dot K_{\rm CMC}\) has **not** been assigned a numeric value because this stage ledger stores only the five outer cells, not the full proper-volume weights and full (K)-RHS profile needed for the denominator. A numerical estimate from outer cells alone would be invented.

That global projection is then used in the lapse solve via `rhs = -kdot / inv` on all interior PDE rows. So the P-minus operation is not just a local (A_a) patch or a local H/M residual change: through the proper-volume projection it can also change the CMC lapse solution across the slice, after which α feeds later (A_a/K) and metric RHS evolution. This is a source-code-derived pathway, not a measured magnitude yet. The exact next offline audit should use an artifact that contains full-grid (w_i) and (L_K[i]) (or reconstruct them from an admitted full checkpoint) to calculate Δ\(\dot K_{\rm CMC}\) exactly; if that state is unavailable, a later minimal read-only recorder should capture those full arrays on the same boundary call. Do not change P-minus or the CMC projection before this value and its effect on the lapse are understood.

### The same boundary assignment also shifts the Hamiltonian constraint directly

The pinned production diagnostic defines the Hamiltonian residual as
[
H=R-(A_a^2+2A_b^2)+\tfrac23 K^2-16\pi\rho,\qquad A_b=-\tfrac12 A_a.
]
On the same primary-predictor slice, (R), (K), and total (ho) are unchanged by the boundary's (A_a[-1])-only assignment. Thus the exact direct Hamiltonian-residual jump is
[
\Delta H=-\tfrac32\left[(A_{a,\mathrm{post}})^2-(A_{a,\mathrm{pre}})^2\right].
]

Using the same accepted historical-source ledger and the cleanly reconstructed pre-/post-boundary primary predictor, the jump is:

| Accepted time | Boundary Δ(A_a[-1]) | (H_{A_a}) term before boundary (-\tfrac32 A_{a,pre}^2) | Exact same-slice Δ(H) |
|---:|---:|---:|---:|
| 46.680 | +2.53865e−5 | −2.39620e−2 | −9.62686e−6 |
| 46.710 | +2.85039e−5 | −2.51860e−2 | −1.10817e−5 |
| 46.740 | +3.17962e−5 | −2.64698e−2 | −1.26730e−5 |
| 46.770 | +3.52480e−5 | −2.78155e−2 | −1.44016e−5 |
| 46.800 | +3.88339e−5 | −2.92252e−2 | −1.62640e−5 |
| 46.830 | +4.25264e−5 | −3.07008e−2 | −1.82547e−5 |
| 46.860 | +4.62923e−5 | −3.22442e−2 | −2.03648e−5 |
| 46.890 | +5.00893e−5 | −3.38573e−2 | −2.25798e−5 |
| 46.920 | +5.38647e−5 | −3.55418e−2 | −2.48786e−5 |
| 46.950 | +5.75540e−5 | −3.72994e−2 | −2.72321e−5 |

At (t=46.95), the P-minus correction makes the Hamiltonian residual (2.7232\times10^{-5}) more negative, which is 0.0730% of the pre-boundary (A_a)-quadratic Hamiltonian term. This formula-based delta is exact for the coded diagnostic constraint at fixed state: its derivation changes no curvature, (K), metric, or matter term. It is a measurable same-slice constraint contribution, not a claim that the absolute (H) norm worsens; the pre-boundary full-H sign/value at that exact stage is not present in this artifact in a way that admits that absolute statement.

Because the boundary only changes (A_a), there is no direct same-slice change in the metric-only conformal connection constraint from this assignment alone. The later evolved (Lambda), (A_a), (K), and metric can of course couple across subsequent primary/final stages. Taken together with the separately derived Δ(M), the existing evidence now quantifies direct (H) and (M) footprint at the outer cell without asserting that P-minus is invalid. The boundary target, full constraint propagation, and corrected-source event causation remain separate open questions.

### What this narrows—and what it still cannot prove

This establishes that the outer (P_-) rule materially corrects the outer (A_a) primary predictor and that the correction grows smoothly near the terminal failure. Because (A_a) feeds (dot a=-2\alpha a A_a) at zero shift and also enters the subsequent primary RHS, this is a concrete pathway for a previous boundary condition to affect future metric updates. It strengthens the case for a **boundary/CMC/geometry compatibility audit** rather than simply clearing the boundary because its spatial residual is near roundoff.

It does **not** establish that the boundary correction is wrong, that it causes the radiation-cone crossing, or that the corrected-rate run has the same correction magnitudes. The rate-fixed trace does not contain paired pre-/post-boundary (A_a) on the same accepted snapshots; the historical ledger's completed-state boundary delta is not independently reconstructable from its current saved final-L2 recorder field, which was overwritten by a later diagnostic invocation. I do not substitute that ambiguous value for a measured final-boundary delta. This report therefore admits the **primary-predictor** correction only, and leaves its transfer to the corrected-rate baseline open.

The next discriminating step is offline-first: reconcile the existing P-minus formula with the production CMC outer normalization and reconstruct the characteristic/constraint budget using the archived source rows. If a matched corrected-rate boundary attribution is later deemed essential, instrument the same pre/post boundary (A_a) directly at the production call site and pin all physics/runtime source blobs unchanged. Do not alter the target, add an incoming-wave law, or launch another trajectory until that analysis supplies a specific preregistered question.

## Cross-workstream CMC target-rank evidence and scope limits

The Library/research branch contains an independent frozen-slice study of the same production `target_kdot` function. Its diagnostic branch is `diag/independent-metric-ricci-audit`; relevant artifacts are [initial target/lapse fixed point, run 38209333687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687), [fixed-target boundary cross-check, run 38293992203](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38293992203), [nine-target rank probe, run 38294955148](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294955148), and [pointwise closure audit, run 38210163221](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221). Source-qualified reports: [fixed-target comparison](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/independent-metric-ricci-audit/outputs/cmc_fixed_target_boundary_crosscheck_initial_2026-10-12.md), [rank audit](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/independent-metric-ricci-audit/docs/CMC_TARGET_CLOSURE_RANK_AUDIT_2026-10-12.md), and [consolidated handoff](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/independent-metric-ricci-audit/outputs/CMC_CMC_TARGET_CLOSURE_CONSOLIDATED_HANDOFF_2026-10-12.md).

### What is established there

- At the archived R80 initial slice (N=160, (R=80), (\Delta r=0.5)), with R160 alpha imposed at the shared outer cell, nine fixed (dot K) trials span (3.46833157\times10^{-7}). The solved target tracks input with slope (0.999999999757); (H(\dot K)=\mathrm{target\_kdot}(\alpha_{solve})-\dot K) varies only (4.3788\times10^{-16}), with finite-difference slope (-2.4269\times10^{-10}). The CMC PDE rows close below (3.71\times10^{-15}).
- The outer lapse-normalization cell is not a PDE row; its pointwise (K)-RHS residual is approximately (-4.5882\times10^{-12}) and it carries 1.86291% of the projection weight on that particular initial grid. The scalar closure (H) matches the weighted mean pointwise residual to (1.14\times10^{-20}). Thus (H) is essentially the weighted normalization-row residual after interior PDE rows have closed, and hardly changes across the trial-(dot K) interval.
- With boundary held at the R160 value, changing only R80's global target from native to the archived R160 value removes the initial shared-core lapse-profile mismatch from RMS (2.5835\times10^{-4}) to (2.45\times10^{-12}). Changing only the boundary to native (alpha=1) creates a separate profile difference. This separates target-value effect from boundary normalization on that frozen initial comparison.

These observations strongly support a **near-rank-deficient / near-redundant discrete target–lapse closure** for that tested slice. They do not prove exact multiplicity, establish that the CMC target is physically wrong, or authorize excluding the outer cell. The prescribed stopping tolerance in the separate Picard map was not reached; it must not be relabeled converged.

### Why this matters to the P-minus trace, but does not quantify it

The `target_kdot` function body is byte-identical between the fixed-point diagnostic branch and the corrected-rate branch's `engine/cmc_gauge.py`; the diagnostic `solve_cmc_lapse` has extra source for alternate outer-alpha controls, and its archived initial data differ from the late radiation trace. The initial CMC audit also uses an R80/R160 pair and D amplitude (10^{-4}), whereas the corrected-rate radiation trace is N=80, (R=80), D amplitude (10^{-10}), at late time. Therefore the initial rank/weight numbers above are **structural cross-checks, not late-time measurements**.

For the historical-source late primary-predictor ledger, a measured P-minus correction of (\Delta A_a=+5.75540\times10^{-5}) at t=46.95 changes the local (K)-RHS by (+2.72321186\times10^{-5}). Because that value enters the all-grid proper-volume projection, the exact global target change is (\Delta\dot K_{CMC}=(w_{-1}/\sum_i w_i)\Delta L_K[-1]). The late ledger lacks the full-grid weights/RHS, so this nonlocal target change and its effect on the solved lapse have not been numerically admitted. The prior CMC audit shows why that missing scalar matters: the target residual may contain little independent selection of (dot K) beyond the outer-row seam.

Next work is therefore a **full-grid, same-slice CMC sensitivity reconstruction using an already archived checkpoint if available**: preserve all state fields, calculate each (w_i), (L_K[i]), the exact target change induced by P-minus, and the resulting lapse profile for fixed versus recomputed target. If no suitable full checkpoint exists, write the minimal read-only recorder specification before any matched candidate is run. No change to (P_-), CMC target average, lapse normalization, or radiation cone is justified yet.

## Full-grid checkpoint recovery: quantify the global CMC target channel at t≈60

The full canonical run_campaign artifact contains twelve native N=120 checkpoints (t=5 through t=60), not only the summary ledger. This enables a whole-grid proper-volume weight calculation without another evolution.

- Canonical production run: [38025979654](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38025979654), artifact [11659369195](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38025979654/artifacts/11659369195). Run status is failure because target tau=22.1 was not reached; artifact upload itself succeeded.
- ZIP SHA-256: 142d92de51bd89a632c9c53d898a389301248d73d171c082e18c0380618eaf24
- N120/ledgers.json SHA-256: d88734d45a6be81e3382aff7b810beed2df8ac63e3f1f80b8d9c328300a369fd
- N120/checkpoint_t00060.000000.npz SHA-256: 81b598c635b58ebc6465dd3525a54a53336e6f654fb8f04d51adb96127e5805c
- Matching boundary trace: [run 38028306184](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184), artifact [11661332148](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184/artifacts/11661332148), trace.json SHA-256 f418a023c44dcd51802696e6468e3e09691d01ea44ecbb5a9463881be6b14e76.

The full-run commit f4a18444899bf78242a40fc2492a37a5d65187ba and boundary-trace baseline a7170c12ed08b4fdb3c52ca98f02e106a3f39b82 have identical Git blob IDs for the relevant production sources (production_kernel.py, cmc_gauge.py, matter_system.py, matter_rhs.py, valencia.py, v55_matter.py, and v55_pirk_adapter.py); the same pinned vendor is used. The terminal/checkpoint states are not bit-identical: corresponding outer values differ by roughly 0.4–1.3e-9, consistent with a small cross-run/runtime trajectory difference. Accordingly, the following CMC perturbation is a near-matched artifact reconstruction, not a bitwise same-state counterfactual.

At canonical checkpoint t=60.00000000000151, tau=21.60324970604194, the native full-grid CMC proper-volume weights are w_i = r_i^2 * sqrt(a_i) * b_i / X_i^3 * dr. Their total is 954906.55541281728; the outer-cell weight is 26653.673999705257, or 2.79123374414% of the slice. The recorded CMC target is Kdot_CMC=0.06692678496752487.

The matching boundary trace's completed-state P-minus row labelled t=59.97000000000151 (the step completing into the t≈60 checkpoint) has Aa_pre[-1]=0.046751434211626593; its post-boundary Aa is 0.046754262454086246. The canonical checkpoint's actual post-boundary value is 0.046754261846389554. Using the canonical checkpoint's post value and the trace's pre value to minimize cross-run mismatch gives local delta_LK[-1] = 1.5*alpha[-1]*(Aa_post^2-Aa_pre^2) = +3.9659993506294643e-7, and delta_Kdot_CMC = (w_outer/sum(w))*delta_LK = +1.107003121672167e-8. That is 1.65405094867e-7 of the recorded K rate. Using the trace's own pre/post pair with the canonical whole-grid weights yields +1.10724103871e-8, so the two artifact pairings agree at 2.4e-12 absolute for this scalar.

This is the first quantified whole-grid CMC-target effect of the P-minus change. At this t≈60 checkpoint, its direct contribution to the proper-volume average is extremely small relative to Kdot_CMC; the evidence does not support claiming this global scalar jump alone explains the radiation-cone failure. It still does not determine the solved-lapse response, the correction at the exact t≈61.0 terminal state, or the total effect of prior boundary history. The boundary correction grows steeply late; applying t=60 weights to the terminal delta_Aa would be an unverified cross-time extrapolation and is not done.

## Frozen-grid CMC lapse response from the recovered t=60 checkpoint

I reconstructed the linear CMC lapse response to the near-matched boundary-induced delta-Kdot-CMC on the full native N=120 canonical checkpoint, without re-running the evolution. The calculation uses the checkpoint's all-grid geometry and scalar/conservative matter fields, the pinned fourth-order cell-centred derivative rows, and the same CMC operator:

- L[delta-alpha] = delta-alpha'' + c delta-alpha' - m delta-alpha.
- c = -0.5 D_r(a)/a + D_r(b)/b - D_r(X)/X + 2/r.
- m = {1.5 Aa^2 + K^2/3 + 4*pi*(rho + pr + 2*pt)} / gamma^rr.
- Frozen outer normalization: delta-alpha[-1] = 0.

The reconstructed total-density maximum is 1.2550018116046338e-4, matching the campaign's t=60 native ledger value to the last displayed digit. The boundary changes only Aa[-1], so the interior PDE-row coefficients m[i] for i < N-1 are unchanged; the target scalar changes the interior RHS by -delta-Kdot-CMC / gamma^rr[i].

For the cross-artifact-qualified delta-Kdot-CMC = +1.107003121672167e-8 at t≈60, solving the frozen linear response gives:

| Measure | Frozen-operator response |
|---|---:|
| max absolute delta-alpha | 7.031980962362896e-7 |
| location of max | r = 102.5 |
| RMS absolute delta-alpha, full grid | 3.049057233153554e-7 |
| boundary delta-alpha[-1] | 0 (normalization fixed) |
| max absolute delta-alpha / absolute checkpoint alpha | 1.54446e-7 |
| max absolute lapse on checkpoint | 4.564922859146704 |

Using the trace's own (rather than canonical post-boundary) Aa pair changes the inferred delta-Kdot by about 0.0215% and the max response to 7.03349e-7, so this cross-pair uncertainty is immaterial to the scale at this checkpoint.

The matrix-system solve is a frozen-coefficient derivative, not a second full nonlinear production step. It does not include previous-step boundary history or the direct Aa-to-a-dot pathway. As a consistency guard, the saved checkpoint lapse itself has a maximum residual around 2.1e-5 under the reconstructed CMC operator and recorded Kdot, with the stored alpha profile about 1.06e-3 from the exact same-slice linear solve. The campaign diagnostic records the last-interior residual near 1.85e-5; the full interior maximum is not that boundary-only diagnostic. This stage mismatch is why I report only the difference between two linear solves on the same operator, not that the reconstructed solve reproduces the entire recorded lapse bit-for-bit.

Conclusion limited to t≈60: the global CMC target scalar channel from the measured P-minus Aa-change is very small in the solved lapse, max about 7e-7. That makes this scalar channel alone an unlikely explanation of the large late outer lapse/geometry structure. It does not exonerate the boundary: direct Aa-driven metric evolution, same-slice H/M changes, outer-face radiation flux, boundary history, and cross-gauge coupling remain separate possible pathways. Do not transfer this t=60 sensitivity to the t≈61 failure without a matched state.

## What this does and does not establish

**Established for the captured frozen configuration:**

1. The old high-resolution early failure from checking a dt=1 surrogate was real software behavior but false evidence of a production-sized matter update failure.
2. After rate extraction was corrected and fixed dt=0.03 used, the engine accepted 1,566 steps before a distinct radiation-cone crossing at the first predictor CMC solve in outer cell i=79.
3. In that exact step, flux and geometric source raise the conservative cone margin when evaluated with the accepted metric; the explicit predictor metric then makes it negative.
4. The raw a and X update raises X²/a enough to account for the stage crossing; the centre-regularity projection has zero change at that cell.
5. The trace reconstructs predictor conservative variables, and its ordered cone-margin and inverse-metric budgets report zero residual.

**Not established:**

- That a specific RHS term is solely incorrect or that a physical interaction law is missing.
- That earlier outer-boundary influence is absent.
- That a clamp, relaxed cone, accepted-metric recovery, or projection modification is justified.
- That the full machine passes centre-series regularity, joint H/M/C_lambda constraint propagation, Misner–Sharp conservation, resolution convergence, dynamic D-to-matter transfer, local/global-Q closure, or production admission.
- That a bounce, turnaround, or cycle exists.

## Next gate

Use this captured failing-step record rather than another long trajectory to perform a bounded geometry/source sensitivity attribution on the same accepted snapshot: verify the evolution-equation contributions to the a and X rates; check source and constraint budgets; and compare the resulting X²/a change to the radiation cone margin with a preregistered criterion. Keep geometry, radiation transport, boundary closure, CMC solve, and physical interpretation as separate ledger columns. No production change follows until a matched comparison distinguishes the defective layer. Keep the exact E ≥ |S| guard.

No change to `main`, production physics, or the canonical GEAR archive was made for this diagnostic report.


## Cross-resolution centre probe: unit-step RHS-extraction guard confirmed — 12 October 2026

**Status:** Diagnostic artifact admitted for call-path attribution only. The counterfactual run exited nonzero for two N=160 cases; this is not four failed physical comparisons and is not evidence for or against the proposed radiation four-force.

- Workflow: [38316415773](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38316415773), completed with expected per-case rejection records and successful artifact upload.
- Diagnostic branch/commit: `diag/radiation-four-force-counterfactual-2026-10-12`, `1f1c49d722807082bf7a1ba112afdab7d5726843`.
- Artifact: [11777115667 — radiation-four-force-counterfactual](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38316415773/artifacts/11777115667).
- Artifact ZIP SHA-256, independently checked: `2af6bb96012655b88ea57f1610819be35c5fa94d2da5b481561ff27fd28abd5d`.
- JSON SHA-256, independently checked against the artifact manifest: `0fcaa70934c6b830da3689e61d44eeaedca18c03fdd69c390ceeee8a7c678e0d`.
- Focused null-identity tests: 3 passed; full existing suite: 25 passed; the four-case runner preserved outputs and returned nonzero for the two N=160 cases.

### What actually passed

The N=80 null and candidate one-step cases both complete at dt = 0.015. Their aggregate outputs are nearly identical. Candidate-minus-null at t = 0.015 includes changes of approximately:
- H_eff: -4.5e-17;
- max admitted Hamiltonian residual: -4.0e-14;
- max admitted momentum residual: -4.0e-18;
- centre radiation margin: +1.01e-8;
- admitted-region radiation margin: +1.29e-10.

This is only the observed response of one short N=80 step under the diagnostic candidate. It does not identify a physical coupling law, establish long-time stability, or admit beta_rad = -0.04 for production. The Misner–Sharp current remains explicitly unavailable/NaN in this first-pass harness.

### Why both N=160 cases reject before the ordinary step

Both N=160 cases begin with an admissible centre primitive at rest. During RHS extraction, the current main-based `species_rhs()` calls `evolve_species(..., dt=1.0)` and subtracts the original conservative state to approximate the semi-discrete rate. The radiation cone guard runs inside that artificial unit update. It rejects at cell 0 before `V55ProductionKernel.step()` has applied its actual dt = 0.0075.

The captured unit-update ratios are:
- beta_rad = 0: |S|/E = 1.3585922116020037;
- beta_rad = -0.04: |S|/E = 1.3584029299128977.

The same structured capture records the actual initial centre-face flux as the existing zero vector and the first interior face radial-momentum HLL flux as about +2.3399337053e-10. The centre flux tendency and geometric-source split is, to the captured precision, the same for the two coefficients: flux momentum tendency -9.3597348210e-10, geometric momentum source +7.4876150659e-10, and energy source -7.9588193414e-13. The small difference in rejected unit-step ratios does not isolate a centre transport change.

The artifact separately forms the diagnostic first-stage candidate U + dt L(U), keeping the pre-step metric fixed. It finds:
- beta_rad = 0: |S|/E = 0.010131364229693135;
- beta_rad = -0.04: |S|/E = 0.010131353702168193.

Both satisfy the existing cone inequality in this limited calculation. **This is not the complete PIRK predictor/final state**, but it directly shows why a unit-sized wrapper update must not be mistaken for the actual timestep-scaled state.

### Cross-agent reconciliation: existing PR #8 is the isolated correction

The same API-contract defect is already isolated in draft [PR #8 — Separate matter RHS extraction from finite-time updates](https://github.com/miahpenn/0-Cosmos-engine/pull/8). Its source-only refactor calls the shared conservative rate calculation directly instead of emulating a dt=1 finite update, while preserving incoming-state recovery checks and admissibility validation on real finite-time updates. Reported checks include [full contract CI 38307704487](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307704487) and the matched coupled one-step A/B [38307235876](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307235876), which completed both CMC arms at N=40, 60, 80, 120, 160, 240 and 320. PR #8 remains draft/unmerged; its passing bounded tests do not admit the full production trajectory.

Therefore do not copy a second version of that fix into PR #9 and do not treat the PR #9 red runner as evidence against beta_rad. The next clean discriminator for the four-force question is a matched four-case one-step test against the already-reviewed rate-extraction implementation, with the candidate still diagnostic-only and PR #9 still draft. Before doing that, reconcile the two branch diffs so the experiment changes only the candidate coefficient plus the approved RHS API correction. Keep actual finite-time predictor/final cone checks; no guard relaxation, clipping, floors, damping or physical promotion.

### Keep the two radiation failures separate

This central N=160 failure is the uncorrected unit-step rate-extraction artifact. It is distinct from the corrected-rate N=80 trajectory that reaches about t = 46.98 and then crosses the radiation cone at the outermost cell during the explicit geometry predictor. PR #8 addresses the former contract defect; it has not repaired or explained the latter late-time metric/matter-stage synchronization failure. Both are software/numerical admission issues until independent evidence establishes more. No bounce, turnaround, missing-physics law, or production readiness conclusion follows from this centre probe.

No changes to `main`, the canonical GEAR archive, or production physics were made by this report addendum.


## Rate-corrected four-case beta_rad one-step comparison — 12 October 2026

**Status:** Diagnostic-only four-case comparison passed. This is one timestep per case, not physical-law validation, a long trajectory, or production admission.

- Workflow: [38317295748](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317295748), success.
- Exact runner commit: `83486b03c0789f571c76e7a5a4287b6afaf13330` on isolated branch `diag/radiation-four-force-ratefix-ab-2026-10-12`.
- Base repair source: PR #10, head `730d6788bff2b53cf84d0504842f95979953f1bd` at branch creation. PR #10 stays draft/unmerged.
- Focused RHS/admissibility repair regressions: 4 passed; four-force null-identity tests: 3 passed; full tests suite: 29 passed.
- Pinned vendor verified: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Artifact: [11777006845 — radiation-four-force-ratefix-ab](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317295748/artifacts/11777006845).
- Artifact ZIP SHA-256: `9d29529a8096b6a3b46b8704ad0901bb2a3e637ad6f5373adbb4767bab227e34`.
- JSON SHA-256, agrees with the artifact's manifest: `3e2a151895160a2e82caa8c8225f80af9b8721a53e5b7343869dc08eeaebcddb`.

### Four actual finite-time cases

| Resolution | beta_rad | dt | Result | Max post-step |S|/E | Min (1-|S|/E) |
|---:|---:|---:|---|---:|---:|
| 80 | 0 | 0.015 | passed | 0.0099086223695 | 0.990091377631 |
| 80 | -0.04 | 0.015 | passed | 0.0099086123118 | 0.990091387688 |
| 160 | 0 | 0.0075 | passed | 0.0099068411682 | 0.990093158832 |
| 160 | -0.04 | 0.0075 | passed | 0.0099068361387 | 0.990093163861 |

All 80 / 160 radiation cells are finite and admissible in the post-step profile. The maximum ratio occurs at the centre; the outer-cell ratio is much smaller. Thus the corrected RHS contract allows both N=160 arms to complete their real timestep; the prior N=160 rejections in [38316415773](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38316415773) were the artificial unit-step extraction guard, not a physical-step rejection.

### Matched candidate-minus-null response

The per-cell conserved profiles give these maximum absolute differences after the one step:

| Witness | N=80 | N=160 |
|---|---:|---:|
| Eulerian radiation E | 1.84241e-14 | 9.21454e-15 |
| (|S|) | 9.33393e-17 | 4.66639e-17 |
| (|S|/E) | 1.00577e-8 | 5.02955e-9 |
| Densitized conservative energy (U_E) | 2.87746e-11 | 1.44783e-11 |
| Densitized radial momentum (U_r) | 1.82071e-15 | 4.55131e-16 |

The global summary changes are also tiny: candidate-minus-null (Delta H_{m eff}) is (-4.51	imes10^{-17}) at N=80 and (-1.30	imes10^{-17}) at N=160; (Deltamax|H|) is about (-3.68	imes10^{-14}) / (-1.95	imes10^{-14}); (Deltamax|M|) is about (-5.11	imes10^{-17}) / (-2.55	imes10^{-16}). These are small one-step counterfactual differences, not an admitted physical effect. The four-force is currently a one-sided diagnostic source in radiation; there is no matching opposite scalar four-force in this candidate, so it is not a demonstrated covariantly closed exchange model. Do not promote (eta_{m rad}=-0.04) from this result.

**Misner–Sharp current/balance remains unavailable** in this harness (reported as NaN, not counted as passing). This run does not address the distinct late-time N=80 outer-cell metric-predictor cone crossing near (t=47).

### Independent follow-up lead: unused radiation rest-current channel

The new post-step profile records a nonzero radiation `rest` slot after evolution despite `initialize_radiation()` initializing it to zero: at the outermost cell it is about (-1.00857	imes10^{-11}) for N=80 and (-5.01262	imes10^{-12}) for N=160. The candidate-minus-null difference is only roundoff-scale. In `engine/matter_system.py`, `_radiation_conserved()` and `_radiation_flux()` explicitly return zero in their first component, but `_hll_flux()` and the outer face currently use the generic `_valencia_conserved()` / `_valencia_flux()` channel for every species. That lets the generic fluid rest-current flux populate the radiation `rest` slot.

The radiation primitive recovery and Einstein stress-energy projection use the radiation energy/momentum variables and ignore this slot, so this observation does not establish a direct effect on (|S|/E) or the late cone failure. It is a separate conservation/representation contract to audit with a regression test before proposing a small species-dispatch correction. No production source is changed here.


## Radiation rest-current dispatch mismatch — bounded same-state A/B, 12 October 2026

**Status:** Diagnostic A/B admitted; minimal representation-repair candidate passed dedicated CI and is tracked as draft stacked [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11). This is an internal representation invariant, not an admitted cause of the late radiation-cone crossing.

### Diagnostic provenance

- Audit run: [38317797540](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317797540) — success.
- Diagnostic branch: `diag/radiation-rest-current-channel-audit-2026-10-12`; source capture commit `09aa20e37d4f1232815cede89dcfbe2514ad6980`.
- Artifact: [11777312270 — radiation-rest-current-channel-audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317797540/artifacts/11777312270).
- Artifact ZIP SHA-256 (locally verified): `e34c99e6b575222f0f8ce5c134f9ac2bdb9a958ece0dc235f593a8a63e38e090`.
- Embedded JSON SHA-256 (matches its SHA256 file): `566b51387ca3272fa76716b19d5e1128450ef02c6ae784ee851f2759697a785b`.
- Pinned vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Test gates: 1 radiation-specific helper-contract test passed; 4 existing radiation RHS/admissibility tests passed; full test suite 27 passed; four diagnostic one-step arms passed.

### Exact isolation result

The one-step diagnostic compares generic all-species HLL dispatch with a diagnostic-process-only dispatch to the existing radiation-specific helpers, at N=80 ((dt=0.015)) and N=160 ((dt=0.0075)).

| Quantity | N=80 generic | N=80 radiation-specific | N=160 generic | N=160 radiation-specific |
|---|---:|---:|---:|---:|
| Max (|mathrm{rest}|) after step | (1.0085731566\times10^{-11}) | 0 exactly | (5.0126168003\times10^{-12}) | 0 exactly |
| Max (|S|/E) after step | 0.0099086223695 | identical exactly | 0.0099068411682 | identical exactly |
| Min (E-|S|) | (8.7823129289\times10^{-9}) | identical exactly | (8.7838426787\times10^{-9}) | identical exactly |

Across both resolutions, the **entire** post-step radiation (E), (S), (E-|S|), (|S|/E), and (1-|S|/E) profiles are bit-identical between the two arms. All stored geometry fields ((a,b,X,\alpha,A_a,K,\Lambda)) and the H, M, C_Lambda summary witnesses are also bit-identical. All real radiation state checks pass; neither arm creates an inadmissible E/S state.

### Source-level interpretation

The production file initializes radiation's `ConservedSpecies.rest` array as exact zero and provides `_radiation_conserved()` and `_radiation_flux()` helpers whose first component is zero. But `_hll_flux()` and the outer-face closure used generic `_valencia_conserved()` / `_valencia_flux()` for all species. That gives radiation a generic fluid rest-current flux which should not be part of the radiation E/S state contract.

The A/B shows that this first component is passive for the tested one-step evolution: removing it changes no E/S value, cone witness, geometry, or constraint witness. This makes it a clean, separately repairable state-representation bug candidate. It does **not** explain or repair the separate late-time N=80 cone crossing, which is driven by the near-cone predicted matter state evaluated against the explicit a/X predictor metric.

### Minimal repair candidate; admission pending

I created separate branch [`fix/radiation-rest-current-channel-2026-10-12`](https://github.com/miahpenn/0-Cosmos-engine/tree/fix/radiation-rest-current-channel-2026-10-12), based on PR #10's unmerged radiation RHS/admissibility repair candidate. It changes only species dispatch in `engine/matter_system.py`: radiation HLL conserved/flux entries and the radiation outer-face flux now use the already-existing radiation-specific helpers. No E/S formula, source, metric stage, time integrator, radiation cone gate, boundary law, or physical coupling is changed. New tests assert the helper identities and exact-zero rest slot after real N=80/N=160 kernel steps.

Dedicated gate [38317961532](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317961532) passed at tested workflow/code head `b9185aa160adaceb5917a9162892e466945ccbd2` (3 focused rest-channel tests, 4 radiation RHS/admissibility tests and 29 full-suite tests passed). The reviewable candidate is draft stacked [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11), based on draft PR #10. `main` remains `2035acc1210db18da0cededea18750b7ed813260`; PR #10, PR #9 and the authoritative GEAR archive are unchanged/unmerged. No long trajectory was launched.


### Repair-candidate gate completed — zero radiation rest-current regression, 12 October 2026

The isolated implementation candidate has now passed its dedicated CI; update the earlier “gate is running” status above accordingly.

- Branch: [`fix/radiation-rest-current-channel-2026-10-12`](https://github.com/miahpenn/0-Cosmos-engine/tree/fix/radiation-rest-current-channel-2026-10-12), stacked on PR #10's unmerged guard-repair baseline `5a5a9657e5d0e45c899addd5851ab2bd8db1cfe0`.
- Tested workflow/head: [38317961532](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317961532), success at `b9185aa160adaceb5917a9162892e466945ccbd2`.
- Source dispatch edit: `1e41d4fa4fc63d4f14bc135670865d4e44f48e11`.
- Focused zero-rest/helper tests: 3 passed; radiation RHS/admissibility gate tests: 4 passed; full repository suite: 29 passed. The exact vendor pin remains `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Repair handoff: [radiation_rest_current_contract_repair_handoff_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/fix/radiation-rest-current-channel-2026-10-12/docs/radiation_rest_current_contract_repair_handoff_2026-10-12.md).

The candidate dispatches radiation HLL conserved/flux vectors to existing `_radiation_conserved()` / `_radiation_flux()`, and the outer face through `_radiation_flux()`. All non-radiation species retain the generic Valencia path. The regression confirms actual N=80 and N=160 production steps preserve an exact-zero radiation rest slot and still pass the original (E\ge|S|) checks. The earlier matched A/B (run 38317797540) proves these dispatches leave the whole radiation E/S profile, all metric geometry, and H/M/C_Lambda diagnostics bit-identical over one step.

This is a bounded representation-invariant repair, not a correction to the physical radiation equations or a cause established for the late-time cone event. It remains stacked on draft PR #10 and has **not** been merged to main. The corrected-rate late-time outer-cell predictor crossing, full Misner–Sharp balance, convergence, and full production admission remain open.


### Checkpoint surface audit for the rest slot

On current `main`, `write_checkpoint_npz()` serializes the `radiation_rest` array alongside the radiation energy/momentum fields. Repository code search found this as a writer field; the radiation-specific primitive recovery and stress-energy projection use energy/momentum and do not use this rest slot. Thus the generic-path drift is observable in saved checkpoints even though the matched one-step A/B shows no change to geometry or the E/S cone witnesses. The candidate's exact-zero invariant prevents that spurious channel from being persisted by new candidate runs. This remains a data-contract cleanup and not evidence that the late-time cone crossing originates here.

Review candidate: [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11), stacked on draft [PR #10](https://github.com/miahpenn/0-Cosmos-engine/pull/10). No merge is performed.


## PIRK Lambda predictor stale-vs-sequential A/B — 12 October 2026

**Status:** Completed bounded diagnostic-only comparison; one ordinary step per arm at N=80 and N=160. A discrete stage-order mismatch is confirmed, but the observed full-step state / radiation-cone changes are negligible at this initialization and timestep. Do not label this the cause of the late-time cone crossing.

- Stage-order test run: [38318119471](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318119471), success.
- Bounded A/B run: [38318309196](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196), success at audit commit `ce9a49560b04ccbe33d5b51ac86219a5c4a3c2f9`.
- Artifact [11777443806 — pirk-lambda-predictor-stage-order-ab](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196/artifacts/11777443806). ZIP SHA-256 `1ec986b7cbbe97ba23184c2d448aba26199257f6ae9d6c4182490c81af3867fe`; embedded report SHA-256 `283c5b337a90670c9a397c9584d19447497017cc44a640ca96bc7303a03186ef`, matching `report.sha256`.
- Pinned geometry/PIRK vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`; focused stage-order test 1 passed and full suite 36 passed.
- Full report: [PIRK_LAMBDA_PREDICTOR_STAGE_ORDER_AB_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/pirk-lambda-predictor-stage-order-2026-10-12/docs/PIRK_LAMBDA_PREDICTOR_STAGE_ORDER_AB_2026-10-12.md).

The stale arm substitutes the Lambda-L2 vector computed on explicit U1 before Aa/K are updated. The sequential arm recalculates Lambda-L2 on stage1 after predictor Aa/K assignments, before the boundary and next CMC lapse solve. Both arms started from bit-identical states.

| Witness | N=80 | N=160 |
|---|---:|---:|
| Max stale-vs-updated Lambda-L2 array difference | 5.3217372642e-8 | 1.3313374252e-8 |
| Predictor Lambda delta, observed | 3.9913029479e-10 | 4.9925153421e-11 |
| Completed-state max Lambda array delta | 6.9388939039e-17 | 2.1684043450e-18 |
| Completed-state max Aa array delta | 1.6329877488e-12 | 1.9526410357e-13 |
| Full-domain H max-absolute delta | -2.43e-17 | 0 exactly |
| Full-domain M max-absolute delta | -4.25e-12 | -1.06e-12 |
| Full-domain C-Lambda max-absolute delta | -6.94e-17 | -2.17e-18 |

Both legacy and candidate arms completed and remained radiation-admissible. Full-domain max |S|/E was 0.005013540746 at N=80 and 0.005015096012 at N=160 in both arms; the minimum cone margins also match exactly (2.2067153165e-9 and 5.5163762656e-10). All other stored matter/scalar and geometry-array differences are zero except the tiny Aa/Lambda/derived lapse deltas shown above.

**Interpretation:** the sequential stage order is supported by the pinned operator contract and now has a passing regression. In this bounded test, however, the Lambda difference almost disappears by the completed state; the cone witnesses are unchanged. It therefore does not explain or fix the separate accepted-t≈46.98 outer-cell N=80 radiation crossing, nor justify physical changes or promotion. Main remains `2035acc1210db18da0cededea18750b7ed813260`; this remains a diagnostic branch only.

### Rest-current stacked PR revalidation after PR #10 advanced — 12 October 2026

The rest-current draft branch was rebased after PR #10 advanced from `5a5a9657e5d0e45c899addd5851ab2bd8db1cfe0` to `6f6fd9d81af739cacce16fcd428bdfeac996e343`. The rebased code head is `ec57b533ba5f8ee1b679ae854c8f89494d66f666`; its dedicated gate [38318631134](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318631134) completed successfully. It verified the pinned vendor, compiled source/tests, passed the radiation rest-current tests, passed the radiation RHS/admissibility tests, and passed the full test suite. The later head `ba73f4de91f67cf11c01d6642d556888a3fef2e3` changed only the handoff documentation; source and tests remain those exercised by the green gate.

Draft [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11) is now mergeable against draft PR #10's current head, and remains unmerged. Production main is still `2035acc1210db18da0cededea18750b7ed813260`. This closes the rest-current representation-gate/rebase issue only; it does not change the separate radiation-cone physics/numerical diagnosis.

## Pause checkpoint — latest independent diagnostic result (12 October 2026)

The separate PIRK stage-order audit completed successfully and has been saved in the diagnostic branch:
- [A/B workflow 38318309196](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196), artifact [11777443806](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318309196/artifacts/11777443806).
- ZIP SHA-256 `1ec986b7cbbe97ba23184c2d448aba26199257f6ae9d6c4182490c81af3867fe`; JSON SHA-256 `283c5b337a90670c9a397c9584d19447497017cc44a640ca96bc7303a03186ef`.
- Pinned-vendor regression run [38318119471](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318119471) passed 1 focused test and 36 total tests.
- The raw Lambda-L2 stage value differs, but the predictor Lambda difference is only 3.99e-10 (N=80) / 4.99e-11 (N=160), and the completed-state Lambda difference falls to roundoff. Radiation E/S and cone witnesses are unchanged in this bounded one-step A/B.

See [stage-order audit](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/pirk-lambda-predictor-stage-order-2026-10-12/docs/PIRK_LAMBDA_PREDICTOR_STAGE_ORDER_AUDIT_2026-10-12.md) and [detailed A/B report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/pirk-lambda-predictor-stage-order-2026-10-12/docs/PIRK_LAMBDA_PREDICTOR_STAGE_ORDER_AB_2026-10-12.md). This confirms a discrete stage-order mismatch but does not explain the accepted-t≈46.98 cone crossing. The radiation rest-current repair remains a separate stacked draft PR #11; neither diagnostic has been merged to main. All production admission criteria remain open as listed in the continuation handoff.
