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

### What this narrows—and what it still cannot prove

This establishes that the outer (P_-) rule materially corrects the outer (A_a) primary predictor and that the correction grows smoothly near the terminal failure. Because (A_a) feeds (dot a=-2\alpha a A_a) at zero shift and also enters the subsequent primary RHS, this is a concrete pathway for a previous boundary condition to affect future metric updates. It strengthens the case for a **boundary/CMC/geometry compatibility audit** rather than simply clearing the boundary because its spatial residual is near roundoff.

It does **not** establish that the boundary correction is wrong, that it causes the radiation-cone crossing, or that the corrected-rate run has the same correction magnitudes. The rate-fixed trace does not contain paired pre-/post-boundary (A_a) on the same accepted snapshots; the historical ledger's completed-state boundary delta is not independently reconstructable from its current saved final-L2 recorder field, which was overwritten by a later diagnostic invocation. I do not substitute that ambiguous value for a measured final-boundary delta. This report therefore admits the **primary-predictor** correction only, and leaves its transfer to the corrected-rate baseline open.

The next discriminating step is offline-first: reconcile the existing P-minus formula with the production CMC outer normalization and reconstruct the characteristic/constraint budget using the archived source rows. If a matched corrected-rate boundary attribution is later deemed essential, instrument the same pre/post boundary (A_a) directly at the production call site and pin all physics/runtime source blobs unchanged. Do not alter the target, add an incoming-wave law, or launch another trajectory until that analysis supplies a specific preregistered question.

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
