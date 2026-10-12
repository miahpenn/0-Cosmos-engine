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
