# Radiation predictor failure after matter-RHS rate-extraction correction

**Date:** 12 October 2026  
**Classification:** Diagnostic capture; production admission remains open.  
**Decision:** The former high-resolution early-step radiation RHS false failure was a unit-increment guard artifact. A *separate, real late-time predictor-metric cone crossing remains* after the rate-extraction correction. Do not clamp or relax the radiation invariant. No production/main/GEAR merge is authorized by this evidence.

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
- The fail-closed workflow proved that the adapted diagnostic file was the only allowed `engine/` difference relative to the corrected rate-fix source; the vendor pin matched. Python 3.12.15, NumPy 2.5.4. Regression suite, compile, trace capture, hashes, and artifact upload all succeeded.
- Configuration: N=80, Rmax=80, dr=1, CFL=0.03, fixed dt=0.03, D amplitude (10^{-10}), radiation ON, scalar profile amplitude 0.01 and width 7. Target t=50; no admissibility retries, and accepted-metric recovery disabled.

### Independent raw-geometry stage ledger

- Workflow: [38309949383](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309949383)
- Artifact: [radiation-predictor-metric-stage-ledger, artifact 11773874630](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38309949383/artifacts/11773874630)
- ZIP SHA-256 (verified locally): `356f067a583ec3596c3b15e8b5d17a18058cc7bd7a2bbc7fd379b6b1e2acbe93`
- Stage-ledger JSON SHA-256: `d2e9d507f8300598c4dd5ac45b4bcd3f77d205adc388e460e6a75ae614a214c9`
- This ledger uses a *different frozen historical source baseline*, `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`, rather than the corrected rate-fix source. Its exact values are therefore a corroborating stage-order result, not a matched numerical A/B. The artifact validates 100 ordered ledger records and reports the same failing time/cell.

## Corrected-rate trace: observed event

The workflow captured 1,566 accepted steps and 66 same-cell radiation predictor budgets with no rejected-step retries or instrumentation errors. The trace itself did **not** reach its target time.

- Last accepted coordinate time: (t=46.98000000000102).
- Central proper time there: (	au=22.016046893100416).
- Next explicit predictor attempt: (t=47.01000000000102), (dt=0.03).
- First failure stage: `explicit_predictor_first_CMC_lapse`, entering `V55ProductionKernel._solve_lapse`.
- Cell: (i=79), (r=79.5), the outermost cell.
- Failure: exact radiation cone test (E\geq |S|) fails. The logged physical values are (E=3.9783663187660695\times10^{-8}), (|S|=3.9785943919623527\times10^{-8}), and (|S|/E=1.0000573283549097). Densitized cone margin is (C=U_E-\sqrt{\gamma^{rr}}|U_r|=-1.699741681811858\times10^{-8}).

This is a *real finite-dt predictor failure* on the corrected-rate source, not the former unit-(dt=1) surrogate check. The exception is raised as the first predictor CMC lapse attempts to invert the predictor radiation state, before a successful first CMC solve for that predictor. As usual, a successful workflow status means the expected failure was captured and checksums/artifact upload succeeded; it is not an evolution pass.

## Same-cell cone-margin budget

The recorded budget uses the same cell and separates (1) accepted state + accepted metric, (2) flux-only matter increment evaluated using accepted metric, (3) geometric-source increment evaluated using accepted metric, and (4) that updated conservative state assessed on the predictor metric.

For cell 79 in the failing attempt:

| Quantity | Recorded value |
|---|---:|
| Accepted cone margin (C_0) | (+4.729090141987707\times10^{-7}) |
| After flux-only increment, accepted metric | (+1.484032098292450\times10^{-6}) |
| After flux + geometric source, accepted metric | (+2.219759745727783\times10^{-6}) |
| Metric-change contribution to (C) | (-2.236757162545902\times10^{-6}) |
| Final predictor cone margin | (-1.699741681811858\times10^{-8}) |
| Budget closure error | (0) |

In state-increment terms at this same cell, the conservative energy flux changes (U_E) by (-4.1564015995\times10^{-7}), and radial momentum by (-1.1796291609\times10^{-6}). The geometric source then adds (+5.7872482608\times10^{-6}) to (U_E) and (+4.1765310727\times10^{-6}) to (U_r). Taken together, the flux and source leave a **positive** margin when evaluated on the accepted metric. With those same predictor conservative values, the metric contribution is just large enough to turn (C) negative. The algebraic stage budget closes exactly under its defined order, and the production predictor reconstruction residuals for (U_E) and (U_r) are zero.

This sharply localizes the *immediate trigger*: the predictor spatial metric, not a unit-time radiation RHS wrapper, creates the final crossing in this step. It does not identify a unique underlying term in the geometry RHS as the ultimate defect, and it is not a proof of a missing physical law.

## Geometry update decomposition

The companion raw-geometry ledger at the same point gives the rate-fixed trace's geometry values:

| Cell-79 metric quantity | Value |
|---|---:|
| Accepted (a) | 0.6089966857200283 |
| Explicit RHS (dot a) | -0.1921657198214381 |
| (a) after Euler update | 0.6032317141253851 |
| Accepted (X) | 0.9438736828908237 |
| Explicit RHS (dot X) | +0.08873930174627159 |
| (X) after Euler update | 0.9465358619432119 |
| (gamma^{rr}=X^2/a), accepted | 1.462893887838753 |
| (gamma^{rr}), after (a) update | 1.476874488513230 |
| (gamma^{rr}), after (X) update | 1.485217234050057 |
| Change from (a) update | +0.013980600674477 |
| Change from (X) update | +0.008342745536827 |
| Change due to centre-regularity projection at this cell | 0 |
| Decomposition closure error | 0 |

The explicit decrease in (a), followed by increase in (X), both increase (gamma^{rr}=X^2/a) at this cell. Together they account for (+0.022323346211304) in (gamma^{rr}). The centre-regularity projection leaves (a,b,X) unchanged at this outer cell, so it is not the direct cause of this particular crossing. This statement is limited to the recorded failing cell and step; it does not exonerate all boundary influence or all regularity effects earlier in the evolution.

## Approach to the crossing

The last three fixed-(dt) attempts show the positive predictor margin on the updated metric shrinking toward zero:

| Accepted time | (C) on accepted metric | (C) on predictor metric |
|---:|---:|---:|
| 46.92000000000102 | (+1.3010870\times10^{-6}) | (+8.7787474\times10^{-7}) |
| 46.95000000000102 | (+9.0393351\times10^{-7}) | (+4.4739351\times10^{-7}) |
| 46.98000000000102 | (+4.7290901\times10^{-7}) | (-1.6997417\times10^{-8}) |

On the failing step the conservative source and flux lift the margin to (2.22\times10^{-6}) on the accepted metric, but the predictor metric removes (2.24\times10^{-6}). That confirms a real, rapidly growing metric sensitivity at the outer edge in this frozen configuration, but does not say whether it is due to geometric RHS discretization, constraint propagation, boundary-conditioned geometry, gauge/CMC coupling, or a deeper physical model requirement.

## Cross-check against the older stage ledger

The earlier frozen source ledger also captured (t=47.01), cell 79, with (gamma^{rr}) increasing from 1.4628939 to 1.4853572 across the raw Euler geometry update and a negative cone margin after the centre-regularity stage. In that old-source run, the margin remained negative after projection ((C=-1.6900205\times10^{-8}), ratio 1.0000569941); the raw metric already crossed the cone and projection only improved the margin slightly. The rate-fixed trace has a slightly different geometry/matter trajectory and should be treated as the authoritative numerical capture for its own source baseline. This cross-check reinforces the ordering—raw predictor metric crossing precedes the lapse rejection—but does not warrant combining the two rows as if they were a paired experiment.

## What this does and does not establish

**Established for the captured frozen configuration:**
1. The old high-resolution early failure caused by checking a unit-(dt) surrogate was real as software behavior but false as a claim about an actual production-sized matter increment.
2. After rate extraction was corrected and fixed (Delta t=0.03) used, the engine still admits 1,566 steps then encounters a distinct radiation-cone crossing at the first predictor CMC solve, outer cell (i=79).
3. In that exact step, flux and geometric source raise the conservative cone margin when held to the accepted metric; the explicit predictor metric then makes it negative.
4. The raw (a,X) update raises (X^2/a) enough to account for the stage crossing; the centre-regularity projection has zero delta at that cell.
5. The trace math reconstructs the predictor conservative variables and closes its ordered cone-margin and inverse-metric budgets at zero reported residual.

**Not established:**
- That a specific RHS term is solely incorrect or a physical interaction law is missing.
- That outer-boundary influence is absent earlier in the evolution.
- That a clamp, relaxed radiation cone, accepted-metric recovery, or projection modification is justified.
- That the full machine has passed constraint propagation, centre-series regularity, Misner–Sharp conservation, resolution convergence, dynamic (D\to\) matter transfer or local/global-(Q) closure, or any full-trajectory/production-admission criterion.
- That a bounce, turnaround, or cycle exists.

## Next gate

Use the captured failing-step record, not another long trajectory, to formulate a **bounded geometry/source sensitivity attribution** on the same accepted snapshot: verify which evolution equations contribute to (dot a) and (dot X), check their source and constraint budgets, and compare the resulting (X^2/a) variation to the cone margin using a predeclared criterion. Keep geometry, radiation transport, boundary closure, CMC solve, and physical interpretation as separate ledger columns. No production change follows until an admissible, matched comparison distinguishes the actual defective layer. Keep the exact (E\ge |S|) guard.

No change to `main`, production physics, or the canonical GEAR archive was made for this diagnostic report.
