# Issue #4 — zero-step full-Bianchi evaluator implementation log

Date: 2026-11-04  
Status: research tooling implementation; accepted-state formula cross-check complete; committed evaluator/tests not yet executed in CI. No scientific admission.

## 1. 0* archive-shape audit

The three accepted-state archives from workflow run `40236402303` were downloaded and independently checked:

| N | ZIP SHA-256 | NPZ SHA-256 | NPZ arrays |
|---:|---|---|---:|
| 60 | `41d646090e05c08684806c5e39d57851f119cf0f38aa90f1d26ab5f72f87d1ce` | `988ccd3c71ac6733bb47d013b32f8ceded9d132d2d21f7483c6ee4953993601b` | 101 |
| 120 | `9f5c2504f28290ff540c96f93e1fb569f5b960b7a2b3fc312dca9362c8e2eacf` | `5751a16949b763ad5e96a4afac9db4fb1d92d55a0c6e7ddaabc4cc6b25652e0b` | 101 |
| 240 | `e8a9cd26fbe3dbfb0bd62fc0606237c387a41efd14ade550df50da9842b38bea` | `29a4bc1aeca718eff0dbc0de068809f8279aab94cb7023b004dea545589bf814` | 101 |

All internal SHA256SUMS entries and all 101 per-array dtype/shape/content hashes pass at each resolution. All arrays are finite. The five scientific-admission flags remain false in every report.

The 0* structural inventory confirms four distinct stored families: frozen input state, semi-discrete tangent, pointwise total-stress divergence, and finite-volume total-stress divergence. It finds no separately stored full-Bianchi predictor-term arrays. Those terms must be reconstructed; the archived residual table is not itself a predictor.

## 2. Two state digests preserved

The archive's upstream accepted-state identity and the tangent helper's canonical state digest are different fields and remain separate:

| N | Upstream archive state digest | Exact tangent-helper state digest |
|---:|---|---|
| 60 | `2ae2a5704b2b87d16ef9239c8a22ef261ed9108d971d3a5c78a3d1133aa9e9e4` | `20c07ec4b7885b59a1e0476d5bab448850518dbc8439ad6db4adaf36fd5a4a33` |
| 120 | `0266c7fafc5e9a292a9a376f1d91fa45bb8dd103b6f8c21b543c25d0930daacb` | `032a5d079b76dfe7bf50b08765f03b9f2dcc730ba1a33329517c153a27bf70cb` |
| 240 | `770d36f3db5b60da87f52beb213f8af4731e17021f95c1519ecbdbe759b41cce` | `dbedf4382bfa58179f87c0d91022574c87f92865c10d125ee66a36181234b5f1` |

The fixture loader reproduces the pinned helper's exact byte convention: ordered field names, `str(shape)`, and contiguous float64 bytes. It also validates the upstream archive identity against the captured-state digest. It does not collapse the two hashes into one.

## 3. Canonical evaluator added

New research-only modules:

- `tools/issue4_accepted_state_fixture.py` — exact archive/member/array checksum validation, two state digests, and 0* operand-family inventory.
- `tools/issue4_full_bianchi_evaluator.py` — one accepted-state evaluator and machine-readable CLI.
- `tools/issue4_bianchi_reference_regression.py` — signed reference witnesses and epsilon-sweep regression gate.

The evaluator verifies both source checkouts, all accepted-anchor source blobs, helper blobs, and vendor pin before loading a state. It imports the reviewed total-stress helper only after the provenance gate. It constructs a digest-guarded frozen-state view from archived arrays; the evaluator only reads that view and checks its canonical digest before and after. It does not construct or call a production kernel RHS.

The evaluator arms in-memory refusal tripwires for `V55ProductionKernel.step` and the exact vendor `damped_moving_puncture_pirk2_step` entry point. The accepted fixture must remain byte-identical under the pinned tangent-helper digest before and after evaluation.

The predictor is transcribed as the single locked expression:

```
shift transport
+ alpha K M_total
- (D_r alpha) H_total / 3
+ alpha (D_r H_total) / 6
+ D_j(alpha W^j_r)
- 8 pi alpha div_T_r_total
```

Metric-only Ricci is evaluated by the pinned vendor with a disposable geometry copy whose Lambda is set to the exact metric-derived Lambda. The evolved-connection response is reconstructed separately. The parity/variable distinction is explicit: `chi=X^2`, `chi_r=-D_1(X)/(2X)`, odd-parity `D_1(C^r)`, and even-parity derivatives of the metric and `F_C`. Pointwise and finite-volume stress-divergence channels are never merged or selected by residual size.

The six-sector report preserves both facts: the raw sector-divergence sum does not reconstruct the total directly, while explicit time-density component accounting does; the unexplained remainder remains separately reported.

## 4. Independent accepted-state formula cross-check

A separate local calculation from the downloaded arrays, pinned vendor operator formulas, and explicit six-sector primitive projections reproduces the recorded signed residuals and epsilon deltas. This is an offline formula cross-check, not a claim that the committed evaluator or its unit suite has run.

| N | Pointwise max, r<=10 | Pointwise annulus max, 5<=r<=10 | FV max, r<=10 | FV annulus max, 5<=r<=10 |
|---:|---:|---:|---:|---:|
| 60 | 1.84058546e-5 | 2.04981956e-6 | 9.64043704e-6 | 1.80888634e-6 |
| 120 | 1.52718464e-6 | 1.64026193e-7 | 1.28267547e-6 | 1.35473916e-7 |
| 240 | 3.98965051e-7 | 1.06674508e-8 | 4.07712412e-8 | 8.98793346e-9 |

The signed cell 0/1/2 witnesses reproduce the recorded signs and values. The tangent epsilon-sweep maxima over r<=10 are:

| N | Tangent change | Pointwise predictor change |
|---:|---:|---:|
| 60 | 1.38236e-13 | 7.10592e-15 |
| 120 | 1.92785e-13 | 1.93322e-14 |
| 240 | 5.59889e-13 | 4.35874e-14 |

The archived signed convention is consistent with `measured tangent - predictor`. The regression gate checks the full signed cell witnesses under one global orientation and the independent pointwise/FV maxima; it does not choose an operator by which residual is smaller.

## 5. Test and execution status

Tests were added for:

- exact two-layer provenance and per-pin refusal;
- exact helper state-digest byte convention;
- archive refusal on wrong digests and the 0* operand inventory;
- signed accepted-state regression witnesses and global sign consistency;
- zero-step static call contract and parity/`chi` distinctions.

The tests are committed but have **not yet been executed** in CI. The committed evaluator itself has also not yet been run against the archives. The offline formula cross-check above is separate evidence, not a substitute for running the committed implementation.

No manufactured or t=0 state adapter has been evaluated. No production step, initializer change, vendor change, workflow change, or admission-state change occurred.

## 6. CI side effect and boundary correction

A draft PR (#428) was opened to expose the repository checks. That unexpectedly triggered the existing `lightweight_validation.yml` workflow, whose standard job includes a separate one-step production-kernel smoke test. The PR was closed when this workflow scope was identified; GitHub continued the already-running job.

The workflow completed successfully:
- complete repository suite: 981 passed, 14 subtests passed;
- curved-slice EM discriminator: 1 passed;
- control-only smoke: passed;
- pre-existing production smoke: one fresh test state at resolution 40, (r_{max}=40), (dt=0.03), advanced to (t=0.03), reporting `PRODUCTION_STEP_SMOKE_PASS`.

That single smoke step was outside the authorized zero-step observer scope. It was not invoked by the observer, did not use any accepted archive, and did not change production source, initializer, vendor, or any persistent project state. Nevertheless, triggering it was my mistake: I should have inspected the full PR workflow before opening the draft. PR #428 remains closed. No further workflow/CI trigger will be initiated without first checking its complete step list and obtaining authorization where a production step is included.

The passing suite corresponds to PR head `631ded29906081436d2b1b49291c2d065b01cbfe`. A subsequent provenance-manifest-only change at the current branch head has not been rerun. The committed accepted-state evaluator itself has not yet been executed against the archives; the separate formula reconstruction in §4 remains the numerical cross-check.
