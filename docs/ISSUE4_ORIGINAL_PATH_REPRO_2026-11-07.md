# Issue #4 original reference comparison reproduction — 2026-11-07

**PR:** [#433](https://github.com/miahpenn/0-Cosmos-engine/pull/433)  
**Instrumentation branch:** `issue4-original-path-delta-instrumentation`  
**Unmodified base:** `8965a1985c99d644137321016dc5f172d5695213`  
**Instrumented head:** `afc0845c9ff4b7582ab0815fefec3be806367884`

## Purpose and scope

Re-run the original research comparison from recovered tip `8965a198` while preserving its original generator entry points, shared helper, committed references, package pins, and numeric fields. The only comparison-path change is to emit Decimal signed deltas after the original generators write their temporary outputs, and preserve those outputs/report as an Actions artifact.

The original generator commands are unchanged:

- `python tests/data/issue4_affine_root/generate_reference.py 30 "$RUNNER_TEMP/root.json"`
- `python tests/data/issue4_implicit_derivative/generate_derivative_reference.py 30 "$RUNNER_TEMP/deriv.json"`

The same workflow package pins remain Python 3.12.15, NumPy 2.5.4, pytest 9.1.1, and mpmath 1.3.0. The helper and generator files are inherited unchanged from base commit `8965a198`.

## Result: original reference mismatch did not reproduce

[Instrumented run 42411237181](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42411237181) completed. The research diagnostic step passed and uploaded artifact:

- Name: `issue4-original-path-deltas-run-42411237181-attempt-1`
- Artifact ID: `13589737358`
- SHA-256: `644f1eafbc0e0c756ccb23a86333f14c6881fb64a42100ce211c5f6265b52f69`
- Contents: `root.json`, `deriv.json`, `original-path-delta-report.json`
- Expiry: 2027-02-05T20:28:45Z

Decimal report:

- Affine root: 20/20 exact matches; mismatch count 0; maximum absolute delta 0.
- `cond_inf`: committed and regenerated `624310584.1`; delta 0.
- Derivative at `2e-5`: 20/20 exact matches; mismatch count 0; maximum absolute delta 0.
- Derivative at `1e-5`: 20/20 exact matches; mismatch count 0; maximum absolute delta 0.

The report's original string-equality predicate is also true for root, condition number, and both derivative vectors. This is a stronger direct reproduction of the original comparison path than the PR #432 report because it runs the original generator scripts and original inline numeric fields inherited from `8965a198`; only reporting is expanded.

**Interpretation:** run 42407237209's reported false comparison is not reproducible under this re-run, despite matching the recorded package pins and preserving the original generator/helper source. This does not identify the cause of the earlier false result. Its generated JSON was not preserved, so there are no historical deltas to compare. A transient environment/numerical/order effect remains possible but unproven.

## Separate test-witness failures: not the same predicate

The lightweight validation failures in run 42408547071 are not failures of the high-precision committed-reference comparison. They are float64 test assertions in `tests/test_issue4_cmc_lapse_tangent.py`:

1. `test_self_consistent_cmc_lapse_fixed_point_and_implicit_direction_two_epsilons`: implicit directional derivative assertion fails.
2. `test_fixed_point_root_matches_independent_affine_map_and_reports_conditioning`: float64 affine-root comparison against `fixed.alpha_in` fails.

These are independent witnesses with different calculations and tolerances; an exact match of the high-precision JSON fixtures does not establish that these test assertions pass.

Importantly, the same two test witnesses also fail in run 42411237181. That workflow's visible-witness step is configured `continue-on-error: true`, so the job reports the step as completed and the job succeeds even though the test command reports 2 failed, 16 passed. The blocking contract suite separately passes (7 tests). Thus the run's overall success must not be described as all witness tests passing.

The two failure classes are therefore separated:

- **Reference generator comparison:** earlier false in 42407237209; exact match in 42411237181; historical cause unresolved.
- **Float64 numerical witnesses:** failed in 42408547071 and again in 42411237181; remain open and reproducible in the visible-witness command.

## What this establishes and what remains open

Established:

- Original reference generator entry points on the exact recovered source base regenerate all committed numeric fields exactly in run 42411237181.
- The missing historical filename `mp_fixture_port.py` is provenance text only in the committed derivative metadata; the actual original generators import `mp_fixture_helper.py`. The re-run used that actual shared helper and did not synthesize the missing filename.
- The float64 witness failures are a separate predicate from the extended-precision fixture string comparison.

Unresolved:

- Why run 42407237209's same original comparison returned false when the instrumented reproduction returns true. Its regenerated JSON was not retained, so direct per-cell historical comparison is impossible.
- Why the two float64 witness assertions fail, including their numerical conditioning/sensitivity. Do not infer a cause from the high-precision fixture result.

## Guardrails

No committed reference, tolerance, production equation, solver, initializer, vendor source, or scientific-admission flag was changed. No fit or replacement was performed. This is diagnostic-only evidence for synthetic research fixtures and numerical witnesses, not a production-physics admission.


## Recovery-lane disposition — frozen 2026-11-07

The original reference-reproduction lane is **closed**. The only preserved original-path reproduction available is run [#42411237181](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42411237181), which regenerated every committed high-precision numeric field exactly under the original generator entry points and comparison predicate.

The earlier false comparison in run [#42407237209](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42407237209) is **unrecoverable from the retained evidence**: that run did not preserve its generated JSON, so no historical per-field deltas, output hashes, or exact numerical discrepancy can now be recovered. The controlled reproduction did not reproduce it. This closes the recovery task; it does not assert a cause or invalidate the earlier recorded false result. No further attempts to recreate that particular runner event are planned.

The two strict float64 CMC witnesses remain **unchanged and non-blocking**. They failed reproducibly in runs [#42408547071](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42408547071) and [#42411237181](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42411237181), while the high-precision fixture comparison passed. Their measured assertion discrepancies are recorded in the workflow/test logs; the tests are retained as visible witnesses and are not reclassified as passing. No tolerance, fixture, production code, or test assertion is changed by this disposition. This record does not claim that the exact numerical floor has been derived from first principles; it freezes the observed float64 witness behavior without using it to block the reference-recovery lane.

**Disposition:** reference recovery complete; historical mismatch cause unrecoverable; float64 witnesses retained unchanged and non-blocking; Issue #4 scientific work proceeds on the separate derivation/evaluator lane. All scientific-admission, acceptance-oracle, physical-exchange, and conservation-closure flags remain false.
