# Issue #4 reference-regeneration artifact plan — 2026-11-07

**Branch:** `issue4-reference-repro-artifacts`  
**Base:** `8965a1985c99d644137321016dc5f172d5695213` (`issue4-integrated-test-diagnostics`)  
**Draft PR:** [#432](https://github.com/miahpenn/0-Cosmos-engine/pull/432)  
**Scope:** synthetic research fixtures only; no production-physics claim.

## Why this follow-up exists

Exact-head workflow run [42407237209](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42407237209) used CPython 3.12.15, NumPy 2.5.4, pytest 9.1.1, and mpmath 1.3.0. Its research job regenerated both reference files and compared only these numeric fields:

- Affine-root: `root` and `cond_inf`.
- Implicit derivative: each `derivative[epsilon]` vector for the committed epsilon list.

Both comparisons returned false. The workflow log did not emit the regenerated JSON or numerical deltas, so the cause could not be established from that run.

## Provenance facts established from repository source

The committed affine-root reference records:

- Python 3.12.3
- NumPy 2.5.4
- mpmath 1.4.1
- map description: port of `_fixed_point_map_output -> physical_k_rhs -> solve_cmc_discrete`

The committed implicit-derivative reference records:

- Python 3.12.3
- NumPy 2.5.4
- mpmath 1.3.0
- `repo_commit=f8d5a65ea0e9bb1064a55f6b0691f860cfaa748a`
- map description names `mp_fixture_port.py`

Both current reference generators import `mp_fixture_helper.py` from `tests/data/issue4_research/`. The repository directory listing on the published diagnostic branch contains `mp_fixture_helper.py`, not `mp_fixture_port.py`. The helper is shared by both generators and calculates `HELPER_SHA256` at import time. The affine-root generator writes that helper hash into newly generated provenance; the committed affine-root JSON does not contain it. The derivative generator does the same; its committed JSON also lacks it.

These are concrete provenance differences, not an explanation for the earlier numeric mismatch. The exact numerical comparison predicate does not include metadata.

## Diagnostic implementation in this branch

`tests/data/issue4_research/compare_regenerated_references.py`:

- Reads committed JSON without modifying it.
- Reads regenerated JSON from runner-temporary paths.
- Uses decimal arithmetic on the serialized decimal strings.
- Reports signed per-cell deltas, exact mismatch counts, maximum absolute delta and its index.
- Reports the `cond_inf` difference independently.
- Reports metadata differences in a separate section, explicitly outside the CI numerical predicate.
- Captures Python/platform/package versions and SHA-256 for generators, committed fixtures, helper, historical test fixture, and tool module.
- Records `mp_fixture_port.py` as missing if absent at its expected path; it does not synthesize or substitute that file.

The workflow writes both regenerated JSON files and the comparison report to the runner temporary directory and uploads them as a 90-day Actions artifact under a run-specific name. Upload is guarded with `if: always()` so numerical mismatch does not discard the evidence.

## Execution result — run 42410428031

Manual dispatch on branch `issue4-reference-repro-artifacts`, exact head `ffa96b745dd06ec08fd0d46cdd2966c82e5e0406`.

[Workflow run 42410428031](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42410428031) completed successfully. All three jobs passed:

- `independent-cmc`: success. Blocking contract tests, visible non-blocking witnesses, and zero-production-call contract step all completed successfully.
- `accepted-state-crosswalk`: success. Frozen accepted-state six-sector and pinned-K RHS crosswalk completed.
- `research-diagnostics`: success. Shared extended-precision port validation passed; regenerated-reference comparison and artifact preservation passed.

Artifact:

- Name: `issue4-reference-regeneration-run-42410428031-attempt-1`
- Artifact ID: `13589236895`
- SHA-256: `8d32d31352e2fd02f616514ecbf387f585068e303328def0f87e3dd961befdcf`
- Expires: 2027-02-05T20:18:36Z
- Contents: `root.json`, `deriv.json`, `comparison-report.json`

Field-by-field numeric findings from `comparison-report.json`:

- Affine-root vector: 20/20 exact matches; mismatch count 0; maximum absolute delta 0.
- Affine-root `cond_inf`: committed and regenerated values both `624310584.1`; signed delta 0.
- Implicit derivative at epsilon `1e-5`: 20/20 exact matches; mismatch count 0; maximum absolute delta 0.
- Implicit derivative at epsilon `2e-5`: 20/20 exact matches; mismatch count 0; maximum absolute delta 0.

The regenerated environment was CPython 3.12.15, NumPy 2.5.4, mpmath 1.3.0, Linux x86_64. The helper SHA-256 was `9ff5ca94d34c0ac397b22e1006007e67c59d1d9f47a662d952cada0eab546381`. The expected historical path `tests/data/issue4_research/mp_fixture_port.py` remains absent; no substitute was created.

**Interpretation:** the earlier reference mismatches did not reproduce in this instrumented run. This run demonstrates exact numerical agreement for the committed fields under the recorded current generator/helper and runner environment. It does not establish why run 42407237209 reported mismatches, nor does it retroactively invalidate that result. The separate lightweight validation run 42408547071 also reported two numerical test failures (affine-root and implicit derivative); those failures remain a conflicting observation to reconcile. No tolerance or reference data was changed to obtain this pass.

## Remaining follow-up

Reconcile why runs 42407237209 and 42408547071 reported numerical mismatches while run 42410428031 regenerated exact matches. Compare the exact checked-out commits, workflow definitions, generator/test entry points, and artifacts/logs for each run before drawing a causal conclusion. Do not overwrite committed references or relax tolerances.

## Guardrails

- No reference fixture is overwritten or regenerated in place.
- No production equations, solver, initializer, vendor code, tolerances, or scientific-admission flags are changed.
- No fitting, tolerance relaxation, or metadata-driven replacement of numeric references.
- Treat the data as synthetic research fixtures, not accepted production physics.
