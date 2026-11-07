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

These are concrete provenance differences, not yet an explanation for the derivative's numeric mismatch. The exact numerical comparison predicate does not include metadata.

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

The workflow now writes both regenerated JSON files and the comparison report to the runner temporary directory and uploads them as a 90-day Actions artifact under a run-specific name. Upload is guarded with `if: always()` so numerical mismatch does not discard the evidence.

## Execution status

This branch's workflow configuration and comparison script have been committed. The new artifact capture has **not yet been observed running**. The previous run 42407237209 predates these changes and contains no comparison artifact from this new script. Do not invent or infer the missing numeric deltas.

Once the branch workflow runs, download and preserve:

1. `root.json` — regenerated affine-root output.
2. `deriv.json` — regenerated implicit-derivative output.
3. `comparison-report.json` — field-by-field numeric deltas, metadata differences, and provenance hashes.

Then append the run ID, artifact name, artifact digest (if exposed), and the exact report findings to this document or a dated result note. Do not replace either committed reference.

## Guardrails

- No reference fixture is overwritten or regenerated in place.
- No production equations, solver, initializer, vendor code, tolerances, or scientific-admission flags are changed.
- No fitting, tolerance relaxation, or metadata-driven replacement of numeric references.
- Treat the data as synthetic research fixtures, not accepted production physics.
