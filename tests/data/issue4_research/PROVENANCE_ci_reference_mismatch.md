# Issue #4 provenance note: research reference comparison, local reproduction

Scope: this note records local reproduction facts only. Besides this note and the diagnostic
workflow commit described below, this branch also changes the two reference generators
(`tests/data/issue4_affine_root/generate_reference.py`,
`tests/data/issue4_implicit_derivative/generate_derivative_reference.py`) and both reference files
(`tests/data/issue4_affine_root/issue4_affine_root_reference.json`,
`tests/data/issue4_implicit_derivative/issue4_implicit_derivative_reference.json`), in separate
commits. It changes no production code, test fixture, tolerance, or `docs/STATUS.md`. Every value
below was either computed in the local run described here or read from the file named next to it.

## Commits and patch

- Base commit: `8965a1985c99d644137321016dc5f172d5695213` (`8965a198`).
- Patch file: `issue4_ci_research_diagnostic.patch`.
  - SHA-256: `5e7e27e4fd077241593c317d4accdd2c0d4e124e2b170620530a58d450ae3038`
  - Patch commit, read from the patch header line `From c178bffe8c938049dd0e355790f2409119e86cdd`:
    `c178bffe`. This object is not present in the repository. Applying the patch locally with
    `git am` produced local commit `1f724c0c94ad83cc4fd795d0c8826aa1ef73c2ca` (different committer),
    which was not pushed.
- Fixed patch commit: `4406b4e9a46441726df416f4aae01133b146c993` (parent `8965a198`; the commit
  immediately below this note on this branch). It amends `1f724c0c`: the diagnostic step hashes only
  numeric inputs and records the non-numeric input `base.d1` by name, type, qualname and
  address-free repr; the commit message states 22 numeric (float64) hashed inputs.
  - Exported as `issue4_ci_research_diagnostic_fixed.patch` (`git format-patch -1 4406b4e9`),
    SHA-256: `c147e8f810f089d95df736d73fa9927a369460e3b3c60e9f67f8c42ba7f63fab`

## Environment

Recorded output of
`python -c "import platform,numpy,mpmath;print(platform.python_version(),numpy.__version__,mpmath.__version__)"`:

```
3.12.15 2.5.4 1.3.0
```

- Python 3.12.15 (CPython, installed with uv 0.13.7)
- numpy 2.5.4
- pytest 9.1.1 (`python -c "import pytest;print(pytest.__version__)"`)
- mpmath 1.3.0
- Platform string written in the diagnostic JSON: `Linux-6.18.44-fc-v310-x86_64-with-glibc2.39`, CPU model `Intel(R) Xeon(R) Processor @ 2.30GHz`.

## Commands, in order

Path variables used below:

```
S=/tmp/claude-0/-home-user/3f8452aa-83fe-5834-bdb1-93689a21c872/scratchpad
PATCH=/root/.claude/uploads/3f8452aa-83fe-5834-bdb1-93689a21c872/21e646b8-issue4_ci_research_diagnostic.patch
```

Environment setup:

```
uv self update                      # failed: "Self-update is only available for uv binaries installed via the standalone installation scripts."
uv python install 3.12.15           # failed with uv 0.11.32: "No download found for request: cpython-3.12.15-linux-x86_64-gnu"
python3 -m pip install --user --upgrade uv     # installed uv 0.13.7
uv python install 3.12.15           # "Installed Python 3.12.15"
uv venv -q --python 3.12.15 $S/py312
source $S/py312/bin/activate
uv pip install -q numpy==2.5.4 pytest==9.1.1 mpmath==1.3.0
```

Test worktree (unpatched base) and patched worktree:

```
git worktree add -q --detach $S/wt-8965a198 8965a198
git worktree add -q --detach $S/wt-patch 8965a198
cd $S/wt-patch && git -c user.name=local -c user.email=local@invalid am $PATCH
```

The two inline scripts of the `research-diagnostics` job were extracted verbatim from the patched
`.github/workflows/issue4-cmc-lapse-tangent.yml` (heredoc bodies, 10-space YAML indent removed) to
`$S/research-diagnostic-run/compare.py` (step "Regenerate references to temporary files and compare
with committed data") and `$S/research-diagnostic-run/diagnostic_step.py` (step "Diagnostic
fingerprints and numerical deltas").

Run 1 (with `RUNNER_TEMP=$S/research-diagnostic-run`, from `$S/wt-patch`):

```
python -m pytest -q -p no:cacheprovider tests/test_issue4_research_port_validation.py
python tests/data/issue4_affine_root/generate_reference.py 30 "$RUNNER_TEMP/root.json"
python tests/data/issue4_implicit_derivative/generate_derivative_reference.py 30 "$RUNNER_TEMP/deriv.json"
python $S/research-diagnostic-run/compare.py
python $S/research-diagnostic-run/diagnostic_step.py
```

Run 2 (full repeat, `PYTHONDONTWRITEBYTECODE=1`, `OUT=$S/research-diagnostic-run2`):

```
sha256sum $PATCH
head -1 $PATCH
python --version
python -c "import platform,numpy,mpmath;print(platform.python_version(),numpy.__version__,mpmath.__version__)"
python -c "import pytest;print(pytest.__version__)"
git -C $S/wt-8965a198 rev-parse HEAD
git -C $S/wt-8965a198 status --porcelain
cd $S/wt-8965a198 && python -m pytest -q -p no:cacheprovider tests/test_issue4_cmc_affine_root_contract.py tests/test_issue4_implicit_derivative_contract.py
cd $S/wt-8965a198 && python -m pytest -q -p no:cacheprovider -rf tests/test_issue4_cmc_lapse_tangent.py
git -C $S/wt-patch log --format='%H %s' -2
git -C $S/wt-patch status --porcelain
cd $S/wt-patch && RUNNER_TEMP=$OUT python -m pytest -q -p no:cacheprovider tests/test_issue4_research_port_validation.py
cd $S/wt-patch && python tests/data/issue4_affine_root/generate_reference.py 30 $OUT/root.json
cd $S/wt-patch && python tests/data/issue4_implicit_derivative/generate_derivative_reference.py 30 $OUT/deriv.json
cd $S/wt-patch && RUNNER_TEMP=$OUT python $S/research-diagnostic-run/compare.py
cd $S/wt-patch && RUNNER_TEMP=$OUT python $S/research-diagnostic-run/diagnostic_step.py
sha256sum $S/research-diagnostic-run/research_diagnostic.json $OUT/research_diagnostic.json
cmp $S/research-diagnostic-run/research_diagnostic.json $OUT/research_diagnostic.json
```

Both worktrees reported an empty `git status --porcelain` after the runs.

## Test results (base `8965a198`, both runs identical)

- Contract tests (`tests/test_issue4_cmc_affine_root_contract.py`,
  `tests/test_issue4_implicit_derivative_contract.py`): `7 passed`.
- Original witness file (`tests/test_issue4_cmc_lapse_tangent.py`): `2 failed, 16 passed`.
  Failing tests:
  - `tests/test_issue4_cmc_lapse_tangent.py::test_self_consistent_cmc_lapse_fixed_point_and_implicit_direction_two_epsilons`
  - `tests/test_issue4_cmc_lapse_tangent.py::test_fixed_point_root_matches_independent_affine_map_and_reports_conditioning`
- Research port validation (`tests/test_issue4_research_port_validation.py`, patched worktree): `2 passed`.
- Regeneration comparison step output (both runs):
  `affine-root reference identical: True | derivative reference identical: True`

## Diagnostic JSON (local)

- Paths:
  - Run 1: `/tmp/claude-0/-home-user/3f8452aa-83fe-5834-bdb1-93689a21c872/scratchpad/research-diagnostic-run/research_diagnostic.json`
  - Run 2: `/tmp/claude-0/-home-user/3f8452aa-83fe-5834-bdb1-93689a21c872/scratchpad/research-diagnostic-run2/research_diagnostic.json`
- Original patch (`1f724c0c`): the file hash did not reproduce between runs 1 and 2; the only
  differing field was the `base.d1` fingerprint, because `base.d1` is a Python function hashed as a
  0-d object array whose bytes are a memory pointer.
- Fixed patch (`4406b4e9`), file SHA-256 `4d3873cd67d1158e0b37751b9c9c82fa1f0d55297906dea59a61ed1d597eefdf`,
  identical in two runs from scripts extracted from the committed workflow
  (`git show 4406b4e9:.github/workflows/issue4-cmc-lapse-tangent.yml`), each run as
  `generate_reference.py 30`, `generate_derivative_reference.py 30`, `compare.py`,
  `diagnostic_step.py` from `$S/wt-patch` (HEAD `4406b4e9`) with `RUNNER_TEMP` set to the
  directory of each output file:
  - `/tmp/claude-0/-home-user/3f8452aa-83fe-5834-bdb1-93689a21c872/scratchpad/research-diagnostic-fixed-run1/research_diagnostic.json`
  - `/tmp/claude-0/-home-user/3f8452aa-83fe-5834-bdb1-93689a21c872/scratchpad/research-diagnostic-fixed-run2/research_diagnostic.json`
  - `input_fingerprint_counts` in both: `{"total": 23, "numeric_hashed": 22, "non_numeric_recorded": 1}`.
- Fields identical in all four runs:
  - `exact_equality_result_reported_by_existing_step`:
    `{"affine_root_identical": true, "derivative_identical": true}`
  - `fixture_source_sha256`: `838544416e36bf0bd0f1a96eaf6ac2b17fb4f59cf49da00755bc347714cbfb35`
  - `provenance.helper_sha256_regenerated_root`: `9ff5ca94d34c0ac397b22e1006007e67c59d1d9f47a662d952cada0eab546381`
  - `provenance.helper_sha256_regenerated_derivative`: `9ff5ca94d34c0ac397b22e1006007e67c59d1d9f47a662d952cada0eab546381`
  - `affine_root.max_abs_root_difference`: `0.0`; derivative `max_abs_difference` at epsilons
    `2e-05` and `1e-05`: `0.0`
- The regenerated `root.json` and `deriv.json` were byte-identical across the two runs
  (SHA-256 `fba1167a348cb780ec5873d1f5f88377bc55d07a6273ffa8008ceae006082053` and
  `e288f20550a43cb691ecd5e854a6287280c4ff0a374e7aa1180c6fd5700de630`).

## mpmath recorded in committed references (at `8965a198`)

- `tests/data/issue4_affine_root/issue4_affine_root_reference.json`, `provenance`: `"mpmath": "1.4.1"`
- `tests/data/issue4_implicit_derivative/issue4_implicit_derivative_reference.json`, `provenance`: `"mpmath": "1.3.0"`

## CI artifact comparison (workflow run 42497223250)

- Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42497223250 (push of
  `20ab5cb88deea2cd75fbb9458497ec29f1c7e3cf`). Job `research-diagnostics`: step "Regenerate references
  to temporary files and compare with committed data" concluded `failure`; the diagnostic and upload
  steps concluded `success`.
- Artifact `issue4-research-diagnostic-42497223250-1` (id 13624181680). Downloaded zip SHA-256
  `7ba133fb95df96e6ecf0fca2413ef373ce355506fe63dbcc55d952232cc5ef57`, equal to the digest GitHub reports.
  Extracted `research_diagnostic.json` SHA-256
  `0b7f4e90828a2bbc6dd6ce72bb2c008c57c3cf44ea6427650e166dd4a1625d96`, saved at
  `/tmp/claude-0/-home-user/3f8452aa-83fe-5834-bdb1-93689a21c872/scratchpad/ci-artifact-42497223250/extract/research_diagnostic.json`.
- Compared field by field with local fixed-run 1 (above):
  - `cpu_model`: CI `AMD EPYC 7763 64-Core Processor`; local `Intel(R) Xeon(R) Processor @ 2.30GHz`.
  - `platform`: CI `Linux-7.0.0-1017-azure-x86_64-with-glibc2.39`; `runner_image` CI
    `ubuntu24 20261101.360.1`.
  - Identical: `python` 3.12.15, `numpy` 2.5.4, `mpmath` 1.3.0, `python_implementation`, `machine`,
    `fixture_source_sha256`, both helper SHA-256 values, all four mpmath provenance values,
    `input_fingerprint_counts`.
  - Input fingerprints: 15 of 23 identical; 8 differ: `base.rho`, `base.pr`, `base.pt`,
    `rates.a_dot`, `rates.alpha_in_dot`, `rates.b_dot`, `rates.pr_dot`, `rates.rho_dot`.
  - CI `exact_equality_result_reported_by_existing_step`:
    `{"affine_root_identical": false, "derivative_identical": false}`.
  - CI `affine_root`: 19 of 20 components changed by string; `max_abs_root_difference`
    `2.9336762861e-20`; `cond_inf` identical (`624310584.1`).
  - CI `implicit_derivative`: 19 of 20 components changed at each epsilon; `max_abs_difference`
    `4.19435617432e-21` (eps 2e-05) and `1.7326678467376e-19` (eps 1e-05).

## AVX-512 test (local)

- Local NumPy (`np.show_runtime()`) reports SIMD extensions found: `X86_V3`, `X86_V4`, `AVX512_ICL`,
  `AVX512_SPR`.
- With `NPY_DISABLE_CPU_FEATURES="X86_V4 AVX512_ICL AVX512_SPR"`, `np.show_runtime()` reports found
  `X86_V3` and not found `X86_V4`, `AVX512_ICL`, `AVX512_SPR`. Under that setting, from `$S/wt-patch`
  (HEAD `4406b4e9`), the same generator, compare and diagnostic commands as the fixed runs above
  (`RUNNER_TEMP=$S/research-diagnostic-noavx512`) gave:
  - Compare step: `affine-root reference identical: False | derivative reference identical: False`,
    exit 1.
  - All 23 input fingerprints identical to the CI artifact.
  - `exact_equality_result_reported_by_existing_step` identical to CI (`false`, `false`).
  - Every `per_component_abs_difference` value (root, and derivative at both epsilons) identical to CI.
  - Diagnostic JSON SHA-256 `cf05eb783ddd81de966888931762e3bcb127d650646b88b77941d31314d340c2`.
  - Contract tests `7 passed`; witness file `2 failed, 16 passed` (base `8965a198`).

## Open

- (a) CI on CPython 3.12.15 mismatches the committed references; the cause is unknown. Source:
  `docs/ISSUE4_DIAGNOSTIC_HANDOFF_2026-11-07.md` on branch
  `docs/issue4-diagnostic-handoff-2026-11-07` (commit `144d5b53`), which records workflow run
  42407237209 printing `affine-root reference identical: False` and
  `derivative reference identical: False`. Not reproduced locally: the local run above printed
  `True` for both.
- (b) mpmath provenance discrepancy: the root reference records mpmath 1.4.1; local regeneration
  used mpmath 1.3.0 and gave identical strings.
- (c) Answered: CI fingerprints and environment fields compared with local values; see
  "CI artifact comparison" and "AVX-512 test" above.
