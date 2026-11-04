# Issue #4 — zero-step accepted-state observer execution handoff

Date: 2026-11-04  
Research branch: `research/issue4-zero-step-bianchi-observer-2026-11-04`  
Workflow implementation commit: `6b8422d1c94822cca44ee1cccf3d93581ab2c21e`  
Status: workflow prepared on the research branch; **not yet dispatched**.

## Purpose and scope

This handoff records the execution wiring for the already-implemented accepted-state full-Bianchi evaluator. It is not a new derivation, a physics review, or a three-class comparison.

The next execution is limited to the archived accepted state at (t=0.60), resolutions (N=60,120,240). It runs the focused observer contract test, then the committed evaluator against the exact original archives. No production step, RHS evaluation, initializer, trajectory, production/vendor edit, or scientific-admission change is authorized by this workflow.

The evaluator's current implementation covers the accepted-state class only. Manufactured regular profiles and frozen production (t=0) are not implemented by this evaluator and must not be described as evaluated.

## Three distinct checkouts

The three roots are not interchangeable:

| Checkout | Exact identity | Runtime purpose |
|---|---|---|
| `runner/` | The dispatched research-branch SHA; observer file blobs are checked explicitly | Supplies `issue4_full_bianchi_evaluator.py`, fixture loader, reference regression and focused tests |
| `diagnostic-source/` | `4aec88d0e3ff028e0a65a8ac36b30948625df792` | Supplies the reviewed total-stress helper, engine modules and vendor checkout |
| `helper-baseline/` | `7f47669a1c0edfafef9fc5f27b2db80840d3f035` | Separately pinned protected helper baseline; not a replacement runtime engine |

The vendor submodule is pinned to `d6052d605673ce9d82cdc99b0df79855fa2ea215` in both source trees. The workflow checks both submodule HEADs. The evaluator also checks the accepted source commit, protected helper baseline, vendor commit, required protected source blobs and reviewed helper blobs before loading an archive.

The runner's current protected blob checks are:

- `tools/issue4_full_bianchi_evaluator.py`: `bdd7574807def8a612e5d39c00c97629d0511ae5`
- `tools/issue4_accepted_state_fixture.py`: `c9ec1bbd4373e3f96cd4dac7dcb9c91c1014b2ba`
- `tools/issue4_bianchi_reference_regression.py`: `fff8882ff09b3f9f1d1f86c4b331873d87a7d729`
- `tests/test_issue4_full_bianchi_evaluator_contract.py`: `cc438c281e056a5ec2b9d6b5a0fc9006fa4a7e6f`

## Accepted-state archives

All three original artifacts are retained from workflow run `40236402303`. The artifact ZIP bytes must be downloaded intact because `load_accepted_fixture` validates the outer ZIP SHA before opening its members. Do not use a download path that extracts and then recreates the ZIP; recompression changes the bytes and fails the archive identity check.

| Resolution | Artifact ID | Artifact name | Required ZIP SHA-256 |
|---:|---:|---|---|
| 60 | `12615136565` | `issue4-exact-profile-momentum-tangent-N60` | `41d646090e05c08684806c5e39d57851f119cf0f38aa90f1d26ab5f72f87d1ce` |
| 120 | `12615555594` | `issue4-exact-profile-momentum-tangent-N120` | `9f5c2504f28290ff540c96f93e1fb569f5b960b7a2b3fc312dca9362c8e2eacf` |
| 240 | `12615550675` | `issue4-exact-profile-momentum-tangent-N240` | `e8a9cd26fbe3dbfb0bd62fc0606237c387a41efd14ade550df50da9842b38bea` |

The archive loader separately validates each archive's internal member hashes, array count, shapes, dtypes, per-array SHA-256 values, accepted source/vendor/helper pins, accepted time, upstream state digest and tangent-helper state digest. The expected upstream state digests are:

- N=60: `2ae2a5704b2b87d16ef9239c8a22ef261ed9108d971d3a5c78a3d1133aa9e9e4`
- N=120: `0266c7fafc5e9a292a9a376f1d91fa45bb8dd103b6f8c21b543c25d0930daacb`
- N=240: `770d36f3db5b60da87f52beb213f8af4731e17021f95c1519ecbdbe759b41cce`

The artifacts were listed as unexpired on 2026-11-04; their recorded expiration is 2027-01-23.

## Exact evaluator interface

From the runner checkout, one resolution is evaluated as:

```bash
python -m tools.issue4_full_bianchi_evaluator \
  --archive ../archives/issue4-exact-profile-momentum-tangent-N60.zip \
  --resolution 60 \
  --diagnostic-source ../diagnostic-source \
  --helper-baseline ../helper-baseline \
  --output ../results/accepted-N60.json
```

Repeat for 120 and 240 with the matching archive and resolution. The workflow already contains these three invocations.

The helper import is intentional: the evaluator loads `audit_total_stress_energy_divergence.py` from the accepted diagnostic source, sets `COSMOS_PINNED_SOURCE_ROOT` to the protected helper-baseline checkout during import, and restores the prior environment afterward. The helper imports runtime engine modules from its own diagnostic-source root. Do not collapse these roots or reverse their roles.

## Workflow dispatch and safety

The registered workflow path is the existing `.github/workflows/issue4-independent-total-stress-energy-divergence.yml`, which already exists on the default branch and already declares `workflow_dispatch`. GitHub requires a manually dispatched workflow to be registered on the default branch, but permits selecting a non-default branch for the run.

To avoid creating a new unregistered workflow or opening a PR that could activate unrelated checks, this research-branch version adds a branch-gated job:

- Manually dispatch this existing workflow with the branch set to `research/issue4-zero-step-bianchi-observer-2026-11-04`.
- On that exact manual-dispatch ref, the older independent-divergence matrix job is skipped and the accepted-state zero-step observer job runs.
- For other events/refs, the existing job retains its previous behavior.
- Do not open a PR merely to make this job dispatchable. The default-branch workflow already provides registration; the selected branch supplies the branch-specific workflow content.

Permissions are read-only: `contents: read` and `actions: read`. Checkouts use `persist-credentials: false`. The workflow downloads the original artifact ZIPs through the GitHub Actions artifact API, verifies their pinned SHA-256 values, captures runner/CPU/BLAS/package provenance, runs the focused contract tests, evaluates only the three frozen accepted archives, and uploads archives/reports/environment data even on failure.

The recorded Python/dependency environment follows the earlier diagnostic-runner precedent documented in project Slack: Python 3.12.15, NumPy 2.5.4, SciPy 1.18.1, SymPy 1.14.0, mpmath 1.3.0, pytest 9.1.1, plus threadpoolctl. The actual CPU model, flags, BLAS dispatch and installed package versions are recorded from the new runner; historical runner hardware is not assumed.

The workflow asserts for every resolution:

- status is `diagnostic_evaluation_complete`;
- the exact upstream accepted-state digest matches;
- all three source identities and the vendor pin match;
- `production_steps_run == 0`;
- `state_advanced == false`;
- `zero_step_tripwire == armed_and_not_triggered`.

It does not contain a production smoke, production step, RHS call, initializer, evolution loop, or admission-flag update. The scientific-admission and conservation-closure flags remain false.

## Current boundary

This workflow is prepared on the research branch but has not yet been dispatched. A successful run will establish the committed evaluator's execution against the accepted-state archives under its exact pins. It will not complete the original three-class comparison, establish a missing physical term, or change scientific admission. The next class extensions remain separate work under the same observer contract.
