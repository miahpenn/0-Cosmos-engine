# Issue #4 — zero-step accepted-state observer execution handoff

Date: 2026-11-04  
Research branch: `research/issue4-zero-step-bianchi-observer-2026-11-04`  
Workflow implementation commit: `6b8422d1c94822cca44ee1cccf3d93581ab2c21e`  
Status: first dispatch refused at N=60 due to an optional matter-schema mismatch; research-only adapter correction committed, awaiting rerun.

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

The first dispatch was attempted and refused at N=60 for the adapter schema mismatch documented below. The corrected evaluator and matching workflow blob pin are now committed on the research branch; a fresh manual dispatch is required to execute the corrected SHA. A successful run will establish the committed evaluator's execution against the accepted-state archives under its exact pins. It will not complete the original three-class comparison, establish a missing physical term, or change scientific admission. The next class extensions remain separate work under the same observer contract.


## First dispatch result and schema correction (2026-11-04)

The first manual dispatch was run at research-branch SHA
`c91bc01c551287f6faf47f71b392123b7c4cb52d` (workflow run
`41932751640`). Provenance gates, archive downloads/checksums, environment
capture, and focused contract tests passed. The N=60 evaluator then refused
before producing a Bianchi result:

```text
AttributeError: 'types.SimpleNamespace' object has no attribute 'charged_dust'
production_steps_run: 0
state_advanced: false
```

The uploaded artifact was inspected directly. All three accepted archives
(N=60/120/240) contain the same 101-array schema: frozen geometry, six scalar
arrays, and three fluid species (dark matter, baryons, radiation), plus tangent
and total-stress-divergence witnesses. No archived array path or NPZ key refers
to electromagnetic fields, charge, Maxwell, or charged dust. The six-sector
tangent witnesses are S, D, COSMOS-phi, dark matter, baryons, and radiation.

The pinned production source distinguishes two optional structures:

- `V55MatterState.charged_dust` is an optional additional pressureless fluid
  carrier. The canonical fluid projection accesses this slot and treats
  `None` as absent.
- `total_matter_projection(..., electromagnetic_projection=...)` accepts a
  separate, already consumer-normalized EM stress/current projection
  (`rho, pr, pt, j`). Passing no projection means no EM contribution is
  included; raw electric fields or charge normalization are not inferred.

Therefore the refusal is an adapter schema mismatch: the research fixture
constructed a `SimpleNamespace` with the three archived species but omitted
the optional `charged_dust` attribute expected by the pinned canonical
projection. The evaluator adapter now explicitly supplies
`charged_dust=None`; it does not create a synthetic species or an EM source.
The workflow's evaluator blob pin was updated to the corrected blob
`4fa0e6d8b23304e06938073db8c2d34281991f24`.

This accepted-state archive set is **not an EM-on dataset**. That does not
mean the codebase lacks an EM formulation: the pinned source has a separate
radial-Maxwell state and normalized electromagnetic projection path. It means
the current accepted-state observer cannot claim to have measured or closed
that contribution. The present run must be described as the archived
six-sector accepted-state comparison, with EM-on accounting remaining outside
its scope.

The correction has not yet been rerun. The first dispatch stopped at N=60
under `set -e`, so N=120 and N=240 were not evaluated. No production step,
RHS, initializer, source/vendor edit, or scientific-admission change occurred.


## Corrected accepted-state dispatch — successful (2026-11-04)

This section supersedes the pending-rerun statements above. The corrected
workflow was manually dispatched successfully:

- Workflow run: [41935352237](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/41935352237)
- Run conclusion: `success`; event: `workflow_dispatch`; attempt: 1
- Branch: `research/issue4-zero-step-bianchi-observer-2026-11-04`
- Runner commit: `e15afa167190e02c5ae554edf1e1d222109a2ba7`
- Corrected evaluator blob: `4fa0e6d8b23304e06938073db8c2d34281991f24`
- Diagnostic source: `4aec88d0e3ff028e0a65a8ac36b30948625df792`
- Helper baseline: `7f47669a1c0edfafef9fc5f27b2db80840d3f035`
- Vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`

The `accepted-state-zero-step-observer` job completed successfully. All
provenance checks passed, focused contract tests passed (6/6), and the
evaluator completed for N=60, 120, and 240. The secondary
`bounded-independent-divergence` job was skipped by the branch-gated
workflow; it did not fail and is not part of this accepted-state observer
result.

Artifact:
- ID: `13374196058`
- Name: `issue4-zero-step-accepted-state-observer-41935352237`
- SHA-256: `31180befe6a50b9fa4fc954f98feecaddf3fa187f4c3e4ac83bdfa4e2d3f46cf`
- Contains the three original accepted archives, environment provenance,
  and `accepted-N60.json`, `accepted-N120.json`, `accepted-N240.json`.

### Verified invariants and cross-checks

For every resolution:
- `status = diagnostic_evaluation_complete`;
- accepted state class is `accepted_t_0p60`;
- original upstream state digest matches the pinned expected digest;
- accepted diagnostic source, helper baseline, and vendor pins match;
- `production_steps_run = 0`;
- `state_advanced = false`;
- input-state digest before and after evaluation is identical;
- `zero_step_tripwire = armed_and_not_triggered`;
- accepted-reference regression passed;
- canonical matter projection cross-checks for `rho`, `j`, `pr`, and `pt` passed;
- metric-only Ricci / evolved-connection offset cross-checks passed (maximum
  reported gaps are below 1e-18).

### Residual and sector-accounting observations

The maximum absolute tangent-minus-predictor residual over r<=10 decreases
with resolution in both independently retained divergence channels:

| N | Pointwise max | Finite-volume max |
|---:|---:|---:|
| 60 | 1.8406e-5 | 9.6404e-6 |
| 120 | 1.5272e-6 | 1.2827e-6 |
| 240 | 3.9897e-7 | 4.0771e-8 |

The corresponding maximum annular (5<=r<=10) residuals also decrease:

| N | Pointwise annular max | Finite-volume annular max |
|---:|---:|---:|
| 60 | 2.0498e-6 | 1.8089e-6 |
| 120 | 1.6403e-7 | 1.3547e-7 |
| 240 | 1.0667e-8 | 8.9879e-9 |

These are resolution-trend observations for the archived accepted-state
comparison, not proof that all centre-cell behavior or the continuum
constraint-propagation identity is closed.

The six-sector ledger names dark matter, baryons, radiation, S, D, and
COSMOS-phi. Its raw six-sector sum does not by itself reconstruct the total
(`raw_sum_passed=false`). After the explicit time-density component
accounting, the corrected reconstruction passes in all four channels at both
epsilon values; the largest unexplained gap across the three resolutions is
of order 1e-20. Preserve both facts in downstream reporting: do not describe
the raw sum as closed, and do not omit the successful component correction.

### Scope boundary and next step

The artifact's structural inventory records that the full-Bianchi predictor
terms are reconstructed from the frozen state under the pinned operators;
they are not present as pre-stored named arrays. The archive contains the
six-sector accepted-state data only. The codebase's separate radial-Maxwell
and normalized EM projection pathway remains present, but these accepted
archives do not carry an EM-on witness/projection. This run therefore does
not test or close an EM contribution.

This successful run verifies the accepted-state observer execution under
the exact pins. It does not evaluate manufactured regular profiles or frozen
production t=0, and it does not change scientific-admission or
conservation-closure flags. Next authorized sequence: extend the same
verified observer contract to the manufactured regular-profile class and
frozen production t=0 class, then compare all three classes under one
contract. Keep production steps at zero and retain all existing research-only
guardrails.


## Readiness audit for the two remaining comparison classes (2026-11-04)

After the accepted-state run succeeded, the existing repository and available
Library handoffs were checked for reusable manufactured and frozen-production
t=0 inputs. This is an input/provenance audit only; no new evaluator or
physics calculation was run.

### Manufactured regular-profile class: existing Ricci oracle is not the Bianchi fixture

The existing
`tools/audit_manufactured_regular_ricci_oracle.py` and
`docs/ISSUE4_MANUFACTURED_REGULAR_RICCI_ORIGIN_2026-10-16.md` provide a
useful analytic geometry control. Its regular-even profile is
`B(r)=1+b2*r^2+b4*r^4+b6*r^6`, with the explicit parity-breaking arm
retained as a negative control. However, this is a Ricci/operator oracle,
not a full matter-coupled Bianchi fixture:

- it does not archive the six matter/scalar stress sectors, total-stress
  divergence witnesses, or an independently supplied M tangent required by
  the accepted-state evaluator;
- its current implementation calls `V55ProductionKernel.initialize()`
  to obtain a base grid/state. That is outside the zero-step observer's
  strict no-initializer contract, even though the diagnostic itself performs
  no time evolution;
- its report explicitly disclaims a Brown/Bianchi propagation proof.

Therefore do not route this existing Ricci oracle through the accepted-state
evaluator or describe it as the manufactured Bianchi class. Reuse its analytic
profile and negative-control mathematics only as a possible input definition,
after the full Bianchi witness contract is independently supplied.

### Frozen production t=0: historical tangent artifact is not pin-compatible

A Library handoff records a historical t=0 discrete-constraint tangent
artifact: workflow `38290633971`, artifact `11766405316`, with N=60/80/120
baseline/candidate states and symmetric directional tangents. Its production
source pin is `6d8b3f2790644591326df57c55908863f994cb0c`, not the currently
accepted diagnostic source pin
`4aec88d0e3ff028e0a65a8ac36b30948625df792`. It also uses a different
resolution set and an older artifact schema. It is useful historical context,
but cannot be silently admitted as the current pinned t=0 class.

The current 2026-10-24 accepted-state artifacts verified in the Library and
the 2026-11-04 accepted-state observer are at t=0.60, not t=0. No matching
current-pin t=0 archive with the complete frozen state, total-stress
divergence witnesses, and tangent arrays was located in the available
handoffs.

### Required next implementation boundary

The next change must be a fixture/schema extension, not an invocation of the
production initializer or RHS inside the observer:

1. Identify or create a source-pinned, immutable manufactured full-state
   fixture that includes the complete frozen geometry, all six named sectors,
   pointwise and finite-volume total-stress divergence witnesses, and signed
   tangent witnesses for both epsilon channels. Keep the analytic regularity
   and parity-breaking controls explicit.
2. Identify a saved production t=0 state archive under the current accepted
   diagnostic source/vendor/helper identities, or produce a new archive in a
   separately authorized capture workflow. The observer itself must only
   load and verify it; it must not call the initializer, RHS, or evolution.
3. Extend the fixture loader with explicit schema/class dispatch and
   fail-closed class-specific provenance/digests. Keep the existing accepted
   archive loader and regression path unchanged.
4. Run the existing accepted-state regression unchanged, then run the
   manufactured and t=0 classes with the same predictor implementation,
   signed term names, epsilon channels, centre/annulus/fixed-radius summaries,
   state immutability digest, and zero-step tripwires.
5. Only then produce the three-class comparison. Do not mix historical
   source pins, infer absent arrays, synthesize matter/EM contributions,
   compute missing tangent arrays inside the observer, or relax any gates.

Disposition: the accepted-state class is verified; the other two classes are
not yet evaluable under the same contract from the currently located inputs.
This is a concrete fixture/provenance prerequisite, not a physics failure.
Scientific admission, full propagation closure, and total-stress-energy
conservation flags remain false.
