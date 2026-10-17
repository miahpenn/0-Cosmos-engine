# Qualified auto-merge policy

This workflow is a merge helper, not scientific admission and not a substitute for source review. It is intentionally fail-closed.

## Required checks on the exact current PR head

Every open, non-draft, same-repository PR targeting `main` must have these newest GitHub Actions check runs completed successfully on the PR's current head SHA:

- `regression-suite`
- `Tests and one-step production smoke`
- `integrated-admission`

A PR touching any of the following must also pass `joined-energy-ledger` on that exact head:

- `tools/audit_em_total_energy_budget.py`
- `tools/recompute_fixed_r_balance_artifact.py`
- `tests/test_research_em_total_energy_budget_contract.py`
- `tests/test_research_fixed_r_independent_recompute.py`
- `tests/test_research_em_adjoint_face_work_v1.py`
- `research/em_adjoint_face_work_v1.py`
- `research/reconcile_em_face_work_artifact_v1.py`
- `.github/workflows/em_total_energy_budget.yml`

The privileged `actions/github-script` dependency is pinned to the verified full commit SHA for v7.1.0. Trusted reviewer authorization for workflow/config edits and current-head reviewer veto state are rechecked on every qualification and immediately before merge. The workflow selects the most recently started check run for each exact check name; attempt number is only a tie-breaker after run timestamp, so an older successful retry cannot outrank a newer failure. It rejects missing, queued, running, skipped, cancelled, timed-out, or failed required checks and requires the GitHub Actions check-run app. A result from an earlier SHA never qualifies a newer SHA. The selected check run and the corresponding completed pull-request workflow run must share the same GitHub check-suite ID. The workflow run must also have an explicit association to this PR number, head SHA, base ref `main`, and current base SHA. If either the check-suite or PR/base association is missing or mismatched, qualification fails closed; timestamps are not accepted as a substitute. Passing a head check that ran before main advanced is not sufficient.

## Research-provenance-first fast path

Established research-backed implementation, test, and diagnostic changes are merge-eligible by default when their provenance is traceable and the applicable exact-head implementation/integration checks pass. Do not require a repeated independent scientific review merely to revalidate an unchanged derivation already supported by the project's validated research/archive. The fast path still requires an open non-draft same-repository PR, the current `main` base, clean mergeability, exact-head and exact-base check evidence, applicable specialized gates (including the joined-energy ledger where relevant), unchanged head/base during qualification, and no trusted non-author human collaborator's current-head `CHANGES_REQUESTED` review.

The workflow retains its checks for a trusted human review on the exact current head when the PR changes `.github/workflows/**`, because workflow/config edits can weaken future validation and merge safeguards. Approval requires a non-author human account whose repository collaborator permission is verified as `write`, `maintain`, or `admin`. Bot accounts, non-collaborators, read-only/triage-only accounts, and reviews whose permission cannot be verified do not qualify; permission lookup errors fail closed. Only a trusted non-author human collaborator's current-head `CHANGES_REQUESTED` review blocks ordinary changes or a workflow/config PR. Untrusted reviews do not independently approve or veto promotion; their feedback remains available for maintainer consideration.

## Native geometry source-pinned gate

A PR that adds or changes files under `engine/native_geometry/**`, or changes `.github/workflows/tracefree-ricci-proper-radial.yml`, must also pass the `tracefree-ricci-proper-radial` check on the exact current head and current `main` base. That workflow verifies the pinned existing production file blobs and vendor SHA; it permits only the explicitly listed native candidate files when they are absent from the pinned production source. Any other engine/vendor/dependency diff fails closed. This prevents the native migration from being confused with a production-runtime switch and prevents that scope guard from being silently bypassed.

## Joined-ledger trigger parity

The `ledgerPaths` set used by the merge helper must match all paths in `.github/workflows/em_total_energy_budget.yml` exactly. A regression test compares the sets; any new path added to the joined-ledger workflow must also be protected by the merge-helper gate before this policy can pass.
 
## Scope and special gates

The fast path does not waive implementation, provenance, regression, integration, or scope gates, and it does not change scientific-admission criteria. A passing software suite is not scientific admission. Truly provisional source-law/index/normalization choices, unexplained mismatches, missing provenance, or production-physics changes outside the declared task scope must remain blocked by their applicable project gates rather than being inferred as accepted from a green CI result.

Workflow/config changes remain approval-gated as described above. For ordinary changes, the merge helper uses exact-head tests and workflow-to-check-suite/PR/base association, checks source/dependency-sensitive ledger paths, and fails closed on stale base, missing/failed checks, unknown mergeability, head movement, or a trusted current-head change request.

## Freshness and merge safety

Before merging, the workflow verifies that the PR is open, not a draft, from this repository (fork PRs are manual), conflict-free, and based on the current `main` SHA. All candidate qualification/merge runs share one concurrency group to reduce races between competing helper runs. It rechecks the head/base, required checks, and reviews immediately before calling GitHub's squash-merge API with the expected head SHA. PR events, completion of the relevant CI workflows, and a 10-minute safety sweep trigger reevaluation. A failed, missing, stale, or ambiguous qualification never causes a merge.

The automation does not alter equations, widen tolerances, interpret a green smoke test as scientific validation, or mark a numerical model as admitted. It is not a substitute for verified effective branch protection: manual/direct writes can bypass it, and GitHub's expected-head merge check cannot atomically lock the base SHA. Keep governance status separate. Changes to this policy or to workflow definitions remain review-sensitive.
