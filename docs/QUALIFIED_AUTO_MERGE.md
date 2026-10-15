# Qualified auto-merge policy

This workflow is a merge helper, not scientific admission and not a substitute for review. It must fail closed whenever exact-head qualification cannot be established.

## Exact-head checks required for every eligible PR

Every open, non-draft, same-repository PR targeting `main` must have the latest GitHub Actions check runs completed successfully on its current head SHA:

- `regression-suite`
- `Tests and one-step production smoke`
- `integrated-admission`

A PR touching any of the following must also pass `joined-energy-ledger` on that same head:

- `tools/audit_em_total_energy_budget.py`
- `tools/recompute_fixed_r_balance_artifact.py`
- `tests/test_research_em_total_energy_budget_contract.py`
- `tests/test_research_fixed_r_independent_recompute.py`
- `.github/workflows/em_total_energy_budget.yml`

A check is not qualified merely because a check-run object is green. The workflow also requires a matching completed pull-request workflow run for the expected workflow, with an explicit association to this PR number, head SHA, base ref `main`, and the exact current base SHA. If GitHub omits that association, the result is ambiguous and the PR stays blocked. Timestamp comparisons are not a substitute for an exact PR/base association.

## Documentation-only no-formal-review fast path

Only changes wholly limited to documentation text qualify for the no-formal-review fast path: files under `docs/` with `.md`, `.rst`, `.txt`, or `.adoc` extensions, plus the explicit documentation files `README.md`, `CONTRIBUTING.md`, `CHANGELOG.md`, `SECURITY.md`, `CODE_OF_CONDUCT.md`, and `.github/pull_request_template.md`.

The fast path still requires all applicable exact-head checks, a current `main` base, clean mergeability, and no current-head `CHANGES_REQUESTED` review. Research code, research/diagnostic tests, audit/recompute tools, engine/vendor code, dependency changes, workflow changes, and all other paths require an `APPROVED` review from someone other than the PR author, submitted against the exact current head SHA. Current-head `CHANGES_REQUESTED` reviews block both paths.

## Freshness and merge safety

Before merge, verify that the PR is open, not draft, from this repository, conflict-free, and based on the live `main` SHA. Recheck head, base, required checks and reviews immediately before the squash merge, and pass the expected head SHA to GitHub's merge API so an intervening PR-head update rejects the merge. The automation does not claim to replace repository branch protection and cannot close a main-branch race by itself.

The workflow uses `pull_request_target` with write permissions only for API metadata and the qualified merge request. It must not check out or execute pull-request code. A green smoke test or unit suite does not imply scientific admission; scientific readiness remains a separate, explicit gate.
