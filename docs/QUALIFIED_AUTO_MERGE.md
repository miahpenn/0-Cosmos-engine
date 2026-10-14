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
- `.github/workflows/em_total_energy_budget.yml`

The workflow checks the latest attempt for each exact check name, rejects missing, queued, running, skipped, cancelled, timed-out, or failed required checks, and requires the GitHub Actions check-run app. A result from an earlier SHA never qualifies a newer SHA.

## Fast path for diagnostic-only work

A PR is eligible for the diagnostic fast path only when **every** changed file is in `docs/**`, `research/**`, `tests/test_research_*.py`, `tests/test_diagnostic_*.py`, `tools/audit_*.py`, or `tools/recompute_*_artifact.py`. These PRs do not need a formal approval, but still need all applicable exact-head checks and no current-head `CHANGES_REQUESTED` review.

## Other changes

Any PR outside that allowlist requires an `APPROVED` review from a reviewer other than the author, submitted against the exact current head SHA. Current-head `CHANGES_REQUESTED` reviews block both paths. Workflow changes, production engine/vendor changes, dependency changes, and general tools/tests therefore do not bypass review.

## Freshness and merge safety

Before merging, the workflow verifies that the PR is open, not a draft, from this repository (fork PRs are manual), conflict-free, and based on the current `main` SHA. It rechecks the head/base, required checks, and reviews immediately before calling GitHub's squash-merge API with the expected head SHA. PR events, completion of the relevant CI workflows, and a 10-minute safety sweep trigger reevaluation. A failed, missing, stale, or ambiguous qualification never causes a merge.

The automation does not alter equations, widen tolerances, interpret a green smoke test as scientific validation, or mark a numerical model as admitted. Changes to this policy or to workflow definitions remain review-sensitive.
