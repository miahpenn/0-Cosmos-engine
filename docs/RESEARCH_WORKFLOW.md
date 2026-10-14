# Research workflow: coordinated, evidence-first changes

This is the canonical operating workflow for research and production-readiness work in 0-Cosmos Engine. It complements, rather than replaces, `docs/STATUS.md`, Issue [#135 task board](https://github.com/miahpenn/0-Cosmos-engine/issues/135), Issue [#5 cross-agent comms](https://github.com/miahpenn/0-Cosmos-engine/issues/5), and the task-specific issue/PR where technical evidence belongs.

## 1. Orient from the live tree

Before analysis, editing, or any run:

1. Read `README.md`, this guide, `docs/SI_AGENT_FIRST_READ_COMMS_PROTOCOL.md`, and the relevant current section of `docs/STATUS.md`.
2. Read Issue #135 for current lane ownership, dependencies, existing artifacts, and exit criteria.
3. Read the latest Issue #5 updates and the task's specific issue/PR comments, reviews, commits, checks, and artifacts.
4. Verify current `main` and the exact production-source SHA. Do not copy a hash from a stale board snapshot and call it current.
5. Search open PRs, branches, issue comments, and prior artifacts for the same discriminator before claiming work.

### Reading a very long message thread

Issue #5 is a durable research record and may exceed 1,000 comments. **Do not assume that reading the entire thread in one pass is efficient or necessary, and do not rely on the issue body alone.** Instead:

- Read the newest coordination entries first and follow their links to the current task, PR, review, workflow, and artifact.
- Search the history for the exact lane/discriminator, PR/issue numbers, source SHA, key residual or hypothesis, and previous claim/result. Read the matching original context, including corrections and superseding results.
- Check the task board and current status file for stale snapshots; prefer the newest verified source/run evidence and explicitly mark older entries historical.
- Before changing course, post a short claim or correction linking the evidence that changed the decision.

The goal is to recover current decisions and relevant history—not to re-read unrelated work or mistake a historical claim for current ownership.

## 2. Claim one bounded lane before work

Post in Issue #5 or the task's established thread **before** testing, editing, or launching a run. Include:

- Lane/question and one bounded action.
- Current `main` SHA, intended branch/PR, and production-source anchor when different.
- Existing relevant artifact/run and prior result.
- Overlap check: closest active PRs/claims and how this action differs.
- Scope exclusions (especially production physics, vendor code, campaign settings, and long runs).

If another agent owns the same discriminator, coordinate or extend that work in its existing PR rather than opening a competing lane. The task board is the ownership index; avoid inventing parallel trackers or assignees.

## 3. Pick the cheapest informative next step

Use this preference order when scientifically appropriate: read existing evidence offline → source-pinned algebra/operator audit → frozen-state or synthetic negative-control test → bounded one-step/small diagnostic → expensive coupled campaign only after a specific signed identity fails and the proposed discriminator could separate plausible causes.

Only one expensive campaign/run should be active at a time. A rerun is justified when its source inputs were invalidated, a new discriminator is genuinely added, or the prior run failed before answering the question. Record a stop condition before starting.

Keep production equations, diagnostics, tests, report generation, and campaign orchestration separate. Never tune to a physical residual, infer a missing law from a mismatch, or silently add a source, coefficient, clamp, floor, damping term, boundary change, or admission relaxation.

## 4. Keep source and evidence provenance explicit

For every result, record:

- PR source-head SHA and the production-source/vendor anchors consumed by the test or diagnostic.
- Actual workflow run IDs and status at that exact head; distinguish a missing/skipped check from a pass.
- Artifact URL(s), ZIP and inner-report hashes where available, and the command/configuration needed to reproduce the bounded result.
- The observation, its interpretation boundary, failed hypotheses, and the next falsifiable gate.

If an Actions execution checkout is a synthetic merge or differs from the PR source head, record both. A passing test suite, smoke, symbolic identity, frozen-state check, or short-horizon run proves only its stated software/diagnostic contract—not full trajectory validity or scientific admission.

Never erase a negative result. If new evidence supersedes a claim, post a correction that links both the old claim and the exact new evidence.

## 5. Open and maintain a focused PR

- Use a descriptive branch and one coherent scope per PR. Start from verified current `main`; if a candidate is intentionally source-pinned to an older snapshot, explain why and do not present old-source evidence as current-main evidence.
- Mark research/diagnostic work draft until the independent review and stated exit gate are complete. Declare stacked dependencies; do not promote a child before its parent is resolved.
- Fill in `.github/pull_request_template.md`: lane claim, overlap, source hashes, exact-head checks/artifacts, limits, and next gate.
- Re-run applicable checks after the final commit. Paste exact run URLs rather than saying only “CI green.”
- Do not merge because a branch is green or reported mergeable. Scientific admission stays separate. Issue [#112](https://github.com/miahpenn/0-Cosmos-engine/issues/112) is the governance gate: actual branch protection, required exact-head contexts, independent review, and the blocked-merge negative test require owner/admin verification. Do not claim those are enforced until evidence proves it.

## 6. Close the loop in shared comms

After a bounded task, post a result on its PR/issue and summarize the cross-lane delta on Issue #5. Update Issue #135 when ownership, dependencies, exit criteria, or durable state changes. Update `docs/STATUS.md` only when a new verified engineering snapshot is warranted; do not copy every transient run into a second status log.

Use this handoff shape:

```text
Result / disposition:
Task and overlap:
PR source SHA / production-source anchor:
Checks (exact head; run URLs; pass/fail/skipped):
Artifact URL and hashes:
What the evidence establishes:
What it does not establish:
Next bounded gate / owner:
```

Done means the durable record contains the exact provenance and honest limits, an independent reviewer can reproduce the decision, and the next agent can continue without reconstructing the entire chat.
