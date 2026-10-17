# SI-agent first-read protocol: shared communications are part of the repository

**Required for every fresh SI agent before analysis, coding, or launching a run.**

This is a coordinated team effort. The repository's communication record is operational research context, not optional background. Skipping it causes duplicated experiments, conflicting branches, lost findings, and avoidable work.

## First-read checklist

1. Read the current `README.md`, this document, [the canonical research workflow](RESEARCH_WORKFLOW.md), and the relevant current sections of `docs/STATUS.md`.
2. Read the [cross-agent message board, Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5), beginning with the newest coordination entries and following their links. Because this thread may exceed 1,000 comments, use the scalable history-reading method below; do not rely on the issue body alone or assume the entire thread must be loaded in one pass.
3. Read the [team task board, Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135) for current ownership, dependencies, existing artifacts, and exit criteria. Read [Issue #4 — Momentum-constraint propagation](https://github.com/miahpenn/0-Cosmos-engine/issues/4) when touching constraints, geometry, centre regularity, or related diagnostics.
4. Inspect relevant pull requests and their **conversation comments, review threads, commits, checks, and linked Actions artifacts**. Search for the same lane/discriminator and inspect current source SHAs before deciding there is no overlap.
5. Verify current `main` and the exact source SHA used by any reported test. A result from a stale or frozen branch is not automatically evidence about current production source.
6. Before starting work, post a coordination claim on Issue #5 (or the task's established shared thread): the bounded action, branch/PR and source SHA, existing evidence, overlap check, and scope exclusions. If another agent is already doing that exact test or edit, coordinate rather than duplicate it.

## Research provenance and no-duplicate-review policy

**Owner direction: apply this policy across research, implementation, documentation, and agent handoffs.**

- Treat the project's existing independently validated research and archive as the source of truth for claims already established there. Record the specific source/provenance once; do not require another independent scientific review or repeat a computation solely to re-prove an already-supported result.
- Research-backed findings, reports, computations, tests, and implementation changes are **merge-eligible by default** when the relevant source/provenance is traceable and applicable implementation/integration checks pass. Prefer automatic merge under that condition. Do not add a redundant scientific-review gate merely because an old task template or stale PR description asks for one.
- Keep checks that serve a distinct engineering purpose: implementation/build correctness, applicable exact-head CI, regression and interface behavior, source/dependency pin integrity, input/state safety, and stated scope limits. Reuse existing evidence when the inputs and relevant sources are unchanged. Rerun only the affected check when a relevant source, dependency, input, configuration, interface, or behavior has materially changed.
- “Correct computations do not need validation” means established scientific results and identical computations should not be independently revalidated as a ritual. A new code implementation can still receive focused checks to detect transcription, sign/index, parity/centre, integration, or regression defects; those are implementation checks, not a new scientific review of the established result.
- Escalate only a genuinely new or unresolved scientific claim, an unexplained discrepancy from the established result, missing provenance, failed implementation checks, a material change in physical assumptions/production behavior, or an explicit admission/safety/security gate. Vendor output is a comparator, not the scientific authority.
- This workflow policy does **not** change scientific-admission criteria, authorize unproven physical-law changes, or imply model admission. Preserve `scientific_admission=false` wherever currently required and do not treat an archive citation or green CI alone as scientific admission.

If an older issue, PR body, checklist, or handoff calls for duplicate independent review of a result already established in the archive, flag the conflict and align the task wording with this policy instead of starting another review cycle.

## Keep the message board alive

- Put durable cross-agent findings, blockers, hypotheses, test links, artifact links, exact source SHAs, and next-step recommendations in the shared issue/PR discussions. Do not leave important context only in a private chat, transient run log, or unlinked branch.
- After a run, report **what was tested, what source it used, what passed/failed, what the evidence does and does not establish, and the next gate**. Link the workflow run and artifacts.
- Separate observation from interpretation. A unit/synthetic test, lightweight CI pass, or source-only algebraic check is not proof of full-trajectory or physical validity.
- Preserve negative results and failed hypotheses; they prevent the team from repeating work.
- Keep physics changes separate from diagnostics, test harnesses, and campaign orchestration changes. Do not promote a hypothesis to production physics without the required admission evidence and review.
- Never silently overwrite another agent's branch or assume another agent's task is abandoned. Coordinate first.
- If a finding changes, post a correction with the affected commit/run and explain what is superseded.

## Message freshness and invalidation contract

**Treat freshness as a property of the claim and its dependencies, not just the age of the comment.** A recent comment can still be wrong if it names an obsolete branch head; an older numerical result can remain useful historical evidence if its full provenance is intact.

### Operational state expires quickly

- Treat snapshots of open PRs, branch heads, check/run status, active ownership, queued work, and “what should run next” as **stale after 5 minutes**. Before acting on one, refresh the corresponding live GitHub state. If it cannot be refreshed, mark it `UNKNOWN—REVALIDATION REQUIRED`; do not infer that a task is still running, finished, abandoned, or safe to duplicate.
- Immediately before any expensive test, campaign, artifact generation, or source edit, run a **preflight no more than 2 minutes old**. Confirm the live `main` SHA; the exact PR/branch head; required check and workflow state; relevant source/vendor/dependency pins; inputs and configuration; relevant existing artifacts; and any active overlapping claim on Issue #5.
- If preflight is older than 2 minutes at launch, or any relevant ref/configuration changes during preflight, repeat it. If a dependency moves between preflight and launch, stop and reconcile rather than spending resources on the old plan.
- A claim's expiry is **not** permission to duplicate it. Expiry means ownership is unknown and must be checked with the owner/thread before starting overlapping work.

### Numerical and test evidence is pinned, not “fresh forever”

A result is reusable only for the source and setup it actually tested. Every result must name, as applicable:

- exact source commit SHA (and PR head plus base SHA for PR validation);
- vendor/submodule and relevant dependency revisions;
- command/workflow revision, configuration/input fingerprint, resolution and run parameters;
- workflow/run URL and artifact URL plus SHA-256 for durable generated artifacts;
- observed result, interpretation, known limits, and next gate.

A result is **invalidated for the intended decision** if any relevant source, vendor/dependency pin, command/workflow, input, configuration, or run parameter differs—or cannot be verified. Revalidate or rerun only the affected discriminator. Do not rerun an expensive campaign by default when a small source/status check will establish whether the existing evidence still applies.

- “CI passed” is incomplete without the exact tested head and run links. Passing on a previous head does not qualify the new head.
- A result pinned to a research branch is not automatically evidence for `main`; a board snapshot is not proof of live state.
- Do not change or re-label historical findings to look current. Preserve them and add a correction/superseding result that points to the original.
- A timestamp alone never supersedes stronger provenance. If a newer note conflicts with a result pinned to a different SHA, stop and reconcile the underlying refs and artifacts before acting.

### Message types and required fields

Prefix coordination comments with one of these labels so agents can scan the newest stream:

- `[CLAIM]` — bounded task claimed, overlap checked, no expensive work started yet.
- `[RUNNING]` — current action and exact source/config being used.
- `[RESULT]` — completed action with exact evidence and limits.
- `[CORRECTION]` — corrects or invalidates an earlier claim/result; links the superseded comment.
- `[BLOCKER]` — a dependency or decision needed before safe progress.
- `[HANDOFF]` — ownership, current state, next bounded action, and acceptance gate.

Include the applicable fields below in the comment body. Use UTC ISO-8601 timestamps; never use relative phrases such as “just now” or “current” without an exact ref.

```text
[type] short summary
as_of_utc: YYYY-MM-DDTHH:MM:SSZ
lane:
owner_or_agent:
status: ACTIVE | BLOCKED | COMPLETE | UNKNOWN | SUPERSEDED
live_main_sha:
work_ref: branch=...; pr=#...; head_sha=...; base_sha=...
inputs_config_and_dependency_pins:
action_or_finding:
overlap_check:
evidence: run/PR/issue/artifact URLs; artifact sha256 where applicable
limits_or_invalidation:
revalidate_by_utc:
next_action_and_gate:
```

For a claim, fill in the planned action, overlap check, current SHA and next preflight deadline. For a result, replace the plan with observed outcomes, run links, artifact hashes and what the evidence does **not** establish. Mark non-applicable fields `N/A`; never omit a required provenance value silently.

When correcting information, use `[CORRECTION]`, link the affected claim/comment and run, state exactly what is invalidated and what remains valid, and provide the replacement source/run SHA. Do not delete the audit trail. When ownership or run state is unclear, use `UNKNOWN` and check directly instead of guessing.

## Chat and system-of-record boundary

Issue #5 is the canonical **durable cross-agent log**; Issue #135 is the task/ownership index; PR discussions and Actions runs are authoritative for their own review/check/evidence state. Keep these roles distinct and cross-link them. GitHub comments are asynchronous and do not guarantee immediate push delivery; do not describe them as a websocket or real-time chat transport. If the team later selects and connects a dedicated chat backend, use it for fast coordination and notifications, but keep decision-critical source pins, results, corrections and handoffs linked to the repository system of record. Never post credentials, access tokens or secrets into either channel.

## Scalable history reading for long threads

Issue #5 is the durable board, but with more than 1,000 comments, a full chronological dump is expensive and can obscure the live edge. Read it as an indexed research record:

- Start from the latest coordination entries and follow links to the live PR, exact-head workflow, review, and artifact.
- Search prior comments for the exact discriminator, issue/PR numbers, branch/source SHA, and the specific residual or hypothesis. Read the matching claim, result, corrections, and superseding result.
- Cross-check current Issue #135 and `docs/STATUS.md`; board snapshots can lag behind live GitHub state.
- If two recent notes conflict, report the conflicting SHAs and evidence and resolve the conflict before running anything.

This changes the reading strategy, not the requirement to review relevant context. When in doubt, ask in the existing shared thread rather than opening a duplicate lane.

## Current communication entry points

- [Issue #5 — Cross-agent coordination / shared message board](https://github.com/miahpenn/0-Cosmos-engine/issues/5)
- [Issue #4 — Momentum-constraint propagation](https://github.com/miahpenn/0-Cosmos-engine/issues/4)
- [Pull requests and their full discussions](https://github.com/miahpenn/0-Cosmos-engine/pulls)
- [GitHub Actions runs and artifacts](https://github.com/miahpenn/0-Cosmos-engine/actions)

**Definition of done for agent onboarding:** the agent has read the latest shared communications, identified the exact source/run relevant to its task, and posted its intended work or findings where the rest of the team can see them. Before expensive work begins, the agent has also passed the 2-minute source/status preflight above.
