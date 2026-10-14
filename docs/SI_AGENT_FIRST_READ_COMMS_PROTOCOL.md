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

## Keep the message board alive

- Put durable cross-agent findings, blockers, hypotheses, test links, artifact links, exact source SHAs, and next-step recommendations in the shared issue/PR discussions. Do not leave important context only in a private chat, transient run log, or unlinked branch.
- After a run, report **what was tested, what source it used, what passed/failed, what the evidence does and does not establish, and the next gate**. Link the workflow run and artifacts.
- Separate observation from interpretation. A unit/synthetic test, lightweight CI pass, or source-only algebraic check is not proof of full-trajectory or physical validity.
- Preserve negative results and failed hypotheses; they prevent the team from repeating work.
- Keep physics changes separate from diagnostics, test harnesses, and campaign orchestration changes. Do not promote a hypothesis to production physics without the required admission evidence and review.
- Never silently overwrite another agent's branch or assume another agent's task is abandoned. Coordinate first.
- If a finding changes, post a correction with the affected commit/run and explain what is superseded.

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

**Definition of done for agent onboarding:** the agent has read the latest shared communications, identified the exact source/run relevant to its task, and posted its intended work or findings where the rest of the team can see them.
