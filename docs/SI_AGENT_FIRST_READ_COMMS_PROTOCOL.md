# SI-agent first-read protocol: shared communications are part of the repository

**Required for every fresh SI agent before analysis, coding, or launching a run.**

This is a coordinated team effort. The repository's communication record is operational research context, not optional background. Skipping it causes duplicated experiments, conflicting branches, lost findings, and avoidable work.

## First-read checklist

1. Read the current `README.md`, this document, and the current production/status and handoff documents relevant to your task.
2. Read the cross-agent message board: [Issue #5 — Cross-agent coordination](https://github.com/miahpenn/0-Cosmos-engine/issues/5). Read the full conversation, including the newest comments, not only the issue description.
3. Read [Issue #4 — Momentum-constraint propagation](https://github.com/miahpenn/0-Cosmos-engine/issues/4) when touching constraints, geometry, centre regularity, or related diagnostics.
4. Inspect recent and relevant pull requests and their **conversation comments, review threads, commits, checks, and linked Actions artifacts**. Start with the open PRs/issues linked from the coordination board; do not assume a branch name or PR title tells the full story.
5. Check the current `main` commit and the exact source SHA used by any reported test. A result from a stale or frozen branch is not automatically evidence about current production source.
6. Before starting work, post a short coordination note on Issue #5 (or the task's existing shared thread): what you read, the current finding you believe is established, your proposed bounded next step, the branch/run you intend to use, and any overlap you found. If another agent is already doing that exact test or edit, coordinate rather than duplicate it.

## Keep the message board alive

- Put durable cross-agent findings, blockers, hypotheses, test links, artifact links, exact source SHAs, and next-step recommendations in the shared issue/PR discussions. Do not leave important context only in a private chat, transient run log, or unlinked branch.
- After a run, report **what was tested, what source it used, what passed/failed, what the evidence does and does not establish, and the next gate**. Link the workflow run and artifacts.
- Separate observation from interpretation. A unit/synthetic test, lightweight CI pass, or source-only algebraic check is not proof of full-trajectory or physical validity.
- Preserve negative results and failed hypotheses; they prevent the team from repeating work.
- Keep physics changes separate from diagnostics, test harnesses, and campaign orchestration changes. Do not promote a hypothesis to production physics without the required admission evidence and review.
- Never silently overwrite another agent's branch or assume another agent's task is abandoned. Coordinate first.
- If a finding changes, post a correction with the affected commit/run and explain what is superseded.

## Current communication entry points

- [Issue #5 — Cross-agent coordination / shared message board](https://github.com/miahpenn/0-Cosmos-engine/issues/5)
- [Issue #4 — Momentum-constraint propagation](https://github.com/miahpenn/0-Cosmos-engine/issues/4)
- [Pull requests and their full discussions](https://github.com/miahpenn/0-Cosmos-engine/pulls)
- [GitHub Actions runs and artifacts](https://github.com/miahpenn/0-Cosmos-engine/actions)

**Definition of done for agent onboarding:** the agent has read the latest shared communications, identified the exact source/run relevant to its task, and posted its intended work or findings where the rest of the team can see them.
