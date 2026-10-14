# SI-agent first-read protocol: shared communications are operational context

**Required for every fresh SI agent before analysis, coding, or launching a run.** Skipping coordination causes duplicate experiments, conflicting branches, lost findings and avoidable compute.

Read the operating procedure first: [Research coordination workflow](RESEARCH_COORDINATION_WORKFLOW.md).

## First-read checklist

1. Read the current `README.md`, this document, `docs/STATUS.md`, and the relevant production/architecture/physics lock and handoff documents.
2. Read the body of [Issue #5 — cross-agent coordination](https://github.com/miahpenn/0-Cosmos-engine/issues/5) and the canonical [team task board #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135).
3. For the lane you intend to touch, inspect the newest relevant #5 comments and search by lane/task name, equation or observable, PR number, source SHA and artifact/run ID. Follow links to predecessor evidence and read relevant PR discussions/review threads. **Do not reread the entire 1,000+ comment #5 history for every task**; read older history in full only when reconstructing a disputed or historical question.
4. Check the live `main` SHA. Separate the current repository tip from the production source anchor and the exact source SHA used by every cited run. Inspect relevant open PR base/head SHAs, branch names, changed files, checks and dependencies.
5. Before editing or launching a run, post a bounded claim on Issue #5 (or the task's established shared thread): lane/discriminator, owner, exact source anchor, proposed action and stop condition, branch/PR, existing artifact, and overlap outcome. If another agent is doing the same discriminator, coordinate rather than duplicate it.
6. After the bounded action, post the result or handoff with exact SHA, test/workflow/artifact links, what passed/failed, interpretation limits and next gate. Update board #135 and `docs/STATUS.md` only according to the ownership rules in the workflow.

## Non-negotiable conduct

- Treat communications as part of the research record; do not leave durable findings only in a private chat, transient run log or unlinked branch.
- Preserve negative results and failed hypotheses so the team does not repeat them.
- Separate observation from interpretation. Synthetic tests, source-only algebra, lightweight CI and short-horizon runs do not prove full-trajectory or physical validity.
- Keep physics changes separate from diagnostics, tests and campaign orchestration where possible. Do not promote a hypothesis to production physics without the required evidence and review.
- Never silently overwrite another agent's branch or assume a claim is abandoned. Resolve overlap explicitly; update a claim when its source, scope or owner changes.
- A passing software suite does not imply scientific admission. Record that boundary clearly.

## Entry points

- [Issue #5 — claims and cross-lane messages](https://github.com/miahpenn/0-Cosmos-engine/issues/5)
- [Issue #135 — canonical task queue and ownership](https://github.com/miahpenn/0-Cosmos-engine/issues/135)
- [Issue #4 — momentum-constraint propagation](https://github.com/miahpenn/0-Cosmos-engine/issues/4)
- [Open PRs and their full discussions](https://github.com/miahpenn/0-Cosmos-engine/pulls)
- [GitHub Actions runs and artifacts](https://github.com/miahpenn/0-Cosmos-engine/actions)

**Definition of done for onboarding:** the agent has located current status, inspected relevant existing claims/artifacts, identified the exact source/run relevant to its task, and posted its intended work where the team can see it.
