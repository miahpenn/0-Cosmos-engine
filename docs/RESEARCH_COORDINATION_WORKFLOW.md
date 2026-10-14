# Research coordination workflow

This is the operating procedure for parallel agents working in `0-Cosmos-engine`. It strengthens the existing Issue #5 claims, board #135, and `docs/STATUS.md`; it does not replace them or authorize a physics change.

## One source of truth per kind of state

| Record | Owns | Do not use it for |
|---|---|---|
| [`docs/STATUS.md`](STATUS.md) | Durable, verified current implementation/status facts, exact source anchors, evidence limits, scientific-admission disposition | Live ownership, every transient run, or an undated "current" PR count |
| [Team board #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135) | Prioritized queue, one owner/next action per task, blockers, completed/superseded disposition | Repeating full diagnostic reports |
| [Cross-agent Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5) | Claims, overlap checks, cross-lane decisions, handoffs, corrections and links to evidence | A second copy of the whole task board |
| Each task PR | Scope, exact changed files, code/docs, exact-head checks, review, artifacts, interpretation boundary | Claims that aren't reflected on #5/#135 when cross-agent coordination is needed |
| GEAR/archive materials | Historical blueprint and unresolved premises | Permission to invent missing coupling, laws, normalization or matter assignments |

PR counts and check states change quickly. Query GitHub live before making merge/priority decisions; date-stamp any count copied into a snapshot.

## 1. Start with a bounded read, not a full-thread reread

Before editing or launching a run:

1. Check the live `main` SHA and `docs/STATUS.md`; identify the production-source anchor separately if it differs from the current documentation tip.
2. Read the relevant row on board #135: owner, previous evidence, open gate and next action.
3. Read the Issue #5 body and newest relevant comments. Search by task/lane, PR number, equation/observable and artifact/run ID; follow links to the specific predecessor discussion.
4. Inspect relevant open PRs, including their base/head SHAs, changed files, discussion, review threads, checks and artifacts. Similar titles are insufficient; identical branch/head is a collision signal, not proof of duplicate intent.
5. Read the full historical thread only when reconstructing a disputed result or premise. With a very large comment history, repeatedly rereading from the beginning adds cost without improving the current claim check.

Fresh agents must still read `docs/SI_AGENT_FIRST_READ_COMMS_PROTOCOL.md`. The efficient rule is **latest relevant context + linked predecessor evidence**, not ignoring communications and not rereading every unrelated comment.

## 2. Claim before work

Post a claim on Issue #5 (or the task's established shared thread) before coding or expensive work. Make one claim about one bounded discriminator. Include:

- **Lane / discriminator:** the precise question and measured observable.
- **Owner:** one accountable agent; identify collaborators/reviewer separately.
- **Scope:** proposed minimal action and files/layers allowed to change; explicit no-touch boundaries.
- **Source:** current `main` SHA, production source anchor and relevant prior artifact/run.
- **Branch / PR:** planned branch; existing PR if continuing one.
- **Overlap check:** relevant board row, issue comments, open PRs and matching changed-file/branch/head findings.
- **Discriminator:** hypothesis, control/counterfactual and result that would falsify or distinguish it.
- **Output / stop condition:** expected artifact or patch, budget/horizon and the exact condition to stop or hand off.
- **Claim state:** `ACTIVE`, `BLOCKED`, `COMPLETE`, `RELEASED` or `SUPERSEDED`.

A claim is ownership of a *specific discriminator*, not a whole physics lane. Independent tasks may run concurrently when their observables, inputs and edit surfaces do not conflict. If multiple tasks must touch the same production seam, designate one integration owner and agree the order first.

## 3. Resolve overlap before spending compute

Check both scientific and code overlap. Two tasks overlap when they answer the same bounded question with materially the same state/source, operator/stage pairing, control and observable—even if their branch names differ.

When a likely collision appears:

1. Pause new expensive runs and shared-file edits.
2. Compare exact purposes, base/head SHAs, changed-file sets, existing artifacts, and already completed checks.
3. Choose one of three outcomes in the claim thread: **join** the existing task, **consolidate** into one surviving PR, or **orthogonalize** the question so the outputs differ in a stated, useful way.
4. Link the decision from #135. If a PR is duplicate/superseded, close it with a clear link to the survivor; preserve its unique evidence first. Never force-push over another agent's work or silently assume an idle claim is abandoned.

Identical head SHAs (for example, two PRs sharing a branch) need immediate purpose/diff review. They can be separate administrative proposals, so do not close one solely on the SHA.

## 4. Minimize run and review cost

- Diagnostics before physics edits; demonstrate the failing discrete identity before changing production equations.
- Run the cheapest discriminating check first. Preflight source pins, inputs, JSON/evidence contracts, expected outputs and stop conditions before a costly campaign.
- One expensive run at a time for a shared discriminator. Do not repeat an unchanged run/artifact unless a reviewer identifies a specific gap.
- Never launch a long campaign merely because a short run passed. State what new question the longer horizon resolves and the gate that decides continuation.
- Keep diagnostics, production physics, test harnesses and campaign orchestration separated where possible. No fitted coefficient, arbitrary damping, floors/clamps, inferred matter law, silent physics, or altered GEAR premise to make a residual look better.
- Prefer one PR per logical change. Base new work on live `main` unless stacking is necessary; a stacked PR must name its parent PR, dependency and merge order. Refresh the base and rerun applicable checks before promotion.

## 5. PR lifecycle and exact-head evidence

A PR moves through these states: **draft/research → review-ready → mergeable → merged**. It can instead be **blocked**, **released**, or **superseded**. A green workflow is not a review and is not scientific admission.

Before declaring review-ready, record:

- Exact PR source-head SHA and, if applicable, the separate synthetic merge/execution SHA.
- Checks relevant to the changed paths, with run links and result at that exact head.
- Artifact links and hashes for diagnostic outputs; source/vendor pins where used.
- Test/campaign controls, actual horizon/resolution and key signed quantities—not only pass/fail or norms.
- What the result establishes, what it does **not** establish, remaining uncertainty, and the next gate.
- Independent review required by the task, plus any stale-base/dependency issue.

Any substantive head change invalidates a claim that the prior head's checks qualify the new head. Re-run applicable checks or explicitly identify the unverified gap. Keep research PRs draft and research-only until the relevant independent scientific and admission gates—not just software CI—are satisfied.

## 6. Keep records short, linked and current

- **Issue #5:** one concise claim before work; one durable result/handoff after the bounded action; add interim notes only for a blocker, material scope/source change or time-sensitive coordination decision. Link artifacts instead of pasting their full contents.
- **Board #135:** update owner, state, exact next action, PR/run links and disposition when work starts, blocks, completes or is superseded. Do not leave completed work looking active.
- **`docs/STATUS.md`:** update only when verified implementation/status, source baseline, evidence boundary or admission disposition changes. Put new current facts at the top, date-stamp snapshots, and clearly mark older sections historical/superseded. A documentation-only coordination change with no new scientific result does not justify rewriting the scientific status.
- **Task PR:** owns detailed patch and validation record. Use links from #5/#135 rather than cloning full reports into every surface.
- **Correction:** if an interpretation or source SHA changes, identify the superseded claim, state why, and link the corrected evidence. Preserve negative results.

During active work, refresh the claim after a meaningful commit or run, when the source/base changes, and before handing off. If no visible update has occurred for roughly 30 minutes, post the current state/blocker or explicitly release it; silence alone never transfers ownership.

## 7. Completion checklist

Before marking a task complete, confirm:

- [ ] The claim has an explicit final state and overlap outcome.
- [ ] The relevant PR/branch is linked, with exact current SHAs and checks.
- [ ] Artifacts are linked and hash/provenance recorded where applicable.
- [ ] Interpretation limits and the next scientific/software gate are explicit.
- [ ] Board #135 reflects the true owner/status/next action.
- [ ] `docs/STATUS.md` is updated if—and only if—a durable verified status fact changed.
- [ ] Duplicate/stacked/superseded work has a clear survivor/dependency and no unique evidence was lost.

**Default rule:** one claim, one owner, one bounded discriminator, one canonical evidence artifact, one updated queue row. Parallelize independent work; coordinate shared seams.
