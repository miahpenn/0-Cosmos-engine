## Summary
<!-- What changed, why it is useful, and the intended outcome. -->

## Scope
- [ ] Documentation / coordination
- [ ] Test / evidence contract
- [ ] Diagnostic / report only
- [ ] Campaign orchestration
- [ ] Production implementation (explain why this lane is authorized)

**Files / behavior changed:**
<!-- List the bounded change and explicitly state what did not change. -->

## Task claim and overlap
- **Canonical lane / issue:** <!-- Issue #5 claim + Issue #135 task or task-specific thread -->
- **Bounded action and stop condition:**
- **Base `main` SHA:**
- **PR source-head SHA:**
- **Production-source / vendor anchor(s), if different:**
- **Overlap check:** <!-- Closest open PRs, claims, and artifacts; explain why this is not duplicate work. -->
- **Dependencies / stacked PRs:**

## Exact-head evidence
| Check / workflow | Run URL | Result at PR head |
|---|---|---|
| Lightweight regression / smoke | | |
| Strict evidence/JSON contract | | |
| Integrated retry/admission | | |
| Additional diagnostic (if applicable) | | |

**Artifact URL(s) and hashes:**  
**Negative results or failed attempts preserved:** <!-- link them; write N/A only if none -->

## Interpretation boundary
**Established by this PR:**

**Not established / open admission gates:**

## Review and disposition
- [ ] Final commit SHA is recorded and applicable checks were run on that exact head.
- [ ] Skipped or missing workflows are identified honestly, not described as passes.
- [ ] Source pins and vendor/dependency boundary were checked where applicable.
- [ ] No unapproved physics, normalization, source law, numeric fitting, floor/clamp, damping, boundary, or campaign-policy change is hidden in the diff.
- [ ] Independent review / derivation gate is complete, or this PR remains draft.
- [ ] Issue #5 and Issue #135 are updated with the result and next gate.
- [ ] Scientific admission is described separately from software/diagnostic success.

**Requested review / next gate:**
