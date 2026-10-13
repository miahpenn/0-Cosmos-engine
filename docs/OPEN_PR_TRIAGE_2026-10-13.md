# Open Pull Request Triage — 2026-10-13

**Repository:** `miahpenn/0-Cosmos-engine`  
**Snapshot:** PR #111 merged at `17961bbbbac135132e6ecd0c149d4481f865d38f`. PRs #40, #108, #114, #117, #119, and #123 are closed unmerged as superseded/overlapping. The live GitHub pulls API reports **29 open PRs** (13 Oct 2026 UTC), including this triage PR #115; #120 is the current-main successor to #114, and #122 is a docs-only evidence reconciliation. Refresh before future promotion decisions.  
**Purpose:** prioritize review effort and keep software promotion separate from scientific admission.

> This is a first-pass triage, not a claim that every active PR diff, exact-head workflow state, and review thread has been re-audited in this pass. Recheck each PR's current head, base, CI, reviews, and changed files immediately before any promotion decision. A green workflow proves only the scope it actually tests.

## Promotion policy proposed

1. **Test-only / documentation changes:** promote only after scope and test quality are reviewed, all applicable workflows pass on the exact current head, and overlap/staleness is resolved. Documented technical review is appropriate; an author cannot count their own review as independent.
2. **Production numerical or physics changes:** require a distinct reviewer to inspect derivation, conventions, units, conservation/source reciprocity, affected integration paths, and regressions, followed by exact-head software checks and appropriately scoped full-machine evidence.
3. **Scientific claims/admission:** keep separate from merge readiness. Neither green CI, a short smoke run, nor publication alone establishes long-horizon scientific admission.
4. **No silent promotion:** do not merge a stale/diverged or overlapping candidate until its disposition is explicit. Do not launch an expensive campaign until diagnostics and the question being tested are clear.

## Triage groups

### A. Near-term, bounded review candidates

| PR | Candidate | Initial disposition / next action |
|---|---|---|
| [#110](https://github.com/miahpenn/0-Cosmos-engine/pull/110) | EM zero-field source-projection parity preflight | Exact head `b66452f` has passed all three required workflows (315 tests each); the expanded test checks initialization, both actual RHS-entry stages, and accepted endpoint with radiation on/off. Now needs distinct review; still not EM evolution. |
| [#113](https://github.com/miahpenn/0-Cosmos-engine/pull/113) | Symbolic radial-electric geometry/CMC consumer gains | Exact head `5f77fe4` passed all three required workflows. The test uses consumer-unit increments and fail-closed metric weights; it does not derive EM normalization or evolve Maxwell fields. Needs distinct review. |
| [#107](https://github.com/miahpenn/0-Cosmos-engine/pull/107) | Flat/curved plus frozen-production-metric Gauss references | Current-main-based head `106ecf7` has all three exact-head workflows green (**363 tests each**). Cell charges integrate `4*pi*rho_e*sqrt(gamma)` on the frozen production metric and are checked against an independent high-order reference; the exactly flat metric is cross-checked against the closed-form Coulomb oracle. Tiny-domain source-edge and overflowing `R_areal**2` cases fail closed. Still research-only; review interpolation/extrapolation and quadrature conventions before promotion. |
| [#101](https://github.com/miahpenn/0-Cosmos-engine/pull/101) | Research-only spherical EM scope/gates documentation | Documentation review: verify units/conventions, symmetry limits, consumer inventory, and admission ladder. No EM implementation is authorized by this document alone. |

### B. Physics/source changes — hold for equation review

| PR | Candidate | Initial disposition / next action |
|---|---|---|
| [#116](https://github.com/miahpenn/0-Cosmos-engine/pull/116) | CMC proper-volume input/grid guard | Latest current-main head `3d2c81d` passes all three exact-head workflows (**331 tests each**): [lightweight](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38421064599), [strict JSON](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38421064302), [campaign/retry](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38421064804). Validates the uniform half-cell grid, including off-centre and reordered centres, plus no-floor proper-volume weights and finite projected RHS. Production utility only; keep draft for distinct maintainer review. |
| [#120](https://github.com/miahpenn/0-Cosmos-engine/pull/120) | Scalar–DM exchange on current-main Valencia state | Exact-main head `046ba85` passes all four required workflows (direct RHS extraction, lightweight/one-step, strict JSON, integrated campaign/retry). The source-to-production RHS derivation is conditional on the archived variable-mass law; it does not derive why 0* microphysics selects that law or the exact `beta=-0.04`. Keep draft for distinct maintainer equation review. This is not an EM charge/current implementation. |
| [#40](https://github.com/miahpenn/0-Cosmos-engine/pull/40) | Earlier covariant scalar–DM source correction plus stacked workflow/evidence changes | **Closed unmerged as superseded on 13 Oct.** The current-main comparison found its mixed-component source mapping incompatible with active Eulerian `sqrt(gamma)*E`. Its campaign workflows and several tests were already byte-identical to main; the remaining validator/test versions were stale, and its A/B source test depends on the old convention. PR #120 is the current-main successor source candidate, still draft pending independent equation review. Preserve #40 as historical evidence; do not cherry-pick its old source formula. |

### C. Refresh, duplicate/overlap, or stale-base review first

| PR | Candidate | Initial disposition / next action |
|---|---|---|
| [#103](https://github.com/miahpenn/0-Cosmos-engine/pull/103) | Production-status documentation refresh | Content refreshed in head `16bab1a` against the PR #111 merge and current open-candidate status. Mergeable metadata is true; review the docs-only scope and the note that no workflow was reported for the document commit. |

### D. Numerical diagnostics and scientific-admission evidence — cross-review before prioritizing more runs

These PRs are evidence/diagnostic work rather than automatic production fixes. Review each report against its pinned source, input artifact, metric/stencil conventions, uncertainty, and explicit admission boundary. Prefer a synthesis/cross-review of existing artifacts before launching another trajectory.

| PR | Topic from title | Initial disposition / next action |
|---|---|---|
| [#122](https://github.com/miahpenn/0-Cosmos-engine/pull/122) | Radiation/constraint evidence reconciliation | Current-main docs-only synthesis of existing radiation replay, CMC, constraint, centre-series and Misner–Sharp artifacts; no new trajectory or causal fix. Cross-review source pins and preserve the stated “NOT ADMITTED” boundary. |
| [#64](https://github.com/miahpenn/0-Cosmos-engine/pull/64) | Radiation outer-metric kinematic identity | Cross-check against current radiation continuation and metric-update evidence. |
| [#71](https://github.com/miahpenn/0-Cosmos-engine/pull/71) | Historical constraint onset across resolutions | Check provenance, common-time comparisons, and whether it resolves an actionable cause. |
| [#89](https://github.com/miahpenn/0-Cosmos-engine/pull/89) | Native Misner–Sharp mass-gradient identity | Compare the independent identity and source/geometry conventions with worldtube balance work. |
| [#84](https://github.com/miahpenn/0-Cosmos-engine/pull/84) | Origin-only Ricci split across PIRK stages | Cross-review centre parity and stage-input matching. |
| [#79](https://github.com/miahpenn/0-Cosmos-engine/pull/79) | Origin half-cell quadrature vs interval weights | Reconcile with initialization quadrature and fixed-worldtube balance findings. |
| [#39](https://github.com/miahpenn/0-Cosmos-engine/pull/39) | PIRK L2 centre-stage terms | Compare with later centre-series and source/operator audits; check for superseded diagnostic logic. |
| [#18](https://github.com/miahpenn/0-Cosmos-engine/pull/18) | Initial radial quadrature and CMC centre row | Candidate changes require matched centre, constraint, and fixed-physical-worldtube evidence; short-time improvement alone is insufficient. |
| [#74](https://github.com/miahpenn/0-Cosmos-engine/pull/74) | Initialized centre Ricci vs proper-radial geometry | Reconcile with parity-corrected centre diagnostics and current source pins. |
| [#35](https://github.com/miahpenn/0-Cosmos-engine/pull/35) | Fixed-worldtube Misner–Sharp quadrature balance | High scientific relevance. Confirm fixed **areal** radius, source/work terms, parity convention, and time integration before drawing conclusions. |
| [#70](https://github.com/miahpenn/0-Cosmos-engine/pull/70) | Joint production CMC K target/lapse solve | Diagnostic candidate only unless the governing operator/target convention is independently established. No gauge change by implication. |
| [#73](https://github.com/miahpenn/0-Cosmos-engine/pull/73) | Current-source initialized centre series | Use as a source-pinned evidence input; cross-check against #74/#84/#18. |
| [#78](https://github.com/miahpenn/0-Cosmos-engine/pull/78) | PIRK L2 metric cross-check on current main | Check independence and source parity relative to #77. |
| [#77](https://github.com/miahpenn/0-Cosmos-engine/pull/77) | Independent metric curvature at PIRK L2 stages | Cross-review formula independence and exact stage/profile provenance. |
| [#72](https://github.com/miahpenn/0-Cosmos-engine/pull/72) | Current-source Hamiltonian constraint tangent | Combine with momentum tangent and accepted-step evidence; tangent closure is not long-horizon convergence. |
| [#68](https://github.com/miahpenn/0-Cosmos-engine/pull/68) | Current-source momentum-constraint tangent | Pair with #72 and actual accepted-step/refinement evidence. |
| [#66](https://github.com/miahpenn/0-Cosmos-engine/pull/66) | H/M constraints on archived radiation replay states | Artifact audit only; retain the recorded failure to reach target and avoid claiming cause from terminal snapshots alone. |
| [#65](https://github.com/miahpenn/0-Cosmos-engine/pull/65) | Frozen terminal-state CMC row/operator consistency | Cross-review with #70 and source-matched K-RHS target; do not change gauge/boundary from a counterfactual alone. |
| [#61](https://github.com/miahpenn/0-Cosmos-engine/pull/61) | Retry-safe RHS and independent L3 capture | Valuable evidence/provenance audit. Radiation remains fail-closed before target; do not label this as a trajectory repair or admission. |

### E. Legacy / experimental track — keep isolated from production promotion

| PR | Candidate | Initial disposition / next action |
|---|---|---|
| [#1](https://github.com/miahpenn/0-Cosmos-engine/pull/1) | Strong-field representation and central-clock diagnostic gates | Experimental branch. Keep isolated; require explicit current-base review and the full strong-field admission ladder before considering production changes. |
| [#2](https://github.com/miahpenn/0-Cosmos-engine/pull/2) | N=320 checkpoint reconstruction audit | Historical reproducibility audit with documented host-dependent frozen-tolerance failures. Preserve the preregistered tolerance; restore missing source documents and resolve/report the failures before review completion. |

### F. Queue metadata (not a production promotion)

| PR | Candidate | Initial disposition / next action |
|---|---|---|
| [#115](https://github.com/miahpenn/0-Cosmos-engine/pull/115) | This dated open-PR triage document | Repository-maintenance only. The live GitHub pulls API reports 29 open PRs after #40/#108/#114/#117/#119/#123 were closed unmerged; this count includes #115 itself. This document is a prioritization aid, not a sign-off or a substitute for exact-head review. Refresh the count before future decisions. |

## Recommended work order

1. **Refresh the live inventory** and verify exact current head/base, draft state, checks, and review threads for the 29 metadata-confirmed open PRs. The exact list was refreshed from GitHub's open-pulls API after closing overlapping PRs #119/#123 and opening the newer #120/#122 candidates; refresh before future decisions.
2. **PR #111 is now merged** at `17961bbbbac135132e6ecd0c149d4481f865d38f`; do not continue treating it as open.
3. **Cross-review #110 and #113 next.** Their latest exact-head software gates are green; keep review claims bounded to the source-projection/consumer contracts, not EM physics implementation.
4. **Request a distinct equation review for #120.** PR #114 was closed unmerged as its current-main successor opened; PR #40 is also closed because its active-source mapping is stale. Do not reuse the old A/B witness without re-deriving it on the current Eulerian state.
5. **Cross-review existing diagnostic evidence**, prioritizing #35 fixed-physical-worldtube balance, #61 radiation replay evidence, #65/#70 CMC target/operator consistency, and #66/#68/#72 constraint evidence. Decide the next expensive run only after this synthesis identifies a falsifiable unresolved question.
6. **Cross-review #122's synthesis** against its hash-pinned replay/constraint artifacts; preserve the distinctions between observed failures and unproven causal links, and do not launch another trajectory without a falsifiable unresolved question.
7. **Review the refreshed #103 status PR** against current `main`, preserving the explicit no-scientific-admission statement.
7. Keep scientific admission **false/open** until the independent long-horizon, multi-resolution, centre-regularity, radiation-continuation, source-closure, and fixed-areal-radius balance gates actually pass.

## Scope and honesty notes

- The live open-pulls API returned 29 open PRs on 13 Oct 2026 after PRs #40, #108, #114, #117, #119 and #123 were closed unmerged. New candidates #120 and #122 are included. Counts are point-in-time snapshots and should be refreshed.
- Grouping is an initial prioritization based on PR titles, descriptions, and already-recorded coordination evidence. It is not a substitute for reviewing every diff.
- No PR is merged, closed, approved, or scientifically admitted by this triage document.
