# 0-Cosmos Engine — Engineering Audit and Handoff
**Audit date:** 2026-10-10  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Scope:** Current covariant source-closure / production CMC diagnostic branch. This is a provenance and engineering-status record, not a claim of physical validation.

## 1. Snapshot at audit start

- **Working branch:** `physics/spatial-beta-covariant-source-closure`
- **Audited head:** `96f2b3a6aecd4d0bf56d446ca1b9cc387d50de58`
- **Head commit:** “Add production CMC lapse residual gate”
- **Main head observed:** `48f06d09a5b88b831ac86f6bdd3d9ecc94f092f6`
- No changes in this investigation have been merged to `main`.
- The branch is experimental and must remain isolated until its gates and full-machine checks are complete.

## 2. CI evidence tied to exact commits

| Gate | Commit | GitHub Actions run | Result |
|---|---|---|---|
| Covariant scalar/DM source-pair correction | `69f8e8b1e39a879e0976f0d36a1b435756067450` | [Run 560](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38062843643) | Pass |
| Homogeneous COSMOS-limit regression | `67732a32a4189ede2c3b841f3110fd5e247389a7` | [Run 562](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38063212493) | Pass |
| Stage-synchronization gate (after correction) | `89469cb7856db173348c14f73d887a97a7cc8a14` | [Run 565](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38063717162) | Pass |
| Production CMC-lapse residual gate | `96f2b3a6aecd4d0bf56d446ca1b9cc387d50de58` | [Run 566](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38063885173) | Pass; **114 tests passed in 31.71 s** |

Run 566's `tests` job passed. The workflow's `center_repair_ab_clock` and `persist_d0_on_artifact` jobs were skipped by workflow conditions; they are not passes and do not indicate that any physical trajectory was run.

Two earlier stage-gate CI attempts failed and remain part of the audit trail:
- [Run 563](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38063399106) — initial stage synchronization gate incorrectly assumed CMC solve idempotence.
- [Run 564](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38063591005) — follow-up stage-gate assertion still mishandled the auxiliary (B)-field reset.
- The gate was corrected to compare each stage with its own recorded solve input/output and account for the expected (B) reset. The corrected gate passed in Run 565. Keep failed attempts visible; do not overwrite or relabel them as successes.

## 3. What is currently covered

### Covariant exchange-source audit
The active scalar convention is
`∂t φ = αΠ + βʳ∂rφ`.
For the archived exchange law `∇μ T_DM^{μν}=+βρ_DM ∇^νφ`, the paired scalar source must include the lapse factor `+α β_DM ρ_DM` in the scalar-`Π` equation. The matter source terms are derived from the same covariant law:
- rest-mass/current source: `q_rest = −α√γ βρ W(Π + vʳ φ_r)`;
- energy source: `q_e = −√γ(q_t − βʳ q_r)`, with `q_t=βρ∂tφ`, `q_r=βρ∂rφ`, reducing to `−√γ αβρΠ`;
- momentum source: `q_s = α√γ βρ φ_r`.

The regression uses non-unit lapse, nonzero shift, scalar gradient, and nonzero radial DM velocity to exercise the live source assembly and the scalar/DM energy-pair cancellation. These checks establish consistency of the tested source terms; they do not establish long-run physical correctness by themselves.

### Homogeneous-limit regression
A homogeneous FLRW slice is compared against `cosmos_rhs` with matched `a,H,φ,Π,ρ_DM,ρ_b,ρ_r`. The exact spherical branch does not import archived Bianchi-I shear as an isotropic fluid. The test checks selected scalar/matter derivatives, the BSSN trace/expansion relation, and Hamiltonian/momentum constraints. This is a deterministic limit regression, not a long-time convergence proof.

### Stage-synchronization regression
The gate checks two RHS evaluations per step; the predictor as an explicit Euler update from RHS0; each stage's lapse against that stage's own CMC solve; accepted scalar/matter values against the Heun combination; and input-state non-mutation. The auxiliary (B) field is intentionally reset by the solve and is handled accordingly in the corrected test.

### Production CMC-lapse residual regression
The latest test calls the production kernel's `_solve_lapse`, inserts the solved lapse into a copied geometry, recomputes the production discrete (K)-RHS using the vendor L2 and matter-coupled L3 operators, and checks the interior residual against the projected target, positivity/finiteness of lapse, and `α[-1]=1`. It passed in Run 566 at the exact audited head.

## 4. Repository and governance audit findings

- GitHub currently reports `main` as **not protected** (`protected: false`). The feature branch also reports `protected: false`.
- The branch-protection API request returned HTTP 403 (“Resource not accessible by integration”), so this connection could not inspect detailed protection settings or change them. Treat the observed `protected: false` metadata as the finding; do not claim that required checks are enforced.
- Two open draft PRs were observed:
  - [PR #1 — Strong-field representation and central-clock diagnostic gates](https://github.com/miahpenn/0-Cosmos-engine/pull/1), head `0star-central-clock` → `main`.
  - [PR #2 — N320 checkpoint reconstruction audit](https://github.com/miahpenn/0-Cosmos-engine/pull/2), head `n320-reconstruction-audit-review-2026-10-09` → `0star-central-clock`.
- The current covariant-source branch is not represented by either of those PRs in the observed open-PR list.
- The README and `docs/STATUS.md` still contain older project-state language (README points to `0star-central-clock`; STATUS is dated 2026-10-06). This audit is deliberately additive and does not rewrite that historical material. A later documentation pass should update the current-status pointer without erasing the dated record.
- The current branch's CI workflow is the repository's regression gate, but green unit/regression tests are not a substitute for run-artifact admission, hash matching, convergence, or a full physical trajectory.

## 5. Next steps — in order

1. Confirm CI for this audit-document commit, and record the exact resulting commit SHA and run URL here or in a dated addendum.
2. Keep the covariant-source branch isolated; do not merge to `main`.
3. Continue diagnostic gates only after reading the exact CI result. Fix failures at the layer that produced them; no arbitrary tolerances, hand tuning, damping, floors/clamps, or physical-law changes to make tests pass.
4. Audit the source-pair equations and test coverage against the archived covariant convention; ensure the tests fail if any required lapse/shift/velocity factor is removed.
5. Review the production CMC solve and its residual gate at more than one deterministic state/resolution before authorizing an expensive trajectory. Any new tolerance must be derived from operator precision/discretization, not selected to conceal residual.
6. Before any trajectory is admitted: pin branch + commit SHA + run URL + artifact SHA-256; verify the completed run's commit equals the expected commit; preserve failed attempts; record configuration, grid/domain/CFL, completion status, checkpoint count, and relevant residual/constraint ledgers.
7. Separately reproduce the late radiation outer-cell realizability failure on the current patched head when authorized. Keep that campaign separate from the source-closure audit; no clipping, damping, speculative radiation changes, or boundary-law edits.
8. Have a repository administrator review enabling branch protection / required CI checks for `main`. This connection could not make that change.

## 6. Scientific interpretation lock

No bounce, turnaround, completed cycle, or validated cosmology is established by these code-level gates. The preferred current interpretation remains open to the machine's results; do not impose a bounce/reset/branch flip or call a single crossing a completed cycle. The governing project rule remains:

**Run the full machine → expose the actual defect → repair only the affected layer → rerun.**

Keep diagnostic engineering, physical-law changes, and numerical campaigns separately identifiable in commits and artifacts.
