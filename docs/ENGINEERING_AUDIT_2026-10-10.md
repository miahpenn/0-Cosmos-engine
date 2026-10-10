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


## 7. Addendum — verified follow-on gates (2026-10-10)

**Current experimental head at this update:** `1ad4b6219bfd6899d5f064cd90f19f6a2dee32c3`  
**Branch:** `physics/spatial-beta-covariant-source-closure`  
**Observed `main` head:** `48f06d09a5b88b831ac86f6bdd3d9ecc94f092f6` (unchanged; no merge)

The audit snapshot above accurately records the starting point. The following later checks extend it; they are listed separately to preserve chronology and the failed attempt record.

| Gate | Exact commit | Run | Result |
|---|---|---|---|
| Audit document CI | `72d6ef48082fc2c96af014f8b9aa70c3b4df25d5` | [#567](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38066744796) | Pass |
| CMC residual at three deterministic initial configurations | `76999d2c9a4d62c88262444e357825ef9d97ffd8` | [#568](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38067318842) | 116 passed |
| CMC residual on actual predictor and accepted slices | `7af42c00fb4b68a4888bebbc7fd3f4deb176f2a9` | [#569](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38067488971) | 117 passed |
| Reject invalid metric in proper-volume projection | `373427dfba0307efb8c1051d16eaba556a6664c8` | [#570](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38067741630) | 121 passed |
| First shared-helper refactor attempt | `8dc13270053d6aba943653571f49c32b6ffd0d5f` | [#571](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38067950792) | **Failed** — stale `expected_shape` reference raised `NameError` in `target_kdot` |
| CMC helper reference correction and full solver preflight check | `1ad4b6219bfd6899d5f064cd90f19f6a2dee32c3` | [#572](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068034703) | **124 passed in 60.02 s** |

Run #571 is intentionally preserved as a failed development gate; #572 is the verified correction, not a replacement of the historical failure. In all these CI runs, the two optional long-campaign jobs (`center_repair_ab_clock` and `persist_d0_on_artifact`) were skipped by their explicit workflow conditions. No long physical trajectory was launched by these checks.

### CMC metric-admissibility finding

The earlier proper-volume expression used `sqrt(max(a, 0))`. That could conceal a negative radial metric coefficient while forming CMC gauge weights. The current branch removes that clamp and uses a shared fail-fast validation path for both `target_kdot` and `solve_cmc_lapse`. It rejects mismatched shapes, non-finite metric data, non-positive `r/a/b/X`, and invalid grid spacing before metric conversion, matter projection, or K-RHS evaluation. Regressions confirm negative `a/b/X` values are rejected before the respective downstream operators are called and non-finite metric data are rejected before K-RHS evaluation.

This is an admissibility guard, not a geometry repair, projection, new physical term, or change in the CMC target. The same discrete residual threshold (`1e-10`) remains unchanged. The result is code-level consistency evidence only; it is not evidence of long-run stability or physical validity.

### Repository status and next gate

- `main` remains at `48f06d09a5b88b831ac86f6bdd3d9ecc94f092f6`; all work above remains isolated.
- The repository still reports `main` and the experimental branch as unprotected. The connected GitHub integration could not inspect or change branch protection (HTTP 403); a repository administrator must review required-check enforcement.
- The source-exchange, homogeneous-limit, stage-synchronization, multi-configuration CMC, actual predictor/final CMC, and invalid-metric rejection gates are green on their listed commits.
- Still **not established**: long-horizon constraint/convergence behavior on this exact patched head, radiation outer-boundary failure reproduction on this head, reciprocal local-COSMOS coupling, a D-to-matter law, global shear closure, turnaround, bounce, or a completed cycle.
- Next: continue the remaining deterministic diagnostic gates, then reproduce the known late radiation realizability/boundary failure on an explicitly pinned current commit. Keep that campaign separate from this source/gauge code audit, preserve all artifacts and hashes, and do not introduce clipping, damping, or speculative boundary physics.

## 8. Addendum — patched-branch radiation predictor trace (2026-10-10)

The late outer radiation failure has now been reproduced on the patched covariant-source/CMC source branch, using an explicitly pinned diagnostic-only workflow. See the full gate report: [Radiation predictor failure trace](RADIATION_PREDICTOR_TRACE_2026-10-10.md).

- [Trace run #7](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938182) completed successfully as a diagnostic capture and uploaded artifact [11676486591](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938182/artifacts/11676486591). [Ordinary CI #576](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938171) passed 127 tests in 59.58 s. The physical trajectory itself stopped at the recorded admissibility failure rather than reaching the target time.
- Run SHA / trigger: `c8a9967126ba07cd12e4da9c86b7066a703d3419`; source baseline before instrumentation commits: `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`; instrumentation: `527b44931b059f4d118a86c898b91d79ca4d96f2`; workflow: `ed0bc6f61103f2439556b27ad588116a149f4add`; vendor submodule: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Configuration: N=80, r_max=80, dr=1, CFL=dt=0.03, D=1e-10, radiation ON, target t=50. The last accepted state is t=46.98, center proper time τ=22.016115393070738, 1,566 steps. The predictor attempting t=47.01 fails at `explicit_predictor_first_CMC_lapse` at outer cell i=79 / r=79.5 with `|S|/E=1.00005699412463422`.
- Ordered same-cell budget: `C` at cell 79 is +4.6971e-7 on the accepted metric, +1.4812e-6 after flux update on that metric, and +2.2171e-6 after sources on that metric. Evaluating the predictor conservative state on the predictor metric contributes −2.2340e-6 and yields `C=−1.6900e-8`. Thus the metric evaluation is the decisive crossing in this budget; flux/source contributions at this cell on the failing step improve the margin under the stated ordering.
- The accepted outer-cell ratio had risen from ~0.99254 at t=46.86 to 0.99839 at t=46.98. Between accepted and predictor metrics, `a` falls 0.9445%, `X` rises 0.2821%, and `gamma_rr_inv=X²/a` rises 1.5239% while boundary `alpha=1`, `beta=0` are unchanged.
- Trace counters: 394 stage snapshots, 66 late-window predictor budgets, 1,566 progress samples, zero instrumentation errors, zero source-capture misses, zero budget closure errors, and zero predictor conservative-variable reconstruction errors. ZIP digest is `bef03745ad4ea1ab62d8a42de9b957c57210571d933cf2728c96a6196b970c43`; trace JSON hash is `642a82b8b97ddc849b08f70b4379672b08e81f842d3ae8bdeaaa1d5374afd8a3`.
- Bounded conclusion: the immediate violation is caused by evaluating the explicit predictor conservative radiation state on its predictor metric in this discretization. This does not prove the metric evolution law is wrong or establish the upstream cause. No physical source, damping, clamp, floor, or boundary law was changed.

### Next gate

Add a diagnostic-only outer-worldtube ledger for the explicit `a`, `b`, and `X` RHS contributions at cells 75–79; reconstruct `Δgamma_rr_inv`; then audit the matching characteristic-boundary variables before any controlled resolution/CFL comparison. No stabilization edit is authorized by this result.