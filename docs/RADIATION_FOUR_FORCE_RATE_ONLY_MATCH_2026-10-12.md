# Radiation Four-Force Matched Rate-Only Comparison — 2026-10-12

## Status and disposition

**Bounded diagnostic gate passed. No physics promotion.**

The matched four-case experiment now uses direct semi-discrete matter RHS extraction rather than the former `evolve_species(dt=1.0)` surrogate. The production kernel also checks radiation admissibility at its real initial, explicit-predictor, primary-predictor, and final states. The optional radiation coefficient remains diagnostic-only and defaults to `beta_rad=0.0`.

- Workflow: [38317343830](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317343830) — success.
- Tested code/workflow head: `278172b3f9e7d2f054909e69647ccb652b3e4593`.
- Branch: `diag/radiation-four-force-rate-only-match-2026-10-12`.
- Vendor gitlink verified: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Focused RHS/admissibility/four-force identity tests: **9 passed**.
- Full repository suite: **31 passed**.
- Matched one-step cases: **4/4 completed**, no partial captures or stage admissibility rejections.
- Artifact: [11776942119 — radiation-four-force-rate-only-match](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317343830/artifacts/11776942119).
- ZIP SHA-256: `dfa7e72bfd9232fa7b3d055f668098b2a7069b3d2d35192cee25131d71f3027f`.
- Internal JSON SHA-256: `bc4a1ab7f6146e157dd5514256a208b24194cce5540d0a63e1175a53e156e8cc`.

The first attempt (run [38317282828](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317282828)) passed both test tiers but did not launch the script because direct script execution omitted the repository root from Python's import path. The workflow invocation was corrected to set `PYTHONPATH=.`; the successful run above used that correction. This was a harness-only failure, not a numerical case result.

## What was combined — and what was not

The diagnostic branch starts from the existing four-force instrumentation commit and combines three already-isolated software responsibilities:

1. The rate-only `_species_rhs_components` structure from draft [PR #8](https://github.com/miahpenn/0-Cosmos-engine/pull/8), so a derivative query no longer performs a fictitious unit-time update.
2. The existing optional `beta_rad` covector-source path from draft [PR #9](https://github.com/miahpenn/0-Cosmos-engine/pull/9), with default zero unchanged.
3. Real-state radiation admissibility checks inspired by draft [PR #10](https://github.com/miahpenn/0-Cosmos-engine/pull/10), evaluated against each stage's metric at initial, explicit predictor, primary predictor, and final state; direct finite-time `evolve_species` calls remain guarded as well.

The diagnostic exception handler no longer manufactures a `dt`-scaled estimate by subtracting the artificial `dt=1` rejected state. Any future exception is instead recorded as a rejection from the real production call path. No floors, clipping, projections, retries, damping, boundary edits, or production-physics changes were introduced.

No changes were made to `main`, the canonical GEAR archive, or the source branches of PRs #8–#10. This branch is not a request to merge any of those PRs.

## Fixed experiment

Exactly one ordinary production step per case:

| Resolution | `beta_rad` | `dt` | Result |
|---:|---:|---:|---|
| 80 | 0.0 | 0.0150 | completed |
| 80 | -0.04 | 0.0150 | completed |
| 160 | 0.0 | 0.0075 | completed |
| 160 | -0.04 | 0.0075 | completed |

All four cases passed the initial, explicit-predictor, primary-predictor, and final radiation-cone checks. The guard remains the existing (E\ge0), (|S|\le E) condition; no state was projected into the cone.

## Matched endpoint comparison

Values below are from each run's saved initial and completed-step diagnostic snapshots. Differences are candidate (`beta_rad=-0.04`) minus null (`beta_rad=0`) at the same resolution.

| Observable | N=80, null | N=80, candidate | N=160, null | N=160, candidate |
|---|---:|---:|---:|---:|
| Completed (t) | 0.015 | 0.015 | 0.0075 | 0.0075 |
| Completed proper time (	au) | 0.015193016568303245 | 0.015193016568303474 | 0.007597378861711842 | 0.007597378861710647 |
| (H_{\rm eff}) | 0.005668672088434478 | 0.005668672088434433 | 0.005668820081124816 | 0.005668820081124803 |
| Max admitted (|H|) | 6.717284995089558e-5 | 6.717284991405613e-5 | 8.402310785286197e-5 | 8.40231078333936e-5 |
| Max admitted (|M|) | 4.206728234575102e-6 | 4.206728234524002e-6 | 4.1778562716321506e-6 | 4.177856271377314e-6 |
| Max (|C_\Lambda|) | 6.392653018548987e-8 | 6.392653018657407e-8 | 3.175797589727129e-8 | 3.1757975898572335e-8 |
| Centre radiation margin | 0.990091377630515 | 0.9900913876881988 | 0.9900931588317909 | 0.9900931638613395 |
| Admitted-region radiation margin | 0.9998741295870003 | 0.9998741297155405 | 0.9999374621794956 | 0.9999374622114284 |

The endpoint differences in (H_{\rm eff}), (|H|), (|M|), and (|C_\Lambda|) are tiny at both resolutions. The centre radiation-margin increases are approximately (1.006\times10^{-8}) at N=80 and (5.030\times10^{-9}) at N=160; admitted-region margin increases are approximately (1.285\times10^{-10}) and (3.29\times10^{-11}), respectively. Proper-time endpoint differences are around (10^{-16})–(10^{-15}).

This is a useful **small one-step response** to the candidate source on these frozen initial profiles. It neither establishes that the four-force is physically correct nor demonstrates a meaningful trajectory-level benefit. The fact that the counterfactual is reachable and changes the RHS is only a software-path test, not a derivation of a radiation coupling law.

## Explicit limitation: Misner–Sharp witness unavailable here

The snapshot fields `MS_current_max` and `MS_current_rms` are serialized as `nan` for all four cases. Do **not** interpret this run as a Misner–Sharp balance pass. A preregistered tolerance and a valid same-stage multi-step balance witness remain required.

## Reconciled defect status

The former N=160 rejection with (|S|/E\approx1.36) occurred during derivative extraction because the old wrapper tested (U+L(U)) at an artificial unit step. The rate-only implementation removes that false test, and all N=160 one-step cases here pass the real stage guards. This settles only the RHS API-contract defect.

It does not resolve the separate corrected-rate long-run event: after 1,566 accepted steps near (t\approx46.98), the explicit geometry predictor produces an outer-cell ((i=79)) radiation-cone crossing. See the [central trace and offline attribution](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/RADIATION_PREDICTOR_TRACE_AFTER_RATE_FIX_2026-10-12.md). That failure still needs a same-step geometry/source attribution; this short N=80/N=160 initial-data test cannot explain or negate it.

## Next gate

1. Have an independent reviewer check the combined diff for algebraic equivalence of the rate helper and exact preservation of the default `beta_rad=0` path.
2. Keep this comparison diagnostic-only. Do not merge PR #9 or promote `beta_rad=-0.04` from one-step data.
3. Resume from the already-saved accepted-step/failed-predictor artifact for the distinct outer-cell crossing; preregister the geometry/source contribution checks and keep actual-stage cone admission fail-closed.
4. Keep the centre-series, joint Hamiltonian/momentum/connection constraint propagation, Misner–Sharp balance with a frozen tolerance, resolution convergence, D-to-matter, and local/global-Q gates separate and open until the evidence closes them together.

No long trajectory was launched for this result.
