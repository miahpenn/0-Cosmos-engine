# Issue #4 split-derivative work — recovery and provenance audit

**Recorded:** 2026-11-06  
**Status:** documentation-only recovery record; original Nov 5–6 commits and bundle were created locally, but the original patch bytes and bundle are not currently available in the active recovery environment  
**Scientific status:** diagnostic-only; `scientific_admission=false`; no production or solver change

## Purpose

Preserve the Nov 5–6 affine-root / implicit-derivative diagnostic work without rewriting history or presenting reconstructed material as original source. This document distinguishes (1) artifacts whose contents have been recovered, (2) the integrated contents described by the contemporaneous Slack handoff, and (3) original Git objects and files that remain unavailable.

This is a recovery/provenance audit, not a new scientific review, not a test implementation, and not evidence of scientific admission.

## Coordination and source trail

- Parent workstream: [Issue #4 — discrete momentum-constraint propagation at the central stencil](https://github.com/miahpenn/0-Cosmos-engine/issues/4).
- Coordination record: [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5), including the recovery claim posted 2026-11-06.
- Nov 6 Slack handoff, part 1: https://0-cosmos-engine.slack.com/archives/C0CA75Z6B1N/p1793948780413219
- Nov 6 Slack handoff, part 2: https://0-cosmos-engine.slack.com/archives/C0CA75Z6B1N/p1793948786306169
- Reported local branch: `issue4-integrated-test-diagnostics`, based on `f8d5a65`.
- Reported two-commit sequence:
  1. `2937ef96` — test-only integrated diagnostic patch.
  2. `8965a198` — CI workflow separation and branch trigger.
- The handoff describes both commits as local-only, unpushed, unmerged, and diagnostic-only. Searches of the connected GitHub repository did not find that branch or either commit. The abbreviated SHAs above are recorded exactly as supplied; no full SHA is inferred.
- A subsequent local bundle operation was reported in the Nov 6 conversation: `/mnt/user-data/outputs/issue4-integrated-test-diagnostics.bundle`, 16,133 bytes, SHA-256 `8544465262bfb4806f14dd18b477387fa6391c930928ef165c0eeab3683d6e1e`. `git bundle verify` passed; its only prerequisite was base `f8d5a65`, and the bundle contained commits `2937ef96` and `8965a198` with `8965a198` as branch head. This is a recorded prior-session verification, not a fresh verification in the current environment.
- The subsequent Slack file-upload attempt returned file ID `F0CFA70KUE9` but failed; Slack reported the file as deleted/no data. No successful Slack file link or upload confirmation was recorded. The handoff messages themselves are available at the two links above, but they are not the bundle attachment.

## Recovered material

### 1. Affine-root contract patch — recovered

Full patch identifier supplied: `484c79ad286003f4fa51e850eaa116c8e2ca0c24`.

The recovered patch is a **root-contract patch only**. It adds three files:

- `tests/data/issue4_affine_root/generate_reference.py`
- `tests/data/issue4_affine_root/issue4_affine_root_reference.json`
- `tests/test_issue4_cmc_affine_root_contract.py`

Declared behavior and scope:

- Synthetic `_stage_state` fixture only.
- Solver residual gate: `2e-12`.
- Root accuracy is bounded using the extended-precision reference and a conditioning-derived estimate involving `sqrt(n) * gate / sigma_min(I-J)`.
- The original strict assertions remain unchanged as historical failure witnesses.
- No production equations, solver implementation, vendor code, tolerance policy, or scientific-admission state is changed.

**Preserved defect:** in the pasted `generate_reference.py`, `compare()` has a `return dict(...)` expression whose closing parenthesis for `dict()` is missing. The generator therefore cannot execute in the recovered pasted form. The test reads the committed JSON and does not import that generator, so this defect does not by itself establish that the test's reference data is unreadable. Do not silently repair the recovered patch.

**Provenance distinction:** the recovered JSON records mpmath `1.4.1`, Python `3.12.3`, and NumPy `2.5.4`. This is distinct from the integrated Slack handoff's stated research-job environment (mpmath `1.3.0`, Python `3.12.13` used locally as the closest available version). The exact pinned CPython `3.12.15` run was still outstanding in the handoff. These environment records must not be collapsed into one claim.

### 2. High-precision implicit-derivative proposal — recovered separately

Recovered as a proposal, not as the missing integrated commit.

It concerns the synthetic Issue #4 fixture and the eight-field direction `(a,b,X,Aa,K,rho,pr,pt)`. The proposal reports:

- The implicit tangent residual check (`<2e-10`) passes.
- The original finite-difference comparison of solver-returned roots at epsilons `2e-5` and `1e-5` fails at the original comparison settings (`rtol=4e-5`, `atol=2e-8`).
- The proposal attributes the discrepancy to float64 root-solve error amplified by division by epsilon, and compares the tangent against high-precision exact-root finite differences instead.
- Reported high-precision tangent-versus-exact-root-FD maximum difference is approximately `7.913e-13` at both primary epsilons; reported tangent residual is approximately `1.6e-17`.
- The proposal keeps the original strict test as a visible witness and calls for pinned-environment regeneration and an explicit owner decision on the diagnostic.

These are proposal-reported values, not independently rerun or promoted by this documentation change.

**Do not conflate two tolerance formulations:** the proposal describes comparison under the original `rtol=4e-5`, `atol=2e-8` settings. The later integrated Slack handoff describes the integrated derivative contract as requiring **both** absolute maximum error `<=2e-8` **and** normwise relative error `<=1e-7`, and mentions a committed negative control rejecting a `1e-6` relative error. Preserve these as distinct proposal and handoff versions until the original integrated patch is recovered.

### 3. October 25 fixed-R ledger correction — recovered separately, not part of this sequence

A separate pasted documentation diff concerns:

- `docs/research/CONSERVATION_LEDGER_SOURCE_CROSSWALK_2026-10-23.md`
- `docs/research/FIXED_R_SHELL_LEDGER_OPERAND_RECOVERY_2026-10-25.md`

It records operand-recovery/provenance corrections and clarifies the interpretation of the fixed-coordinate shell-ledger sweep. It is not one of the Nov 5–6 Issue #4 split-derivative commits and must not be listed as part of the integrated patch.

## Integrated Nov 6 handoff — described contents, original patch missing

The Slack handoff identifies the first reported commit, `2937ef96`, as a broader test-only integration. In addition to the root-contract work, it describes:

- `tests/test_issue4_implicit_derivative_contract.py`.
- A shared research helper (exact path/name not recovered).
- A research port guard comparing the float64 tool map with an mpmath port.
- A detection-floor record (exact path/name not recovered).
- Regeneration of both reference datasets into temporary files and comparison against committed data.
- A blocking derivative contract with tangent residual `2e-10`, exact-root central-difference checks at epsilons `2e-5` and `1e-5`, absolute maximum error `<=2e-8`, and normwise relative error `<=1e-7).
- Visible non-blocking historical witnesses in `tests/test_issue4_cmc_lapse_tangent.py`, reportedly byte-identical to base `f8d5a65`.
- A non-blocking research-diagnostics job; the existing accepted-state crosswalk job reportedly remained unchanged.

The handoff reports local results under Python `3.12.13`, NumPy `2.5.4`, pytest `9.1.1`, and mpmath `1.3.0`: blocking contract 7 passed; historical/original tests 2 failed and 16 passed (described as expected); port guard 2 passed; regenerated references string-identical. The exact CPython `3.12.15` CI run had not yet been completed. These are handoff-reported results, not results independently reproduced in this recovery task.

## Missing original artifacts

| Artifact | Recovery status |
|---|---|
| `484c79ad286003f4fa51e850eaa116c8e2ca0c24` root-contract patch | Full pasted contents recovered; generator syntax defect preserved above |
| High-precision derivative proposal | Recovered separately; not the integrated commit |
| Oct 25 fixed-R ledger correction | Recovered separately; unrelated to this sequence |
| `2937ef96` integrated test-only patch | Exact original patch bytes and full commit object not recovered |
| `8965a198` workflow-separation patch | Exact original patch bytes and full commit object not recovered |
| `issue4_integrated_test_only.patch` | Not recovered as a file; only handoff description is available |
| `issue4_ci_separation_workflow.patch` | Not recovered as a file |
| Original two-commit Git bundle | Created and `git bundle verify` passed in the prior local session; exact size and SHA-256 recorded above. Slack upload failed (`F0CFA70KUE9`, deleted/no data). Bundle bytes are not present in the current runtime, so a fresh byte-level verification or transfer is still outstanding. |
| Shared helper, port-guard, detection-floor and derivative-reference exact paths/content | Not fully recovered |
| Exact pinned CPython `3.12.15` workflow result | Outstanding in the handoff; no run is established by this audit |

The handoff describes the second commit, `8965a198`, as adding workflow separation plus a branch trigger to `.github/workflows/issue4-cmc-lapse-tangent.yml`. This document does not recreate that workflow or claim that it was pushed or run.

## Recovery rules and next actions

1. Preserve the recovered `484c79ad` patch and derivative proposal as separate artifacts. Do not rewrite either to make them appear to match the later integrated handoff.
2. Recover the already-created bundle bytes from the prior local output path or a user-provided copy. The recorded prior-session path is `/mnt/user-data/outputs/issue4-integrated-test-diagnostics.bundle`; it does not exist in the current runtime. Do not treat the failed Slack file ID as a usable attachment.
3. When the bytes are available, recompute SHA-256, run `git bundle verify`, confirm the sole prerequisite is `f8d5a65`, and inspect the two commit objects and their file inventory/diffs against this ledger.
4. Separately recover the original patch files if available; the bundle may be sufficient to recover the commit trees, but the exact named patch files were not confirmed as uploaded.
5. Keep the generator syntax defect and all environment/version differences visible. Any repair or regeneration must be a separately identified change with its own provenance. Do not create substitute commits under the missing SHAs or label a reconstructed implementation as recovered original history.

## Scope and admission

This record changes documentation only. It does not add or modify tests, solver code, workflows, production equations, initialization, vendor code, tolerances, scientific-admission logic, or accepted-state evidence. It launches no numerical run and does not establish constraint propagation, conservation closure, or scientific admission. Preserve `scientific_admission=false` and all existing open gates.
