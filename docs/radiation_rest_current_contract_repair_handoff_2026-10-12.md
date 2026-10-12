# Radiation rest-current representation repair — candidate handoff

**Date:** 12 October 2026  
**Status:** isolated repair candidate; not merged, not production-admitted.  
**Branch:** `fix/radiation-rest-current-channel-2026-10-12`  
**Current stacked base:** draft PR #10 head `6f6fd9d81af739cacce16fcd428bdfeac996e343`. The original implementation/first passing repair gate used PR #10 base commit `5a5a9657e5d0e45c899addd5851ab2bd8db1cfe0`; this branch has since been rebased onto the latest PR #10 head.

## Admitted diagnostic evidence

Diagnostic run [38317797540](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317797540) completed successfully on isolated audit branch `diag/radiation-rest-current-channel-audit-2026-10-12`. The workflow verified the pinned vendor, compiled the audit/tests, passed 1 radiation-helper test, 4 RHS admissibility tests, and the 27-test full repository suite, then completed four one-step arms (generic vs radiation-specific dispatch at N=80/160).

- Artifact [11777312270 — radiation-rest-current-channel-audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317797540/artifacts/11777312270).
- Artifact ZIP SHA-256: `e34c99e6b575222f0f8ce5c134f9ac2bdb9a958ece0dc235f593a8a63e38e090`.
- JSON SHA-256: `566b51387ca3272fa76716b19d5e1128450ef02c6ae784ee851f2759697a785b`, verified against the archive's own SHA-256 file.
- Vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

The generic path produced post-step outer rest values `-1.0085731565667868e-11` at N=80 and `-5.012616800265632e-12` at N=160, with nonzero rest values at every cell. A diagnostic-only dispatch through existing `_radiation_flux` / `_radiation_conserved` helpers produced exactly zero rest values. Across both resolutions, every radiation E/S value, every derived cone observable, every geometry field and H/M/C_Lambda summary was **bit-identical** between arms.

This isolates the extra component as passive for these tested steps: radiation recovery/projection uses the energy/momentum channels, which did not change. It supports fixing the channel wiring, but does not establish that this issue caused the separate late-time outer-cell cone crossing.

## Minimal candidate change

On this branch only, `engine/matter_system.py` dispatches `_hll_flux` to the existing radiation-specific conserved/flux helpers when `species is Species.RADIATION`, and dispatches the outer face through the radiation-specific flux helper as well. Other species continue to use generic Valencia helpers. No equation for E or S, source term, cone test, time integration, boundary state, or geometry operator changes.

The added regression `tests/test_radiation_rest_current_contract.py` checks:
1. existing radiation-specific helper first components are exact zero and the energy/momentum components match generic helpers exactly;
2. the real production kernel step preserves an exact-zero radiation rest slot at N=80 and N=160 and still passes the existing E≥|S| state validator.

## What must pass before review

The dedicated repair gate verifies the vendor pin, helper contracts, actual N=80/N=160 one-step invariant, existing radiation admissibility tests, and the full repository test suite. Any test failure keeps the candidate unpromoted. The diagnostic A/B is evidence for only this representation contract, not a production physics or full-run admission.

## Guardrails

- No changes to `main`, the canonical GEAR archive, PR #10, or the separate beta_rad diagnostic branch.
- No clamp/floor, source, damping, limiter, radiation inequality change, or new coupling law.
- The late-time corrected-rate N=80 cone crossing remains unresolved and must not be attributed to this rest channel based on these data.
- Misner–Sharp balance remains untested by this one-step rest-current audit.


## Dedicated repair CI result — 12 October 2026

The candidate gate is complete and successful:
- Workflow: [38317961532](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317961532) — success at tested head `b9185aa160adaceb5917a9162892e466945ccbd2`.
- Source dispatch edit: `1e41d4fa4fc63d4f14bc135670865d4e44f48e11`.
- Focused rest-current contract tests: **3 passed**.
- Radiation RHS / stage-admissibility tests: **4 passed**.
- Full repository suite: **29 passed**.
- Pinned vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215` verified.
- No evolution artifact is produced by this regression gate; the audited one-step evidence remains in [artifact 11777312270](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317797540/artifacts/11777312270).

The tested source makes two narrow species-dispatch changes in `engine/matter_system.py`: radiation HLL face fluxes use `_radiation_flux` and `_radiation_conserved`, and the outer closure uses `_radiation_flux`. It does not change E/S formulas: the helper-contract test proves the existing radiation helper E/S components equal the generic helpers exactly, while the actual-step regression proves radiation's rest-current remains zero on both resolutions.

This candidate remains **stacked on draft PR #10**, not promoted to main. The diagnostic A/B and passing candidate gate establish the radiation rest-channel contract only. The late-time corrected-rate outer cone crossing and Misner–Sharp conservation remain open.


## Review disposition

The candidate gate passed at tested head `b9185aa160adaceb5917a9162892e466945ccbd2`: 3 focused rest-current tests, 4 radiation RHS/admissibility tests and 29 full-suite tests passed. The code candidate is now draft stacked [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11), with base [PR #10](https://github.com/miahpenn/0-Cosmos-engine/pull/10). No merge was made. The later report-only branch commit does not change the source or test files covered by that run.

The main-branch checkpoint writer currently serializes `radiation_rest`; the separate A/B confirms the generic path can put the spurious rest drift into that saved field. Radiation-specific recovery and stress-energy projection continue to use E/S. Fixing the species dispatch therefore preserves an existing state representation contract in runtime state and checkpoint output without changing E/S dynamics.


## Stacked-base revalidation status — 12 October 2026

PR #10 advanced after the original passing PR #11 repair gate, briefly making PR #11 non-mergeable. The rest-current branch has been re-applied onto latest PR #10 head `6f6fd9d81af739cacce16fcd428bdfeac996e343` at commit `ec57b533ba5f8ee1b679ae854c8f89494d66f666`. GitHub now reports [draft PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11) mergeable against that base. The new stacked gate [38318631134](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318631134) is queued/running for the rebased head; the prior pass [38317961532](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317961532) proves the pre-rebase source candidate only. Do not call the rebased stack fully revalidated until run 38318631134 succeeds.

No main or physics change was made. The minimal source diff and existing regressions remain the same; main is still `2035acc1210db18da0cededea18750b7ed813260`.
