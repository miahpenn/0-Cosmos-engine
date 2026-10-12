# Radiation rest-current representation repair — candidate handoff

**Date:** 12 October 2026  
**Status:** isolated repair candidate; not merged, not production-admitted.  
**Branch:** `fix/radiation-rest-current-channel-2026-10-12`  
**Current stacked base:** draft PR #10 head `6f6fd9d81af739cacce16fcd428bdfeac996e343`. The rest-current code was reapplied on this parent and revalidated below.

## Admitted diagnostic evidence

Diagnostic run [38317797540](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317797540) completed successfully on isolated audit branch `diag/radiation-rest-current-channel-audit-2026-10-12`. It verified the pinned vendor, compiled the audit/tests, passed 1 radiation-helper test, 4 RHS-admissibility tests, and **27 tests from the `tests/` subtree** before completing four one-step arms (generic vs radiation-specific dispatch at N=80/160). The `python -m pytest -q tests` command is not the whole-repository suite.

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
2. the actual radiation HLL dispatch keeps the first flux component exactly zero for moving left/right states, while the generic dispatch produces a nonzero first component;
3. the real production kernel step preserves an exact-zero radiation rest slot at N=80 and N=160 and still passes the existing E≥|S| state validator.

## Regression and CI evidence

### Previous gate (subtree only)

The original candidate gate [38317961532](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317961532) passed at code head `b9185aa160adaceb5917a9162892e466945ccbd2`:
- focused rest-current contract tests: **3 passed**;
- radiation RHS/admissibility tests: **4 passed**;
- `python -m pytest -q tests`: **29 passed in the `tests/` subtree**;
- pinned vendor verified.

The rebase gate [38318631134](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318631134) also passed at rebased code head `ec57b533ba5f8ee1b679ae854c8f89494d66f666`, using the same focused tests and `tests/`-subtree command. No evolution artifact is produced by this repair CI; the one-step A/B evidence remains in artifact 11777312270.

### Repository-wide validation

The workflow was upgraded from `python -m pytest -q tests` to repository-wide `python -m pytest -q`.

- Workflow commit: `85e34da02c97cb6829ece387ea917dcd6895e9f8`.
- Run [38318740928](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318740928): **3 rest-current tests passed, 4 RHS/admissibility tests passed, 86 repository-wide tests passed**.
- A direct moving-state HLL-dispatch regression was added in test-only commit `d100d27fe94efdecbd759a53fecd5ba2560229c5`. It asserts the radiation-specific HLL branch returns exact-zero rest flux for nonzero opposing face velocities, while the generic species path returns a nonzero first-component flux.
- Run [38318805927](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318805927) passed on that source/test head: **4 focused rest-current tests passed, 4 radiation RHS/admissibility tests passed, 87 repository-wide tests passed**. Vendor pin and module/test compilation also passed.
- The latest branch head is documentation-only commit `6c843a7f79663a348b585b234266422d1c9b17a4`; source and test files are unchanged from tested head `d100d27fe94efdecbd759a53fecd5ba2560229c5`.

## Guardrails and interpretation

- No clamp/floor, source, damping, limiter, radiation inequality change, or new coupling law.
- The rest-current patch preserves an existing representation invariant and does not change radiation E/S dynamics at captured one-step precision.
- The late-time corrected-rate N=80 cone crossing remains unresolved and must not be attributed to this rest channel based on these data.
- Misner–Sharp balance remains open.

## Disposition

- Draft PR #11: https://github.com/miahpenn/0-Cosmos-engine/pull/11
- Parent draft PR #10: https://github.com/miahpenn/0-Cosmos-engine/pull/10
- PR #11 is currently mergeable against PR #10's latest head but remains draft and unmerged.
- Main remains `2035acc1210db18da0cededea18750b7ed813260`; no main or physics change was made.
- Independent code review remains outstanding.
