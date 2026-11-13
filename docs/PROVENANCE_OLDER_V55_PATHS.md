# Older V5.5 paths: what is live, what is provenance only (2026-11-08)

Rule: everything that is part of production, or that produces readiness evidence, runs on the
current production implementation and is kept correct. Code kept only to reproduce earlier results
is frozen by commit, not maintained as a second live implementation.

## Provenance-only modules

| Module | Why it is not live |
|---|---|
| `engine/reference_pirk_unified.py` | Older copy of the V5.5 equations and initializer. `make_initial` maps the archive operating point into the wrong variables (φ ← a = 33.8983, Π ← φ = 0.179055; φ̇ = 0.00341913 unused), integrates the areal reconstruction with a plain trapezoid from r = 0 (the origin defect #439 fixed in production), and uses an outer-region CMC target that production replaced with the proper-volume projection. |
| `engine/true_cmc_reference.py` | CMC lapse/step for that stack. Its centre row was repaired in #442, but it still runs on the stack above. |

No live module imports them (enforced for the two drivers by
`tests/test_research_observables_consistency.py`). Do not use them for new results.

## Re-pointed to production

`engine/coupled_engine.py` and `engine/bidirectional_production.py` keep their function names but are
now thin drivers over `V55ProductionKernel` (same initial data, CMC gauge, matter projection,
worldtube extraction, proper time, cycle ledger and handoffs as production). Evidence:
`results/v55-drivers-on-production-kernel-2026-11-08/`.

## Reproducing earlier results

| Commit | State |
|---|---|
| `afc4d96f101cd0ce5960ce798746fd6b4ee127da` | Before the 2026-11-08 centre repairs (#436 CMC centre row, #439 initializer origin quadrature, #442 reference CMC centre row). |
| `41ee0cb4b8a641c889834d5062829f7d31fdb245` | After those repairs, before the drivers were re-pointed: the older V5.5 drivers still evolve `reference_pirk_unified`. |

Ledger, cycle or handoff numbers produced by `coupled_engine` / `bidirectional_production` at or
before `41ee0cb4` were computed on the mis-mapped initial data and should be treated as unverified
until re-run on current `main`.

## Other older paths (live, unchanged here)

The reduced two-lane machine (`engine/cosmos.py`, `local.py`, `interface.py`, `diagnostics.py`,
`run_controls.py`, run in CI by `controls.yml` and `lightweight_validation.yml`) and
`research/gear_archive_machine_v1.py` do not use the V5.5 reference stack. The GEAR-63/65 barrier in
`engine/local.py` is conditional physics (see `docs/PHYSICS_LOCK.md`), an owner decision rather than a
numerical defect. **Retired to provenance on Nov 11 2026** (DECISION-map-row-5): the GEAR-67 quartic is the one arrest
closure; `engine/local.py` is kept, unchanged in behaviour, for archive comparisons and the control runner only.
