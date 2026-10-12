# Radiation Rest-Current Patch — Compatibility Check Against Latest RHS Contract

**Date:** 2026-10-12  
**Disposition:** latest-base compatibility gate passed; original stacked PR remains draft and has not been rebased/merged by this check.

## Provenance

- Base at test time: PR #10 current head `d1338a8473aa8bc2474019f99b13001edde221a6`.
- Compatibility branch: `diag/review-rest-current-on-latest-contract-2026-10-12`.
- Tested workflow/code head: `2ffe88951e8d226e5231d0bb5e9de074a29e8085`.
- CI [38318625839](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318625839): **success**.
- Vendor submodule verified at `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Radiation rest-current contract: **3 passed**.
- Existing radiation RHS/admissibility gate: **4 passed**.
- Full repository tests: **29 passed**.

## Exact scope

A clean branch was created directly from the latest PR #10 head. Only the two source hunks from PR #11 were then applied to the current `engine/matter_system.py`:
1. For `Species.RADIATION`, `_hll_flux` uses the existing `_radiation_flux` and `_radiation_conserved` helpers; all other species still use the generic Valencia helpers.
2. The radiation outer-face closure uses `_radiation_flux`; all other species keep the generic outer flux.

The latest PR #10 finite-update/admissibility contract was preserved without replacement. The exact same rest-current regression file from PR #11 was used. A dedicated workflow ran helper/source compilation, the rest-current regression, the current PR #10 RHS/admissibility regression, and the full suite. There were no physics coefficients, E/S formulas, source terms, geometry equations, boundary conditions, timesteps, or cone limits changed.

## Interpretation

This gate confirms the small PR #11 source patch composes with the **current** PR #10 code and tests cleanly. It does not repeat the diagnostic A/B (the existing report already shows the rest channel becomes exactly zero while all tested E/S, cone, geometry, and H/M/C\(_\Lambda\) profiles remain bit-identical). It also does not explain or repair the distinct long-time outer radiation-cone crossing, and it does not establish Misner–Sharp balance.

## Branch hygiene

Original [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11) still points to `fix/radiation-rhs-contract-...` history that branched before the latest PR #10 changes; GitHub compare reported it six commits ahead and two behind its intended base. This diagnostic branch does not force-update or merge the shared PR. Before PR #11 can be reviewed against its intended stack, its owner should rebase/recreate the PR head on the latest PR #10 head and confirm the current gate. Until then, treat the diagnostic branch only as a compatibility proof, not as evidence the stale PR stack is merge-ready.
