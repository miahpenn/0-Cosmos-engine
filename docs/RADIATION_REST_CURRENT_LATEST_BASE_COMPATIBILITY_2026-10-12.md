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

## Branch hygiene and follow-up

At the time of the initial check, original [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11) had branched from an earlier PR #10 head; compare reported it six commits ahead and two behind. I did **not** force-update the shared PR branch. The PR owner subsequently rebased it onto the newer PR #10 head, which resolved the divergence:

- PR #10 current head at rebase: `6f6fd9d81af739cacce16fcd428bdfeac996e343`.
- PR #11 rebase/source commit: `ec57b533ba5f8ee1b679ae854c8f89494d66f666`.
- PR #11 own rebase gate: [workflow 38318631134](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318631134) — success: 3 rest-current tests, 4 radiation RHS/admissibility tests, 29 full-suite tests; pinned vendor verified.
- PR #11 current head: `5d4cddd665a0ab697723147dd9aa1680507def61`; compare against base is now **3 commits ahead, 0 behind**, and PR metadata reports mergeable/clean.
- The two commits after tested source head `ec57b533` change only the documentation handoff; source and tests are unchanged from the successful gate.

This compatibility branch remains a separate duplicate confirmation on base `d1338a8473aa8bc2474019f99b13001edde221a6`. The newer PR-owned gate is the stronger current-stack evidence. PR #11 is still draft: green regressions are necessary, not sufficient for merging or production admission.
