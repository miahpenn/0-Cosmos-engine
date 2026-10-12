# Rate-corrected radiation four-force one-step comparison — handoff

**Date:** 12 October 2026  
**Branch:** `diag/radiation-four-force-ratefix-ab-2026-10-12`  
**Base:** PR #10 head `730d6788bff2b53cf84d0504842f95979953f1bd`  
**Purpose:** Run the same diagnostic-only N=80/160 × beta_rad=0/-0.04 matrix using the corrected radiation RHS/admissibility API, so the candidate is not judged by a unit-step wrapper artifact.

## Why this branch exists

PR #9's earlier run [38316415773](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38316415773) passed both N=80 one-step cases but rejected both N=160 cases while `species_rhs()` attempted `evolve_species(dt=1.0)` for derivative extraction. That routine's ordinary finite-time guard checked the artificial unit-step intermediate rather than the real caller timestep. The captured dt-scaled first-stage estimate had |S|/E about 0.01013, but that estimate was not the full PIRK state and could not answer the physical question.

Draft PR #10 isolates the API-contract repair without a radiation four-force: RHS extraction disables only that algebraic intermediate check, direct finite-time update validation remains enabled by default, and the production kernel checks the actual initial, explicit-predictor, primary-predictor, and completed-final states against each stage metric. Its current repair gates report 4 focused tests and 26 suite tests passing, with a full baseline N=160 step.

This branch starts from the PR #10 repair candidate and adds the optional PR #9 beta_rad source only as a diagnostic coefficient. It does not merge or modify PR #9 or PR #10.

## Fixed case matrix

- N=80, beta_rad=0.0
- N=80, beta_rad=-0.04
- N=160, beta_rad=0.0
- N=160, beta_rad=-0.04

Default initial data, r_max=40, radiation ON, CFL=0.03, one ordinary kernel.step per case. Full focused tests and the full repository suite run before the four-case comparison. The vendor pin is frozen to `d6052d605673ce9d82cdc99b0df79855fa2ea215`. A JSON file and SHA-256 are uploaded even if any case is rejected.

## Interpretation rules

A rejected case must be diagnosed from its actual stage and metric, not relabelled as the previous unit-step artifact. A passing one-step result is not evidence that beta_rad=-0.04 is a physically correct radiation interaction law; it is an observed counterfactual response from one step. The Misner–Sharp current/balance field is not assembled by this first-pass harness and must remain explicitly unavailable/NaN.

No projection, clipping, floors, damping, retry, relaxed cone, altered flux/source/recovery, or production promotion is authorized. The existing radiation admissibility inequality is preserved. PR #9, PR #10, PR #8 and the campaign-admission drafts remain open/draft; nothing is merged to main by this experiment.

## Result record

Workflow and artifact IDs, exact commit, pass/fail status for each case, JSON SHA-256, and the distinction between source-level and executed evidence will be appended here after the new Actions run is available.


## Result — full radiation profiles and matched-arm deltas

- Workflow [38317295748](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317295748): **success**; exact runner commit `83486b03c0789f571c76e7a5a4287b6afaf13330`.
- Regression results: four radiation RHS guard-repair tests passed, three null-identity tests passed, and all 29 tests in the repository suite passed.
- Artifact [11777006845](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38317295748/artifacts/11777006845); ZIP SHA-256 `9d29529a8096b6a3b46b8704ad0901bb2a3e637ad6f5373adbb4767bab227e34`; JSON SHA-256 `3e2a151895160a2e82caa8c8225f80af9b8721a53e5b7343869dc08eeaebcddb`, matches the manifest.

All four actual finite-time cases passed and all post-step radiation cells are finite/admissible. Maximum (|S|/E): N=80 null/candidate `0.0099086223695 / 0.0099086123118`; N=160 null/candidate `0.0099068411682 / 0.0099068361387`. Candidate-minus-null maximum full-profile changes: Eulerian E `1.84e-14 / 9.21e-15`; (|S|) `9.33e-17 / 4.67e-17); (|S|/E) `1.01e-8 / 5.03e-9`; densitized energy `2.88e-11 / 1.45e-11); densitized momentum `1.82e-15 / 4.55e-16`, for N=80 / N=160 respectively.

The very small response is a one-step counterfactual only. This `beta_rad` candidate adds a radiation source without a matching scalar source, so it is not yet a covariantly closed exchange model. Misner–Sharp current/balance remains NaN/unavailable in the harness; no production or full-machine admission follows.

### New separate audit lead: radiation rest-current slot

The post-step profile shows the initially-zero radiation `rest` slot drifting to `-1.00857e-11` (N=80 outer cell) and `-5.01262e-12` (N=160 outer cell). `initialize_radiation`, `_radiation_conserved`, and `_radiation_flux` indicate this first component is intended to be zero for radiation, but `_hll_flux` and the outer-face fill currently use the generic fluid `_valencia_conserved` / `_valencia_flux` for all species. Radiation energy/momentum recovery ignores this rest slot, and candidate-minus-null differences are only roundoff-scale; therefore do not attribute the late cone failure to it. Add an isolated regression/audit for this invariant before changing the transport code.
