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
