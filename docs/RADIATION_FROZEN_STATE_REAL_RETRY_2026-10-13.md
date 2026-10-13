# Real Adaptive Radiation Retry from the Archived Failed Production State

**Date:** 13 October 2026 UTC  
**Classification:** Bounded end-to-end software/one-step gate; not a long-horizon physical admission.  
**Candidate branch:** `fix/radiation-adaptive-step-on-cone-2026-10-12`  
**Exact tested candidate head:** `651dc124484e7c056b784033b69bedba6da0b5a5`  
**Source baseline:** `cdcec8282bff0443f2fe4674b15929822ff38736` (stacked radiation RHS/rest-current candidate).  
**Pinned numerical vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Workflow:** [38325032498](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38325032498) — success.  
**Artifact:** [11781015455 — radiation-frozen-state-real-retry](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38325032498/artifacts/11781015455)  
**Artifact ZIP SHA-256:** `efa04c7274ba0d68a2683c4f874c7cea4aac2f6d9b9f1a183d4f672779a993f2`  
**Embedded report SHA-256:** `3e3911c6574ce5de67fc9914f37ab5bbcb56ab412a8d9058e2188bdc2337c763`  
**Accepted-state NPZ SHA-256:** `ccd1fa0fbabe65c51c3dd7ac2c9227d6d25aa680c169eef1afd3fca4326294cf`  
**Accepted-state metadata SHA-256:** `0e6251ec58b843d5d5f0dfc8911afea0dbc5da27b409f34d3435628c61676858`

## Executive result

**The adaptive retry now succeeds on the actual archived failing production state.** This closes the previous gap where retry was covered only by a synthetic fake kernel.

The workflow recovered the exact accepted numerical arrays immediately preceding the production failure at \(t=47.130000000001026\), verified the state-file SHA-256 and exact source/vendor metadata, then invoked the real production kernel through the candidate's `advance_with_radiation_admissibility_retries` helper with requested \(dt=0.03\).

The actual attempt sequence was:
1. \(dt=0.03\): rejected by the typed `RadiationAdmissibilityError` at `stage='explicit-predictor'`, cell \(i=79\), \(r=79.5\), ratio \(|S|/E=1.0004012501659771\).
2. \(dt=0.015\): accepted by the actual production kernel, reaching \(t=47.145000000001026\), with full-grid maximum \(|S|/E=0.9996001341484262\), and minimum densitized cone margin \(C=2.038846200348128\times10^{-9}>0\).

Retry count is **one**; accepted step size is exactly \(0.015\). A numerical-array digest before and after the helper is identical (`6ec904c223640d1b43d84946b54b909032705705727bc395dd1b87c774ccedf9`), verifying the failed attempt did not mutate the saved accepted input arrays. The helper changed only timestep and did not clip or project matter, change a boundary condition, relax the cone, or alter any physical source.

## Gate results

- Exact source/vendor scope check passed; only the allowed retry candidate/test/harness paths differ from the frozen stacked baseline.
- Focused adaptive + radiation admissibility/rest-current tests: **16 passed**.
- Full repository test suite: **95 passed**.
- Archived accepted-state NPZ SHA-256 checked against expected digest.
- Actual first rejection matched the stage, cell, and known ratio within a fail-closed numerical tolerance.
- Second production step was accepted, its coordinate-time increment matched the accepted \(dt\), and the final full-grid radiation state satisfied \(C\geq0\).
- Artifact uploaded with exact provenance and report checksum.

Separately, the underflow logic no longer enforces an unexplained \(8\mathrm{ULP}\) minimum. It stops halving only when the candidate step no longer advances the stored floating-point coordinate time, and a regression verifies that \(dt=10^{-16}\) can be accepted at \(t=0\) when it genuinely advances time. This change passed the same pinned full-suite gate above.

## Full-state witnesses before and after the accepted retry

| Witness | Archived accepted state at \(t=47.13\) | Actual accepted retry at \(t=47.145\) |
|---|---:|---:|
| Maximum full-grid \(|S|/E\) | 0.9988576965550366 | 0.9996001341484262 |
| Maximum Hamiltonian residual | 0.084934723579255 | 0.08505320892648258 |
| Maximum momentum residual | 0.018343263859204154 | 0.01840765888000545 |
| Maximum connection residual | 0.05829447235243363 | 0.05841266581648077 |
| Maximum determinant-constraint residual | \(1.63892784490\times10^{-5}\) | \(1.64444345989\times10^{-5}\) |
| Misner–Sharp mass at \(R_\sigma\) | 4.501726366496262 | 4.591636438124435 |
| Surface areal radius \(R_\sigma\) | 10.494223796016815 | 10.49079806567956 |
| \(\chi_\sigma\) | 0.1420563437564698 | 0.12463543585956849 |
| \(H_{\rm eff}\) | -0.09091019278724598 | -0.09180154849558375 |
| Outer CMC residual maximum | \(2.36009193898\times10^{-4}\) | \(2.33779117741\times10^{-4}\) |
| Last-interior CMC residual | \(2.25537742968\times10^{-5}\) | \(1.81489931586\times10^{-5}\) |

These witnesses are included to stop the retry from being judged solely on the radiation cone. They show that the one-step state advances and remains radiation-admissible, while full sourced constraints and CMC residuals remain nonzero. The values of \(M_{\rm MS},H_{\rm eff}\), etc. changing in one successful step are not independently a conservation balance and do not imply convergence.

The archived NPZ contains the numerical field arrays but not the prior Python `CycleLedger` or history list objects. The audit intentionally uses fresh ledger/history objects; those objects do not enter the step's radiation predictor or local cone condition. This is documented in the embedded report.

## Interpretation and integration disposition

**Established:** the retry wrapper works against the real failing kernel state; it catches the intended typed radiation-admissibility rejection, retries from unchanged numerical arrays, accepts the first smaller step that passes, and records the attempt details. Unrelated errors are still outside the retry catch type and propagate normally in the tested helper contract.

**Not established:** a one-retry success from one saved state is not evidence that a full adaptive campaign reaches the desired horizon, that subsequent steps remain stable, that H/M/\(C_\Lambda\) propagate within an admitted tolerance, or that Misner–Sharp balance/convergence gates pass.

The candidate remains unmerged. Before integrating with campaign admission, preserve the separate strict-campaign-admission PR #3 plan for the shared `engine/run_production.py` surface; do not overwrite or independently merge the overlapping branches. Next safe gate is a bounded runner-level test that confirms a real retry is written into the production retry ledger on a controlled short run, followed only later by one preregistered multi-step test. No physics or boundary-law promotion is authorized by this artifact.

No production `main`, pinned vendor, or canonical GEAR archive was modified by this audit.
