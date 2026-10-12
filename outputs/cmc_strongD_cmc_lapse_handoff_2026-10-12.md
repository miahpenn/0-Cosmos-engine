# CMC strong-D lapse / global-target investigation — continuation handoff

**Prepared:** 12 October 2026 (UTC)  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Active diagnostic branch:** `diag/independent-metric-ricci-audit`  
**Branch head at preparation:** `d6142ac81a4bf56905de0fae5298f8fc540f7d98`  
**Production/main:** `main` verified at `7604bd6cc87c48c9dc8e24f2facbfe1458364ca5`; untouched by this diagnostic work.  
**Purpose:** Preserve verified numerical findings, record the interpolation audit and its repaired workflow, and give the next session an exact continuation path without mistaking a diagnostic counterfactual for a physical result.

## Executive state

The (0.18\%) comparison between the saved-time integral of the central-lapse difference and the R80–R160 final proper-time difference is a **bookkeeping/quadrature check**, not causal attribution. Since (\tau=\int \alpha_c\,dt), the proper-time difference is by definition the integral of the central-lapse-history difference; close agreement checks sampling/integration consistency. The same-slice lapse re-solve correction is negligible at the relevant scale.

The global CMC target (\dot K) is a credible domain-sensitive gauge mechanism: the target is a proper-volume projection of the full-domain K RHS. R80 and R160 can therefore obtain different target histories even when a shared local region initially matches. The frozen-state target swap demonstrates that this global target changes lapse profiles, but does not reproduce the observed center-lapse/proper-time history; the evolving feedback test remains the key experiment.

**New completed result:** the R160 target-history interpolation audit now passes after a workflow dependency correction. Its artifact shows a smooth resolved target curve and a small held-out interpolation error over a deliberately doubled interval. This is empirical evidence only, not a rigorous bound at every evolution stage.

**Still pending:** full evolving R80 trajectory using interpolated archived R160 target history, run [38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450). At the last status check it remained in progress and had not yet uploaded its output artifact. No conclusion about target-history causality may be drawn until that artifact is inspected.

## 1. Core state and numerical ledger

### 1.1 Baseline comparison

- R80 boundary-matched trajectory: [run 38162318255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255), archived profile artifact `11710858125`; source commit `35072830f9e4d7029bb10243a0520b9a273a17dd`.
- R160 reference: archived profile artifact `11694568292`; exact source ZIP SHA-256 `04949c7e0332c495385f6afef231b54270225892b1d124f2c19973c7d5b4872d`.
- R80 profile artifact `11706797173`; exact source ZIP SHA-256 `8b118d1f81e01f2a5d291d160b8bac42abc6a11051d8453efd9c19d0826d6d8b`.
- Settings relevant to comparison: R80 (N=160,r_{max}=80,dr=0.5); R160 reference (N=320,r_{max}=160); strong-D trajectory with initial (D) amplitude (10^{-4}), scalar amplitude (0.01), width (7), radiation on, CFL (0.0075), through (t=12).
- There are 49 shared saved-time records from (t=0) to (12), spaced approximately a quarter time unit apart.

Observed final accumulated proper-time difference:
\[
\tau_{80}-\tau_{160}=+0.01062109.
\]
The trapezoidal integral of the *actual saved central-lapse history difference* is (+0.01063961). Difference between these is (1.8512\times10^{-5}), about (0.18\%) of the accumulated gap. This supports consistency of the saved-time ledger; it does not establish which code path caused the lapse histories to separate.

Largest saved central-lapse difference: approximately (0.00474876) near (t=6.49875), not at the final saved slice.

### 1.2 Frozen-slice re-solve correction

- Strong-D same-slice Jacobian/residual run: [38205542397](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205542397).
- Full saved-time correction-series run: [38205709871](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205709871), artifact `11727266816`.

At (t=12), the input-slice discrete CMC residual maximum was (2.3133\times10^{-10}) (relative maximum (2.2690\times10^{-10})); same-target solve relative residual closed to about (4.04\times10^{-15}). Re-solving the lapse on the frozen state changed the center by only (4.0921\times10^{-9}), with RMS correction (9.4410\times10^{-8}).

Across 49 snapshots, the trapezoidal integral of the frozen-slice central-lapse re-solve correction was (2.6818\times10^{-8}) signed and (2.6860\times10^{-8}) absolute, only (2.53\times10^{-6}) of the observed (0.01062109) proper-time gap. This route is therefore much too small to be the direct explanation. It remains a frozen-state diagnostic, not a re-evolution.

### 1.3 Frozen-state global-target counterfactual

- Run [38205934746](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205934746), artifact `11727417247`.
- Script: [cmc_domain_target_counterfactual.py](../blob/8b2bafa5ec192e31fb70906ff99e9c14f4696930/diagnostics/cmc_domain_target_counterfactual.py).

Native R80 and archived R160 full-domain targets:
- At (t=0): R80 (dot K\approx6.19884\times10^{-5}), R160 (dot K\approx6.23352\times10^{-5}). Supplying the R160 target on the initial R80 slice reproduces the R160 lapse over the shared core to about (10^{-12}).
- At (t=12): R80 (dot K\approx2.45708\times10^{-4}), R160 (dot K\approx7.50581\times10^{-5}).

On frozen slices, swapping to R160's target significantly reduces shared-core radial lapse-profile RMS error at later times, but overshoots the central lapse. The signed integral relative to the R160 center history is about (-0.03777), opposite sign and larger than the observed (+0.01062) proper-time gap. Interpretation: the target affects the lapse profile, but a frozen-state target transplant does not reproduce the center history; evolving geometry/matter feedback must be included.

## 2. Interpolation audit — resolved first-run setup failure

### 2.1 Failed initial attempt was not a numerical finding

- Initial audit run [38206744317](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206744317) failed before evaluating any target samples.
- The actual exception was `ModuleNotFoundError: No module named 'numpy'` at script import. The workflow omitted installation of the audit's NumPy dependency in the commit that launched it.
- This was an environment/setup failure, not a failed interpolation result.

### 2.2 Dependency fix and verified rerun

- Fix commit [`d6142ac81a4bf56905de0fae5298f8fc540f7d98`](https://github.com/miahpenn/0-Cosmos-engine/commit/d6142ac81a4bf56905de0fae5298f8fc540f7d98), message: `[cmc-kdot-interpolation-audit] install NumPy before audit`.
- Successful audit run [38206848868](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868).
- Artifact [`cmc-kdot-interpolation-audit`, ID `11727926565`](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868/artifacts/11727926565), ZIP SHA-256 `993340f756973e2a0251da03269ea37d38dd1591a0cffa45e9540148c471c39b`.
- Audit source: exact R160 archived profile ZIP above; internal JSON reports profile SHA-256 `e0178a5e14ba3806fcdd11f9a1c912abffef0f640da468f06037f088d1bc47fe`.

### 2.3 Measured held-out interpolation errors

The audit retains every other interior saved target and estimates the withheld target from adjacent retained samples. The test interval is therefore about twice the archive's ordinary snapshot spacing: roughly (0.499)–(0.503) instead of roughly (0.25).

| Quantity | Measured value |
|---|---:|
| Saved target samples | 49 |
| Time span | (0\) to (12) |
| Native snapshot spacing | (0.2475)–(0.25125) |
| Held-out points | 24 |
| Held-out interpolation span | (0.49875)–(0.50250) |
| Target range | (5.94546\times10^{-5}) to (7.50581\times10^{-5}) |
| Maximum absolute error | (5.29437\times10^{-8}) |
| RMS absolute error | (1.82802\times10^{-8}) |
| Median absolute error | (4.65641\times10^{-9}) |
| 95th percentile absolute error | (4.82885\times10^{-8}) |
| Maximum relative error at held-out samples | (8.8267\times10^{-4}) ((0.0883\%)) |
| RMS relative error at held-out samples | (3.0437\times10^{-4}) ((0.0304\%)) |
| Largest error location | near (t=6.25125) |

This is encouraging evidence that the stored (\dot K) history is smooth at the saved-sample scale. It is **not** a rigorous error bound for the linear interpolation at each internal PIRK/evolution stage. It cannot rule out faster between-snapshot variation, and it does not translate directly into a bound on accumulated (\tau), because the gauge/evolution response can be nonlinear. Treat it as a measured interpolation-risk diagnostic, not certification that interpolation has zero influence.

## 3. Full evolving target-history counterfactual — outstanding

Run [38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450), launched on diagnostic branch commit `51045b8e0dfa8a00d167f7df78da9f07cfa772c5).

It runs one full R80 trajectory with:
- archived R80 outer-cell matter-boundary history;
- archived R160 outer-alpha history;
- archived R160 full-domain (\dot K) target history, linearly interpolated in time and prescribed at every CMC lapse solve;
- the same production evolution equations and positivity gate; no production source edits.

The workflow verifies both input archive ZIP checksums before running. At the last status check, it was still in the evolution step and had not uploaded an artifact. Do not infer the final result from the in-progress state.

Once it completes, download artifact `cmc-r80-r160-kdot-counterfactual` and inspect `report.json`, `comparison.json`, `campaign_progress.json`, and the profile JSONL. Required admission checks:
1. Evolution reaches (t=12); 49 records including the initial slice; no nonfinite or nonpositive geometry/lapse.
2. All target calls are finite; maximum prescribed-versus-returned target delivery error is effectively zero; target call times cover all lapse solves rather than only checkpoints.
3. The report's failure field is null and the comparison uses the exact baseline R80/R160 archives verified by SHA.
4. Compare the target-forced final (\tau), central lapse history, and full shared-core lapse-profile errors with both native R80 and R160. Report the absolute gaps and fraction closed; do not describe a partial closure as full attribution.
5. Note the counterfactual's limitations: (\dot K), outer alpha, and outer matter-face histories are imported from archived trajectories. This is a gauge sensitivity/causal intervention on the R80 evolution, not a self-consistent R160 spacetime solved on an R80 domain.

Interpretation:
- If target forcing closes most of the proper-time and central-lapse history gap, it is evidence that the prescribed global target history causally drives the divergence within this counterfactual. It still does not establish that the archived R160 target is uniquely physically correct on R80.
- If it does not close the gap, global-target mismatch alone is insufficient. The next investigation should target evolving geometry/matter feedback and the shared boundary forcing, keeping those layers separate.
- If it fails a positivity or evolution gate, preserve the partial progress and failure time/profile; do not tune or soften the gate.

## 4. Next-work order

1. **Finish and admit run 38206154450.** Capture its exact conclusion, artifact ID/hash, target-call delivery statistics, failure field, final time and (\tau), and checkpoint comparisons.
2. **Reconcile interpolation risk with the trajectory result.** The held-out test is small at resolved scale but is not a stage-level bound. If target forcing appears causally important, the next single expensive diagnostic can compare the native quarter-unit target history with a deliberately coarsened every-other-sample history, using the same initial state and same boundary histories, to directly measure outcome sensitivity to interpolation resolution. Do not start this while the current full run is still active.
3. If interpolation sensitivity is small, attribute only to target-history sensitivity to the degree directly shown by the trajectory. Otherwise report the sensitivity interval and do not overclaim.
4. If target history alone does not close the discrepancy, isolate the evolving geometry/matter response next; avoid changing several layers at once.
5. Keep all work diagnostic-only. Do not alter the production equations, add fitted feedback, clamp/floor/damp, or modify `main`.

## 5. Provenance and guardrails

- This handoff is a diagnostic record, not an admission that the engine has passed the broader (0^*\) criterion or established a physical bounce/no-bounce interpretation.
- The (0.18\%) proper-time/integrated-lapse agreement is a consistency check; it is not causal evidence.
- The tiny same-slice correction is quantitatively insufficient as the direct source of the proper-time gap.
- The global full-domain CMC target is a concrete plausible route for domain-size dependence; frozen target swap establishes influence on profiles but does not close the observed center-history gap.
- The interpolation audit is now successful with a setup fix; it measured rather than guaranteed interpolation accuracy.
- Full evolving target-swap result is still pending.
- The verified `main` head remains `7604bd6cc87c48c9dc8e24f2facbfe1458364ca5`; do not promote diagnostic findings to production without a separate review and exact successful CI.
