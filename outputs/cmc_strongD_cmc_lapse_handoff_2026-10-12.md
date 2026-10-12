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

**Completed:** full evolving R80 trajectory using interpolated archived R160 target history, [run 38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450), artifact [11727976800](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450/artifacts/11727976800). It reached t=12 in 3,200 accepted steps with 49 records, null failure, finite prescribed targets, and zero measured target-delivery error across 12,801 lapse solves. It closed 99.9941% of the final proper-time gap in the boundary-matched comparison. The next test is sensitivity to coarser target-history sampling.

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

## 3. Full evolving target-history counterfactual — completed

- Run [38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450): successful.
- Artifact [cmc-r80-r160-kdot-counterfactual, ID 11727976800](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450/artifacts/11727976800).
- Artifact ZIP SHA-256: 63835669c695f9577d0416a2aa1fb3740aef9d4a90f36aaf9b9dbb94a5783006.
- report.json SHA-256: 454515cf736fcbe75e210c06285259ac138aca7b730863107aeae402b7697994.
- comparison.json SHA-256: 2cd4146a6e04b9f07a6e6cc665ab9d058c23a1140e8761e7d21c87efc696bb43.
- Profiles JSONL SHA-256: bbb35f6377fc0e399ecc9f7b66e3557ebce084f5aaf9e5d7aaa81a868ae9bd41.

The trajectory reached t=12 in 3,200 accepted steps and emitted 49 records including t=0. Final tau was 6.9067799130064635; final center lapse was 0.017277267256200452; failure was null. The positivity/finiteness assertion passed at each accepted step. Both source artifacts were checksum-verified. All 12,801 CMC lapse solves received a finite prescribed target, spanning stage times 0 through 12; maximum prescribed-versus-returned target difference was exactly zero. The matter-boundary provider ran 6,400 times. Production physics/kernel source was unchanged.

### 3.1 Proper-time closure

The previous boundary-matched R80 run with native R80 target history [38162318255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255) had final tau 6.917400377028174 and final gap versus R160 of +0.010621093774165757. Its implied R160 reference tau is 6.906779283254008.

With the archived R160 global Kdot history prescribed at every lapse solve:
- Target-forced R80 final tau: 6.9067799130064635.
- Remaining final gap versus R160: +6.297524555165523e-7.
- Final proper-time gap closed: 99.9940707%.
- Maximum saved-slice absolute proper-time difference decreased from 0.0121899232135 in the native-target boundary-matched trajectory to 0.0000121914168: 99.89999% reduction.

At t=12, the maximum absolute lapse difference on the shared core versus R160 decreased from 0.0521882686 in the native-target boundary-matched run to 0.00000209140736 in the target-forced run: 99.99599% reduction. The all-saved-slice lapse-profile RMS error decreased from 0.0302839287 to 0.000000782941. The report also shows major reductions in shared-core discrepancies in geometric K, Aa, a, b, X and Lambda. Some matter differences remain near the finite outer region.

**Interpretation:** within this boundary-matched R80 counterfactual, replacing the native R80 global target history with archived R160 target history closes essentially all of the proper-time and lapse-profile gap. This is strong causal evidence that target history drives the observed divergence in this controlled setup. It is not proof that the R160 target is uniquely physical or self-consistent on R80: the target, outer-alpha boundary and matter-face histories are imported from archived trajectories.

The completed interpolation audit [38206848868](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868) measured maximum held-out Kdot interpolation error 5.29437e-8 and RMS error 1.82802e-8 for an approximately 0.5-time-unit interval. This is small at the measured resolved scale but is not a rigorous stage-level bound. A coarsened-history full evolution is the next test of its actual impact on the coupled trajectory.

## 4. Next-work order

1. The full-resolution target-forced run is admitted; its result and hashes are recorded above.
2. Run one coarsened target-history counterfactual: retain every other archived R160 Kdot sample (while retaining the final t=12 endpoint), then linearly interpolate between retained samples. Keep initial data, R160 outer-alpha boundary history, R80 matter-face history, equations, grid, CFL, positivity gate and saved comparison times unchanged.
3. Compare coarsened-target tau, center-lapse history and lapse-profile errors against both the full-resolution target-forced trajectory and the native-target boundary-matched R80 trajectory. If close, it supports low sensitivity to measured interpolation resolution. If materially different, report the sensitivity range and do not assign the full 99.994% closure to the target alone until resolved.
4. Keep all work diagnostic-only. Do not alter production equations, add fitted feedback, clamp/floor/damp, or modify main. Do not start another expensive run until the coarsened-target run completes.

## 5. Provenance and guardrails

- This handoff is a diagnostic record, not an admission that the engine has passed the broader (0^*\) criterion or established a physical bounce/no-bounce interpretation.
- The (0.18\%) proper-time/integrated-lapse agreement is a consistency check; it is not causal evidence.
- The tiny same-slice correction is quantitatively insufficient as the direct source of the proper-time gap.
- The full evolving counterfactual strongly implicates the global full-domain CMC target history: final proper-time gap closure is 99.9941% and maximum saved-slice gap reduction is 99.90%, conditional on the archived boundary inputs.
- The interpolation audit is now successful with a setup fix; it measured rather than guaranteed interpolation accuracy.
- The next run is a coarsened target-history sensitivity test; do not change production physics before reviewing it.
- The verified `main` head remains `7604bd6cc87c48c9dc8e24f2facbfe1458364ca5`; do not promote diagnostic findings to production without a separate review and exact successful CI.


## 8. Follow-up started 12 October: target-history stride-2 sensitivity

After the full target-forced result was admitted, the diagnostic script was extended to support an optional integer target sampling stride. The existing default remains stride 1, preserving the full archived target history. The new diagnostic run uses stride 2, retaining every other archived R160 Kdot row and preserving the final endpoint at t=12; the linear interpolation span is therefore about 0.5 time units.

- Script update commit: [2657c5d4f13ddee30ed0c446a99bcb60720ebf0f](https://github.com/miahpenn/0-Cosmos-engine/commit/2657c5d4f13ddee30ed0c446a99bcb60720ebf0f).
- New workflow commit: [7d53a3d78d628aeec1af88b392a43136bd087142](https://github.com/miahpenn/0-Cosmos-engine/commit/7d53a3d78d628aeec1af88b392a43136bd087142).
- Stride-2 full-evolution workflow: [run 38207554072](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207554072).
- Expected artifact name: `cmc-r80-r160-kdot-coarse-counterfactual`.

At the last status check the workflow was in progress and had not yet been reviewed. The next task is to inspect the artifact and compare:
1. final \\(\\tau\\) and final gap to R160;
2. maximum saved-slice absolute \\(\\tau\\) difference;
3. shared-core center-lapse history and t=12 profile errors;
4. geometry differences and finite/positive gate;
5. 12,801 target deliveries, all finite, with zero prescribed-versus-returned target mismatch;
6. report metadata confirming 25 retained target samples with the t=0 and t=12 endpoints preserved.

The full-stride-1 reference is artifact `11727976800` from [run 38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450); its ZIP SHA-256 is `63835669c695f9577d0416a2aa1fb3740aef9d4a90f36aaf9b9dbb94a5783006`. The stride-2 run is the sole new expensive numerical trajectory active for this investigation. Do not launch another trajectory until it completes.

