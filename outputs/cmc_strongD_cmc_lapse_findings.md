# CMC lapse-map and strong-D proper-time investigation

**Status:** diagnostic-only full target-forced trajectory completed; coarsened target-history sensitivity run in progress.  
**Branches:** weak-field Jacobian study on `diag/bond-sweep-self-consistent`; strong-D studies on `diag/independent-metric-ricci-audit`. Production sources were not edited; `main` was not changed by this work (last observed head `7604bd6cc87c48c9dc8e24f2facbfe1458364ca5`).

## 1. Weak-field fixed-point study

- Fixed-target / Picard test: [run 38204727992](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38204727992)
- Lapse-map Jacobian and input-slice residual: [run 38205208414](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205208414)
- Script: [cmc_lapse_jacobian_study.py](../blob/85b32b8b6b9abc0a46e563fbc4639b4a76c4291d/diagnostics/cmc_lapse_jacobian_study.py)

The fixed-target solve removes the repeated weak-field offset, but ordinary repeated substitution did not contract. Finite differences on four radial lapse modes found a small map response rather than an identity-like translation. The sampled (J_F=J_G-I) had no near-zero singular value in that tested subspace. This is local evidence only, not proof of a global lapse root.

## 2. Strong-D single-slice audit

- Archived snapshot: [boundary-matched R80 run 38162318255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255), artifact `11710858125`, source commit `35072830f9e4d7029bb10243a0520b9a273a17dd`.
- Strong-D Jacobian/residual study: [run 38205542397](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205542397).
- Script: [cmc_strongD_snapshot_jacobian.py](../blob/2ad1593dec8e7b5b62aa411aa35249b1c1acd3f5/diagnostics/cmc_strongD_snapshot_jacobian.py).

At (t=12), (N=160), (r_{max}=80), (D(0)=0.0686763):
- The input-slice discrete CMC residual maximum was (2.3133\times10^{-10}) (relative max (2.2690\times10^{-10})).
- The lapse re-solve correction was (9.4410\times10^{-8}) RMS; at the center it was only (4.0921\times10^{-9}).
- The same-target elliptic solve closed to (4.04\times10^{-15}) relative max residual.
- On the four tested radial directions, (J_G)'s largest eigenvalue was about (0.1012), and (J_F)'s smallest singular value was (0.7321). No near-degenerate fixed-point residual direction was found in that subspace.

The immediate lapse re-solve mismatch at (t=12) is not large enough at the center to explain the accumulated proper-time discrepancy by direct integration.

## 3. Full saved-time lapse correction series

- Run: [38205709871](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205709871)
- Artifact: `11727266816`, `cmc-strongD-snapshot-series`.
- Script: [cmc_strongD_snapshot_series_audit.py](../blob/b1776c1f5af08eae1e00a89a979f8fce4f7254fc/diagnostics/cmc_strongD_snapshot_series_audit.py).

Across 49 snapshots from (t=0) to (12), the trapezoidal integral of the frozen-slice central lapse re-solve correction was (2.6818\times10^{-8}) signed, (2.6860\times10^{-8}) absolute. This is (2.53\times10^{-6}) of the matched R80–R160 proper-time discrepancy (0.01062109). Every same-target solve closed to relative max below (5\times10^{-15}). This is a scale test on frozen states, not a re-evolution.

## 4. Direct R80–R160 central lapse history

The R80 matched-boundary and R160 reference captures share 49 saved times. The integral of the actual stored central lapse difference is (0.01063961), compared with the observed final accumulated τ difference (0.01062109); the saved-time quadrature differs by (1.8512\times10^{-5}) (about 0.18%). The largest center-lapse difference is (0.00474876) near (t=6.49875), not at the final slice.

This points to divergence of the evolving lapse histories as the proximate numerical source of the accumulated τ difference—not the tiny frozen-slice re-solve correction.

## 5. Global CMC target counterfactual on frozen states

- Run: [38205934746](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205934746)
- Artifact: `11727417247`, `cmc-domain-target-counterfactual`.
- Script: [cmc_domain_target_counterfactual.py](../blob/8b2bafa5ec192e31fb70906ff99e9c14f4696930/diagnostics/cmc_domain_target_counterfactual.py).

R80's native target (\dot K) and R160's full-domain target are nearly equal initially but diverge strongly by (t=12):
- At (t=0): R80 (6.19884\times10^{-5}), R160 (6.23352\times10^{-5}). Transferring only the R160 target reproduces the initial R160 lapse on the shared core to about (10^{-12}).
- At (t=12): R80 (2.45708\times10^{-4}), R160 (7.50581\times10^{-5}).
- On frozen slices the target swap substantially reduces core radial lapse-profile RMS error at later times, but the center correction overshoots: its signed integral relative to R160 is (-0.03777), not the observed (+0.01062) proper-time difference.

Therefore the global target contributes to lapse differences, but transplanting the target onto frozen R80 states does not reproduce the center-lapse history. Its feedback through the evolved geometry/matter must be tested.

## 6. Full evolving target-history counterfactual — completed

[Run 38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450) completed successfully. Artifact [11727976800](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450/artifacts/11727976800), ZIP SHA-256 `63835669c695f9577d0416a2aa1fb3740aef9d4a90f36aaf9b9dbb94a5783006`.

The run reached t=12 using 3,200 accepted steps and 49 saved records; failure was null. All 12,801 lapse solves received finite prescribed targets, the maximum prescribed-versus-returned target difference was exactly zero, and the positive-finite geometry gate passed at each accepted step. Production equations/source were not edited.

Against the earlier boundary-matched R80 trajectory using its native R80 target (run [38162318255](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255)), final τ changed from `6.917400377028174` to `6.9067799130064635`; the final gap versus R160 dropped from `+0.010621093774165757` to `+6.297524555165523e-7`, closing `99.9940707%` of the final gap. Maximum saved-slice absolute τ gap declined from `0.0121899232135` to `0.0000121914168` (about `99.90%`). Shared-core central-lapse differences also collapsed: t=12 maximum from about `0.05218827` to `2.0914e-6`, with all-saved-slice lapse RMS dropping from about `0.03028393` to `7.8294e-7`.

**Interpretation:** the full-evolution intervention provides strong causal evidence that the archived R160 global CMC target history drives the R80/R160 lapse and proper-time divergence in this controlled boundary-matched setup. It does not establish that the R160 target is uniquely physically correct or self-consistent on R80, because the R160 target and outer-alpha history and R80 matter-face history are imported archived inputs. This remains diagnostic-only.

[Coarsened target-history sensitivity run 38207554072](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207554072) remains in progress at the last check. It retains every other archived R160 target sample, preserves t=12, and interpolates across approximately 0.5-time-unit intervals. A supplementary offline comparison found that the stride-2 and stride-1 target inputs differ by at most 5.29437e-8 over 3,201 accepted-step times (RMS 1.05494e-8); the maximum occurs near t=6.25125. This is a forcing-input difference, not a bound on the nonlinear trajectory response. The full-resolution reference is run 38206154450/artifact 11727976800. Compare tau, maximum saved-slice tau gap, shared-core lapse differences, geometry differences, target delivery and admission gates. This is the only active expensive numerical trajectory; do not modify production physics or launch another trajectory until it is complete.


## 7. 12 October continuation: interpolation audit and full target-forced trajectory complete

Detailed dated handoff: [CMC strong-D lapse / global-target handoff](../blob/78720e52f61ac5f0aa1e0c6478a4d624c8ab9db7/outputs/cmc_strongD_cmc_lapse_handoff_2026-10-12.md).

The first Kdot interpolation audit attempt, [run 38206744317](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206744317), failed before analysis with `ModuleNotFoundError: No module named 'numpy'`. This was a workflow setup failure, not a numerical finding. Dependency correction: commit [d6142ac81a4bf56905de0fae5298f8fc540f7d98](https://github.com/miahpenn/0-Cosmos-engine/commit/d6142ac81a4bf56905de0fae5298f8fc540f7d98) adds explicit NumPy installation. The corrected run [38206848868](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868) succeeded and uploaded artifact [11727926565](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206848868/artifacts/11727926565), ZIP SHA-256 `993340f756973e2a0251da03269ea37d38dd1591a0cffa45e9540148c471c39b`.

The held-out interpolation test uses 49 archived target samples over t=0–12, holds out 24 points, and predicts each from adjacent retained values over a doubled interval (about 0.499–0.503). It measured maximum absolute error (5.29437\times10^{-8}), RMS absolute error (1.82802\times10^{-8}), p95 absolute error (4.82885\times10^{-8}), maximum relative error 0.0883%, and RMS relative error 0.0304%. The largest absolute error was near t=6.25125. The target series is smooth at resolved sample scale, but this remains an empirical check—not a formal bound at every internal evolution stage and not a direct bound on accumulated proper-time response.

The full evolving R80 run with the prescribed archived R160 target has completed: [run 38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450), artifact [11727976800](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450/artifacts/11727976800). It closes 99.9941% of the final τ gap and reduces the maximum saved-slice τ gap by about 99.90% relative to the native-target boundary-matched R80 trajectory. Hashes and admission details are in the dated handoff.

At this checkpoint, the recorded conclusion remains:
- The 0.18% integral agreement is a proper-time/lapse bookkeeping check, not attribution.
- Frozen-state lapse re-solve correction is only 2.53e-6 of the observed tau gap.
- Frozen-state target swap changes radial lapse profiles but is insufficient by itself to predict central proper-time history.
- The full evolving target-forced counterfactual closes 99.9941% of the final τ gap, strongly implicating the global target history in this controlled setup.
- Coarsened target-history sensitivity is the next decisive test; the interpolation audit measures sample-scale error but does not certify every evolution-stage interpolation.
- `main` was verified at `7604bd6cc87c48c9dc8e24f2facbfe1458364ca5`; production equations remain unchanged by this diagnostic work.
