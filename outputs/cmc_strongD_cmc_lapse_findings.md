# CMC lapse-map and strong-D proper-time investigation

**Status:** diagnostics only; awaiting the full evolving R80 run with the archived R160 global CMC target.  
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

## 6. Next decisive run

[Full R80 trajectory with archived R160 CMC target — run 38206154450](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450) is in progress. It uses the same R80 matter-face and R160 outer-alpha boundary histories, while prescribing the interpolated R160 full-domain (\dot K) history at every CMC lapse solve. This changes no production files and evolves one counterfactual trajectory to determine whether target history alone closes the final τ gap, leaves it unchanged, or triggers a positivity/failure gate.

**Interpretation rule:** the forced-τ result is a gauge-target sensitivity test, not proof that the prescribed target is a self-consistent physical solution on (R=80). Do not change production physics until the full trajectory and admission diagnostics are reviewed.
