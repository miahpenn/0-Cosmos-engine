# Accepted-state propagation candidate crosswalk — 2026-10-25

**Status:** diagnostic only. This is a measured comparison against existing, explicitly provisional propagation candidates; it is not a frozen or scientifically accepted matter-coupled propagation identity.

**Execution:** [workflow run 40156795380](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40156795380), all N=60/120/240 jobs succeeded. **PR:** [#420](https://github.com/miahpenn/0-Cosmos-engine/pull/420), which is the CI execution head for [#419](https://github.com/miahpenn/0-Cosmos-engine/pull/419). The tested observer commit is `affb42f7dd964838cfc402417f933bf17a9d039c`.

## Purpose

Use the same saved, source-pinned accepted state at (t=0.60) for each comparison, rather than comparing initial data to evolved data or regenerating a trajectory. On each frozen accepted state, reuse the established signed momentum tangent and the existing Brown Eq. (11b)/stress-divergence candidate evaluator. Report the centre cells separately from fixed-radius and annular regions.

Three arms are evaluated:

- `production_baseline`: archived production-initializer evolution and production CMC lapse treatment;
- `partial_joint`: the PR #416 partial first-half-cell correction plus copied regular-origin CMC row;
- `full_continuous_profile`: PR #418 full continuous-source radial metric-profile shadow plus copied regular-origin CMC row.

The latter two use the matching copied regular-origin CMC solve for each directional perturbation. The measurement uses two symmetric directional epsilons ((2\times10^{-5}) and (10^{-5})); their outputs agree closely. No evolution steps are advanced by this crosswalk.

## Results

### 1. The centre-specific geometric improvement survives the propagation comparison

The direct covariant-(M_r) candidate residual at the first two cell-centre samples is:

| Resolution | Arm | cell 0 | cell 1 |
|---:|---|---:|---:|
| 60 | production baseline | (-4.6040\times10^{-4}) | (-5.0878\times10^{-4}) |
| 60 | partial joint | (+1.1269\times10^{-4}) | (+9.8363\times10^{-5}) |
| 60 | full profile | (+8.7893\times10^{-8}) | (+8.9752\times10^{-6}) |
| 120 | production baseline | (-7.5501\times10^{-3}) | (-9.8622\times10^{-4}) |
| 120 | partial joint | (+8.2874\times10^{-4}) | (+1.4131\times10^{-4}) |
| 120 | full profile | (-9.3173\times10^{-7}) | (+4.3685\times10^{-7}) |
| 240 | production baseline | (+3.2485\times10^{-2}) | (-4.0490\times10^{-3}) |
| 240 | partial joint | (-3.9274\times10^{-3}) | (+4.0300\times10^{-4}) |
| 240 | full profile | (+3.0083\times10^{-9}) | (+6.6256\times10^{-9}) |

The full-profile arm drives this **candidate's central residual** markedly lower, especially at N=240. This is consistent with the improved (A_a/r^2) near-centre consistency reported in PR #418. It is a local, candidate-dependent finding—not a proof of the complete propagation identity.

### 2. The complete candidate does not close globally

Relative RMS mismatch of the direct covariant-(M_r) candidate:

| Resolution | Arm | Full domain | (r\ge5) |
|---:|---|---:|---:|
| 60 | production baseline | 0.7182 | 0.9193 |
| 60 | partial joint | 0.9068 | 1.3752 |
| 60 | full profile | 0.7970 | 0.8338 |
| 120 | production baseline | 0.7814 | 0.3393 |
| 120 | partial joint | 0.7407 | 0.3615 |
| 120 | full profile | 0.9565 | 0.8419 |
| 240 | production baseline | 1.9567 | 0.1004 |
| 240 | partial joint | 1.1106 | 0.1004 |
| 240 | full profile | 0.6566 | 0.8178 |

The full-profile candidate is not globally convergent to zero. At N=120 and N=240 its (r\ge5) relative mismatch is worse than the production baseline. This separates the near-centre geometric improvement from the still-unclosed propagation law and argues against using this candidate as an admission oracle.

### 3. Matter-source mapping remains essentially unchanged across the arms

The current evaluator compares the actual matter contribution to the total momentum-constraint tangent,
[
\dot M_{\rm geom}^{\rm matter}-8\pi\dot J,
]
with the provisional stress-divergence source,
[
-8\pi\alpha Q_r.
]

The source-mapping relative RMS residual is:

| Resolution | Production baseline | Partial joint | Full profile |
|---:|---:|---:|---:|
| 60 | 0.999108 | 0.999091 | 0.999091 |
| 120 | 0.999940 | 0.999940 | 0.999939 |
| 240 | 0.999996 | 0.999996 | 0.999996 |

At cell 0/cell 1, the full-profile source-mapping residuals are (+2.2331\times10^{-5},+5.8090\times10^{-5}) (N=60), (+1.1415\times10^{-5},+3.3032\times10^{-5}) (N=120), and (+5.7316\times10^{-6},+1.7037\times10^{-5}) (N=240). The corresponding independently reconstructed (Q_r) values are substantially smaller—around (10^{-9}), (10^{-10}), and (10^{-11}), respectively, at the centre. The mismatch is nearly the full scale of the actual matter-source composite.

Crucially, this source-mapping residual is almost identical among the baseline, partial and full-profile arms at each resolution. The metric-profile change therefore does **not** resolve this source-mapping gap. The source sign/index/normalization and scalar/sector exchange budget remain provisional and incomplete; this is not a verdict that the production equation is wrong.

### 4. The algebraic chain rule is checked, but does not certify physics

The exact discrete derivative of the existing (M_{\rm geom}) definition agrees with its measured tangent to approximately (10^{-13}) absolute or better across the saved states and epsilon choices. This verifies the observer's implementation bookkeeping. It is not a continuum constraint-propagation identity and cannot replace the independent source-aware oracle.

## Provenance and scope

- Production/source anchor: `4aec88d0e3ff028e0a65a8ac36b30948625df792`.
- Vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Each input artifact manifest was checked inside the workflow; the protected-source diff was empty; all frozen-state observers preserved the input digest; all output profiles were finite.
- Directional epsilon choices: (2\times10^{-5}) and (10^{-5}); their reported candidate metrics agree closely.
- The candidate keeps Brown Eq. (11b) as a structural comparison and uses provisional stress-divergence/source mapping. Its source index and normalization have **not** been admitted as the complete matter-coupled propagation law.
- No production initializer/equation, engine, vendor or dependency file is changed by this work.
- `full_matter_coupled_propagation_identity_closed=false`
- `matter_source_sign_and_normalization_reviewed=false`
- `reference_metric_connection_coupling_convention_frozen=false`
- `total_stress_energy_sector_conservation_closed=false`
- `scientific_admission=false`

## Decision

Keep the full continuous-profile arm as a useful **geometric diagnostic/control**, and the partial half-cell arm as a negative control—not a proposed production correction. The next scientific task is not another metric-profile sweep: it is to resolve the frozen source mapping/index/normalization and complete the total sector conservation/exchange identity using the existing archive-backed convention map and exact accepted-state ledger. No production repair is justified from this candidate comparison alone.
