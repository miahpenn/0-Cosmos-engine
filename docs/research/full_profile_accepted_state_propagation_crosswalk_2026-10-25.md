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


## Addendum — source-bridge and fixed-R sampling re-evaluation (2026-10-25)

This addendum supersedes any reading of the source-map relative-RMS number or the earlier linear-sampled shell mismatch as a standalone failed conservation identity. The original measurements remain valid for their declared observables, but the later same-state crosswalks reveal what they do and do not mean.

### 1. Brown matter-source bridge closes as signed bookkeeping

The follow-up [Brown source bridge report](brown_matter_source_bridge_absolute_residual_2026-10-25.md), from the same SHA-verified accepted t=0.60 states and source/vendor pins, defines

- (S_{act}): actual geometric matter injection plus current-tangent composite;
- (S_Q=-8\pi\alpha Q_r): provisional stress-divergence source;
- (ΔS=S_{act}-S_Q);
- (B_{tot}): Brown Eq. 11b candidate evaluated with matter-inclusive constraints;
- (B_{geom}): the same Brown construction with geometric-only constraints;
- (D_{tot}=\dot M_{measured}-(B_{tot}+S_Q));
- (D_{geom}=\dot M_{measured}-(B_{geom}+S_{act})).

The signed relation is
[
(B_{tot}-B_{geom})-ΔS=D_{geom}-D_{tot}.
]
On the full-profile arm this reassembly closes to maximum absolute (5.08\times10^{-21}), (7.46\times10^{-21}), and (1.06\times10^{-20}) at N=60/120/240. Separately, the Brown total-minus-geometric source difference matches (ΔS) to max absolute (3.625\times10^{-10}), (3.009\times10^{-11}), (2.015\times10^{-12}). At fixed r=2, both are approximately (5.81\times10^{-5}) across refinement.

**Interpretation:** the nearly resolution-invariant (S_{act}-S_Q) term is almost exactly balanced by the matter-dependent difference between Brown evaluated on total versus geometric-only constraints. Do not interpret (ΔS) in isolation as proof of a missing physical source or a production-equation defect.

### 2. Absolute full-candidate residual converges despite poor relative RMS

The direct candidate (D_{tot}) has full-profile global max absolute residual
[
8.975\times10^{-6},quad9.317\times10^{-7},quad4.133\times10^{-8}
]
at N=60/120/240. At fixed r=2 it is (8.975\times10^{-6},6.122\times10^{-7},4.030\times10^{-8}), observed orders 3.874 and 3.926; at r=8 it is (-1.194\times10^{-6},-7.328\times10^{-8},-4.591\times10^{-9}), orders 4.028 and 3.992. The normalized relative RMS of the earlier candidate remains a conditioning diagnostic and should not be the only convergence criterion when the reference derivative becomes very small.

This is substantial numerical evidence that the full-profile Brown-plus-provisional-(S_Q) candidate has a decreasing absolute remainder on the tested accepted states. It is not scientific admission: Brown Eq. 11b is still a structural scaffold here, and the exact matter-coupled source sign/index/normalization, reference-metric connection convention, regular-origin limit, and total sector conservation law have not all been independently frozen and admitted.

### 3. The fixed-R shell mismatch is also interpolation-sensitive

The separate [Misner–Sharp probe addendum](bssn_areal_misner_sharp_current_probe_2026-10-25.md) records the same accepted-state shell observer after consistent cubic fixed-(R) sampling. At dt/2, maximum absolute shell residual over R=1–10 is:

| N | Cell-centred cubic | Face-reconstructed cubic | Face-reconstructed linear |
|---:|---:|---:|---:|
| 60 | 4.543e-5 | 4.233e-5 | 2.923e-5 |
| 120 | 3.219e-6 | 2.638e-6 | 1.061e-5 |
| 240 | 2.078e-7 | 1.561e-7 | 5.622e-6 |

The face-cubic maximum shows observed orders 4.00 and 4.08; its N=240 face-minus-cell cubic-query discrepancy is only 5.18e-8, while the linear-query discrepancy is much larger. Together with the centered fixed-(R) mass/current result and the exact target-map manufactured control, this makes the previously large linear-sampled face/cell gap **not sampling-independent**. It does not select cubic interpolation as production authority.

### 4. Current disposition and next gate

The scalar spatial-force crosswalk separately finds an approximately fourth-order annular gap between the pinned expanded D1/D2 RHS form and the cell-centred divergence-form control. Taken together, the scalar spatial operator, fixed-(R) mass/current, and cubic shell witnesses are consistent with finite-grid interpolation/operator remainders in the tested annulus. None alone proves the full continuum law.

The next work should therefore stop repeating interpolation variants and return to the archive-backed full matter-coupled propagation identity: explicitly freeze the covariant versus contravariant (M_r/M^r) mapping, normalized Hamiltonian-gradient term, reference-metric connection derivative/coupling, matter-source contribution, and sector-exchange conventions on the exact accepted states. Keep the origin and fixed-radius annulus as separate gates.

No production equations, initializers, vendor source, dependencies, source pins, or accepted states were changed. `full_matter_coupled_propagation_identity_closed=false`, `total_stress_energy_sector_conservation_closed=false`, and `scientific_admission=false`.
