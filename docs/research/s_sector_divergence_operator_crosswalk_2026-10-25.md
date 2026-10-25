# Same-state S-sector RHS-to-divergence crosswalk — 2026-10-25

**Status:** research-only diagnostic. No evolution was rerun; no production, initializer, vendor, dependency, or physics convention was changed. Scientific admission remains false.

## Provenance

- Accepted full-continuous-profile states at (t=0.60): [workflow run 40156226086](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40156226086), resolutions (N=60,120,240).
- Sector-resolved observer and signed profiles: [workflow run 40207821733](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40207821733), all three resolution jobs passed.
- Frozen production source: `4aec88d0e3ff028e0a65a8ac36b30948625df792`.
- Pinned vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- The accepted-state and diagnostic `report.json` / `profiles.npz` SHA-256 manifests for all six artifacts were checked before analysis.

This note responds to Lumen's request to compare the S-only pointwise and finite-volume (FV) stress divergence against the scalar RHS's spatial Laplacian contribution, including the origin measure.

## 1. What is being compared

At the frozen source pin, the production S-sector RHS uses the cell-centred derivative form

```text
L_S[S]=(X^2/a)[D2(S)+
(D1(b)/b-D1(a)/(2a)+2/r-D1(X)/X)D1(S)].
```

Source inspection confirms that this is the accepted update's actual spatial form: `engine/scalar_system.py` evaluates it with parity-aware cell-centred D1/D2 derivatives, and `engine/production_kernel.py` evaluates that RHS in the evolution stages. The S update is not implemented as an HLL flux or a telescoping face-flux difference. Therefore, the expanded derivative-form operator is authoritative for describing the accepted semidiscrete S update; the FV stress-divergence calculation remains an independent diagnostic and must not be substituted as though it were the update itself. See [scalar RHS](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/engine/scalar_system.py#L53-L93) and [production-stage call site](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/engine/production_kernel.py#L286-L302).

The requested isolated spatial contribution is `L_S[S] * D1(S)`. In the scalar momentum RHS, the corresponding Laplacian term is `alpha * L_S[S]`; a separate lapse-gradient term ( (X^2/a)D_1alpha D1(S)), curvature, potential, shift, and (for the coupled φ field) conditional exchange terms are separate RHS operands.

The sector stress-divergence observer instead evaluates the full (Q_r^S=nabla_mu(T_S^mu_r)) from the time-density/current derivative, radial pressure flux, and metric-source terms. Its FV form reconstructs combined face fluxes and integrates sources and the measure per cell. Therefore, **`Q_r^S` is not the same quantity as the isolated S-sector Laplacian contribution `L_S[S] * D1(S)`**. Their raw ratio is a useful crosswalk witness, but a mismatch does not by itself indicate a missing source or a bad scalar RHS.

## 2. The pinned RHS term reconstruction closes

The two signed current-tangent pieces from the scalar RHS,
`scalar_S_current_rate_from_Pdot_lapse_laplacian_second_derivative` and
`scalar_S_current_rate_from_Pdot_lapse_metric_gradient`, were combined and divided by `8*pi*alpha`. Their mapped result matches `L_S[S] * D1(S)` to floating-point roundoff:

| Resolution | Maximum absolute mapping gap | RMS mapping gap |
|---:|---:|---:|
| N=60 | 1.06e-22 | 2.11e-23 |
| N=120 | 7.94e-23 | 1.26e-23 |
| N=240 | 1.06e-22 | 9.80e-24 |

This confirms that the diagnostic reconstructs the pinned RHS spatial contribution correctly. It is an implementation/source-budget identity, not an independent proof of the full scalar stress-energy conservation law.

## 3. Origin-cell measure

The FV covariant-divergence observer uses the spacetime measure
```text
sqrt(-g)=alpha * sqrt(gamma),
```
not the spatial measure (sqrt(gamma)) alone. It reconstructs `sqrt(-g)` to faces with the regular `r^2 * g(r^2)` form and computes each measure by Simpson integration,
```text
V_0=dr/6[w(0)+4w(dr/2)+w(dr)],  w=sqrt(-g), w(0)=0.
```
It does **not** use the centre value times Δr. The familiar `(4/3) * dr * w(dr/2)` is an excellent near-origin approximation, not the exact Simpson value.

| Resolution | Δr | Observer (V_0=integralalpha * sqrt(gamma) dr) (Simpson) | Centre value × Δr | 4/3 centre approximation | Relative Simpson / approximation difference |
|---:|---:|---:|---:|---:|---:|
| 60 | 1.333333 | 0.7974845181 | 0.5980217458 | 0.7973623278 | 1.53e-4 |
| 120 | 0.666667 | 0.09966964221 | 0.07474927498 | 0.09966569998 | 3.96e-5 |
| 240 | 0.333333 | 0.01245818990 | 0.00934354931 | 0.01245806575 | 9.97e-6 |

For clarity, the same Simpson reconstruction using the *spatial* measure (sqrt(gamma)) alone gives (0.7984648512), (0.09979184239), and (0.01247345586) for (N=60,120,240). Those are not the measure values used by this covariant spacetime-divergence observer.

## 4. Direct S-only comparison requested by Lumen

Cell 0 is **not** a fixed-radius sequence: its coordinate radius is `r=dr/2`: 0.6667, 0.3333, 0.1667 for the three resolutions. Do not infer a central convergence order from these three cell-0 samples.

| Resolution | `L_S[S] * D1(S)`, cell 0 | Pointwise `Q_r^S`, cell 0 | Pointwise / isolated RHS term | FV `Q_r^S`, cell 0 | FV / isolated RHS term |
|---:|---:|---:|---:|---:|---:|
| N=60 | 1.98443e-7 | 7.92656e-7 | 3.9944 | 6.72404e-9 | 0.033884 |
| N=120 | 1.01414e-7 | 1.13000e-7 | 1.1143 | 9.04460e-10 | 0.008919 |
| N=240 | 5.09515e-8 | 1.45989e-8 | 0.28653 | 1.13385e-10 | 0.002225 |

The cell-0 pointwise result is itself the signed sum of large cancelling terms, not the isolated RHS Laplacian. The recorded components are:

| Resolution | Time-density/current | Radial pressure-flux derivative | Metric source | Sum: pointwise Q_r^S |
|---:|---:|---:|---:|---:|
| N=60 | -7.45445e-7 | -3.49372e-5 | +3.64753e-5 | +7.92656e-7 |
| N=120 | -3.78191e-7 | -7.32343e-5 | +7.37255e-5 | +1.13000e-7 |
| N=240 | -1.89775e-7 | -1.47658e-4 | +1.47862e-4 | +1.45989e-8 |

On the fixed-coordinate annulus `2 <= r <= 8`, the RMS of the isolated spatial contribution stays close to 3.9e-7, whereas the full S-sector divergence decreases strongly:

| Resolution | RMS `L_S[S] * D1(S)` | Pointwise `Q_r^S` RMS | Pointwise ratio RMS | FV `Q_r^S` RMS | FV ratio RMS |
|---:|---:|---:|---:|---:|---:|
| N=60 | 3.96049e-7 | 8.20593e-8 | 0.16921 | 6.45932e-9 | 0.037790 |
| N=120 | 3.90021e-7 | 3.30293e-9 | 0.006548 | 4.30698e-10 | 0.002724 |
| N=240 | 3.89688e-7 | 2.17498e-10 | 0.0004336 | 2.71842e-11 | 0.0001752 |

Fixed-coordinate samples make the distinction especially clear:

| (r) | (N) | Isolated RHS term `L_S[S] * D1(S)` | Pointwise `Q_r^S` | FV `Q_r^S` |
|---:|---:|---:|---:|---:|
| 2 | 60 | 4.91731e-7 | 1.79071e-7 | 4.51275e-9 |
| 2 | 120 | 4.85866e-7 | 1.29268e-8 | 2.82897e-10 |
| 2 | 240 | 4.91780e-7 | 7.98422e-10 | 1.71462e-11 |
| 5 | 60 | 3.73132e-7 | -8.07318e-9 | 7.90633e-9 |
| 5 | 120 | 3.74392e-7 | -7.69929e-10 | 5.38652e-10 |
| 5 | 240 | 3.74244e-7 | -4.91024e-11 | 3.40465e-11 |
| 8 | 60 | 3.66119e-8 | -1.94645e-9 | 3.88405e-9 |
| 8 | 120 | 2.96932e-8 | -8.63048e-11 | 2.40370e-10 |
| 8 | 240 | 2.79714e-8 | -4.70553e-12 | 1.49531e-11 |

The ratios are not a constant factor across the annulus or resolutions. At fixed `r=2`, both pointwise and FV `Q_r^S` decrease by about fourth order under refinement, consistent with the previously recorded stress-divergence crosswalk. The evidence points away from a simple `4*pi/8*pi`-type normalization error. It does **not** authorize choosing FV over pointwise solely because the FV residual is smaller.

## 5. Why sector-summed divergence differs from the expanded total candidate

A follow-up source trace resolves the apparent S-sector/total mismatch. Two separate diagnostic observers are using continuum-equivalent, but discretely different, realizations of the same covariant radial divergence:

- The sector-resolved observer in [audit_total_stress_energy_divergence.py](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/tools/audit_total_stress_energy_divergence.py) takes the radial-flux derivative as the cell-centred derivative of the *combined* flux, (D_1(\alpha\sqrt{\gamma}\,p_r)/(\alpha\sqrt{\gamma})), then adds the metric-source terms. Each sector uses the same operation.
- The expanded candidate in [audit_hybrid_origin_integral_patch_newton.py](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/tools/audit_hybrid_origin_integral_patch_newton.py) evaluates the continuum-expanded expression, (D_1(p_r)+(D_1\alpha/\alpha)(\rho+p_r)+2(D_1R/R)(p_r-p_t)), alongside the temporal tangent and metric-time terms.

In the continuum these forms are equivalent. On the frozen grid, (D_1) does not obey an exact product rule, so the combined-flux total and expanded total candidate need not match pointwise. This is not an S-only sector-additivity defect: summing the sector-resolved combined-flux divergences reproduces the combined-flux total to roundoff-scale maximum differences (about (1.4\times10^{-17}), (3.8\times10^{-17}), (3.0\times10^{-17}) for N=60/120/240). The signed sector source-composite/residual decomposition is additive to below (8.0\times10^{-16}) in the same runs. The report's strict raw-sum boolean remains false, so these roundoff-scale comparisons are evidence about the observed arrays—not an admission-gate pass.

For the full-continuous-profile arm at the same accepted (t=0.60) states (directional tangent epsilon (2\times10^{-5})), the combined-flux total minus expanded total candidate is:

| Witness | N=60 | N=120 | N=240 |
|---|---:|---:|---:|
| max absolute gap, (r\le3) | (7.95344\times10^{-7}) | (1.12971\times10^{-7}) | (1.45822\times10^{-8}) |
| RMS gap, (2\le r\le10) | (7.89360\times10^{-8}) | (3.10751\times10^{-9}) | (2.03876\times10^{-10}) |
| full-domain RMS gap | (1.05669\times10^{-7}) | (1.09514\times10^{-8}) | (1.02051\times10^{-9}) |
| signed gap at fixed r=2 | (+1.87233\times10^{-7}) | (+1.34634\times10^{-8}) | (+8.33415\times10^{-10}) |
| signed gap at fixed r=5 | (-9.26693\times10^{-9}) | (-8.87058\times10^{-10}) | (-5.66681\times10^{-11}) |
| signed gap at fixed r=8 | (-5.90229\times10^{-9}) | (-3.49240\times10^{-10}) | (-2.13651\times10^{-11}) |
| algebraic gap-decomposition closure max abs | (2.33\times10^{-21}) | (2.13\times10^{-20}) | (2.41\times10^{-20}) |

The recorded decomposition assigns essentially all of the difference to the radial-flux product-rule gap; time-product and lapse/metric-chain gaps are much smaller. At fixed radii (r=2,5,8), the absolute gap decreases at approximately fourth order. Do not infer a fixed-radius order from cell 0: its coordinate radius changes with resolution.

This explains why a sector-summed (Q_r) of order (10^{-7}) on the coarse central cell can coexist with a much smaller expanded total (Q_r): they are not the same finite-difference operator. It does **not** decide which continuum-equivalent stress-divergence discretization is the accepted conservation diagnostic. The accepted S-sector update is still the expanded cell-centred D1/D2 RHS documented above; the combined-flux observer remains a separate conservation diagnostic. The full matter-coupled propagation and physical sector-conservation gates remain open.

## 6. Same-state scalar RHS/current/measure pairing — 2026-10-25

**Exact-head run:** [40212458916](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40212458916), all N=60/120/240 jobs passed. This was post-processing of SHA-verified accepted states at t=0.60; the accepted states were not advanced or mutated. Production source remains 4aec88d0e3ff028e0a65a8ac36b30948625df792, vendor d6052d605673ce9d82cdc99b0df79855fa2ea215, and the protected engine/vendor/dependency diff was empty.

The audit now records the sector-specific RHS rates, the current tangent, and the spatial-measure rate:
d_t ln(sqrt(gamma)) = 0.5*a_dot/a + b_dot/b - 3*X_dot/X.
The current tangent reconstructed from the pinned RHS matches the directional witness, and (j_dot + j*d_t ln(sqrt(gamma)))/alpha matches the observer's time-density contribution. The measured maximum absolute gaps are:

| N | Production-expanded EOM residual max abs | Current-tangent RHS vs directional witness | Reconstructed time-density vs observer |
|---:|---:|---:|---:|
| 60 | 8.67e-19 | 4.47e-17 | 7.37e-17 |
| 120 | 3.47e-18 | 3.51e-16 | 3.55e-16 |
| 240 | 1.73e-18 | 7.40e-16 | 7.13e-16 |

The production-expanded EOM residual is expected to close because it uses the same pinned Pdot RHS; it is an implementation reconstruction check, not independent physics validation. The separate divergence-form EOM residual differs from the expanded production form by the finite-grid spatial product-rule/operator gap.

The continuum identity is nabla_mu(T_S^mu_r) = EOM_S * D_r(S). On the grid, the pointwise stress divergence, volume-averaged FV divergence, and EOM-gradient pairing are deliberately kept as distinct discrete objects. Their observed gaps are:

| N | Pointwise pairing-gap RMS (full domain) | Pointwise max abs | FV pairing-gap RMS (full domain) | FV max abs |
|---:|---:|---:|---:|---:|
| 60 | 1.010e-7 | 7.610e-7 | 3.903e-9 | 2.491e-8 |
| 120 | 1.050e-8 | 1.084e-7 | 3.761e-10 | 3.706e-9 |
| 240 | 9.793e-10 | 1.400e-8 | 3.427e-11 | 4.848e-10 |

Do not assign a formal convergence order to the full-domain maxima: the innermost cell represents a different physical radius at each resolution. At fixed coordinate radius r=2, the signed pointwise gap decreases from +1.76049e-7 to +1.26966e-8 to +7.85394e-10, with observed orders 3.79 and 4.01. The FV gap at the same radius decreases from +1.49114e-9 to +5.27063e-11 to +4.11783e-12. This supports a refinement-decreasing discrete pairing remainder; it does not pick a preferred conservation operator.

The shell-integrated comparison uses the FV observer's Simpson-reconstructed alpha*sqrt(gamma) cell measures and cell-centre radius masks. The domain r<10 integrated gap is:

| N | Pointwise divergence minus EOM-pairing integral | FV divergence minus EOM-pairing integral |
|---:|---:|---:|
| 60 | +1.26115e-6 | +1.45554e-6 |
| 120 | +7.96113e-8 | +9.47507e-8 |
| 240 | +4.86911e-9 | +5.94881e-9 |

Both shell-integrated differences decrease at approximately fourth order over this refinement sequence. These are observer-consistent, volume-weighted crosswalk sums; because the pointwise/EOM operands are cell-centre samples and the FV operand is already volume-averaged, they are not an exact finite-volume conservation identity or a scientific-admission test. The algebraic split of the pointwise pairing gap into its time-density and spatial/metric parts closes to max absolute 2.65e-21, 1.23e-20, and 1.01e-20 at N=60/120/240.

**Artifact provenance:** workflow 40212458916 uploaded N60 artifact 12605564698, N120 artifact 12606088948, and N240 artifact 12605739292. SHA-256 (report.json / profiles.npz): N60 d22faa6fcd80348018054e88bda6d06485482560814f78256c3c98ac2edad19d / d321016747d6d5606458270dfaa0013d8d69b5462fffe470d75401c2992311e4; N120 e6cd0f788802d9a0623eae4d0a5e795bc1d13d27abd519286e7e252ceab713ce / b5b21465a71ef81cb1b7548f9a60039c3fe9bb82c46e36585be3843c5e548096; N240 1545eace2555c8a8ffd95860ec4a98c39676adcf3280674a77706a5fa14355d5 / e17520b8e6169eeccc423f6fa9f9edfbcfcc52541fb994f3ca40bada8508a020.

## 7. Archive-locked DM/COSMOS-phi exchange comparator — same accepted states

The same saved-state sector ledger also compares the radial covariant divergences with the archive-locked exchange comparator (Q_r^{DM}=+\beta\rho_{DM}D_r\phi), (Q_r^{\phi}=-\beta\rho_{DM}D_r\phi), using the normalized DM projection from the same stress observer. This is a **diagnostic comparison**, not acceptance of the exchange law or its normalization.

At fixed coordinate radius (r=2), the DM-sector observed divergence is close to its comparator and the signed difference shrinks across the saved-state refinement sequence. The COSMOS-φ observed divergence is much larger than the comparator signal at these resolutions, although its signed difference also decreases rapidly. The two central points are:

| N | DM expected (Q_r) | DM observed (Q_r) | DM observed minus expected | COSMOS-φ expected (Q_r) | COSMOS-φ observed (Q_r) | COSMOS-φ observed minus expected |
|---:|---:|---:|---:|---:|---:|---:|
| 60 | +3.86998e-14 | +4.31002e-14 | +4.40038e-15 | -3.86998e-14 | -2.42207e-10 | -2.42169e-10 |
| 120 | +3.75918e-14 | +3.86402e-14 | +1.04847e-15 | -3.75918e-14 | -1.76788e-11 | -1.76412e-11 |
| 240 | +3.75661e-14 | +3.76963e-14 | +1.30167e-16 | -3.75661e-14 | -1.12531e-12 | -1.08775e-12 |

For COSMOS-φ, the fixed-r=2 discrepancy drops by factors of about 13.7 then 16.2, approximately fourth-order refinement. Nevertheless, at N=240 it remains about 29 times larger than the expected source itself. The N=60/120/240 values are stable under halving the directional-tangent epsilon from (2\times10^{-5}) to (1\times10^{-5}), so this does not look like a simple finite-difference epsilon artifact. It still does **not** distinguish the small archived exchange signal from the spatial product-rule/operator remainder, nor does it prove that either normalization is wrong.

The DM comparison is encouraging numerical consistency but is too narrow to close the gate: it is one fixed-radius diagnostic on accepted background states. The COSMOS-φ comparison is unresolved at these resolutions because the comparator signal is below the observed discrete remainder. This is consistent with the standing Issue #210 guardrail: no reciprocal source, factor, or density conversion should be selected from a numerically convenient match. The archive-backed exchange normalization and the complete total-sector conservation ledger remain open.

## 8. Per-sector expanded divergence resolves the apparent COSMOS-phi mismatch — 2026-10-25

**Corrected exact-head workflow:** [40213110427](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40213110427), all N=60/120/240 jobs passed. This follow-up builds the continuum-expanded radial-divergence expression for each stress sector from the same saved-state directional RHS endpoints, per-sector projected stress, and geometry-rate arrays. The archive-locked expected-source arrays are joined explicitly from the sector-source ledger; no missing comparator is replaced by zero.

This **refines the interim interpretation in Section 7**. The combined-flux pointwise COSMOS-phi divergence is far from the small archived exchange signal, but the sector-expanded divergence matches that comparator closely. The difference is the discrete product-rule/operator-pairing gap, not evidence on these data of a missing physical source.

At fixed coordinate radius r=2:

| N | COSMOS-phi: pointwise minus expanded Qr | COSMOS-phi: expanded minus expected exchange | DM: expanded minus expected exchange | S: expanded minus zero-source comparator | Expanded-sector sum minus total expanded candidate, max abs |
|---:|---:|---:|---:|---:|---:|
| 60 | -2.42168e-10 | -6.07801e-16 | -8.90094e-15 | -8.40417e-9 | 6.14e-18 |
| 120 | -1.76411e-11 | -6.80868e-17 | +1.15190e-17 | -5.54197e-10 | 2.16e-17 |
| 240 | -1.08770e-12 | -4.72254e-17 | +6.77129e-17 | -3.60805e-11 | 1.66e-17 |

The phi pointwise-minus-expanded gap decreases by factors of approximately 13.7 and 16.2 on successive resolution doublings (about fourth order). At N=240, the expanded phi residual is \(4.72\times10^{-17}\) versus an expected source of about \(-3.76\times10^{-14}\); the N=60/120/240 expanded residuals are about 1.6%, 0.18%, and 0.13% of the expected-source magnitude. The direct projected-current tangent also matches the constraint-derived sector tangent to max absolute \(3.23\times10^{-27}\), \(1.62\times10^{-27}\), and \(1.29\times10^{-26}\) for the phi sector across N=60/120/240.

The DM expanded comparator residual is coarse-grid limited at N=60 but falls to about \(10^{-17}\) at N=120/240. The S-sector expanded residual against its zero-source comparator also decreases at about fourth order at fixed r=2. The sector-expanded sum matches the expanded total candidate to about \(10^{-17}\), while the pointwise sector sum separately matches the combined-flux total to \(1.39\times10^{-17}\), \(3.78\times10^{-17}\), and \(3.03\times10^{-17}\). Thus the two operators each preserve their own sector additivity; the discrepancy is between the operators.

**Limit:** this is an implementation-level consistency check under the provisional archive-locked comparator and the current sector stress normalization. It is not an independent derivation of the historical source law, and it does not close full matter-coupled Brown propagation, total stress-energy conservation, or scientific admission. No evolution was rerun; no production or vendor code was changed.

**Artifact provenance:** exact head 2b723e7aaf129f82ab9270fe28a2587c5f8860b0; artifact IDs N60 12606725843, N120 12606133211, N240 12606611057. All source/vendor guard checks passed. Report/profile SHA-256: N60 cb8ceae8e81dc351150306996f14d90450a7d6c68e0ca6c164f2649b7d8c4bd7 / 47b93a6a859b602bde5515ab8923e7517bd355d059c8703e3397c4bbcf61191f; N120 47e170b42998be47cff4437cedb0a014b1b951513a5b7b3449e1e7055dd497f1 / bc1a3d23b5c7f6b2358c170d4fca87972607a79b36dfc5b6da50f951d931aaae; N240 42683a5e5cb187decdd22bce6d4b118d6dadc0d1fe71780bca7d2c4ba018d4e9 / 595079087e65ff9c47568b0a8bbf7dc395c4462d2c02fa1f124f9f81ede16b3e.

## 9. Disposition and next discriminating check

1. The origin-volume suspicion is resolved: the observer reconstructs regular (r^2g(r^2)) face profiles and Simpson-integrates the α-weighted spacetime measure; it does not use an arithmetic face-weight average or centre value times Δr.
2. The pinned S-sector RHS spatial terms reconstruct exactly from the existing current-tangent ledger. That isolated term is not the full `Q_r^S`, which also contains time-density, radial pressure-flux, and metric-source contributions.
3. The saved-state RHS/current/measure pairing is now instrumented at N=60/120/240, including the exact RHS-rate reconstruction, the spatial-measure log rate, pointwise and FV pairing gaps, and the volume-weighted shell comparison. The current and measure time-density operands close to about 1e-15, while the pointwise/FV pairing remainders decrease at fixed coordinate radii and in shell-weighted comparisons. They remain discrete operator-pairing diagnostics, not an admitted stress-energy conservation identity. Do not select a divergence form based solely on the smaller residual.
4. The sector-expanded comparator now tracks the provisional archive exchange law on these saved states, while the pointwise-vs-expanded product-rule gap is isolated; this remains implementation consistency only. The full matter-coupled Brown momentum-propagation identity, total stress-energy sector conservation, archive provenance/normalization admission, and scientific admission all remain open.

**Disposition:** diagnostic only; production equations and vendor sources unchanged; `scientific_admission=false`.

## 10. Brown $M_r$ residual: signed geometry versus source ledger — 2026-10-25

This pass reuses the corrected accepted-state artifacts from [workflow run 40213110427](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40213110427), resolutions N=60/120/240, accepted time t=0.60, and both recorded directional epsilons (2e-5 and 1e-5). The frozen source/vendor pins and diagnostic-only scope are the same as above; no evolution was advanced.

The published-covariant-$M_r$ measured-minus-predicted profile is compared with the independent geometric remainder `Mgeom_dot_nonmatter_measured - geometric_Brown_vacuum_source`. These two arrays agree to floating-point reconstruction accuracy. `exact_chain_rule_residual` remains at about 1e-13 or smaller. Thus the reported Brown-candidate residual is localized in the geometric/connection propagation comparison, not in the exact mass-decomposition chain rule. The source-mapping remainder `actual_matter_source_composite - stress_divergence_source` is a separate unresolved ledger item; it must not be substituted for the Brown residual or called closed.

### Matched-coordinate samples: full-continuous-profile arm

Values below are linear interpolations of the saved cell-centred profiles to the same coordinate radius, not comparisons at the same cell index.

| N | Brown measured-minus-predicted at r=2 | Brown geometric-source sum at r=2 | Connection-second term | Connection-first term | H-gradient term | Separate source-mapping remainder |
|---:|---:|---:|---:|---:|---:|---:|
| 60 | +8.975e-6 | +7.398e-7 | -6.478e-7 | -4.195e-7 | +1.914e-6 | +5.809e-5 |
| 120 | +6.12e-7 | +5.40e-8 | -6.82e-8 | +1.63e-8 | +1.13e-7 | +5.78e-5 |
| 240 | +4.03e-8 | +4.01e-9 | -6.03e-9 | +3.56e-9 | +6.99e-9 | +5.84e-5 |

The matched-r residual drops by roughly factors 14.7 and 15.2 when N doubles from 60 to 120 to 240, consistent with near-fourth-order behavior over this sample. The result is an implementation-consistency check on this arm, not an admitted propagation identity. In particular, the source-mapping remainder stays around 5.8e-5 and is visibly distinct from the small Brown geometric residual.

### Origin-cell signature: production-baseline arm

The first cell is at r=40/N, so the following rows are an origin-cell diagnostic, **not** a matched-coordinate convergence table.

| N | First-cell radius | Brown measured-minus-predicted | Connection-second term | Connection-first term | Separate source-mapping remainder |
|---:|---:|---:|---:|---:|---:|
| 60 | 0.667 | -4.604e-4 | -1.660e-3 | +1.278e-3 | +2.225e-5 |
| 120 | 0.333 | -7.550e-3 | -4.747e-3 | +2.501e-3 | +1.137e-5 |
| 240 | 0.167 | +3.249e-2 | +1.700e-2 | -2.712e-2 | +5.706e-6 |

The production-baseline origin-cell residual grows and changes sign under refinement, while the two signed connection terms grow and strongly cancel. The source-mapping remainder on those same cells is much smaller than this residual. By contrast, the full-continuous-profile residual stays small at the origin and decreases strongly at matched radii. This flags an origin-localized geometric/connection regularity or endpoint-response discrepancy in the production-baseline arm; it does **not** yet prove which construction causes it, because the arms also have different accepted-state/profile inputs and different lapse-perturbation policies.

Halving the directional epsilon changes the measured-minus-predicted residual by at most 1.3e-12 across these saved profiles. The observed central mismatch is therefore not explained by the finite-difference epsilon.

### Disposition and next discriminating check

1. The Brown $M_r$ residual and the matter-source mapping remainder are separate: the former tracks the non-matter geometric remainder; the latter still blocks total stress-energy conservation admission.
2. The production-baseline first-cell connection response is non-regular under refinement in this observer. Do not suppress it by changing the production equation, dropping a term, or tuning a coefficient.
3. The next diagnostic should hold **one identical accepted base state** fixed and compare its perturbed endpoints with only the lapse-response path changed (production CMC solve versus copied regular-origin CMC row). Record signed changes in lapse and radial-lapse derivative, each connection/H-gradient term, and the Brown residual. This isolates the current arm/lapse-policy confound without another evolution run.
4. Full matter-coupled Brown propagation, source sign/normalization, reference-metric connection convention, total stress-energy sector conservation, and scientific admission remain open.

**Disposition:** diagnostic only; production equations, vendor sources, and dependencies unchanged; `scientific_admission=false`.

## 11. Same-base paired endpoint-lapse test — 2026-10-25

The follow-up holds the **same exact accepted production-baseline state** at each resolution and reuses its identical geometry/scalar/matter RHS displacement. It creates two matched +/−epsilon endpoint pairs, changing only the lapse recomputation path: pinned production CMC solve versus copied regular-origin CMC row. The observer checks the accepted-state digest before/after, checks the RHS displacement digest in both arms, and requires every non-lapse endpoint field to be bit-identical.

This was completed by [workflow run 40214349711](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40214349711) for N=60/120/240. All three matrix jobs succeeded. Both epsilon values (2e-5 and 1e-5) were checked for every resolution.

| N | Shared accepted-state SHA256 | Δ(dot α) at r=2 (regular-origin − production) | Δ(Dᵣ dot α) at r=2 | Max absolute Δ Brown residual over full profile | Max absolute Δ of either connection term |
|---:|---|---:|---:|---:|---:|
| 60 | `46948673bfb78ee03ce17ab6a854823f107f81871b971f2b5d065df4b4ea25d7` | -1.154312e-3 | +1.233369e-3 | 0 | 0 |
| 120 | `88e97b0c57c7177079fc7a75c50d3517034e8780b1e6b127a8d114afc4ee7c1c` | -1.506529e-4 | +6.233551e-5 | 0 | 0 |
| 240 | `fea489afccc1d86cb9a6220c519c94fd3ca5576eda8bca651cd496cf66bdd480` | -1.848800e-5 | +9.602464e-6 | 0 | 0 |

The lapse-response difference is not a null perturbation: at r=2 the directional alpha rate and its radial gradient differ at all three resolutions. Nevertheless, the entire signed profile delta is **exactly zero** for each of these 11 reported quantities at both epsilons: direct covariant Brown residual, total Brown-plus-actual-source residual, first and second connection terms, H-gradient term, raised-KM term, alpha-H term, index-conversion term, provisional matter source, matter-source mapping remainder, and exact chain-rule residual. The endpoint alpha profiles themselves differ; the accepted base state is immutable and all non-lapse endpoint fields match exactly.

**Interpretation:** the choice of *recomputed perturbed-endpoint lapse path alone* cannot explain the observed Brown residual in this observer. This falsifies that specific endpoint-solve explanation; it does not prove the correct continuum identity or exclude a base-state metric/lapse coupling in the full derivation. The production-baseline central residual therefore remains an accepted-base geometric/connection or operator-construction problem. No lapse solver, coefficient, or production equation is changed.

### Reproducibility manifest

| N | Artifact ID | report.json SHA256 | profiles.npz SHA256 |
|---:|---:|---|---|
| 60 | 12606259731 | `8bc2d77b3e68cb459150d5139a0a03d3edd4ce5e84a1ad5094801cf6377da374` | `dbe81ba7f0204d7f7491c4e7eecdb73c74ef059a624cc5409404799a7d944a0e` |
| 120 | 12606503531 | `dc0ee329e992cab7b1e2292e2f16b35153d37aeca8a442ea335d4447b7296df8` | `22826082c288f37a105b1d1ed37ef60075b69a66c28fd0874de978b48dcf6dc1` |
| 240 | 12607355089 | `391df1aaf841583b0b54dcb91373207fba7f73d70536c4ab21a7ed295e393d83` | `e4e22742d6da68dc6d621c4ca1fb4d5b2e984e64f0534b1909de5e978d015d14` |

The input artifact SHA manifests were verified in each run. Frozen production source remains `4aec88d0e3ff028e0a65a8ac36b30948625df792`, vendor remains `d6052d605673ce9d82cdc99b0df79855fa2ea215`, and the protected source/vendor/dependency diff is empty.

### Disposition / next discriminating check

The next pass should decompose the **accepted-base geometric** residual using the Brown terms in one fixed covariant convention and compare each signed contribution to the measured geometric \(\dot M_r^{geom}\) tangent: \(\alpha K M_r\), \(-D_r\alpha H/3\), \(\alpha D_r H/6\), the two reference-connection terms, and the exact discrete chain-rule pieces. Keep the actual matter-source composite separate from that geometric budget. The same-base lapse endpoint test is complete; do not rerun it or tune a lapse parameter.

**Disposition:** diagnostic only; no evolution rerun; production equations, vendor sources, and dependencies unchanged; `scientific_admission=false`.


## 12. Signed geometric Brown term accounting — 2026-10-25

**Status:** diagnostic-only reassembly on the existing SHA-verified accepted states. No evolution was rerun. The observer-only report addition is in [commit 6f0e7e254be83e03c4f81138b7c210450258a99e](https://github.com/miahpenn/0-Cosmos-engine/commit/6f0e7e254be83e03c4f81138b7c210450258a99e). Exact-head workflow [40218549760](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40218549760) passed N=60/120/240.

The pinned source anchor remains `4aec88d0e3ff028e0a65a8ac36b30948625df792`; vendor remains `d6052d605673ce9d82cdc99b0df79855fa2ea215`. The observer verifies the input artifact manifests, frozen-state immutability, finite profiles, and the protected source/vendor/dependency diff. No production source, initializer, vendor, dependency, or physics convention was changed.

### Definition and accounting check

The measured geometric tangent is separated before comparison:

```text
G_measured = measured Mgeom_dot - explicit Aa/K matter-injection contribution
B_geom = shift + alpha*K*Mgeom
         - (D_r alpha)*Hgeom/3 + alpha*(D_r Hgeom)/6
         + reference-connection-second + reference-connection-first
R_geom = G_measured - B_geom
```

On these accepted CMC states the shift-transport term is zero. The report now emits every signed term at native cells 0–2 and at matched coordinate radii r=2, 5, 8, while retaining the exact discrete chain-rule residual and matter-source mapping remainder as separate ledger entries.

The saved arrays reconstruct `B_geom` with maximum absolute sum gap **0** at all three resolutions and both directional epsilons (2e-5 and 1e-5); the independently recorded reconstruction-gap profile also matches with maximum absolute difference **0**. The configured fail-closed bound is (1.14\times10^{-13}). This is a profile serialization/arithmetic-accounting check—the candidate source itself is defined as the same signed sum—not evidence that the full continuum propagation identity is correct.

The discrete chain-rule residual remains small (maximum absolute value over the three arms and two epsilons is between (1.2\times10^{-13}) and (1.2\times10^{-12}), depending on resolution/arm). Halving the directional epsilon changes the tracked Brown-source, non-matter-tangent, residual, and connection-term profiles by at most about (1.3\times10^{-12}). The central discrepancy is therefore not explained by a finite-difference epsilon artifact.

### Production-baseline origin-cell ledger

The origin cell is at (r=40/N); do not interpret these values as a fixed-coordinate convergence sequence. Values below are the actual saved signed terms.

| N | cell-0 r | \(\alpha K M_{geom}\) | \(-D_r\alpha H_{geom}/3\) | \(\alpha D_r H_{geom}/6\) | connection second | connection first | Brown term sum | measured non-matter \(\dot M_{geom}\) | measured minus Brown |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 60 | 0.667 | -2.158e-5 | +1.163e-8 | -1.044e-4 | -1.660e-3 | +1.278e-3 | -5.082e-4 | -9.686e-4 | -4.604e-4 |
| 120 | 0.333 | +1.835e-5 | +5.582e-9 | -2.915e-5 | -4.747e-3 | +2.501e-3 | -2.256e-3 | -9.806e-3 | -7.550e-3 |
| 240 | 0.167 | +1.154e-5 | +1.878e-8 | -5.494e-3 | +1.700e-2 | -2.712e-2 | -1.560e-2 | +1.689e-2 | +3.249e-2 |

The two signed connection terms are individually large and strongly cancel. At N=240 the H-gradient term is also large; the full residual cannot be assigned to one connection term from this table. The measured-minus-Brown residual changes sign and grows under refinement, while the separate source-mapping remainder is smaller on these origin cells. This is evidence of a central geometric/connection or regular-origin construction problem in the production-baseline arm, not identification of a single wrong term.

### Matched-coordinate comparison

All entries are measured-minus-Brown residuals. The coordinate r=2 is a cell centre at each of N=60/120/240; it is not an origin-cell sequence.

| Arm | r=2, N=60 | r=2, N=120 | r=2, N=240 |
|---|---:|---:|---:|
| production baseline | -5.088e-4 | +2.526e-4 | -6.914e-5 |
| partial joint | +9.836e-5 | -3.848e-5 | +1.670e-6 |
| full continuous profile | +8.975e-6 | +6.122e-7 | +4.029e-8 |

The full-continuous-profile arm retains approximately fourth-order absolute decrease at r=2. The production-baseline r=2 sequence remains sign-changing and non-monotonic. At r=5, the baseline residual changes from -4.154e-5 (N=60) to -1.759e-8 (N=120) and +5.462e-9 (N=240); at r=8 it decreases from -1.719e-6 to -7.688e-8 to -5.547e-9. This supports a near-centre localization, while the production-baseline r=2 witness is still unresolved.

The separate source-mapping remainder at r=2 stays around (5.8\times10^{-5}) across arms/resolutions. It must not be equated with `R_geom` or substituted into the geometric Brown budget. The same-base endpoint lapse-path test in Section 11 remains a negative control: despite a nonzero lapse-rate change, all tracked Brown/connection/source/chain-rule profiles have exactly zero delta between lapse paths.

### Reproducibility and disposition

Accepted production-baseline state digests are unchanged from Section 11:

| N | accepted-state SHA-256 |
|---:|---|
| 60 | `46948673bfb78ee03ce17ab6a854823f107f81871b971f2b5d065df4b4ea25d7` |
| 120 | `88e97b0c57c7177079fc7a75c50d3517034e8780b1e6b127a8d114afc4ee7c1c` |
| 240 | `fea489afccc1d86cb9a6220c519c94fd3ca5576eda8bca651cd496cf66bdd480` |

Artifact `report.json / profiles.npz` SHA-256:

| N | report.json | profiles.npz |
|---:|---|---|
| 60 | `b0251c8caf79f74c9b3105f3346de9b97bb4603b4ba5ede039b2e0c02b496b0b` | `dc55a6876c19acd498e6096215896c65faec4e3bed937331e84e610429a85ffd` |
| 120 | `639cfd7d8bb1ef32aa7b3ce5b56401840ddd7cb813874c2c5bfb6fa6f7571b14` | `1a34088e870fe443bd3a00905b0b4302596479c556eaeea52833b4f598fa9cda` |
| 240 | `e35aff78341efedf271f9b981029f21d01fceebb020cbde75e482c7813eb62bf` | `e9b15827c49e2c0aedbc450aeba324967abb5d7e5b5bdd2f50197562c5a1febf` |

**Next discriminator:** do not repeat the lapse-path test or tune/drop individual terms. Freeze the reference-metric connection convention and origin-limit form, then derive/compare the convention-pinned coupled propagation identity against the exact accepted-state tangent, retaining cells 0–2 separately from fixed-r samples. Keep matter-source sign/normalization, total-sector conservation, and scientific admission open.

**Disposition:** diagnostic only; `full_matter_coupled_propagation_identity_closed=false`, `total_stress_energy_sector_conservation_closed=false`, `scientific_admission=false`.

## 13. Brown Hamiltonian-slot convention audit — 2026-10-25

**Status:** research-only post-processing of the already accepted, SHA-verified \(t=0.60\) states. No evolution replay, production edit, source-pin change, or admission change.

### Convention finding

Brown Eq. (11b)'s \(H\) and \(D_r H\) slots use the **metric-only Hamiltonian constraint** for its geometric/vacuum comparator. The repository explicitly distinguishes this from the evolved BSSN Hamiltonian: the existing audit helper **tools/audit_brown_vacuum_momentum_propagation.py** provides a metric-only construction and a Lambda-only discriminator. At fixed metric, the discriminator leaves \(H_{\rm metric}\) unchanged but confirms that the BSSN Hamiltonian responds as

\[
\delta H_{\rm BSSN}=X^2 D_r(\delta\Lambda).
\]

The saved-state geometric Brown ledger had sourced \(H\) from the vendor BSSN constraint inside **tools/audit_hybrid_origin_integral_patch_newton.py**. That is the Lambda-dependent BSSN representation, not the metric-only geometric \(H\) required for Brown's two \(H\)-slots. The ledger's signed term-sum arithmetic remains valid, but its earlier label as a Brown geometric comparator was incomplete: it was a **hybrid convention comparator** in those slots. This is an observer/convention issue, not evidence of a production-equation defect.

### Same-state H-slot counterfactual

Using saved accepted-state \(a,b,X,A_a,K,\alpha\) profiles and the repository's parity-aware fourth-order radial operators, I reassembled only

\[
T_{\alpha H}=-\frac{D_r\alpha}{3}H,\qquad
T_{\nabla H}=\frac{\alpha}{6}D_rH
\]

with metric-only \(H\), retaining the other signed geometric Brown and connection terms unchanged. The corrected residual is

\[
R_{\mathrm{metricH}}=R_{\mathrm{saved}}
+T_{\alpha H}^{\mathrm{saved}}+T_{\nabla H}^{\mathrm{saved}}
-T_{\alpha H}^{\mathrm{metric}}-T_{\nabla H}^{\mathrm{metric}}.
\]

This is a narrow convention counterfactual on the saved arrays, not an independent accepted-state oracle.

| N | Full-profile residual at fixed \(r=2\), saved H slot | Same residual with metric-only H | Observed order for metric-H residual |
|---:|---:|---:|---:|
| 60 | \(+8.9748281\times10^{-6}\) | \(+9.7763622\times10^{-6}\) | — |
| 120 | \(+6.1215367\times10^{-7}\) | \(+6.3049133\times10^{-7}\) | 3.955 |
| 240 | \(+4.0294555\times10^{-8}\) | \(+4.1874593\times10^{-8}\) | 3.912 |

The fixed-radius absolute remainder retains approximately fourth-order decrease after this slot correction, so the H-slot mismatch does **not** explain the full-profile annular remainder. It does materially change some baseline/partial-arm centre terms: at \(N=240\), the production-baseline cell-0 residual changes from \(+3.2485246\times10^{-2}\) to \(+3.3113882\times10^{-2}\), and cell 1 from \(-4.0490374\times10^{-3}\) to \(-4.2418994\times10^{-3}\). The origin discrepancy plainly remains; no single connection term is thereby identified as defective.

### Observer patch and next bounded action

The research-only observer patch is now on this PR branch. The geometric Brown comparator uses the existing metric-only \(H\); the provisional total candidate uses \(H_{\rm metric}-16\pi\rho\), and the BSSN-representation profiles remain separately recorded. The accepted-state driver also requires and reports the new H-slot witness arrays. Changes are confined to **tools/audit_hybrid_origin_integral_patch_newton.py** and **tools/audit_full_profile_propagation_candidate.py**; no production or vendor files were changed.

The same-state validation [workflow run 40224825903](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40224825903) **completed successfully for N=60, 120, and 240**. The workflow verified each input-artifact checksum, kept each frozen accepted-state observer immutable, confirmed the exact source/vendor pins and empty protected engine/vendor/dependency diff, compiled the observer, and evaluated the signed ledger against the saved \(t=0.60\) accepted states. The accepted-state digests match the pre-patch reports exactly; no evolution replay was performed.

The corrected observer's direct covariant-\(M_r\) candidate now reports:

| N | Residual at fixed \(r=2\) | Residual at fixed \(r=8\) | Full-domain RMS |
|---:|---:|---:|---:|
| 60 | \(+9.77672\times10^{-6}\) | \(-1.33948\times10^{-6}\) | \(1.60784\times10^{-6}\) |
| 120 | \(+6.30515\times10^{-7}\) | \(-7.88855\times10^{-8}\) | \(1.39001\times10^{-7}\) |
| 240 | \(+4.18762\times10^{-8}\) | \(-4.83828\times10^{-9}\) | \(1.42940\times10^{-8}\) |

At fixed \(r=2\), the observed orders are 3.955 and 3.912; at fixed \(r=8\), they are approximately 4.086 and 4.027. The corrected geometric Brown-plus-actual-source remainder agrees closely with the direct candidate at these radii. The signed Brown term sum reconstructs exactly (zero arithmetic gap at the stored precision), and the discrete chain-rule observer remains a separate implementation check.

**Correction after the patched signed-bridge check:** the source-mapping difference \(\Delta S=S_{\rm act}-S_Q\) is large and nearly resolution-independent at fixed \(r=2\), but it is **not an isolated missing-source residual**. On the same H-corrected accepted-state profiles it is balanced by the Brown source difference \(B_{\rm tot}-B_{\rm geom}\); see Section 15 below. Its persistence alone does not imply failed conservation. The still-open gate is the archive-pinned interpretation and normalization of the complete sector source identity, alongside the independent origin/reference-connection and total-sector conservation gates—not driving \(\Delta S\) to zero by fitting or cancelling it.

The regular-origin question also remains open: native cell 0 samples are at different radii for the three resolutions, and the corrected direct residual there changes sign (\(+8.19\times10^{-7}\), \(-6.57\times10^{-7}\), \(+1.36\times10^{-7}\)). Those values are **not** a matched-radius convergence sequence.

Next, keep the existing accepted-state ledger and frozen pins. Derive the matter-source mapping from the validated sector equations—including scalar–DM exchange and the repository's normalization/index definitions—rather than fitting the residual; retain the reference-connection/origin-limit gate as a separate check. Keep \(M_r/M^r\) conversion and metric-rate contribution explicit. No duplicate sweeps, production equation changes, or admission changes.

This does not close the full matter-coupled identity. The current flags remain false: full_matter_coupled_propagation_identity_closed=false, reference_metric_connection_coupling_convention_frozen=false, total_stress_energy_sector_conservation_closed=false, and scientific_admission=false.

## 14. Index-conversion scale check — 2026-10-25

**Status:** diagnostic inspection of the pre-H-slot-correction accepted-state ledger. The conversion formula itself is unchanged by the H-slot patch; the residual values below belong to the earlier comparator and are not presented as updated post-patch residuals.

The observer keeps the covariant Brown variable \(M_r\) and the raised-index comparison \(M^r=\gamma^{rr}M_r\) distinct. With \(\gamma_{rr}=a/X^2\), lowering the raised-index rate requires

\[
\partial_t M_r
=\gamma_{rr}\partial_t M^r
+\left(\partial_t\ln\gamma_{rr}\right)M_r,
\qquad
\partial_t\ln\gamma_{rr}=\frac{\dot a}{a}-2\frac{\dot X}{X}.
\]

The existing comparator explicitly adds this metric-rate contribution as \((\partial_t\ln\gamma_{rr})M_r\); the direct covariant-\(M_r\) candidate is computed separately. Thus the term is not silently omitted from the raised-index conversion.

At fixed \(r=2\), the earlier saved full-profile ledger gives:

| N | Metric-rate index-conversion term | Earlier direct covariant-\(M_r\) residual | Conversion-term scale / residual |
|---:|---:|---:|---:|
| 60 | \(+7.1804\times10^{-8}\) | \(+8.9752\times10^{-6}\) | 0.800% |
| 120 | \(+4.9837\times10^{-9}\) | \(+6.1218\times10^{-7}\) | 0.814% |
| 240 | \(+3.3844\times10^{-10}\) | \(+4.0296\times10^{-8}\) | 0.840% |

The ratio is only a **scale comparison**: the conversion term belongs to the raised-index comparator and is not an additive term in the direct covariant residual. It suggests the explicitly recorded metric-rate contribution is not, by magnitude alone at this radius, a dominant explanation for the earlier residual. It does not establish the complete index convention, and the direct-candidate residual will change under the H-slot correction now being validated.

No nonzero-shift conclusion follows from these accepted arms. The index conversion, Brown covariant equation, connection convention, and matter-source mapping remain separately tracked gates; scientific admission stays false.


## 15. H-corrected signed source bridge — 2026-10-25

**Status:** post-processing of the existing accepted \(t=0.60\) profiles from workflow [40224825903](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40224825903). The run passed at N=60/120/240. Input profile hashes and accepted-state digests were unchanged; the source/vendor pins and protected engine/vendor/dependency files were unchanged; no evolution was advanced.

Use the established definitions
\[
\begin{aligned}
B_{\rm tot}&=\text{Brown }M_r\text{ source using the total metric-only }H,\\
B_{\rm geom}&=\text{geometric-only Brown }M_r\text{ source},\\
S_{\rm act}&=\text{actual }M_{\rm geom}\text{ matter injection}-8\pi\dot J_r,\\
S_Q&=-8\pi\alpha Q_r,\qquad \Delta S=S_{\rm act}-S_Q,\\
D_{\rm tot}&=\dot M_{\rm measured}-(B_{\rm tot}+S_Q),\\
D_{\rm geom}&=\dot M_{\rm measured}-(B_{\rm geom}+S_{\rm act}).
\end{aligned}
\]
The bridge to test is
\[
(B_{\rm tot}-B_{\rm geom})-\Delta S=D_{\rm geom}-D_{\rm tot}.
\]

### Patched accepted-state results at fixed \(r=2\)

| N | \(B_{\rm tot}-B_{\rm geom}\) | \(\Delta S\) | Difference \((B_{\rm tot}-B_{\rm geom})-\Delta S\) | Max absolute difference over profile | Bridge arithmetic error, max abs |
|---:|---:|---:|---:|---:|---:|
| 60 | \(+5.80896357\times10^{-5}\) | \(+5.80899982\times10^{-5}\) | \(-3.62482\times10^{-10}\) | \(3.62482\times10^{-10}\) | \(7.62\times10^{-21}\) |
| 120 | \(+5.77976279\times10^{-5}\) | \(+5.77976520\times10^{-5}\) | \(-2.40650\times10^{-11}\) | \(3.00874\times10^{-11}\) | \(7.94\times10^{-21}\) |
| 240 | \(+5.83505119\times10^{-5}\) | \(+5.83505136\times10^{-5}\) | \(-1.60891\times10^{-12}\) | \(2.01493\times10^{-12}\) | \(1.23\times10^{-20}\) |

Thus the total-versus-geometric Brown source difference tracks the matter-source composite difference increasingly closely with refinement. The signed identity \((B_{\rm tot}-B_{\rm geom})-\Delta S=D_{\rm geom}-D_{\rm tot}\) holds to roundoff in the saved arrays. This is a strong algebraic crosswalk of the current observer's terms; it is **not** a standalone derivation of the archive's physical matter source, and it does not close the full matter-coupled propagation identity.

### Stress-divergence and sector accounting controls

At fixed \(r=2\), the corrected observer reports the existing pointwise \(\nabla_\mu T^\mu{}_r\)-style divergence candidate as \(+1.78829\times10^{-7}\), \(+1.29092\times10^{-8}\), and \(+7.97334\times10^{-10}\) for N=60/120/240. The independent finite-volume divergence reports \(+4.49074\times10^{-9}\), \(+2.81380\times10^{-10}\), and \(+1.70483\times10^{-11}\). Their mismatch with the expanded \(Q_r\) candidate decreases strongly at fixed radius, and the pointwise-minus-expanded decomposition closes to roundoff. This remains a comparison of discretizations/conventions, not permission to select whichever form gives the smaller residual.

The component-corrected sector decomposition's unexplained total-divergence gap is at roundoff scale (maximum absolute \(3.12\times10^{-21}\), \(2.99\times10^{-20}\), \(5.47\times10^{-20}\) for N=60/120/240). The raw, uncorrected sector-sum booleans remain false, but the measured pre-correction maximum absolute gaps are only \(1.39\times10^{-17}\), \(3.78\times10^{-17}\), and \(3.03\times10^{-17}\); the saved decomposition identifies the time-density finite-difference-rate difference and accounts for it, leaving the quoted roundoff-scale unexplained gap. This is an accounting/tolerance flag, not evidence by itself of a physical conservation failure. The complete total-stress-energy gate nevertheless remains open because the authoritative sector normalization/source-law interpretation and regular-origin/reference-connection identity have not been admitted. The archive exchange pair remains locked as
\[
\nabla_\mu T_{\rm DM}^{\mu\nu}=+\beta\rho_{\rm DM}\nabla^\nu\phi,\qquad
\nabla_\mu T_\phi^{\mu\nu}=-\beta\rho_{\rm DM}\nabla^\nu\phi,
\]
with opposite signs and net exchange zero. The local \(S/D\) projections and COSMOS-\(\phi\)/fluid projections do not all share one normalization; do not propagate a single \(\kappa\) factor across sectors without the archive's specific definitions.

### Remaining bounded work

1. Preserve the algebraically closed source bridge as a **consistency witness**, not as a fitted source law.
2. Review the specific source normalization and index conversion against the archive definitions for each stress sector—especially the local \(S\) contribution and the reciprocal COSMOS-\(\phi\)/DM exchange—without modifying production equations.
3. Keep the reference-metric connection and regular-origin identity as a separate gate. Native cells 0–2 are not matched-radius samples; no origin closure is inferred from the fourth-order fixed-radius trend.

The current flags remain false: full_matter_coupled_propagation_identity_closed=false, reference_metric_connection_coupling_convention_frozen=false, total_stress_energy_sector_conservation_closed=false, scientific_admission=false.


## 16. Brown reference-connection convention and origin-limit check — 2026-10-25

**External convention anchor:** Brown, [*BSSN in Spherical Symmetry*, arXiv:0705.3845](https://arxiv.org/abs/0705.3845), defines \(G^r=\Gamma^r-g^{bc}\tilde{\Gamma}^{r}{}_{bc}\) in Eq. (10c) and gives the vacuum \(M_r\) propagation law in Eq. (11b). With the repository variables \(\tilde\gamma_{rr}=a\), \(\tilde\gamma_{\theta\theta}=r^2b\), and \(\chi=X^2\), its reference-metric connection is
\[
\tilde\Gamma^r_{\rm metric}
=\frac{D_ra}{2a^2}-\frac{D_rb}{ab}
+\frac{2}{r}\left(\frac1b-\frac1a\right),
\qquad
C^r_\Lambda=\Lambda^r-\tilde\Gamma^r_{\rm metric}.
\]
The two propagation terms are
\[
T_{G''}=\frac{2}{3}\alpha\chi D_r^2C^r_\Lambda,\quad
T_{G'}=\left(\frac23\alpha'\chi-\frac13\alpha\chi'
+\alpha\chi\left[\frac{D_rb}{b}+\frac2r\right]\right)D_rC^r_\Lambda.
\]

### Same-state convention reconstruction

I reconstructed \(C^r_\Lambda\) from the stored accepted-state \(a,b,X,\Lambda\) arrays using the **pinned vendor's actual fourth-order parity-aware derivative stencil**, then independently rebuilt \(T_{G''}\) and \(T_{G'}\). These are same-state post-processing checks; no RHS/evolution replay was used.

| N | \(C^r_\Lambda(r=2)\) | Max \(|C^r_\Lambda|\) on the saved profile | Reconstructed-vs-stored \(C^r_\Lambda\) max difference |
|---:|---:|---:|---:|
| 60 | \(+3.90956\times10^{-6}\) | \(3.90956\times10^{-6}\) | \(0\) |
| 120 | \(+2.78488\times10^{-7}\) | \(2.82346\times10^{-7}\) | \(6.37\times10^{-16}\) |
| 240 | \(+2.07568\times10^{-8}\) | \(2.14683\times10^{-8}\) | \(1.67\times10^{-16}\) |

At the fixed radius \(r=2\), the connection-constraint profile decreases at observed orders approximately 3.81 and 3.75. The two Brown connection terms also reproduce the corresponding saved signed profiles to floating-point precision (maximum difference \(3.31\times10^{-15}\) for \(T_{G''}\), \(5.59\times10^{-17}\) for \(T_{G'}\)). This independently checks the **sign, index, radial \(r^2b\) factor, and coefficient assembly** in this observer against Brown's Eq. (10c)/(11b); it does not prove that the complete matter-coupled law or the continuum origin limit is closed.

For context, at fixed \(r=2\), the N=240 connection contributions are \(T_{G''}=-6.03\times10^{-9}\) and \(T_{G'}=+3.56\times10^{-9}\), both much smaller than the direct full-profile residual \(+4.19\times10^{-8}\). This agrees with the prior diagnostic caution that no one connection term explains that annular residual.

The first-cell values are **not** a matched-radius convergence sequence: the centers are \(r=0.667,0.333,0.167\) for N=60/120/240. They must not be fit as if they were observations at \(r=0\). Moreover, the large N=240 production-baseline first-cell residual remains a different profile arm from the smooth full-continuous-profile state used for the fixed-\(r\) convergence table. Therefore this reconstruction rules out an obvious reference-connection **observer coefficient/sign mismatch**, but it does not resolve the baseline regular-origin discrepancy.

### Decision

No further change to Brown's connection coefficients is indicated by this saved-state reconstruction. Keep the current observer convention, keep the regular-origin gate open, and isolate the baseline first-cell mismatch using a matched-origin/regularity argument on the *same named state family* before drawing a physical conclusion. No production or vendor file changed; all scientific-admission flags remain false.



## 17. Matched connection-profile contrast across the existing origin arms — 2026-10-25

**Method:** I independently reconstructed the reference-metric connection and Brown \(T_{G''},T_{G'}\) terms from the already archived \(t=0.60\) production-baseline, half-cell-integral-only, regular-CMC-row-only, joint-integral-and-CMC-row, and full-continuous-profile state arrays. The derivative stencil is copied algebraically from the pinned vendor's **r0_operators.py** at **d6052d605673ce9d82cdc99b0df79855fa2ea215**: centered fourth-order five-point stencil with parity ghosts and its polynomial-moment outer closure. This is offline post-processing only; no evolution was run.

As a check of implementation fidelity, the reconstructed connection terms reproduce the previously saved signed term arrays:
- N=60: max differences \(2.17\times10^{-19}\) for \(T_{G''}\) and \(3.23\times10^{-27}\) for \(T_{G'}\).
- N=120: \(8.67\times10^{-19}\) and \(4.04\times10^{-28}\).
- N=240: \(1.12\times10^{-14}\) and \(1.27\times10^{-16}\).

### First-cell and full-profile scale

Here \(C_0\) denotes \(C^r_\Lambda\) at the first cell center for each resolution; \(\max|T_{G''}|\) and \(\max|T_{G'}|\) are full-domain profile maxima. These are **separate state arms**, not a same-state perturbation pair.

| Arm | N | \(C_0\) | \(\max |T_{G''}|\) | \(\max |T_{G'}|\) |
|---|---:|---:|---:|---:|
| production baseline | 60 | \(+1.07480\times10^{-3}\) | \(1.65991\times10^{-3}\) | \(1.27766\times10^{-3}\) |
| production baseline | 120 | \(+7.11697\times10^{-4}\) | \(4.74691\times10^{-3}\) | \(2.50145\times10^{-3}\) |
| production baseline | 240 | \(-8.58895\times10^{-4}\) | \(1.70039\times10^{-2}\) | \(2.71199\times10^{-2}\) |
| joint integral + CMC-row arm | 60 | \(-8.91046\times10^{-5}\) | \(1.46478\times10^{-4}\) | \(8.22120\times10^{-5}\) |
| joint integral + CMC-row arm | 120 | \(-6.20906\times10^{-5}\) | \(4.48570\times10^{-4}\) | \(1.17276\times10^{-4}\) |
| joint integral + CMC-row arm | 240 | \(+5.15855\times10^{-5}\) | \(1.20086\times10^{-3}\) | \(1.23446\times10^{-3}\) |
| full continuous-profile arm | 60 | \(+3.12809\times10^{-6}\) | \(2.48456\times10^{-6}\) | \(9.23599\times10^{-6}\) |
| full continuous-profile arm | 120 | \(+1.68878\times10^{-7}\) | \(5.17465\times10^{-7}\) | \(1.98616\times10^{-6}\) |
| full continuous-profile arm | 240 | \(-8.60088\times10^{-10}\) | \(4.71529\times10^{-8}\) | \(7.69101\times10^{-8}\) |

For N=240, the first four baseline \(C^r_\Lambda\) samples are \([-8.59,-2.98,+1.60,-0.897]\times10^{-4}\), changing sign on adjacent cell centers; its connection contributions are \(+1.7004\times10^{-2}\) and \(-2.7120\times10^{-2}\) at cell 0. In the full continuous-profile arm the first four connection-constraint samples are only \([-0.0860,+0.388,+1.049,+1.415]\times10^{-8}\), and the first-cell connection contributions are \(+4.7153\times10^{-8}\) and \(+2.9170\times10^{-8}\).

The joint integral+CMC-row arm reduces the N=240 full-profile \(\max|C^r_\Lambda|\) from \(8.59\times10^{-4}\) to \(5.16\times10^{-5}\) (about \(16.7\times\)), and reduces the maxima of the two connection terms by about \(14.2\times\) and \(22.0\times\). It still does not reach the full-continuous-profile scale. The half-cell-integral-only arm also reduces the baseline's max connection term levels, while the regular-CMC-row-only arm leaves the baseline-scale defect largely present. Those arm-to-arm changes are **localization evidence** consistent with an inner connection/initialization regularity defect; because the state arrays differ between arms, they do not isolate a single algebraic cause by themselves.

At fixed \(r=2\), the baseline \(C^r_\Lambda\) values are \(-2.39837\times10^{-4}\), \(+2.83213\times10^{-5}\), and \(-1.57285\times10^{-6}\) for N=60/120/240, while the large near-origin derivatives/connection terms grow at N=240. Thus a convergent fixed-radius tail can coexist with a badly behaved first-cell neighborhood. The native first-cell locations differ with resolution (\(0.667,0.333,0.167\)); this is not a matched-\(r=0\) sequence and must not be called a continuum origin-limit test.

### Decision

The pinned Brown connection observer formula and saved term arrays are consistent; the large baseline origin residual is co-localized with a sharp, alternating \(C^r_\Lambda\) inner profile and large \(D_r C^r_\Lambda\)/\(D_r^2 C^r_\Lambda\) operands. This makes the connection-constraint regularity at the origin the leading diagnostic target, not a reason to alter Brown's coefficients. The next bounded step should be an archive-pinned regularity/initialization identity for the same baseline state family, using matched-origin extrapolation and separate \(C^r_\Lambda\), \(D_rC^r_\Lambda\), and \(D_r^2C^r_\Lambda\) checks. Do not infer closure from the arm comparison; keep all admission flags false.
