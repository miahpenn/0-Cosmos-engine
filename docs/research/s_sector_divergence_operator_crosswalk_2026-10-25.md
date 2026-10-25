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

## 6. Disposition and next discriminating check

1. The origin-volume suspicion is resolved: the observer reconstructs regular (r^2g(r^2)) face profiles and Simpson-integrates the α-weighted spacetime measure; it does not use an arithmetic face-weight average or centre value times Δr.
2. The pinned S-sector RHS spatial terms reconstruct exactly from the existing current-tangent ledger. That isolated term is not the full `Q_r^S`, which also contains time-density, radial pressure-flux, and metric-source contributions.
3. Keep pointwise and FV divergence results separate and convention-labelled. Next, compare each with a fully reconstructed *scalar equation-of-motion residual* assembled from all signed RHS terms and the time-varying measure/current contributions, then integrate the matching identity over shells. Do not select a divergence form based solely on the smaller residual.
4. The full matter-coupled Brown momentum-propagation identity, total stress-energy sector conservation, the archive-locked φ/DM exchange admission gate, and scientific admission all remain open.

**Disposition:** diagnostic only; production equations and vendor sources unchanged; `scientific_admission=false`.
