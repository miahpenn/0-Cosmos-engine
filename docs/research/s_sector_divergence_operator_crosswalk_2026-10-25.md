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

$
\mathcal L_S[S]=\frac{X^2}{a}\left[D_2S+
\left(\frac{D_1b}{b}-\frac{D_1a}{2a}+\frac{2}{r}-\frac{D_1X}{X}\right)D_1S\right].
$

The requested isolated spatial contribution is (\mathcal L_S[S]D_1S). In the scalar momentum RHS, the corresponding Laplacian term is (\alpha\mathcal L_S[S]); a separate lapse-gradient term ( (X^2/a)D_1\alpha D_1S), curvature, potential, shift, and (for the coupled φ field) conditional exchange terms are separate RHS operands.

The sector stress-divergence observer instead evaluates the full (Q_r^S=\nabla_\mu T_S^\mu{}_r) from the time-density/current derivative, radial pressure flux, and metric-source terms. Its FV form reconstructs combined face fluxes and integrates sources and the measure per cell. Therefore, **(Q_r^S) is not the same quantity as the isolated S-sector Laplacian contribution (\mathcal L_S[S]D_1S)**. Their raw ratio is a useful crosswalk witness, but a mismatch does not by itself indicate a missing source or a bad scalar RHS.

## 2. The pinned RHS term reconstruction closes

The two signed current-tangent pieces from the scalar RHS,
`scalar_S_current_rate_from_Pdot_lapse_laplacian_second_derivative` and
`scalar_S_current_rate_from_Pdot_lapse_metric_gradient`, were combined and divided by (8\pi\alpha). Their mapped result matches (\mathcal L_S[S]D_1S) to floating-point roundoff:

| Resolution | Maximum absolute mapping gap | RMS mapping gap |
|---:|---:|---:|
| (N=60) | (1.06\times10^{-22}) | (2.11\times10^{-23}) |
| (N=120) | (7.94\times10^{-23}) | (1.26\times10^{-23}) |
| (N=240) | (1.06\times10^{-22}) | (9.80\times10^{-24}) |

This confirms that the diagnostic reconstructs the pinned RHS spatial contribution correctly. It is an implementation/source-budget identity, not an independent proof of the full scalar stress-energy conservation law.

## 3. Origin-cell measure

The FV covariant-divergence observer uses the spacetime measure
$
\sqrt{-g}=\alpha\sqrt{\gamma},
$
not the spatial measure (\sqrt{\gamma}) alone. It reconstructs `sqrt(-g)` to faces with the regular (r^2 g(r^2)) form and computes each measure by Simpson integration,
$
V_0=\frac{\Delta r}{6}\left[w(0)+4w(\Delta r/2)+w(\Delta r)\right],\quad w=\sqrt{-g},quad w(0)=0.
$
It does **not** use the centre value times Δr. The familiar (\frac43\Delta r\,w(\Delta r/2)) is an excellent near-origin approximation, not the exact Simpson value.

| Resolution | Δr | Observer (V_0=\int\alpha\sqrt{\gamma}\,dr) (Simpson) | Centre value × Δr | (\frac43) centre approximation | Relative Simpson / approximation difference |
|---:|---:|---:|---:|---:|---:|
| 60 | 1.333333 | 0.7974845181 | 0.5980217458 | 0.7973623278 | (1.53\times10^{-4}) |
| 120 | 0.666667 | 0.09966964221 | 0.07474927498 | 0.09966569998 | (3.96\times10^{-5}) |
| 240 | 0.333333 | 0.01245818990 | 0.00934354931 | 0.01245806575 | (9.97\times10^{-6}) |

For clarity, the same Simpson reconstruction using the *spatial* measure (\sqrt{\gamma}) alone gives (0.7984648512), (0.09979184239), and (0.01247345586) for (N=60,120,240). Those are not the measure values used by this covariant spacetime-divergence observer.

## 4. Direct S-only comparison requested by Lumen

Cell 0 is **not** a fixed-radius sequence: its coordinate radius is (r=\Delta r/2=0.6667,0.3333,0.1667) for the three resolutions. Do not infer a central convergence order from these three cell-0 samples.

| Resolution | (\mathcal L_S[S]D_1S), cell 0 | Pointwise (Q_r^S), cell 0 | Pointwise / isolated RHS term | FV (Q_r^S), cell 0 | FV / isolated RHS term |
|---:|---:|---:|---:|---:|---:|
| (N=60) | (1.98443\times10^{-7}) | (7.92656\times10^{-7}) | 3.9944 | (6.72404\times10^{-9}) | 0.033884 |
| (N=120) | (1.01414\times10^{-7}) | (1.13000\times10^{-7}) | 1.1143 | (9.04460\times10^{-10}) | 0.008919 |
| (N=240) | (5.09515\times10^{-8}) | (1.45989\times10^{-8}) | 0.28653 | (1.13385\times10^{-10}) | 0.002225 |

The cell-0 pointwise result is itself the signed sum of large cancelling terms, not the isolated RHS Laplacian. The recorded components are:

| Resolution | Time-density/current | Radial pressure-flux derivative | Metric source | Sum: pointwise Q_r^S |
|---:|---:|---:|---:|---:|
| N=60 | -7.45445e-7 | -3.49372e-5 | +3.64753e-5 | +7.92656e-7 |
| N=120 | -3.78191e-7 | -7.32343e-5 | +7.37255e-5 | +1.13000e-7 |
| N=240 | -1.89775e-7 | -1.47658e-4 | +1.47862e-4 | +1.45989e-8 |

On the fixed-coordinate annulus (2\le r\le8), the RMS of the isolated spatial contribution stays close to (3.9\times10^{-7}), whereas the full S-sector divergence decreases strongly:

| Resolution | RMS (\mathcal L_S[S]D_1S) | Pointwise (Q_r^S) RMS | Pointwise ratio RMS | FV (Q_r^S) RMS | FV ratio RMS |
|---:|---:|---:|---:|---:|---:|
| (N=60) | (3.96049\times10^{-7}) | (8.20593\times10^{-8}) | 0.16921 | (6.45932\times10^{-9}) | 0.037790 |
| (N=120) | (3.90021\times10^{-7}) | (3.30293\times10^{-9}) | 0.006548 | (4.30698\times10^{-10}) | 0.002724 |
| (N=240) | (3.89688\times10^{-7}) | (2.17498\times10^{-10}) | 0.0004336 | (2.71842\times10^{-11}) | 0.0001752 |

Fixed-coordinate samples make the distinction especially clear:

| (r) | (N) | Isolated RHS term (\mathcal L_S[S]D_1S) | Pointwise (Q_r^S) | FV (Q_r^S) |
|---:|---:|---:|---:|---:|
| 2 | 60 | (4.91731\times10^{-7}) | (1.79071\times10^{-7}) | (4.51275\times10^{-9}) |
| 2 | 120 | (4.85866\times10^{-7}) | (1.29268\times10^{-8}) | (2.82897\times10^{-10}) |
| 2 | 240 | (4.91780\times10^{-7}) | (7.98422\times10^{-10}) | (1.71462\times10^{-11}) |
| 5 | 60 | (3.73132\times10^{-7}) | (-8.07318\times10^{-9}) | (7.90633\times10^{-9}) |
| 5 | 120 | (3.74392\times10^{-7}) | (-7.69929\times10^{-10}) | (5.38652\times10^{-10}) |
| 5 | 240 | (3.74244\times10^{-7}) | (-4.91024\times10^{-11}) | (3.40465\times10^{-11}) |
| 8 | 60 | (3.66119\times10^{-8}) | (-1.94645\times10^{-9}) | (3.88405\times10^{-9}) |
| 8 | 120 | (2.96932\times10^{-8}) | (-8.63048\times10^{-11}) | (2.40370\times10^{-10}) |
| 8 | 240 | (2.79714\times10^{-8}) | (-4.70553\times10^{-12}) | (1.49531\times10^{-11}) |

The ratios are not a constant factor across the annulus or resolutions. At fixed (r=2), both pointwise and FV (Q_r^S) decrease by about fourth order under refinement, consistent with the previously recorded stress-divergence crosswalk. The evidence points away from a simple (4\pi/8\pi)-type normalization error. It does **not** authorize choosing FV over pointwise solely because the FV residual is smaller.

## 5. Disposition and next discriminating check

1. The origin-volume suspicion is resolved: the observer reconstructs regular (r^2g(r^2)) face profiles and Simpson-integrates the α-weighted spacetime measure; it does not use an arithmetic face-weight average or centre value times Δr.
2. The pinned S-sector RHS spatial terms reconstruct exactly from the existing current-tangent ledger. That isolated term is not the full (Q_r^S), which also contains time-density, radial pressure-flux, and metric-source contributions.
3. Keep pointwise and FV divergence results separate and convention-labelled. Next, compare each with a fully reconstructed *scalar equation-of-motion residual* assembled from all signed RHS terms and the time-varying measure/current contributions, then integrate the matching identity over shells. Do not select a divergence form based solely on the smaller residual.
4. The full matter-coupled Brown momentum-propagation identity, total stress-energy sector conservation, the archive-locked φ/DM exchange admission gate, and scientific admission all remain open.

**Disposition:** diagnostic only; production equations and vendor sources unchanged; `scientific_admission=false`.
