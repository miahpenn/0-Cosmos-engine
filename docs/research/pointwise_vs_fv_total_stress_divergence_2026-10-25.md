# Pointwise versus finite-volume total-stress divergence on accepted full-profile states — 2026-10-25

**Scope:** post-processing of the already accepted t=0.60 states; no evolution rerun. This records a comparison between two distinct diagnostics from the existing sector-resolved stress-energy observer. It does not promote either discrete form as the physical authority.

**Source run:** [40157745275](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40157745275), with N=60/120/240 artifacts 12585232569, 12585182574 and 12585666768. Source and vendor pins match the frozen research anchor, protected source diffs are empty, and scientific admission remains false.

## Full continuous-source arm

The pointwise observer differentiates the cell-centred covariant flux and applies the pointwise metric-source terms. The finite-volume observer reconstructs spherical face fluxes and cell-integrates the time-density and metric-source terms before dividing by cell measure. These are not identical discrete quantities; the second is a cell-averaged integrated balance.

At the origin cell, for the first directional epsilon \(\epsilon=2\times10^{-5}\):

| N | Pointwise \(\nabla_\mu T^\mu{}_r\) | FV cell-averaged divergence | Pointwise/FV magnitude ratio |
|---:|---:|---:|---:|
| 60 | +7.91543e-7 | +6.68940e-9 | 118.3 |
| 120 | +1.12839e-7 | +9.00292e-10 | 125.3 |
| 240 | +1.45780e-8 | +1.12866e-10 | 129.2 |

The result is stable to the directional epsilon: changing \(\epsilon\) from \(2\times10^{-5}\) to \(10^{-5}\) changes these origin-cell values only at roundoff-scale relative differences (roughly \(10^{-7}\) or less).

The origin-cell values decrease by factors 7.01 and 7.74 for the pointwise form and 7.43 and 7.98 for the FV form on successive refinements. The FV origin residual therefore scales approximately as \(h^3\) in this range.

On the fixed physical annulus \(2\le r\le10\):

| N | Pointwise divergence RMS | FV divergence RMS |
|---:|---:|---:|
| 60 | 7.48087e-8 | 5.97270e-9 |
| 120 | 2.85680e-9 | 3.78225e-10 |
| 240 | 1.88120e-10 | 2.38704e-11 |

The FV annulus RMS falls by factors 15.79 and 15.84, approximately fourth order under each grid halving. The pointwise annulus RMS also decreases strongly, with a higher-resolution observed order near four.

## Sector attribution and interpretation

The S scalar sector dominates the pointwise origin residual: its contribution is +7.92656e-7, +1.13000e-7 and +1.45989e-8 at N=60/120/240, closely tracking the total. This is consistent with the previously measured radial-flux product-rule gap in the scalar spatial operator. It is not a fluid-HLL issue, and the fluid HLL pairing does not define the scalar S-sector spatial flux.

The FV result is numerically much smaller because it evaluates a cell-integrated balance with reconstructed spherical face fluxes and integrated source terms; the pointwise result evaluates a local derivative expression. The size difference is a reason to keep the representations distinct, not by itself proof that one is the authoritative physical convention.

## Convention gate

The archive-backed unified first-law record establishes the importance of an independent integrated current balance, but its areal-radius formula cannot be transplanted directly into the present isotropic/BSSN radial coordinate without the correct areal-radius and measure transformation. The production source convention, the scalar sector's covariant exchange, and the full matter-coupled momentum propagation identity remain unclosed.

Next: document the coordinate/measure map for the present BSSN state and compare the shell-integrated total-current identity with these cell-integrated residuals. Until that map is frozen from the archive and implementation, keep pointwise and FV forms as separate, convention-labelled diagnostics.

No production equations, initializer, vendor or dependencies changed. \`scientific_admission=false\`.
