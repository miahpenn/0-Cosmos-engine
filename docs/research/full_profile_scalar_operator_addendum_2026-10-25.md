# Full-profile propagation: scalar-current and discrete-operator addendum — 2026-10-25

**Status:** research diagnostic only. No production equation or convention is accepted by this result.

**PR:** [#420](https://github.com/miahpenn/0-Cosmos-engine/pull/420)  
**Exact-head run:** [40157745275](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40157745275), all N=60/120/240 jobs passed.  
**Accepted state:** t = 0.60, source-pinned trajectory artifacts; tangent observer advances no evolution steps.  
**Provenance pins:** source 4aec88d0e3ff028e0a65a8ac36b30948625df792; vendor d6052d605673ce9d82cdc99b0df79855fa2ea215.

## 1. What was measured

The v5 observer joins four existing diagnostics on the same accepted state / same signed tangent endpoints:

- the existing current-source/Brown Eq. (11b) propagation candidate;
- the existing additive sector-resolved total-T divergence observer;
- the existing production fluid HLL/source operator snapshot;
- the existing scalar RHS term budget, used to reconstruct the S-sector current tangent.

The directional epsilons are 2e-5 and 1e-5. The saved accepted-state inputs were SHA-verified; state observers were immutable; the production/vendor/dependency guard passed. The result remains diagnostic-only.

## 2. The discrete pointwise-versus-expanded gap is a radial-flux product-rule gap

For the full-continuous-profile arm, at cells 0/1:

| N | Pointwise total div(T)^r, cell 0 | Pointwise minus expanded Q_r, cell 0 | Radial-flux product-rule gap, cell 0 |
|---:|---:|---:|---:|
| 60 | +7.9154e-7 | +7.9534e-7 | +7.9534e-7 |
| 120 | +1.1284e-7 | +1.1297e-7 | +1.1297e-7 |
| 240 | +1.4578e-8 | +1.4582e-8 | +1.4582e-8 |

Across all arms, resolutions and epsilon values, the decomposition closes to at most about 2.6e-20 absolute. The time-density product term and metric-chain terms are much smaller than the radial-flux product-rule gap. At N=240 full-profile cell 0, the radial-flux term is +1.4582214e-8; the lapse/metric-chain and areal-chain terms are O(1e-16) or below.

This identifies the dominant *discrete representation difference*: applying the cell-centred derivative to the full radial flux is not equivalent on this grid to splitting the derivative into pressure-gradient and metric-prefactor terms. It does not select which representation is authoritative for the physical ledger.

## 3. The pointwise radial divergence is S-sector dominated

Full-profile pointwise sector-S div(T)^r at cell 0 is:

| N | S-sector div(T)^r |
|---:|---:|
| 60 | +7.9266e-7 |
| 120 | +1.1300e-7 |
| 240 | +1.4599e-8 |

At N=240, the S radial-flux and metric-source terms are individually about -1.4766e-4 and +1.4786e-4; the time-density term is about -1.8978e-7. Those large terms cancel, leaving an O(1e-8) remainder. This is consistent with why pointwise and expanded forms separate on the centre cells, but does not prove either form is the accepted conservation convention.

## 4. The S-current rate reconstructs, and it scarcely changes between metric arms

The existing scalar RHS term budget reconstructs the exact signed S-current contribution
\[
\dot M_S=-8\pi\,\partial_t j_S,\qquad j_S=-P_S D_r S,
\]
term by term. Maximum absolute reconstruction gaps across the arms and both epsilons are:

| N | Max absolute S-current reconstruction gap |
|---:|---:|
| 60 | 8.4e-15 |
| 120 | 2.16e-14 |
| 240 | 2.98e-14 |

For N=240, full-profile cell 0, the largest signed contributions are:

| Existing pinned scalar RHS term | Contribution to S-current tangent |
|---|---:|
| Pdot potential derivative | +1.075e-5 |
| Sdot = alpha P_S lapse-momentum term | -7.261e-6 |
| Pdot lapse-metric-gradient term | +8.530e-7 |
| Pdot lapse-Laplacian term | +4.260e-7 |
| Pdot trace-curvature term | -1.432e-7 |
| Actual reconstructed S-current contribution | +4.621e-6 |

The potential-derivative and lapse-momentum terms substantially cancel. The measured cell-0 S-current contribution is almost unchanged between the partial-joint and full-profile arms at N=240: 4.62000e-6 versus 4.62056e-6. Therefore the full-profile improvement in the regularity witness is not explained by a large change in this source-rate budget.

## 5. The existing fluid-HLL comparison is clean, but not the scalar-sector answer

On these accepted states, the existing production HLL observer reports an exactly zero origin face flux for each fluid species and independently contracted geometric momentum source gaps no larger than about 6.7e-24. Fluid HLL versus projected-stress face fluxes are close to roundoff on these near-rest fluid states. The scalar S/D/phi sectors are not fluid HLL face fluxes, so this observation does not explain or close the S-sector divergence discrepancy.

## 6. Related curved regular finite-volume calibration

The existing test at commit [0a8f184](https://github.com/miahpenn/0-Cosmos-engine/commit/0a8f1847a33a3a0a9f3b74fbf7bcf7794397a75d) uses flat Minkowski space represented in a nonlinear areal-radius chart with static isotropic stress. Its exact-head workflows passed, and its assertions require approximately 8x reduction in the first-cell finite-volume error and 16x reduction in the fixed-physical-annulus RMS when dr is halved. This is a valuable operator calibration, but it is not a production evolved-state identity and does not choose pointwise versus finite-volume divergence as the authoritative convention for the S sector.

## 7. Scientific disposition

The evidence now separates three matters:

1. The initial radial metric profile strongly affects the near-centre regularity witness.
2. The difference between pointwise and expanded Q_r is almost entirely a discrete radial-flux product-rule gap, dominated by S-sector stress.
3. The S-current RHS source budget itself reconstructs to tight tolerance and is nearly the same in the partial and full-profile arms.

Next, use the archive-backed convention/source map to determine which integrated/discrete divergence identity governs this ledger, then reconcile that identity with the actual matter-source composite, the reference-metric connection response, scalar spatial-operator terms, and total sector exchange/conservation. Keep pointwise and FV profiles distinct until that convention choice is documented.

No production/initializer/vendor/dependency files were changed. The current source mapping and complete matter-coupled propagation identity remain unclosed, total sector conservation remains unclosed, and scientific_admission remains false.
