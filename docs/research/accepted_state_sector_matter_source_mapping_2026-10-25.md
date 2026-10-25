# Accepted-state sector matter-source mapping — 2026-10-25

**Status:** research diagnostic only. This report does not promote a source law, a discrete divergence convention, or a scientific-admission decision.

**Execution:** [workflow run 40159027760](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40159027760), all N=60/120/240 jobs passed.  
**Review surface:** [draft PR #420](https://github.com/miahpenn/0-Cosmos-engine/pull/420).  
**Accepted state:** t=0.60, source-pinned production-baseline, partial-joint, and full-continuous-profile arms. The observer advances no evolution step.  
**Provenance pins:** source 4aec88d0e3ff028e0a65a8ac36b30948625df792; vendor d6052d605673ce9d82cdc99b0df79855fa2ea215.

## 1. Per-sector additivity closes

For each of S, D, COSMOS phi, dark matter, baryons, and radiation, the observer separately records the existing geometric matter-injection contribution, the signed current tangent, and the pointwise covariant stress divergence. It then compares the sector sum against the existing total ledger.

| Additivity check, maximum absolute gap | N=60 | N=120 | N=240 |
|---|---:|---:|---:|
| Sum of sector geometric injections vs total injection | 1.20e-19 | 8.30e-20 | 2.26e-19 |
| Sum of sector current tangents vs total composite minus injection | 1.54e-16 | 5.43e-16 | 4.16e-16 |
| Sum of sector divergences vs pointwise total divergence | 1.39e-17 | 3.78e-17 | 3.03e-17 |
| Sum of sector mapping residuals vs pointwise-total residual | 3.32e-16 | 7.99e-16 | 7.97e-16 |

These establish that the diagnostic decomposition is additive. The pointwise total is deliberately not equated to the expanded Q_r candidate: its separate discrete radial-flux product-rule gap has maximum absolute size 7.95e-7, 1.13e-7, and 1.46e-8 at N=60/120/240. That difference is retained as its own operator comparison.

## 2. The provisional source-mapping residual is dominated by S

For the full-continuous-profile arm, the first two cell-centred source-mapping residuals are:

| N | S, cell 0 | S, cell 1 | COSMOS phi, cell 0 | Dark matter, cell 0 |
|---:|---:|---:|---:|---:|
| 60 | +4.232e-5 | +6.279e-5 | -2.603e-8 | +1.157e-9 |
| 120 | +1.425e-5 | +3.390e-5 | -3.127e-9 | +1.034e-9 |
| 240 | +6.096e-6 | +1.715e-5 | -6.504e-11 | +1.566e-9 |

The D contribution is around 1e-20 to 1e-21 in these central cells; baryon and radiation pieces are much smaller than S. The residual in the existing provisional relation between geometric/current momentum-constraint sources and the stress-divergence source is therefore S-dominated.

This is not a finding that the physical S sector violates its own conservation equation. It shows that the current candidate source mapping does not yet account for all terms needed to relate the geometric momentum constraint to the matter stress divergence. The normalized Hamiltonian-gradient term, reference-metric connection response, and scalar spatial-operator convention remain open.

## 3. The S mapping residual persists at fixed physical radii

For the same full-profile arm, the S-sector mapping residual at fixed radii is:

| N | r=2 | r=5 | r=8 |
|---:|---:|---:|---:|
| 60 | +6.279e-5 | +6.166e-5 | +2.151e-5 |
| 120 | +5.813e-5 | +6.250e-5 | +2.094e-5 |
| 240 | +5.836e-5 | +6.245e-5 | +2.079e-5 |

The same residual is nearly unchanged at these fixed physical radii in baseline, partial-joint, and full-profile arms. Thus the full metric profile improves near-centre regularity and the one-step connection-constraint rate, but does not resolve this source-mapping residual at finite radius.

At r=5 and N=240, the S-sector geometric matter injection is approximately +8.536e-6, its current tangent is +5.392e-5, and their composite is +6.246e-5. The pointwise stress-divergence source is about +1.233e-9. This difference is much larger than the separately isolated pointwise-versus-expanded Q_r gap at that radius; it requires the complete source-aware constraint-propagation identity, not a relabeling of the divergence array.

## 4. The pointwise S-sector stress divergence itself converges down

The pointwise S-sector covariant divergence, whose archival exchange source is zero for the uncoupled S sector, decreases strongly at fixed physical radii:

| Radius | N=60 | N=120 | N=240 |
|---|---:|---:|---:|
| r=2 | +1.791e-7 | +1.293e-8 | +7.984e-10 |
| r=5 | -8.073e-9 | -7.699e-10 | -4.910e-11 |
| r=8 | -1.946e-9 | -8.630e-11 | -4.706e-12 |

That refinement trend contrasts with the nearly constant S source-mapping residual at r=2, 5, and 8. The measured sector divergence and the provisional geometric/current source balance are different witnesses; they should not be conflated.

## 5. Conditional comparison with the archive-locked DM/phi exchange

The master checkpoint and the explicit locked-sign statement in the V5.5 covariant gate give the exchange convention

\[
\nabla_\mu T_{\rm DM}^{\mu\nu}=+\beta_{\rm DM}\rho_{\rm DM}\nabla^\nu\phi,\qquad
\nabla_\mu T_{\phi}^{\mu\nu}=-\beta_{\rm DM}\rho_{\rm DM}\nabla^\nu\phi.
\]

This comparison uses the production coupling and normalized, independently reconstructed sector stress projections on the same saved accepted state. It remains a diagnostic comparator, not a fitted source term.

**Dark matter:** the measured pointwise Q_r closely tracks the locked exchange source at fixed radii. The measured-minus-expected residuals at r=2 are +4.40e-15, +1.05e-15, +1.30e-16; at r=5 they are -8.04e-15, -6.43e-15, -2.23e-15; at r=8 they are +3.03e-15, +6.39e-16, +1.62e-16, for N=60/120/240 respectively.

**COSMOS phi:** the absolute residual decreases with refinement but is not yet a closure pass over all radii at N<=240. At r=2 it is -2.422e-10, -1.764e-11, -1.088e-12; at r=5 it is +1.956e-11, +1.639e-12, +1.045e-13; at r=8 it is +6.290e-12, +3.421e-13, +2.029e-14.

The archive gate has an earlier sentence that conflicts with its explicit locked-sign statement. This comparison follows the latter, corroborated by the master checkpoint; the contradictory earlier sentence remains visible as archive provenance debt until the canonical record is reconciled. Nothing in this experiment edits that record or changes production physics.

## 6. Disposition

The result separates three facts:

1. The per-sector geometric/current/divergence decomposition is additive to tight numerical tolerance.
2. The S-sector pointwise stress divergence falls rapidly with refinement, while the provisional S geometric/current source-mapping residual remains around 5.8e-5 to 6.25e-5 at fixed r=2 and r=5.
3. The archive-locked dark-matter divergence comparator tracks closely; the COSMOS phi comparator improves with refinement but is not yet quantitatively closed.

The next step is to reconcile the convention-pinned Brown Eq. 11b scaffold with the reference-metric connection response, normalized Hamiltonian-gradient term, scalar spatial operator, and complete sector exchange identity. Do not call the S mapping residual a physical conservation violation or edit a production term until that identity is independently derived and closed.

**Guardrails:** source/vendor pins match; accepted-state artifacts are SHA-verified and immutable during the observer; sector additivity closes; no production, initializer, engine, vendor, or dependency files changed; full matter-coupled momentum propagation remains open; total stress-energy sector conservation remains open; scientific_admission=false.
