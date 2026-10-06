# 0-Cosmos Engine

Canonical writable workspace for the interconnected GEAR/COSMOS research engine.

## Objective

Build one reproducible machine in which the cosmic lane, local Einstein-scalar lane, dynamically generated proper-time interface, stress-energy/flux ledger, and cycle observer evolve together without carrying the full machine in chat.

The governing engineering rule is simple:

**run the complete machine → expose the actual defect → repair only that layer → rerun.**

## Physics locks

- GEAR-03: kappa=0, g=1; omega_S^2=1 and alpha_D^2=1.
- GEAR-63/65 local cubic geometric barrier is retained.
- G=1.
- Corrected COSMOS potential normalization is retained; legacy V0=1 is a control only.
- No fitted feedback coefficient.
- No imposed clock conversion.
- No D-to-matter conversion law.
- No artificial bounce, stop, reset, ejection threshold, or lifetime threshold.
- Reciprocal coupling comes only from the unified stress-energy/geometry system.
- A reduced diagnostic source is never promoted to physical coupling without conservation and resolution evidence.
- The exact spherical branch does not isotropize the archive's Bianchi-I shear.
- lambda_cycle is not invented; only archive-supported phase, phi, H_eff, proper time, and derived e-fold coordinates are carried.

## 0* diagnostic rule

0* is a **quantification lens, not a construction rule**.

Authoritative scorecard:

- 0*_eff iff t and b and (w >= 1)
- 0*_smooth iff t and b and (w > 1)
- 0*_global iff 0*_eff and G

The atomic folds are measured from the evolved system rather than imposed by the gauge:

- t: turnaround, H=0 and dot H<0.
- b: bounce, H=0 and dot H>0.
- w>=1: anti-BKL preservation condition; w>1 buys active shear dilution.
- G: global reset/timescale condition.

For contraction, the archive scorecard uses N=ln(a) rather than treating coordinate time as the invariant contraction clock. Diagnostics must not manufacture rho_tot to pay the criterion; use the prescribed invariant quantities instead.

**Important:** 0* does not prescribe a lapse or boundary condition. Any custom gauge must first be independently justified by the archive equations and then run through the 0* scorecard.

## Current time/gauge investigation

The current production CMC gauge uses the archive-supported outer normalization alpha(R)=1. A diagnostic branch is testing the archived central-clock relation d tau = alpha(0) dt as a **consistent coordinate transformation**.

The key requirement now under test is:

> A lapse rescaling is a coordinate transformation only when the entire evolution/time interface is transformed consistently with the new time parameter.

The first central-lapse-only experiment failed before turnaround, at an inner radiation cell. That result is not being interpreted as a physical failure of central proper time; it exposed that the engine contains time-dependent components whose clock bookkeeping must be transformed together.

The current experiment therefore tests the time-interface layer itself, with no new physical source, coefficient, boundary condition, or bounce rule.

## Architecture

### Physics and stress-energy

- engine/scalar_system.py — S/D/COSMOS scalar evolution and projections.
- engine/cosmos.py — archive homogeneous COSMOS lane, signed-H evolution, beta exchange, radiation, and homogeneous shear semantics.
- engine/local.py — GEAR-03 S/D sector and GEAR-63/65 barrier.
- engine/valencia.py — metric-aware conservative/primitive identities and exact mixed-tensor source contraction.
- engine/matter_system.py — conservative spherical dust/radiation transport and DM four-force.
- engine/v55_matter.py — corrected V5.5 matter bundle and archive operating point.
- engine/stress_energy.py — single total Einstein source assembly.
- engine/spherical_scope.py — explicit exact-spherical matter scope.

### Geometry and evolution

- vendor/bb-palatini-unified-r0 — pinned reference-metric spherical BSSN/PIRK numerical kernel.
- engine/v55_pirk_adapter.py — V5.5 adapter over the pinned geometry kernel.
- engine/v55_initial.py — corrected production initial-data mapping.
- engine/production_kernel.py — synchronized geometry + scalar + conservative-matter state evolution.
- engine/production_contract.py — production capability contract.
- engine/cmc_gauge.py — production CMC gauge and proper-volume Kdot projection.
- engine/true_cmc_reference.py — historical/reference CMC layer; not the production entry point.
- engine/reference_pirk_unified.py — historical archive reference; not the production initializer.

### Coupling and observables

- engine/worldtube.py — invariant areal radius, chi, Misner-Sharp mass, marginal roots, and current residual.
- engine/handoff.py — proper-time/cycle handoff observation.
- engine/cosmology_observables.py — H_eff and derived e-fold coordinate.
- engine/invariant_diagnostics.py — invariant numerical witnesses.
- engine/interface.py — bookkeeping only; no phenomenological source.
- engine/coupled.py — orchestration data structures.
- engine/coupled_engine.py — historical continuous CMC driver; not production.

### Campaign and controls

- engine/campaign.py — long-horizon campaign configuration.
- engine/run_production.py — canonical final campaign runner.
- engine/run_controls.py — cheap control-only harness.
- engine/zero_star_clock_gauge.py — experimental archive-derived central proper-time coordinate; not yet production.
- tests/ — regression suite.
- docs/ — provenance, locks, status, and cycle-coordinate decisions.
- runs/ — generated campaign checkpoints/ledgers.

## Validation status

The machine has passed substantial component and campaign diagnostics, including multi-resolution worldtube-radius tests and boundary-layer isolation. The remaining question is the **time-interface consistency** of the production evolution through turnaround and contraction.

The current experimental branch is deliberately not merged into production.

No completed multi-cycle or physical-bounce claim is made until the full evolution, invariant witnesses, resolution behavior, and 0* scorecard support such a result.