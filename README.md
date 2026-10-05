# 0-Cosmos Engine

Canonical writable workspace for the interconnected GEAR/COSMOS research engine.

## Objective
Build one reproducible machine in which the cosmic lane, local Einstein-scalar
lane, dynamically generated proper-time interface, stress-energy/flux ledger,
and cycle observer can evolve together without carrying the full machine in chat.

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
- tests/ — regression suite; not executed during construction.
- docs/ — provenance, locks, status, and cycle-coordinate decisions.
- runs/ — generated campaign checkpoints/ledgers.

## Scientific status
The production state graph is implemented but NOT yet campaign-validated.
The next and final engineering gate is a complete repository test plus the
long synchronized numerical campaign at multiple resolutions.

No completed multi-cycle or physical-bounce claim is made until that campaign
actually runs and the invariant witnesses support it.