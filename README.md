# 0-Cosmos Engine

Canonical writable workspace for the interconnected GEAR/COSMOS research engine.

## Objective
Build one reproducible machine in which the cosmic lane, local Einstein-scalar lane, dynamically generated proper-time interface, stress-energy/flux ledger, and cycle observer can be evolved and tested without carrying the full machine in chat.

## Physics locks
- GEAR-03: kappa=0, g=1; omega_S^2=1 and alpha_D^2=1.
- GEAR-63/65 local cubic geometric barrier is retained.
- G=1.
- Corrected COSMOS potential normalization is retained; legacy V0=1 is a control only.
- No fitted feedback coefficient.
- No imposed clock conversion.
- No D-to-matter conversion law.
- No artificial bounce, stop, reset, ejection threshold, or lifetime threshold.
- Reciprocal coupling must come from exposed stress-energy, geometric flux, local energy/transport, and COSMOS ledger quantities.
- A reduced diagnostic source is never promoted to a physical source without conservation and resolution tests.
- The exact spherical branch does not isotropize the archive's Bianchi-I shear.

## Architecture

### Physics and stress-energy
- engine/cosmos.py — cosmic background equations and corrected component budget.
- engine/local.py — GEAR-03 S/D sector and GEAR-63/65 barrier.
- engine/unified.py — unified scalar stress-energy and Misner–Sharp identities.
- engine/matter_closures.py — kinematic fluid stress projections.
- engine/matter_evolution.py — conservative/primitive Valencia identities, DM four-force exchange, baryon dust, and radiation EOS.
- engine/spherical_scope.py — explicit spherical-sector scope; shear cannot be silently enabled.

### Geometry and evolution
- engine/reference_pirk_unified.py — archived unified V5.5 spherical reference kernel.
- engine/true_cmc_reference.py — stage-aware CMC reference evolution.
- engine/production_contract.py — required capabilities for the final horizon-penetrating production kernel.

### Coupling and observables
- engine/bidirectional_production.py — derived Misner–Sharp/stress-energy ledger.
- engine/handoff.py — dynamical handoff and H=0 cycle-event observation.
- engine/coupled_engine.py — continuous one-metric coupled driver.
- engine/interface.py — stress-energy/flux bookkeeping.
- engine/coupled.py — orchestration state.
- engine/cycle.py — cycle event observation.
- engine/invariant_diagnostics.py — invariant trapping/current/constraint witness interfaces.

### Controls and final campaign
- engine/diagnostics.py — cheap algebraic/conservation diagnostics.
- tests/ — final controls and regression suite; not executed during construction.
- docs/ — locked equations, status, and decision records.
- runs/ — generated outputs/checkpoints.

## Scientific status
The repository is still under construction. The CMC/reference driver remains a reference layer, not the final production strong-field kernel. The final engine must connect the frozen V5.5 stress-energy to a horizon-penetrating reference-metric BSSN/PIRK implementation, complete DM momentum, baryon momentum, and radiation evolution, then run the full campaign.

No completed multi-cycle or physical-bounce claim is made until that final campaign is passed.
