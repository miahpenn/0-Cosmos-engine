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

## Architecture
- engine/cosmos.py — cosmic background dynamics and corrected component budget.
- engine/local.py — GEAR-03 S/D local sector and GEAR-63/65 barrier.
- engine/interface.py — common-clock, stress-energy, flux, and exchange bookkeeping.
- engine/coupled.py — orchestration of bidirectional evolution.
- engine/cycle.py — cycle coordinates and event detection.
- engine/diagnostics.py — constraints, conservation, convergence, and seam diagnostics.
- tests/ — cheap controls before expensive campaigns.
- docs/ — locked equations, status, and decision records.
- runs/ — generated outputs/checkpoints.

## Scientific status
Archived V5.3/V5.4 results remain validation evidence, not silently promoted production physics. The known final seam is a covariant/conserved reciprocal interface: the earlier moving-boundary current is diagnostic until a joint spacetime/interface treatment closes conservation and radial convergence.
