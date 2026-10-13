# Initial radial quadrature + CMC centre-row candidate — 13 October 2026

## Scope and provenance

This isolated candidate is stacked on the strict campaign/radiation integration source `9f0bd16d403ada22a4f21413048c93b37c9924a7` with vendor submodule `d6052d605673ce9d82cdc99b0df79855fa2ea215`. It changes only the numerical initial-data quadrature in `engine/v55_initial.py` and the first finite-difference row in `engine/cmc_gauge.py`, with focused regressions. It does not change the physical source equations, matter model, parameters, radiation cone, RHS or vendor source.

## Why both layers are included

Exact-source diagnostics separated two issues:

- With the default production initializer, the CMC alternatives regularize the even `alpha(r²)` fit and reduce the maximum `|K_rhs-target|` from about `1.58e-3` to `6.3e-6` at N=160, but the centre `Aa_rhs/r0²` remains about `1.67e-2` and the Hamiltonian witness remains about `4.0e-3`. The row-only change is therefore not a complete repair.
- With the diagnostic high-order radial integral, the same CMC row alternatives greatly suppress the centre `Aa_rhs/r0²` and improve the immediate M/connection witnesses. This indicates the initial constrained radial integration is a second, independent numerical layer.

## Candidate numerical changes

1. Replace the origin half-cell trapezoid with the analytic even-series weights `(401/720, -31/480, 11/1440)`, exact for a smooth even source through the (r^4) term.
2. Use integrated cubic-Lagrange interval weights, including an even-ghost first interval and a one-sided final interval, for cumulative radial integration.
3. Replace `alpha[0]=alpha[1]` in the CMC banded system with the discrete CMC PDE row obtained by even reflection at the origin. Interior rows and the outer `alpha[-1]=1` normalization remain unchanged.

## Admission boundary

This is a candidate, not a production promotion. Focused polynomial tests and the repository suite must pass, followed by a bounded full-machine step and review of centre-series/H/M/connection, radiation admissibility, Misner–Sharp conservation, matched spatial refinement, and the pinned archived-state radiation continuation. A short smoke cannot admit the full trajectory. No long campaign should be launched until the cross-lane gates are reconciled.
