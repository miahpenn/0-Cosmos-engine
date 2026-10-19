# Issue #4 — shared-state momentum / total-stress calculation

Date: 2026-10-19  
Status: implementation committed; exact-head CI execution and numerical results pending.

## Implementation

- PR: [#408](https://github.com/miahpenn/0-Cosmos-engine/pull/408), draft, diagnostic-only.
- Shared-state adapter: [audit_shared_state_momentum_total_stress_divergence.py](https://github.com/miahpenn/0-Cosmos-engine/blob/diagnostic/accepted-state-raw-momentum-witness/tools/audit_shared_state_momentum_total_stress_divergence.py).
- Workflow: [accepted-state-raw-momentum-witness.yml](https://github.com/miahpenn/0-Cosmos-engine/blob/diagnostic/accepted-state-raw-momentum-witness/.github/workflows/accepted-state-raw-momentum-witness.yml).
- Production source pin: `7f47669a1c0edfafef9fc5f27b2db80840d3f035`.
- Vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Initial digest: `a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0`.
- Accepted digest: `c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a`.

The adapter advances the pinned production kernel exactly one step at resolution 60, then constructs one RHS-direction plus/minus pair at epsilon `1e-5`. It passes those same in-memory state objects to the signed momentum-constraint observer and the independent total-stress-energy divergence observer. It emits signed cell-centred and finite-volume divergence arrays plus the signed momentum tangent, and records state digests before/after observation.

A wrong-digest negative control must be rejected. Source/vendor provenance and both expected state digests fail closed.

## Scope and limits

This is the first shared-state calculation, not yet the full matter-coupled momentum-propagation identity. It reports residuals; it does not assert that they close. The existing six-sector observer still excludes charged dust and radial Maxwell and therefore cannot be promoted to an all-sector oracle. Brown's existing propagation audit remains vacuum-only.

No production equations, source/vendor pins, numerical bounds, or scientific-admission rules changed. `scientific_admission=false` and `acceptance_oracle_satisfied=false` remain mandatory. PR #408 stays draft until the exact-head workflow runs, signed results and controls are inspected, and the complete independent propagation RHS is assembled.
