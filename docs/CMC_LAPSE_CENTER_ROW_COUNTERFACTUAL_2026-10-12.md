# CMC lapse centre-row counterfactual — 2026-10-12

**Status:** completed diagnostic-only counterfactual; no production patch authorized.  
**Production source:** \`6d8b3f2790644591326df57c55908863f994cb0c\`  
**Workflow:** https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38280503385  
**Artifact:** \`cmc-lapse-center-row-counterfactual\`, ID \`11760957038\`  
**Downloaded ZIP SHA-256:** \`ae1226c2f0e2da7bd3b82c876e8c9f84d2eb9bd8c22b7a9d2306f9cd7edd61d7\`  
**Contained JSON SHA-256:** \`02fb38a3b574c1d5eb96ec37937f2251e4eed5b63e03b61bf432795abcb6f038\`

## What was compared

Ten frozen states were tested: N=40/60/80/120 at initial time, saved N=40/60/80/120 states at t≈2.16, and saved N=60/80 states near t≈10. Each state compared three lapse solves:

1. Current second-order solve with current centre row.
2. Second-order solve using a cell-centred even ghost.
3. Production-matched fourth-order operator solve over all rows.

The fixed CMC target was held constant for the residual comparison; target refreshes were recorded separately. No trajectory was advanced, and no production source or initial state was modified.

## Findings

- The current lapse solve has a large first-cell Kdot-minus-target residual: roughly \`−1.44e−3\` to \`−1.57e−3\` on initial N=40–120 states, and \`+1.23e−3\` to \`+1.35e−3\` on t≈2.16 states. Near t≈10 it is about \`+1.50e−3\` for N=60/80. The largest residual is at the first cell; outer \`r>=5\` residuals are much smaller.
- The production-matched fourth-order solve reduces the pointwise Kdot residual to floating-point roundoff (\`~1e−15\` or less). **This is not independent validation**: its matrix is assembled from the same discrete operators used by the production K RHS, so this cancellation is expected by construction.
- The second-order solve with a cell-centred even ghost reduces the first-cell residual compared with the current centre row, but does not eliminate it. Initial first-cell residuals range from about \`1.12e−5\` (N=120) to \`8.49e−5\` (N=40); t≈2.16 residuals range from \`−3.98e−7\` (N=120) to about \`−5.87e−5\` (N=40). At t≈10, first-cell residuals are about \`−2.77e−6\` (N=80) and \`−2.63e−5\` (N=60).
- The alternate lapse profiles differ most at the first cell. The production-matched solve changes first-cell alpha by up to \`2.04e−3\` in the tested states; the even-ghost variant changes it by up to \`2.42e−3\`. These are diagnostic counterfactual changes, not approved production values.
- Recomputing the target from each alternate lapse gives close but not identical values to the target from the current lapse. For example, at N=40 initial time the current refresh is \`6.1959551e−5\`, while the production-matched refresh is \`6.1973859e−5\`. This distinction should remain explicit in future comparisons.

## Interpretation and limitations

The data support a real centre-row / discrete-operator compatibility issue in the CMC lapse solve as a contributor to the central Kdot residual. They do **not** establish that it is the sole cause of the evolving Hamiltonian, momentum, and connection residuals. The production-matched solve's tiny residual is partly tautological and must not be presented as an independent convergence result. The even-ghost result is a counterfactual only; no modified-lapse trajectory was run, so it has not passed one-step propagation, short-run H/M/C, or Misner–Sharp balance gates.

## Next action

Use the current centre-row, even-ghost, and production-matched operators in one explicitly matched, frozen-input RHS identity check. Compare the signed CMC Kdot terms and the centre-series rates at identical inputs; report both the fixed-target residual and refreshed-target result. Only if one specific discrete row/term fails an independent identity should a one-step numerical-only counterfactual be considered. Keep production physics, sources, scalar/matter equations, gauge targets, and boundary conditions locked. Do not use damping, floors, clamps, fitted centre corrections, or an invented physical source. Do not launch the t=160 campaign until t=12 constraint and balance admission gates pass.

## Cross-links

- Centre regularity and constraint-propagation handoff: \`docs/CONTINUATION_HANDOFF_CENTRE_REGULARITY_2026-10-12.md\`
- Centre-series curvature RHS lead: \`docs/CENTRE_SERIES_CURVATURE_RHS_LEAD_2026-10-12.md\`
- CMC lapse investigation workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38280503385
