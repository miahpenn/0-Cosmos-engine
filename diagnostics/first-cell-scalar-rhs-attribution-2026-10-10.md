# First-cell scalar RHS attribution diagnostic

Date: 2026-10-10
Branch: \`diag/discrete-consistent-initial-review\`

## Why this diagnostic is next

The candidate initializer's central Hamiltonian residual at t=3 decreases from about -6.54e-5 (N=40) to -1.124e-5 (N=80) to -2.289e-6 (N=160), while its t=0 residual is at the Newton-solver residual scale. The next question is not whether the initializer should be admitted; it is how the residual begins to reappear during the actual production step.

## Scope

- Uses the opt-in discrete-consistent candidate initializer at N=160, r_max=40, S amplitude 0.01, width 7, D amplitude 1e-10, radiation enabled, CFL 0.03.
- Evolves the unchanged production kernel to t=3.
- At t=0, 0.75, 1.5, and 2.25, records instantaneous cell-0 scalar RHS values (S, PS, D, PD, phi, Pi) and conservative fluid RHS values (dark matter, baryons, radiation).
- For the next accepted step at each probe, computes a symmetric Shapley attribution over 11 state groups: metric, lapse/shift/gauge, A_a, K, Lambda, S/PS, D/PD, phi/Pi, dark matter, baryons, and radiation.
- Verifies finiteness, candidate Newton convergence, first-slice Hamiltonian closure, and component-sum attribution closure.

## Guardrails

This is a numerical diagnostic only. The group contributions describe the discrete accepted step and are not proof of unique physical causality. No evolution equation, source term, gauge choice, initializer default, damping, floor, clamp, or fitted coefficient is changed. Nothing is admitted to production from this result alone.

## Expected artifact

\`runs/first-cell-scalar-rhs-diagnostic/report.json\`, uploaded as \`first-cell-scalar-rhs-diagnostic\`.
