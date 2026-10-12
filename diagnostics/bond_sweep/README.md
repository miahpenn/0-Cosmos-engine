# Self-consistent S-D bond sweep (diagnostic only)

This branch is isolated from `main`. No production source file is changed.

## Test
- Grid: N=40, Rmax=40
- CFL: 0.0075
- Final coordinate time: t=5
- Initial data: A=0.01, width=7, D amplitude=1e-4, radiation enabled
- Couplings: g=0, +0.01, -0.01, +0.05, -0.05, +0.10
- Candidate potential: V_bond = g S D (a diagnostic choice, not derived from the archive)

## Self-consistency closure
The in-process wrapper adds the potential gradients to the scalar momentum RHS:
- PS_t <- PS_t - alpha g D
- PD_t <- PD_t - alpha g S

It adds the potential's stress projections in both independent consumers:
- delta rho = +g S D
- delta pr = -g S D
- delta pt = -g S D
- delta j = 0

At the initial slice and t=2,3,4,5, it independently reconstructs the original scalar projection and fails closed unless (a) both consumers have precisely the expected bond contribution and (b) the total geometry/lapse source projection agrees with the constraint-diagnostic projection to absolute max error <= 5e-12. The initial slice is logged explicitly so the baseline constraint offset is visible, rather than hidden by evolved-only samples.

## Interpretation gate
This run tests the internally consistent consequences of the *chosen* bond form. It does not establish that the archive requires this bond or that it is physical. A prior provisional result was reported with a different baseline, but its cited commit/branch is not present in the current repository refs; therefore the old and current H-center numbers are not yet an apples-to-apples comparison. Read Hamiltonian/momentum residuals as evolved diagnostics, not as automatically zero constraints. Do not admit the bond as physics without an independently justified interaction and convergence/closure checks.

Outputs: `bond_sweep_self_consistent.json`, `bond_sweep_self_consistent.csv`.
