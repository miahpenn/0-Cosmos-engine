# Self-consistent S-D bond sweep (diagnostic only)

This branch is isolated from `main`. No production source file is changed.

## Test
- Grid: N=40, Rmax=40
- CFL: 0.0075
- Final coordinate time: t=5
- Initial data: A=0.01, width=7, D amplitude=1e-4, radiation enabled
- Couplings: g=0, +0.01, -0.01, +0.05, -0.05, +0.10
- Candidate potential: V_bond = g S D (a diagnostic choice, not derived from the archive)

## Self-consistency repair
The in-process wrapper adds the potential gradients to the scalar momentum RHS:
- PS_t <- PS_t - alpha g D
- PD_t <- PD_t - alpha g S

It also adds the potential's stress projections everywhere the scalar stress projection is consumed:
- delta rho = +g S D
- delta pr = -g S D
- delta pt = -g S D
- delta j = 0

This includes the total stress assembly used by constraint diagnostics and the scalar projection used by the CMC solve / matter geometry source. The base `scalar_system.py` and production files are untouched.

## Interpretation gate
This run tests the internally consistent consequences of the *chosen* bond form. It does not establish that the archive requires this bond or that it is physical. Compare g=0 against the earlier provisional baseline first; differences may indicate source revision or wrapper effects. Read Hamiltonian/momentum residuals as evolved diagnostics, not as automatically zero constraints. No bond result should be admitted as physics without an independently justified interaction and convergence/closure checks.

Outputs: `bond_sweep_self_consistent.json`, `bond_sweep_self_consistent.csv`.
