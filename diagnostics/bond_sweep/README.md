# Paired provisional/self-consistent S-D bond sweep (diagnostic only)

This branch is isolated from `main`. No production source file is changed.

## Paired test
- Grid: N=40, Rmax=40
- CFL: 0.0075
- Checkpoints: initial slice, t=2,3,4,5
- Initial data: A=0.01, width=7, D amplitude=1e-4, radiation enabled\n- Initial metric: exact-origin cumulative Simpson radial quadrature, followed by the existing discrete-consistent initializer solving the production H constraint at t=0.
- Couplings: g=0, +0.05, -0.05, +0.10
- Candidate potential: V_bond = g S D (diagnostic choice, not derived from the archive)

Each case is independently initialized from the same scalar profiles and numerical parameters. The radial metric estimate uses cumulative Simpson quadrature with the exact source value at r=0, then the existing discrete-consistent initializer solves the unchanged production Hamiltonian constraint on that model's initial slice. The `rhs_only` model omits bond stress from the initial radial source, lapse, and constraints. The `self_consistent` model includes the same candidate stress in the initial radial metric source, lapse, geometry consumer, and constraint-diagnostic consumer. Thus the initial geometry is allowed to differ between stress prescriptions rather than forcing identical metric data that would violate one model's initial constraint. Both H-residual definitions are still computed on every exact saved state.

## Closure checks
The in-process wrapper adds:
- PS_t <- PS_t - alpha g D
- PD_t <- PD_t - alpha g S

In the self-consistent mode, its stress projection is:
- delta rho = +g S D
- delta pr = -g S D
- delta pt = -g S D
- delta j = 0

At the initial slice and t=2,3,4,5, the script reconstructs the unbonded scalar stress independently, verifies the expected bond contribution in both geometry/CMC and constraint-diagnostic consumers, and fails closed if either per-consumer projection error or cross-consumer total-source error exceeds 5e-12 (absolute max norm). Non-finite projection checks also fail closed.

For every state, it calculates both Hamiltonian residual definitions on that exact same slice: with bond stress omitted and with the bond stress included. Paired trajectory deltas are reported separately, so the algebraic stress-definition change cannot be confused with evolution/backreaction changes.

## Interpretation gate
The initial constraint closure is checked separately from the paired evolution. This tests the consequences of a *chosen* bond form; it does not establish that the archive requires this term or that it is physical. The earlier provisional commit/branch is not present in current repository refs, so its values cannot be treated as an apples-to-apples comparison. These same-revision paired runs replace that missing control. Hamiltonian/momentum residuals are evolved diagnostics, not automatically zero constraints.

Outputs: `bond_sweep_mode_pair.json`, `bond_sweep_mode_pair.csv`, `bond_sweep_mode_pair_deltas.json`, `bond_sweep_mode_pair_deltas.csv`.

This paired comparison is run from the isolated diagnostic branch using the `[bond-sweep-self-consistent]` workflow gate.
