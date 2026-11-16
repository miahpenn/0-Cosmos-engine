# Task 2, research: the zero-flux (lattice-wall) closure against production's open edge (research only)

HANDOFF Nov 15 (3), task 2; one-stone scope confirmed ("identical neighbours, in step: one-stone simplification").
The zero-flux closure already exists as research: `research/lattice_cell_v1.LatticeCellKernel` (W1–W5: mirror
ghosts at the wall, no flux, no boundary data, and K̇ from the wall's compatibility). The CMC lapse is solved with zero
gradient at both ends, α = K̇·u, and K̇ is fixed by normalising α = 1 at the wall ("i") or by ⟨Q⟩ ("ii"). **K̇ never
uses the lapse already in the state.** Here it is compared with production's open edge on the merged machine (#503).
Script: `research/zero_flux_closure_v1.py`. Nothing in `engine/` changed.

## Summary

| test | open edge | wall |
|---|---|---|
| **B** background, production default (N 20, r_max 40, t ≤ 200) | lane to 3e-8 in H and 1e-8 in φ; spurious lapse spread grows from 0 to 2e-10 | lane to 3e-8 and 1e-8; lapse uniform to 1e-14 throughout |
| **S** same, lapse seeded 1% off before the first solve | never recovers: H −2% at t = 10, −13% at 100, −27% at 200; the seeded mode grows 0.95% → 4.0% | **identical to the unseeded run**: it keeps no memory of the lapse |
| **Q** the same pocket initial state for both (cell's conformally flat x_c 2.0 pocket, λ 5000, N 120, r_max 24, same dt) | the centre lapse runs away (1.17 → 1.29 at t = 10 → 23.6 at t = 16), core H falls 65% below the wall's, then the outermost cell fails (t = 16.9) | runs to t = 30; core S–D energy steady (5.07–5.13e2); lapse 0.90–0.94 at the centre, 1.002 at the wall |

- **Derived (W4):** with both ends at zero lapse gradient, the lapse equation fixes only the lapse's shape, and K̇ is one
  number fixed by the normalisation. No step's lapse feeds the next, so the open edge's undamped mode
  (KERNEL_DRIFT.md) cannot exist.
- **Read:** B and S confirm it on production's own background. Q shows the mode matters with matter present. A real
  pocket's lapse is non-uniform, and through the open edge that non-uniformity feeds K̇ and runs away.
- **Caution on Q:** the open edge was built for production's areal initial slice, matched to an exterior. Started
  from the cell's conformally flat slice (ψ′ = 0 at the edge), it is outside what it was made for. With its own areal
  initial data, the same pocket ran past t = 18 with a centre lapse of 0.94 (interim run, not kept). So Q is a
  stress test, not a like-for-like race. B and S are like for like.
- **The wall reflects.** The S–D energy outside r = 12 grows 0.06 → 0.31 (× dV, scaled) by t = 30 in the wall run: outgoing
  ringing and D_r stay in the cell. That is the identical neighbours sending theirs in: correct under the one-stone
  scope, and wrong for an isolated set-up (static and isolated runs keep the open edge).

## Toward the engine PR (for Jeremiah's review, not done here)

- Promote the wall as production's edge **for expanding runs**: K̇ from compatibility, mirror and W1b ghosts at the
  wall, no flux, no boundary data. Keep the open edge for static and isolated runs. Label: "identical neighbours, in
  step: one-stone simplification".
- Engineering the research version doesn't have yet: it subclasses the grid and installs its wall flux globally
  (`install_wall_flux` replaces a matter-RHS function module-wide). The engine version needs this as a kernel option
  with no global side effects.
- Initial data: the wall needs a conformally flat slice (`conformal_cell_slice`; K² = 24π⟨ρψ⁵⟩/⟨ψ⁵⟩). The background
  default is already flat (ψ = 1), so production's default needs no change. The pocket fixtures would need the cell's
  slice.
- Chat's follow-up: the new coupled kernel gets an audit of its own with this closure, on a start that isn't flat (a
  pocket fixture).

Unpaid inputs: λ = 5000 one-budget (Q); the one-stone simplification; closure (i) as the clock (the wall's proper
time); N = 20 (B, S) and 120 (Q).
