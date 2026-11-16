# Task 2, engine: the zero-flux (lattice-wall) closure in the production kernel (for Jeremiah's review)

HANDOFF Nov 15 (3), task 2. Decision 2 (confirmed by Jeremiah): expanding runs use the zero-flux closure, labelled
**"identical neighbours, in step: one-stone simplification"**. Static and isolated runs keep the open edge. The research
comparison is `results/zero-flux-closure-2026-11-16/` (#504). The engine change is `engine/lattice_wall.py` plus
dispatch in `engine/production_kernel.py`. Script: `research/zero_flux_wall_engine_checks_v1.py`. Short runs only.

**The wall reflects outgoing waves.** That is correct for identical neighbours, which send their own ringing back in,
and wrong for an isolated set-up. In check 5, the S–D energy outside r = 12 grows from 1.2e-6 to 7.1e-6 by t = 30: the
ringing stays in the cell.

## Checks (two rulers per result)

| # | check | result | |
|---|---|---|---|
| 1 | Open edge bit-identical where kept | 25-step SHA-256 fingerprints (N 60; geometry, scalars, fluids, t, tau, e-folds and six history columns) are **equal to main b119f56** for: the hand-floor reference, the background start on the open edge, and the rescaled pocket on the open edge (`fingerprint_*.txt`). #503's own fingerprint script also reproduces its recorded hash (309da765…). The tracefree audit's `profiles.npz` and `profiles.csv` are byte-identical to main's; `report.json` differs only in the head commit and kernel blob it records. | Read |
| 2 | Background on the wall vs the production COSMOS lane (`engine.cosmos`), r_max 40, t ≤ 200 | **N 20:** max \|dH/H\| 3.1e-8, max \|dφ\| 1.3e-8. **N 40:** 1.9e-8 and 2.9e-8. The two rulers are the mean over r < 20 and the proper-volume mean over the whole cell; they agree to all printed digits (the cell stays uniform). Lapse spread ≤ 2.8e-14. The centre's proper time against t: 2.1e-12 (N 20), 1.3e-12 (N 40). Same level as the open edge (#503: 3.3e-8 / 5.2e-8), and it does not shrink with N, so it is not the edge. | Read |
| 3 | Lapse seeded α = 1 + 0.01(1 − r²/R²) before the first solve, N 20, t ≤ 200 | **Wall: bit-identical to the unseeded run** (every field, every row, and tau). Open edge: H is −2.0% at t = 10 and −27.3% at t = 200; the seeded mode grows from 0.95% to 4.0%. A stronger version, with the seed left in the state and no re-solve: after 50 steps the evolved fields are still bit-identical. Only tau differs, by 6.0e-4, because the first step's ledger quadrature reads the stored lapse (ledger only). | Read; mechanism Derived (W4) |
| 4 | Engine wall vs research `LatticeCellKernel` (closure i), same state, same dt | **Bit-identical** after 100 steps (N 60, r_max 24), for both the background and the rescaled pocket on the conformally flat slice (`research_kernel_compare.txt`). The engine's conformal slice is also bit-identical to research `conformal_cell_slice` with the potential written out by hand. | Read |
| 5 | Pocket fixture on the wall (rescaled pocket, x = 2 by hand, conformally flat slice), r_max 24, t ≤ 30, N 60 and 120 | **Runs to t = 30.** Hamiltonian residual, proper-volume L2: t = 0: 9.8e-9 → 2.5e-9 (×4.0, second order: the slice); t = 10: 8.8e-8 → 9.1e-9 (×9.7); t = 30: 1.2e-7 → 5.4e-9 (×22). Momentum residual, L2: t = 30: 1.4e-7 → 2.2e-8 (×6.3). Last cell \|M\| ≤ 7e-11 throughout. (3/R)\|Aa(wall)\|, the residual an even Aa would leave there, is 1–3e-8 at most snapshots. Lapse 0.9986–1.0004 at the centre, 1 at the wall. Core H falls 5.69e-3 → 5.14e-3. | Read |
| 6 | Full suite | **1097 passed, 2 xfailed** (main: 1087 passed). The production one-step smoke, the strict-retry campaign smoke, `run_controls` and the workflow's fail-closed blob step all pass locally. | Read |

The new audit is `tests/test_lattice_wall_kernel.py::test_audit_coupled_kernel_on_the_wall_non_flat_start_converges`.
It runs the coupled kernel (#503) with the wall on a non-flat start, N 30 / 60 / 120 to t = 1, using maxima over the run:

| N | H residual (L2) | M residual (L2) | M, wall cell | (3/R)\|Aa(wall)\| |
|---|---|---|---|---|
| 30 | 3.9e-8 | 3.3e-9 | 3.8e-10 | 2.2e-8 |
| 60 | 9.8e-9 | 6.8e-10 | 1.9e-10 | 2.2e-8 |
| 120 | 2.5e-9 | 1.6e-10 | 9.6e-11 | 2.2e-8 |

H and M converge at second order. The wall cell converges at **first** order. That is Read, and expected: W1b's ghosts
for a and b are first order in ln a (exact for Aa).

## Derived / Read / Open

- **Derived (W4):** with zero lapse gradient at both ends, the lapse equation fixes only the lapse's shape. K̇ is one
  number set by the clock, never by the stored lapse. So the open edge's undamped lapse mode cannot exist on the wall.
- **Read:** checks 1–6 above.
- **Open:**
  - **The clock.** Clock "i" (α = 1 at the wall: the proper time of the wall shell, where the matter is at rest by
    symmetry) is the default. The lane's own clock is the rest-mass-weighted proper time of its fluids. The two
    coincide for the background. In a pocket cell they differ by the pocket's inner redshift. A rest-mass-weighted
    normalisation is a third option, not built.
  - **The wall-cell momentum residual converges at first order** (W1b on a, b).
  - **Open-edge seam contracts.** Eleven open-edge contract modules run on the open edge
    (`tests/conftest.py::OPEN_EDGE_REFERENCE_MODULES`). For two of them the wall has its own pinned equivalent (stage
    ordering, native a/b/X). For the rest (constraint-consumer normalisation, stage shadows, native-constraint and
    λ-L3 contracts, the puncture-gauge reproduction, EM) it does not yet.

## Unpaid inputs

- The one-stone simplification itself.
- Clock "i".
- H/m_S, the one unpaid number.
- V0 and φ_today (borrowed).
- λ (a ruler, matched today).
- Fixture only: x = 2, chosen by hand.
- N: 20 and 40 (checks 2, 3); 60 and 120 (checks 4, 5).
