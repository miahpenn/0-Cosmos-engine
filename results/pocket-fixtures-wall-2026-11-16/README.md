# The pocket fixtures on the lattice wall by default (engine PR, for Jeremiah's review)

Decision (Jeremiah, Nov 16): "Pocket fixtures on the wall: yes, as a separate later step. Their pinned data changes, so
it should be its own labelled change." This is that change. Script: `research/pocket_fixtures_wall_checks_v1.py`
(`fingerprint`, `init`, `run STRUCT EDGE N T`, `research STRUCT N steps`). Short runs only (t <= 30). Floor, lambda,
coupling (#503) and the wall's evolution (#505) are not touched.

**What changed.** `initial_structure="rescaled_pocket"` and `"inherited_pocket"` now get the lattice wall when no edge is
asked for (`engine.lattice_wall.EXPANDING_STRUCTURES = ("none", "rescaled_pocket", "inherited_pocket")`). The hand-floor
reference and the radial-Maxwell experiment keep the open edge. Inside `lattice_wall.open_edge_reference()` the fixtures
get the open edge too.

**Labels (in code: `lattice_wall.POCKET_FIXTURE_WALL_LABEL`, the docstrings of `engine/lattice_wall.py`,
`engine/production_kernel.py`, `engine/v55_initial.py`):**
- **The fixtures' data changed.** On the wall they start on the conformally flat cell slice
  (`conformal_wall_initial_data`: a = b = 1, Aa = 0, K uniform with K^2 = 24 pi <rho psi^5>/<psi^5>), not on the areal
  slice of `engine/v55_initial.py`. S and D have the same profiles in the coordinate r on both slices, but r means
  something different on each.
- **Why:** the one-stone scope for expanding runs ("identical neighbours, in step: one-stone simplification").
- **The wall reflects outgoing waves.** Each fixture's ringing stays in the cell.
- **The open-edge versions** are still available by name (`outer_edge="open_edge"`). They are the labelled reference and
  are bit-identical to main.
- **Clock:** "i" for evolution (alpha = 1 at the wall). Physical comparisons are read on the fluids' clock:
  `fluids_clock_rate` and tau_fluids = its integral over t.

## Checks

| # | check | result | |
|---|---|---|---|
| 1 | Open-edge fixtures bit-identical to main when asked for by name | 25-step SHA-256 (N 60, r_max 40, CFL 0.03; geometry, scalars, fluids, t, tau, e-folds and six history columns). **Rescaled 548b148d…, inherited e21cf0d6…** are the same on main 19bc0ad (default edge = open edge there, and by name) and on the branch (by name) (`fingerprint_*.txt`). The hand-floor reference (dc76477b…), the background on the open edge (ed7409d6…), the background on its default wall (9e880c9d…) and #505's rescaled hash are also unchanged (`fingerprint_reference_and_background_*.txt`). | Read |
| 2a | Both fixtures, default edge: initial constraints, r_max 24, N 60 → 120 | **Rescaled:** H L2 9.79e-9 → 2.45e-9 (×3.99), max 1.56e-7 → 3.93e-8 (×3.96). **Inherited:** H L2 1.70e-5 → 4.25e-6 (×3.99), max 7.27e-4 → 1.83e-4 (×3.96). Second order, as expected for the slice. M = 0 to round-off (≤ 1e-16) because Aa = 0 and K is uniform. The slice: rescaled K = −1.70604e-2, X(wall) 1.000201; inherited K = −3.2847e-2 / −3.2844e-2, X(wall) 1.1296, alpha(0) 1.301 (`init.txt`). | Read |
| 2b | Short runs on the default edge, r_max 24, CFL 0.06, to t = 30, N 60 and 120 | **Both run to t = 30.** **Rescaled:** the numbers equal #505 check 5 to every printed digit. H L2 at t = 10 / 20 / 30: 8.8e-8 → 9.1e-9 (×9.7), 6.5e-7 → 4.0e-8 (×16), 1.2e-7 → 5.4e-9 (×22). M L2 at t = 30: 1.39e-7 → 2.22e-8 (×6.3). alpha(0) 0.9986–1.0004. **Inherited:** H L2 at t = 5 / 10 / 15 / 20 / 25 / 30: ratios ×6.6 / 7.3 / 7.9 / 8.3 / 8.0 / **3.5** (t = 30: 3.3e-4 → 9.6e-5). M L2 ratios ×7.4 / 5.5 / 2.5 / 2.8 / 4.1 / 3.0 (N 60: 0.9–1.4e-4; N 120: 2.0–3.5e-5). alpha(0) 0.77–1.30. In every run the largest \|H\| is in the **centre** cell, never at the wall (`run_*_wall_N*.txt`). | Read |
| 2c | The wall-cell condition (W1b) | **Rescaled holds:** \|M(last cell)\| ≤ 7.2e-11 (N 60) and ≤ 3.0e-11 (N 120). (3/R)\|Aa(wall)\| ≤ 2.7e-8 at both N (the even-Aa residual, as in #505). **Inherited holds until its outgoing S–D energy reaches the wall, then degrades.** \|M(last)\| ≤ 3e-8 to t = 15 at both N. At t = 20 it is 1.8e-6 / 2.1e-6 (N 60 / 120: not converging). At t = 25 it is 2.1e-4 / 3.5e-5 (×5.8). At t = 30 it is 2.2e-4 / 7.6e-5 (×2.9). (3/R)\|Aa_w\| reaches 1.0e-4 / 9.2e-5. | Read |
| 3 | Fixture on the default edge vs research `LatticeCellKernel("i")`, same state, same dt, N 60, r_max 24, 100 steps | **Bit-identical** (every field and tau) for the rescaled fixture, and also for the inherited one (`research_kernel_compare.txt`). | Read |
| 4 | What changes physically between the edges (r_max 24, N 60, t ≤ 30) | See below. | Read |
| 5 | Full suite | **1109 passed, 2 xfailed** (main 19bc0ad: 1100 passed, 2 xfailed; +9 new tests). | Read |

### Check 4: wall against open edge, the same fixture

Rulers: S–D energy above D's phi-dependent rest, split two ways: at coordinate r = 12, and at areal radius 12. The
whole-cell total is given beside both. Core H = −<K>/3 over r < 12, read at the same t and at the same tau_fluids.

**Rescaled (x = 2, by hand).**
- **At t = 0 the edges differ only through the slice.** Core H is 5.6868e-3 on the wall against 5.6690e-3 on the open
  edge (+0.31%), because on the wall K is fixed by the cell average.
- **At t = 30, the S–D energy is the same on both edges to 0.3% under both rulers.**
  - Coordinate r = 12: inside 1.0540e-4 (wall) vs 1.0541e-4 (open); outside 6.225e-6 vs 6.243e-6.
  - Areal radius 12: inside 9.808e-5 vs 9.807e-5; outside 1.354e-5 vs 1.358e-5.
  - Whole cell: 1.1163e-4 vs 1.1166e-4.
- **Core H at t = 30:** 5.1355e-3 (wall) vs 5.1237e-3 (open), +0.23%. The fluids' clock and t agree to 3e-5 here (tau_fluids = 30.000 on both), so the comparison at the same tau_fluids gives the same numbers.
- **The open edge's centre constraint residual grows. The wall's does not.** On the open edge, H L2 reaches 1.8e-5
  (max 4.7e-3, centre cell) by t = 30. On the wall at the same N it is 1.2e-7 (max 1.8e-5).
- **Read, against #505's README.** The growth of the S–D energy outside r = 12 (1.2e-6 → 6.2e-6 by t = 30) is the
  same on the **open edge**. So at t ≤ 30, r_max 24, it is the pocket's own outward ringing. It is not yet reflection
  from the wall. #505's note ("the ringing stays in the cell") is true of the wall, but this growth does not show it.

**Inherited (amplitude 0.01, x ~ 157, inherited).**
- **The data change is large here.**
  - On the wall the pocket's energy enters the cell average that fixes K. Core H starts at **1.0949e-2, against
    5.6690e-3 on the open edge (×1.93)**. In the one-stone reading, a lattice of such pockets expands almost twice as
    fast as the background.
  - S–D energy at t = 0: 3.53e-2 (wall) vs 4.29e-2 (open), −18%, under both rulers.
  - The fluids' clock starts at rate 1.0145 (wall) and 1.0185 (open). The pocket's inner redshift sits against
    alpha(0) = 1.30 and 1.36.
- **On both edges the pocket disperses by t ~ 20.** Energy inside r = 12 falls to 2% of its start by t = 22.5:
  7.8e-4 (wall) and 1.8e-3 (open).
  - Wall: the whole-cell energy goes 3.54e-2 → 2.81e-2 (t = 22.5) → 2.68e-2 (t = 30). Core H 8.03e-3 at tau_fluids
    22.33.
  - Open edge: 4.30e-2 → 3.55e-2 (t = 22.5). Core H 3.57e-3 at tau_fluids 22.38.
- **The open edge does not reach t = 30.** At t = 25 the whole-cell energy jumps to 4.97e-2 and alpha(0) to 1.167.
  The step after t = 26.544 fails with "dust rest-current must be positive" (`run_inherited_open_edge_N60.txt`).
  This is the reference path, bit-identical to main: an existing behaviour, not a new one.

## Uses of the fixtures on the default edge: each one sorted

| place | before | now |
|---|---|---|
| `tests/test_lattice_wall_kernel.py::test_which_runs_get_which_edge` | asserted the fixtures default to the open edge | **updated**: both fixtures default to the wall; by name and inside `open_edge_reference()` they get the open edge |
| `tests/test_hand_floor_removal_step2.py::test_pocket_fixtures_are_labelled_and_converge` | areal slice, fourth-order check (ratio ≥ 8, N 40 → 80, r_max 40) | **pinned to `outer_edge="open_edge"`** with its reason (it checks the areal slice of `v55_initial.py`). On the wall it failed loudly (ratio 3.8: the conformal slice is second order). The wall's own convergence check is the new module. |
| `tests/test_hand_floor_removal_step2.py::test_default_start_refuses_pocket_arguments` | `inherited_pocket` with omitted vs explicit old defaults; compares S | **switches to the wall, deliberately left there**: it compares only S, which is the same on both edges, so it means the same thing |
| `tests/test_lattice_wall_kernel.py` (other fixture uses), `research/zero_flux_wall_engine_checks_v1.py` | already pass `outer_edge` explicitly | unchanged |
| `tests/conftest.py::OPEN_EDGE_REFERENCE_MODULES`, the hand-floor reference modules, `tools/`, the seven research drivers pinned in #506, the workflows | no fixture use | unchanged (checked by grep: the fixture names occur only in the files above) |
| `engine/coupled_engine.initialize(initial_structure=...)` | forwards to the kernel | gets the wall too (tested) |

New: `tests/test_pocket_fixtures_on_wall.py`, 8 tests. It covers the default and the labelled slice (same S and D,
different X and K); the open edge by name equals the pre-step `open_edge_reference()` default array for array after 3
steps; initial H converges at second order (N 60 / 120); and a short run on the wall stays finite, with zero face flux,
the N/A mismatch, alpha(wall) = 1 and a small wall-cell M.

## CI

Exact-blob authorisations, recomputed with `git hash-object`:
- New `AUTHORIZED_POCKET_FIXTURES_WALL_HEAD_BLOBS`: `engine/production_kernel.py` a481da7e…, `engine/v55_initial.py`
  5640b485…. The `v55_initial.py` change is a docstring only.
- Off-audit `engine/lattice_wall.py` re-authorised to 466ae720….
- The same blobs are in the workflow `.github/workflows/tracefree-ricci-proper-radial.yml`, in
  `tools/audit_tracefree_ricci_proper_radial.py`, and in the pinning test
  `tests/test_tracefree_ricci_proper_radial.py::test_pocket_fixtures_wall_authorization_is_exact_blob`.

## Derived / Read / Open

- **Derived:** none new here. The slice's K is fixed by its solvability condition (#505).
- **Read:** checks 1–5 above.
- **Open:**
  - **The inherited fixture's wall cell.** After its S–D energy reaches the wall (t ≳ 20), the W1b residual rises to
    ~1e-4 and converges poorly (×1–6). H L2 converges at about second order or better until t = 25, then ×3.5 at
    t = 30. The fixture is barely resolved at N 60–120 (x ~ 157). Not run further: short runs only.
  - **The centre cell carries the largest H residual** in every run, on both edges.
  - **The fixtures on the wall start at the cell-average K.** The inherited one starts at ~2× the background H. Whether
    that is the one-stone answer wanted for a fixture, or a reason to keep the inherited fixture on the open edge, is
    for Jeremiah.
  - **The inherited fixture's open-edge reference fails at t ≈ 26.5 at r_max 24.** It is not tested at r_max 40 here.
  - **S's action** beside its energy is not computed here.

## Unpaid inputs

- The one-stone simplification.
- Clock "i".
- H/m_S, the one unpaid number.
- V0 and phi_today (borrowed).
- lambda (a ruler, matched today).
- Fixture only: x = 2 (by hand); the inherited amplitude 0.01; the D seed 1e-10; the cold start's T = 200, gamma = 0.3.
- N: 60 and 120. r_max 24.
