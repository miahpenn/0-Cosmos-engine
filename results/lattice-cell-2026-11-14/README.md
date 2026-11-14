# The lattice-cell wall, its acceptance tests, and the universe of pockets with expansion

HANDOFF to projects, 14 Nov 2026 (4) (GEAR-107, AMENDMENT-GEAR-107, AUDIT-2). RESEARCH ONLY: no production file is
changed. The wall is a subclass of the production grid and kernel (`research/lattice_cell_v1.py`).

Scripts: `lattice_cell_v1.py` (the wall, the lapse closures, the initial slice), `lattice_cell_run_v1.py` (driver),
`cell_background_ode_v1.py` (the cosmic-lane ODE for T1), `lattice_cell_t1_v1.py`, `lattice_cell_tests_v1.py` (T2–T5),
`lattice_cell_t4_clock_v1.py` (T4's second ruler), `lattice_cell_universe_v1.py` (Task 3). Raw outputs: `T1.txt`,
`T2_T3.txt`, `T4.txt`, `T4_second_ruler.txt`, `T5.txt`, `universe_*.txt`, `spike.txt`.

Marks: **Derived** = follows from the equations; **Read** = read off the runs; **Open** = not settled here.

## Unpaid inputs and forced components (every result below)

- λ = 5000 (D* = 0.02). One-budget rule (AMENDMENT-REV20-2 §4): hand floor off, D's ejected vacuum is the floor,
  V0 and present densities × 304.92. β_D = β_S = 0. G_reset not involved.
- Pocket: radius 8, edge width 1, x_c = 2.0 (and 1.5 in Task 3), D seed 1e-10, cold start (D relaxed against S by
  friction, **forced**) or hot start (D at its analytic settled value, unrelaxed).
- Cell radius r_max (12, 16, 24, 36), grid spacing 0.1 (0.075 as a ruler), time step 0.2 dr.
- **Forced (AUDIT §D):** spherical symmetry; a spherical cell standing in for a space-filling one; all neighbours
  identical and in step (that is what the wall *is*); CMC slicing with zero shift; a conformally flat initial slice;
  the pocket and its start placed by hand.

## Task 1. The wall

The outer boundary is the last cell face, treated as a symmetry surface.

| | rule | status |
|---|---|---|
| W1 | Scalars and every even field have zero slope at the wall; ordinary centred stencils run up to it. Ghost values are a quartic through the last cells with zero slope (the plain mirror image is kept as a second ruler). | as specified |
| W1b | The trace-free variables cannot be even. The momentum constraint at the wall gives Aa′ = −3Aa/R there. The quantities with zero slope are r³b^{3/2}Aa/X³ and r³(1 − b^{3/2}). | **Derived** here; not in GEAR-107 |
| W2 | No flux of anything. The dust flux at the wall face is the flux between the last cell and its mirror image (velocity reversed): mass and energy flux vanish identically. | as specified |
| W3 | No boundary data. The vendor's incoming-light overwrite of Aa at the last cell is not applied. Shift zero. | as specified (E2) |
| W4 | Lapse: zero gradient at the wall and the centre. That fixes the lapse's shape and leaves one number, K̇. Closure (ii): K̇ = ⟨Q⟩, the proper-volume average of Q = K_ijK^ij + 4π(ρ + S). Closure (i): α = 1 at the wall. | (ii) is my reading of "K fixed by the cell average"; **for chat to confirm** |
| W5 | No sponge, no layer, no friction. | as specified (E6) |

- **Derived:** the two closures give the same slices. They differ only by relabelling time, dt(i) = α(ii)(wall) dt(ii).
  They are one foliation with two clocks, and agree identically in a uniform cell.
- **Derived:** with zero slope of the conformal factor at both ends, the initial Hamiltonian constraint can be solved
  only for K² = 24π⟨ρψ⁵⟩/⟨ψ⁵⟩. The cell's expansion rate is set by the cell's mean density, pocket included. It is
  not an input.
- **Read:** the areal initial slice used in #475 cannot meet the wall's symmetry (Hamiltonian residual 0.13 at the
  wall). The cell starts from a conformally flat slice instead, with the pocket laid down in proper radius at the
  centre.
- Things that had to be right before the tests passed, each kept as a switch or note in the module: Aa even (momentum
  residual exactly 3Aa/R at the wall); a second-order lapse solve against a fourth-order K equation (sawtooth at the
  wall); piecewise-constant dust in the last cell (±0.5% sawtooth); ψ normalised at the wall (pocket size depended on
  the cell).

## Task 2. Acceptance tests

| test | result | verdict |
|---|---|---|
| T1 no pocket | uniform to 1e-12 through the turnaround; H(t) matches the cosmic-lane ODE to 7e-6 H0 | pass, with one qualification (accuracy) |
| T2 constraints at the wall | wall residuals at or below the quiet interior's for r_max ≥ 24; in the small cells (12, 16) 5–8 times the quiet zone beside the wall, though 6–20 times below the pocket's | pass for r_max ≥ 24; qualified for 12 and 16 |
| T3 ledger and shell equation | ledger closes to 4–9e-4, shell equation to 6e-6–4e-5; flux term 1e-5 of the work | pass |
| T4 same pocket in three cells | cells of radius 24 and 36 agree to 0.4% in S's action; smaller cells differ in proportion to the pocket's share of the cell | pass for r_max ≥ 24; the central spike is not converged in grid spacing |
| T5 closures | same fields at equal wall proper time to 4e-4–9e-4; clocks differ by 0.6% (16), 0.3% (24), 0.1% (36) | pass |

### T1. No pocket (`T1.txt`)

Cells r_max = 12 and 24, both closures, to t = 70 (turnaround at t = 55.0–55.2, after 1.988 e-folds).

- **Read:** every field stays uniform across the cell to round-off: spreads 1e-13–1e-12 (K: 1e-10), |a/b − 1| ≤ 8e-12,
  |Aa| ≤ 1e-13, α − 1 ≤ 1e-13. Through the turnaround.
- **Read:** H(t) against the cosmic-lane ODE: largest difference over the run 7.1e-6 of H0 at time step 0.02, 2.1e-6
  at 0.01 (2.6e-7 and 6.4e-8 at t = 40). Turnaround time agrees to the snapshot spacing.
- **Qualification:** GEAR-76's 2.5e-11 is not reached at these time steps. The difference falls by 3.4–4 when the
  step is halved (second order in time, as the stepping is), and is the same in both cells and both closures.
  Reaching 2.5e-11 would need a step ~500 times smaller.
- Before, on the production edge (PR #477 Task 6), the same no-pocket background at t = 40: dust +2% at the edge, a/b
  off 1 by 1.3%, central lapse 1.8% above the edge's.

### T2. Constraints at the wall (`T2_T3.txt`)

Pocket x_c = 2.0, cold start, t = 0–40, largest residual over the run:

| r_max | wall \|H\| | interior (r > 11) \|H\| | pocket (1–11) \|H\| | centre (r < 1) \|H\| | wall \|M\| | pocket \|M\| |
|---|---|---|---|---|---|---|
| 12 | 5.4e-6 | 7.0e-7 | 3.3e-5 | 3.6e-2 | 2.0e-6 | 2.3e-4 |
| 16 | 1.9e-6 | 3.8e-7 | 3.7e-5 | 2.2e-2 | 1.4e-7 | 1.6e-4 |
| 24 | 4.0e-7 | 3.9e-7 | 3.2e-5 | 2.4e-2 | 2.8e-8 | 1.8e-4 |
| 36 | 8.4e-8 | 4.0e-7 | 4.1e-5 | 2.2e-2 | 6.1e-9 | 1.8e-4 |
| 24, dr 0.075 | 3.0e-7 | 2.2e-7 | 1.9e-5 | 1.5e-2 | 2.2e-8 | 8.5e-5 |

(16πρ at the wall falls 5.9e-2 → 1.7e-3 over the run.)

- **Read:** the wall is never the worst place: its residuals are 6 to 500 times below the pocket's, and fall as the
  cell grows (roughly as 1/r_max⁴ in H).
- **Read, qualification:** against the quiet zone just inside it (r > 11) the wall is equal at r_max = 24 (4.0e-7
  against 3.9e-7), 5 times lower at 36, but 5 times higher at 16 and 8 times higher at 12. Strictly "no larger than
  the interior" holds for r_max ≥ 24 only. In the small cells the pocket's edge is 1–5 from the wall, so the wall is
  not in a quiet zone at all.
- **Read, and a limit:** the worst residuals are at the centre (r < 1), 1.5–3.6e-2, in every cell alike, late in the
  run when the dust spike forms. They fall only by 1.5 from dr = 0.1 to 0.075. The centre is under-resolved there;
  that is not the wall's doing but it bounds what can be said about the spike (T4, Task 3 f).
- Mirror ghosts instead of quartic ones: wall residual twice as large, everything else identical.

### T3. Ledger and shell equation (`T2_T3.txt`)

At the wall, t = 0–40: dM/dt (Misner–Sharp mass inside the wall) against the pressure work plus the flux term, and
(dR/dτ)² against 2M/R + Γ² − 1 with M integrated from the centre.

| r_max | ledger residual / rms dM/dt | flux term / work | shell equation, largest difference / U² | baryon rest mass in the cell, largest change |
|---|---|---|---|---|
| 12 | 8.8e-4 | 7e-5 | 3.9e-5 | 6.7e-16 |
| 16 | 8.1e-4 | 3e-5 | 2.6e-5 | 2.2e-16 |
| 24 | 5.9e-4 | 5e-6 | 1.2e-5 | 2.2e-16 |
| 36 | 4.4e-4 | 1e-7 | 6.0e-6 | 4.4e-16 |

- **Read:** the wall's mass changes by the wall pressure's work and nothing else; the flux term is zero to 1e-4 of the
  work or better. The ledger residual is the 0.5-unit snapshot differencing of dM/dt.
- **Read:** no baryon crosses the wall (rest mass conserved to round-off). Dark-matter rest mass is not conserved in
  the engine by construction (it has the φ source), so it is not a test.
- **Read:** Γ² − 1 at the wall is −2e-3 (r_max 12) to −2.5e-4 (36): the wall shell is slightly bound, less so as the
  cell grows.

### T4. The same pocket in different cells (`T4.txt`, `T4_second_ruler.txt`)

Pocket x_c = 2.0, cold start, t = 0–40. Reference: r_max = 36. Largest difference from the reference over six windows:

| measure | r_max = 24 | 24, dr 0.075 | 16 | 12 |
|---|---|---|---|---|
| pocket's share of the cell's mean density (from H0²) | 1.8% (36: 0.55%) | 1.8% | 5.9% | 13% |
| H0 | +0.7% | +0.7% | +2.7% | +6.3% |
| S's action J (r < 11), ratio to its first window | 0.4% | 0.4% | 1.4% | 3.0% |
| D at rest (r < 8), in units of D* | 0.002 | 0.002 | 0.007 | 0.013 |
| D's oscillation energy, r < 8 (slow + travelling) | +2.5% | +2.3% | +5–11% | +13–26% |
| slow ring / travelling split, filter over the whole cell | 0.90–0.91 / 1.04–1.05 | 0.90–0.93 / 1.02–1.05 | 1.08–1.29 / 1.03–1.12 | 0.65–0.76 / 1.20–1.39 |
| same split, filter over r < 12 only | 1.00–1.05 / 1.01–1.03 | 1.01–1.04 / 1.00–1.02 | 1.00–1.20 / 1.04–1.11 | 1.02–1.43 / 1.12–1.27 |
| dust profile at r = 2, 4, 8 (over its value at r = 11) | 0.4% | 0.6% | 1.8% | 5% |
| dust at the innermost cell | 2.2% | 10% | 12% | 39% |

- **Read (stated tolerance):** with the pocket well inside the cell (r_max ≥ 24, the wall at least 13 from the
  pocket's edge) the pocket is the same to 0.5% in S's action, 0.3% of D* in D at rest, 3% in D's oscillation energy
  and 1% in the dust profile outside the innermost cells. At that tolerance T4 passes.
- **Read (how it changes as the cell shrinks):** every difference grows in step with the pocket's share of the cell's
  mean density, roughly as 1/r_max³: 0.55% → 1.8% → 5.9% → 13%. A smaller cell is a denser universe: it expands faster
  (H0 up 6% at r_max = 12), S's action falls slightly faster (3% lower by t = 36), D carries more oscillation energy
  in the pocket (+20%), and dust gathers more slowly at the centre.
- **Derived:** H0² rises as c/r_max³ with c = 2.2–2.6 across the four cells, as it must if the pocket adds a fixed mass
  to the cell.
- **Two rulers, and what is not invariant:** the slow-ring / travelling partition moves with the filter's domain (slow
  part 1.4e-3 with the whole-cell filter, 0.9e-3 with the r < 12 filter, same fields). Only the range is invariant;
  their sum is. Comparing in windows of the centre's proper time instead of the label (`T4_second_ruler.txt`) gives
  the same picture with ±1–2% scatter from where the 0.5-unit snapshots fall in the windows.
- **Read:** mirror ghosts give the same numbers as quartic ghosts to all printed digits (r_max = 16).
- **Not converged:** the dust at the innermost cells late in the run. At t = 40 the centre stands at 14.2 (r_max 36),
  13.9 (24), 16.2 (24, dr 0.075), 12.7 (16), 9.1 (12) times the value at r = 2. The grid-spacing difference (15%) is
  larger than the cell-size difference between 24 and 36 (2%).
- **Open:** the spike is earlier and stronger here than on the production edge (PR #477 Task 5: 2.4 at t = 38,
  r_max = 36). The two differ in initial slice (conformally flat against areal), in the clock, and in the edge; which
  of these moves the spike has not been isolated.

### T5. The two closures (`T5.txt`)

Pocket x_c = 2.0, cold start, t = 0–40, closures (ii) and (i) in the same cell.

| r_max | α(ii) at the wall − 1, range | label difference at t = 40 | S(0), largest difference at equal wall proper time (fine series) | K | dust centre/wall |
|---|---|---|---|---|---|
| 16 | −1.8e-2 … +2.9e-2 | 0.244 | — | 1.1e-4 | 2.0e-4 |
| 24 | −1.3e-2 … +1.8e-2 | 0.118 | 9e-4 | 6.7e-5 | 1.1e-4 |
| 36 | −7.0e-3 … +8.5e-3 | 0.041 | 4e-4 | 3.0e-5 | 1.4e-4 |

- **Read:** compared at equal wall proper time the two closures give the same fields (as derived: one foliation). At
  equal label they differ by 5–33% in S(0), because S oscillates and the labels drift apart.
- **Read:** the clocks themselves converge as the cell grows: the label difference falls 0.244 → 0.118 → 0.041 of 40.
- The 1–2% differences that `T5.txt` prints for S(0) and D(0) come from interpolating 0.5-unit snapshots of
  oscillating fields; with the fine time series they are 4e-4–9e-4.

## Task 3. The universe of pockets, with expansion

Runs (one-budget rule, λ = 5000, dr = 0.1, closure (ii), t = 0–70; `universe_*.txt`, `spike.txt`,
`spike_resolution.txt`):

| run | x_c | start | r_max | note |
|---|---|---|---|---|
| main | 2.0, 1.5 | cold, hot | 24 | four runs |
| control | 2.0, 1.5 | — | 24 | same pocket, no D |
| cell size | 2.0 | cold, hot | 16, 36 | second ruler: the cell |
| closure | 2.0 | cold | 24 | closure (i): second ruler, the clock |
| resolution | 2.0 | cold, hot | 16 | dr = 0.05 |

Rulers used: time label against wall proper time (and centre proper time for the settling clock); closure (ii)
against (i); cell 16 / 24 / 36; dr 0.1 against 0.05; for D's split, filter over the whole cell against r < 12.

### What limits every result below

- **The background turns round inside the run (Read; known from the one-budget cosmos and T1).** K = 0 at t = 55.0 in
  every run, after 1.75–1.79 e-folds. From t = 55 the cell contracts. "As expansion damps the ring" therefore covers
  t < 55 only; t = 55–70 is contraction.
- **The grid is comoving, the fields are not (Derived).** The areal size of one grid cell at the centre grows from 0.13
  to 0.6 by t ≈ 45. D's and S's Compton lengths at D* are 0.71 and 0.5. From t ≈ 30 the fields' own lengths are one or
  two cells. Bulk read-outs (action, D at rest, energies) agree between dr = 0.1 and 0.05 to 0.1–1%; the centre does
  not (items c, f).
- **Hot starts violate the constraints (Read).** Cold runs: |H| ≤ 3e-4, |M| ≤ 1e-3 outside r < 1. Hot runs: 2e-2 and
  5e-2 through the pocket in t < 30, comparable with the local source 16πρ (6e-2 → 5e-3). Hot results are
  qualitative.
- The wall stays clean in every run: |H| at the wall ≤ 4e-7 (r_max 24), 9e-8 (36); baryon rest mass conserved to 1e-15.

### (a) S's action over the whole cell. Does "D costs matter no action" survive?

J over the whole cell, last window (60–66) over the first (0–6):

| run | one oscillator, fine series: with D / without D | local sum, snapshots: with D / without D | with D ÷ without D |
|---|---|---|---|
| x_c 2.0 cold | 0.970 / 0.987 | 0.793 / 0.810 | 0.983 / 0.979 |
| x_c 2.0 hot | 0.937 / 0.987 | 0.762 / 0.810 | 0.949 / 0.941 |
| x_c 1.5 cold | 0.980 / 0.990 | 0.797 / 0.808 | 0.990 / 0.987 |
| x_c 1.5 hot | 0.982 / 0.990 | 0.796 / 0.808 | 0.992 / 0.985 |

- **Two rulers, and only the range is invariant.** "One oscillator" is √(⟨ΣP_S²⟩⟨ΣS²⟩) from the fine series; "local
  sum" is Σ√(⟨P_S²⟩⟨S²⟩)dV from 0.5-unit snapshots. S's action falls by 2–3% on the first and 20% on the second, with
  or without D. The local sum can only be the smaller (Cauchy–Schwarz); they separate once S's frequency differs from
  place to place. The local sum also moves by 3% with where the snapshots fall in the first window (label against wall
  clock: 0.793 against 0.818; the one-oscillator figure moves by 0.0005). The first window contains D's fall into its
  well, so it is not a clean starting value.
- **Read:** the cost of D, against the control without D, is the same under both rulers and both clocks: 1–2% of S's
  action for cold starts and for x_c = 1.5, 5–6% for the hot start at x_c = 2.0.
- **Read:** in the cold runs the cost is paid while D settles (t < 12). After that the ratio J(with D)/J(without D)
  over the whole cell is flat: 1.074–1.089 for x_c = 2.0 and 1.155–1.165 for x_c = 1.5, from t = 12 to 66, through the
  turnaround. In the hot x_c = 2.0 run it keeps falling, 1.020 → 0.989, mostly after t = 36.
- **Verdict:** "D costs matter no action" survives expansion and the turnaround at the 1% level once D has settled
  from a cold start. It does not survive the settling itself (1–2%), and not the hot start at x_c = 2.0 (5–6%, a run
  with constraint violation; **Open**).
- **Read:** D keeps S's action in the pocket. J(r < 11) with D over without D rises from 1.10 to 1.22 (x_c 2.0 cold),
  1.19 to 1.29 (1.5 cold). Without D the pocket loses 13% (one oscillator) to 28% (local sum) of its action to the
  rest of the cell; with D 3% to 21%.
- Cell size: last-window J is 0.9740 / 0.9697 / 0.9712 at r_max 16 / 24 / 36 (cold), 0.9329 / 0.9366 / 0.9375 (hot).
  Closure (i): 0.9705. dr = 0.05 (r_max 16): 0.9760 against 0.9740.

### (b) The settling clock. Does LOCK appear?

Lead of the pumping peak over S's saddle crossing, per half-cycle; LOCK = |lead| < 0.1 half-period. With no lock the
share expected by chance is about 0.2.

| run | share locked, P1 (energy into D), by thirds of the run | P3 (d/dt of D's kinetic energy) | mean \|lead\| / half-period, P1 |
|---|---|---|---|
| x_c 2.0 cold | 0.08 / 0.08 / 0.14 | 0.31 / 0.15 / 0.14 | 0.35 / 0.36 / 0.39 |
| x_c 2.0 hot | 0.00 / 0.07 / 0.14 | 0.00 / 0.21 / 0.07 | 0.30 / 0.30 / 0.32 |
| x_c 1.5 cold | 0.23 / 0.00 / 0.00 | 0.23 / 0.08 / 0.07 | 0.28 / 0.45 / 0.44 |
| x_c 1.5 hot | 0.14 / 0.21 / 0.00 | 0.14 / 0.36 / 0.07 | 0.23 / 0.28 / 0.37 |

- **Read: no LOCK.** In no run, cell, closure or resolution does the locked share rise above chance as the ring dies
  down, and the mean lead does not shrink. The same holds with the lead measured in the centre's proper time, and for
  the whole-cell sums.
- **Read, a limit:** the clock runs out of signal. S's amplitude on its initial profile falls from 0.7 to 0.02 by
  t = 13 and 0.004–0.006 after t = 40. The late crossings are of a remnant 1% the size of the pocket.
- P2 (energy out of S) shows 0.36–0.54 in the hot runs. P2 contains S·P_S itself, so its peak is tied to S's phase by
  construction once D² is steady. It is not evidence of lock.

### (c) Does the pocket collapse in a self-consistent background?

| run | smallest α(0), at t | trapped surface inside the cosmological horizon | dust centre / r = 2 at t = 64 |
|---|---|---|---|
| x_c 2.0 cold (r_max 16 / 24 / 36) | 0.852 / 0.853 / 0.847, t = 28–33 | none | 160 / 171 / 172 |
| x_c 2.0 cold, closure (i) | 0.847, t = 33.5 | none | 169 |
| x_c 2.0 cold, dr 0.05 (r_max 16) | 0.851, t = 33 | none | 1173 |
| x_c 1.5 cold | 0.898, t = 1 | none | 1.22 |
| x_c 2.0 hot (16 / 24 / 36) | 0.660 / 0.674 / 0.675, t = 13 | none | 1787 / 1748 / 1736 |
| x_c 1.5 hot | 0.755, t = 11.5 | none | 1929 |
| x_c 2.0 hot, dr 0.05 (r_max 16) | 0.001 at t = 65.5; run stops | flagged at t = 65.5, not clean (below) | 10169 |
| no D | 0.854, t = 1 (the opening transient) | none | 1.06 |

- **Read:** the pocket as a whole does not collapse. In every run the lapse outside r < 1 stays within 0.85–1.0 of the
  wall's, and there is no trapped surface inside the cosmological horizon at dr = 0.1. The Kodama clock √|1 − 2m/R|
  goes to zero only at the cosmological horizon.
- **Read:** the cosmological horizon behaves as it should. Spheres beyond areal radius 9.1 are anti-trapped at t = 0
  (the pocket, areal radius ~9, is horizon-sized at the start); the horizon moves out to 118 by t = 48, leaves the cell
  around the turnaround, and returns as a *trapped* region moving inward during contraction (areal 117 at t = 60, 54
  at 64, 29 at 68).
- **Not a clean result (Open):** the hot run at dr = 0.05 ends at t = 65.6 with the central lapse falling 0.72 → 0.001
  in 1.5 time units, during contraction. The innermost cell carries 2m/R ≈ 0.9–1.0 from t ≈ 52 on, with a Hamiltonian
  residual of 1e2 rising to 2e4. The same run at dr = 0.1 reaches 2m/R ≈ 0.4 in its innermost cell and does not
  collapse. This is an unresolved cusp failing, not a measured horizon. Whether the hot core collapses is open.

### (d) The core

Core = r < 3 (comoving); its areal radius grows 3.85 → 19.3 at the turnaround → 18.3.

| run | x = 2λ⟨S²⟩, r < 3: first window → last | x at the centre, largest window value | S's areal half-width, late | D at rest / D*, r < 3, after t = 24 |
|---|---|---|---|---|
| x_c 2.0 cold | 1.17 → 0.066 | 0.74 (t = 24–30) | 2.0–2.3 | 0.954–0.982 |
| x_c 2.0 hot | 1.18 → 0.050 | 1.22 | 2.1–4.1 | 0.957–0.987 |
| x_c 1.5 cold | 0.82 → 0.031 | 0.06 | 3.8–4.4 | 0.971–1.001 |
| x_c 1.5 hot | 0.88 → 0.033 | 4.5 (t = 60–66) | 1.5–3.1 | 0.946–0.989 |

- **Read:** the core does not compress as a whole. Its mean x falls by a factor 18–27 as S spreads and the cell expands.
- **Read:** a small central knot of S does form where the dust spike is (x_c = 2.0 cold; both hot runs): areal
  half-width ~2, which stays fixed in areal size while the cell expands (0.3–0.5 in comoving radius, 3–5 grid cells
  at dr = 0.1). x at the very centre climbs back to 0.7 (cold), 1.2 and 4.5 (hot). At x_c = 1.5 cold there is no knot.
- **Read:** in the same runs D at the centre does not rest at D*. It swings through zero, between about +1.9 D* and
  −1.0 D* (cold; first sign change at t = 23, both resolutions), +3.3 and −2.7 D* (hot; from t = 2). The sign change
  sits at comoving r = 0.05–0.45 (areal 0.3–2.1). At x_c = 1.5 cold D(0) stays within 0.7–1.1 D* and never changes
  sign.
- **D in three parts, r < 3:** at rest 0.95–0.99 D* (above). Slow ring: 2.3e-4 falling to 5e-5 (whole-cell filter),
  1.1e-4 to 2e-5 (r < 12 filter); only the range is invariant. Travelling: 3.5e-3 rising to 7.7e-3 (cold), 5–7e-3
  (hot), the same to 5–10% under both filters. These are sums over a zone whose proper volume grows 120-fold, so the
  travelling energy *density* falls by a factor ~50.
- dr = 0.05 against 0.1 (r_max 16): D at rest within 0.004 D*, travelling energy within 1.5%, centre x within 10%, S's
  half-width 1.9 against 2.35.

### (e) The pocket's mean D² and S's local frequency

- **Read:** the pocket's mean D² (ring included) rises 0.50 → 0.91 D*² by t = 6–8 and to 0.98–1.00 by t = 18. It then
  holds 0.994–1.001 through the expansion and the turnaround.
- **Read:** the ring's amplitude (swing of the pocket mean of D² within a window of 2) is largest in the window
  t = 6–8 in every run, 0.014–0.049 D*². By t = 36 it is 0.001–0.002, by t = 68, 0.0002–0.0013.
- **Read: yes, S's mass turns there, but the two are not independent.** S's frequency ω_S rises at +0.07 to +0.09 per
  unit time up to the ring's peak and is flat after it (−0.005 to +0.005). S's mass from the mean D² rises at +0.054
  to +0.061 before and +0.005 to +0.007 after, settling at 1.998–2.001. The ring peaks as D arrives at D*, which is
  the same moment S's mass stops rising. There is no second, later peak in the ring during the expansion at which to
  test the turn separately.
- Ruler on ω_S: window to window it scatters 1.95–2.13 around 2.0 (a window of 2 against a half-period of 1.6). The
  mass from ⟨D²⟩ is the steadier measure. Whole-cell values agree with the pocket's to 0.001.

### (f) The spike and the dust profile against cell size

Dust density over its wall value at t = 64, x_c = 2.0:

| run | r = 0.05 | r = 1 | r = 2 | r = 4 | r = 8 |
|---|---|---|---|---|---|
| cold, r_max 16 / 24 / 36 | 276 / 320 / 333 | 1.64 / 1.65 / 1.66 | 1.72 / 1.87 / 1.93 | 1.66 / 1.69 / 1.70 | 1.12 / 1.13 / 1.13 |
| cold, closure (i), 24 | 316 | 1.65 | 1.87 | 1.69 | 1.13 |
| cold, r_max 16, dr 0.05 (innermost cell at 0.025: 2033) | 1054 | 1.63 | 1.73 | 1.66 | 1.12 |
| hot, 16 / 24 / 36 | 3511 / 3645 / 3690 | 2.75 / 2.84 / 2.86 | 1.97 / 2.09 / 2.13 | 1.68 / 1.74 / 1.75 | 1.15 / 1.17 / 1.17 |
| hot, r_max 16, dr 0.05 | 9990 | 2.66 | 1.94 | 1.67 | 1.14 |
| x_c 1.5 cold, 24 | 1.68 | 1.42 | 1.38 | 1.59 | 1.11 |
| no D, 24 | 1.46 | 1.38 | 1.39 | 1.37 | 1.12 |

- **Read:** the dust profile from r = 1 outward does not depend on the cell (1–4% between 24 and 36, up to 10% at 16),
  on the closure, or on dr (≤ 3%).
- **Read: the spike is a cusp, not a number.** It lives in r < 0.5 and its central value grows without limit as the
  grid is refined (276 → 2033 at the innermost cell, cold; 3511 → 19766 hot). The dust *mass* in it is converged:
  share of the cell's dust inside r < 0.5 is 1.62e-4 at both resolutions (cold), 9.3e-4 and 9.4e-4 (hot) (r_max 16).
  Inside r < 1 at r_max = 24 it is 2.3 (cold) and 7.9 (hot) times a uniform cell's share. "ρ(0)/ρ(2)" should not be quoted as a result; the mass inside a stated radius can be.
- **Read:** the spike is not the wall's doing: 16 / 24 / 36 agree to 6% (hot) and 20% (cold) at the centre, and the
  closure changes it by 1%.
- **Read:** the spike needs D and goes with D's sign change at the centre. No D: 1.06. x_c = 1.5 cold (D(0) never
  changes sign): 1.22. x_c = 2.0 cold: the spike leaves 1.04 at t = 20, is 1.2 at t = 24; D(0) first changes sign at
  t = 23. Hot runs: D(0) swings from t = 2; spike 20 by t = 16. **Open:** the mechanism (β_D = 0, so dust feels D only
  through gravity).
- This answers part of T4's open item: the cell's spike differs from the production edge's mainly because it is a
  cusp read at one cell, not because of the wall.

### Task 3 in one table

| question | answer | mark |
|---|---|---|
| (a) does "D costs matter no action" survive | yes to 1% after D has settled from a cold start, through expansion and turnaround; the settling costs 1–2%; hot start at x_c 2.0 costs 5–6% | Read; hot case Open |
| (b) does LOCK appear as expansion damps the ring | no; locked share stays at chance; the clock's signal falls to 1% | Read |
| (c) does the pocket collapse | no trapped surface, α(0) ≥ 0.66, at dr 0.1; one hot run at dr 0.05 fails at the centre during contraction, unresolved | Read; centre Open |
| (d) does the core compress | not as a whole (x falls 18–27-fold); a knot of S of areal width ~2 forms at the centre with D swinging through zero there | Read; knot under-resolved |
| (e) does S's mass turn where the ring peaks | yes, at t = 6–8, but that is D's arrival at D*; no independent later peak | Read |
| (f) spike and dust profile against cell size | profile outside r = 1 independent of cell, closure, dr; the spike is a cusp whose mass converges and whose height does not | Read |

Contradictory or negative results kept: no LOCK; the 20% fall of the local-sum action against 2–3% for the one-oscillator
action; the hot start's 5–6% action cost; the dr = 0.05 hot run's failure at the centre; the non-converging spike height.

## For Jeremiah's review: the wall for map row 7 (no production change made)

Map row 7 (`docs/GEAR_TO_CODE_MAP.md`): "Outer boundary not specified by the gears. Production applies the vendor's
incoming-light constraint at a finite radius." GEAR-107 now specifies it. Proposal: replace the production outer edge
by the lattice-cell wall W1–W5 with W1b, as in `research/lattice_cell_v1.py`.

Before = production edge (PR #477, Task 6). After = the wall, same pocket, r_max = 24, dr = 0.1.

| | before (production edge) | after (wall) |
|---|---|---|
| no pocket, t = 40: uniformity | dust +2% at the edge; a/b off 1 by 1.3%; central lapse 1.8% above the edge's | every field uniform to 1e-12; a/b − 1 ≤ 8e-12; α − 1 ≤ 1e-13; to t = 70, through the turnaround |
| no pocket: H(t) against the cosmic-lane ODE | not available (drifts) | 7e-6 H0 at dt = 0.02, second order in dt |
| pocket without D, t = 40: α(0) / α(r = 12) | 1.335 / 1.354 (running away) | 0.993 / 0.998 |
| same: dust at the edge | +17.5% pile-up | −0.1% against r = 18; last cell against its neighbour 1e-6 |
| same: a/b at the last cell | 1.28 and growing | 1.009 (1.015 at t = 70) |
| full pocket, t = 40: α(r = 12) | 1.173 (r_max 24), 1.038 (36) | 0.998 (24), 0.996 (36) |
| constraints at the edge | not zero by construction | \|H\| ≤ 4e-7, \|M\| ≤ 3e-8 (r_max 24) |
| baryon rest mass | 6% mass loss recorded when D reached r_max (backlog 7) | conserved to 1e-15 |
| run length | runaway from t ≈ 24 (r_max 24) | t = 70 completed in every dr = 0.1 run |

What the change would involve: ghost cells by zero-slope quartic for even fields and W1b for the trace-free ones;
mirror flux for dust at the last face; removal of the incoming-light overwrite of Aa at the last cell; the lapse
solved with zero gradient at the wall and K̇ = ⟨Q⟩; a conformally flat initial slice with K² = 24π⟨ρψ⁵⟩/⟨ψ⁵⟩.

To settle before any change:
- Closure (ii), K̇ = ⟨Q⟩, is my reading of GEAR-107's "K fixed by the cell average". Chat to confirm. (Closure (i)
  gives the same slices; the choice is only the clock.)
- The wall is a different physical statement from the old edge: identical neighbours in step, not an isolated pocket
  in an open background. The expansion rate becomes an output (the pocket's share of the mean density is 1.8% at
  r_max = 24, 5.9% at 16). Existing production runs are not reproduced by it; they remain PROVISIONAL (AUDIT §C).
- T2 and T4 pass for r_max ≥ 24 only. Smaller cells are legitimate but are denser universes.
- The initial slice changes (the areal slice of #475 cannot meet the wall).
- T1's accuracy is 7e-6, not GEAR-76's 2.5e-11, at these time steps.
- The centre, not the wall, is now the weakest place (cusp, comoving grid).
