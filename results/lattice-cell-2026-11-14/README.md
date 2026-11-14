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

Runs in progress; this section is filled in when they finish.
