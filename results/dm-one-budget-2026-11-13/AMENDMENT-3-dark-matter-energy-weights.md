# Amendment 3 to the dark-matter one-budget results (PR #475): the dark-matter energy in table (b) was mis-weighted

Filed 14 Nov 2026 with the lattice-cell work (`results/lattice-cell-2026-11-14/`). Found while checking the cell's
zero-flux wall with the baryon rest mass.

## What was wrong

The fluid's conserved densities (`rest`, `energy_t`, `momentum_r`) already carry sqrt(gamma), which includes r².
PR #475's harness summed `energy_t x grid.volumes` (another r² dr) for "the dark matter's energy in the pocket", so
each cell was weighted by an extra r². The quantity that is conserved by the flux form is Σ (density) dr.
(Check: Σ baryon rest x dr is constant to 3e-16 in a closed cell; Σ rest x volumes drifts by 2e-3 in six time units.)

Affected: the "DM energy in" column and the ratio column of table (b) in the PR #475 README, and the same columns in
Amendment 1. Not affected: every dark-matter density statement (ρ_dm profiles, table (a), the spike), which use the
density itself; S's energy and action (proper volumes); Task 4's ring-against-dust table of PR #477 (densities).

## Corrected table (b), x_c = 2.0, β_D = 0 (rerun of Amendment 1, dr = 0.1, r_max = 24, proper gradient)

Per steradian, per unit coordinate time. Dark matter: conserved energy in r < 11 minus the no-matter control.

| window | DM energy in, as filed | corrected | S's energy, full − no-D | ratio DM/S, corrected (Amendment 1 had) |
|---|---|---|---|---|
| 0–4 | +6.7e-3 | +1.3e-4 | +2.7e-3 | 0.05 (2.5) |
| 4–12 | +1.0e-2 | +2.5e-4 | +1.1e-3 | 0.23 (9.4) |
| 12–24 | +4.3e-3 | +1.9e-4 | +1.6e-4 | 1.2 (27) |

## What it means

PR #475's reading "dark matter's energy enters the pocket 2–12× faster than S's energy changes through D" is
**withdrawn**. Corrected, the dark-matter inflow is 30–50× smaller than filed: well below S's exchange with D early
(0.05–0.23 of it) and comparable to it by t = 12–24, when S's exchange has nearly stopped. The excess dark-matter
energy in the pocket (0.015 at the start, 0.020 at t = 24, per steradian) is about a quarter to a third of S's energy
(0.057). All of this remains PROVISIONAL with every expanding run on the production edge (AUDIT §C).
