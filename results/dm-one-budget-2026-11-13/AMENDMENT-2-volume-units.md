# Amendment 2 to the dark-matter one-budget results (PR #475): absolute energies are per steradian

Filed 14 Nov 2026 with the HANDOFF Nov 14 runs (`results/balance-drain-2026-11-14/`).

The gravity harness (`research/dm_one_budget_v1.py`) integrates with the production grid's cell volumes,
`grid.volumes` = ∫ r² dr over each cell, without the factor 4π. Every absolute energy it reports (S's energy, the
dark-matter energies, the rates in table (b)) is therefore per steradian: multiply by 4π for the total.
Ratios are unaffected (every quantity carries the same factor), so all ratio statements of PR #475, of
Amendment 1 and of PR #476 stand.

One statement of PR #475 used an absolute comparison and is corrected here:
- "D's ejected vacuum … releases 1/(2λ) per unit volume where D falls into its well (≈ 0.2 inside R = 8, several
  times E_S)". With the 4π restored, S's energy at the start is 0.50 (x_c = 1.5) and 0.73 (x_c = 2.0), so the
  released vacuum energy inside R = 8 (0.21) is ≈ 0.3–0.4 of S's energy, not several times it.
- Amendment 1's "E_S0 = 5.6852e-02" is per steradian (0.714 in total).
