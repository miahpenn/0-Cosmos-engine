# Reference CMC lapse centre-row repair (2026-11-08)

`engine/true_cmc_reference.cmc_lapse` carried the same centre row that #436 removed from
`engine/cmc_gauge.py`: α₀ − α₁ = 0 with zero right-hand side. A smooth even lapse has
α₁ − α₀ = α″(0)h² there, so the row leaves an O(h²) kink at the origin. This solver is used by
`engine/coupled_engine.py` and `engine/bidirectional_production.py`, not by the V5.5 production kernel.
The repair is the same even-ghost row (α₋₁ = α₀, the grid's parity-ghost convention in
`reference_pirk_unified.Grid._d4`); diff in `centre_row_fix.diff` (SHA-256 `eb7c1dea…2231`).

## Evolution comparison

`tools/true_cmc_reference_centre_row_evolution_comparison.py` (SHA-256 `09ba96e5…add4b`), run from
the fix branch on `main` at `c861ca1` plus the change (`03fb985`, #441, touches no file used here);
Python 3.12.15, numpy 2.5.4, scipy 1.18.1. The old row is loaded from `main`'s file; nothing else
differs. `q.make_initial(g, 0.01, 7.0)`, r_max = 40, dt = 0.06·dr, real `step_true_cmc` steps to
t = 0.60. Output `compare.json` (SHA-256 `851d0a57…c677`).

Cells 0–2 max, old → new:

| | N = 40 | N = 80 | N = 160 |
|---|---|---|---|
| M | 1.54e-4 → 1.68e-4 | 4.98e-4 → 2.75e-4 | 9.65e-4 → **9.12e-5** |
| H | 4.01e-4 → 2.43e-4 | 1.35e-3 → 6.38e-4 | 1.04e-2 → 1.36e-3 |
| C | 1.64e-4 → 1.85e-4 | 4.83e-5 → 4.96e-5 | 5.65e-4 → 6.12e-5 |
| α₁ − α₀ | 0 → 8.4e-6 | 0 → 2.3e-6 | 0 → 6.0e-7 |

- Old row: the centre momentum and Hamiltonian constraints grow with resolution. New row: centre M and C
  at N = 160 are about 10× lower than with the old row; the lapse is smooth and even (α₁ − α₀ ∝ h²).
- Away from the centre: annulus 5 ≤ r ≤ 10 values agree to two digits; final-field differences at
  r ≥ 5 are ≤ 9e-7, 1.3e-7, 1.7e-8 (N = 40, 80, 160), falling with resolution.

Not fixed here: the centre H still grows with resolution on this path (2.4e-4 → 1.4e-3). The reference
initializer `reference_pirk_unified.make_initial` integrates (r B)′ = src with a plain trapezoid from
r = 0 (line 121), the same kind of O(h³) origin constant that #439 removed from `engine/v55_initial.py`.
That module is kept as an audit/reference artifact, so changing it is a separate owner decision.
