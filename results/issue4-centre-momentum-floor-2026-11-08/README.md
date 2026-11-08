# Centre momentum-constraint "floor": temporal vs spatial (2026-11-08)

Follow-up to `results/initializer-origin-quadrature-repair-2026-11-08/`, which left open why the total
momentum constraint at cells 0–2 at t = 0.60 stopped falling (5.0e-8 at N = 120, 5.4e-8 at N = 240) with
dt = 0.06·dr.

Tool `tools/issue4_centre_momentum_floor.py` (SHA-256 `9a2e04bc…3050`), run from `main` at `c861ca1`
(CMC centre row #436 and initializer origin quadrature #439 in place; vendor `d6052d60`; Python 3.12.15,
numpy 2.5.4, scipy 1.18.1). Default initial data, r_max = 40, real `kernel.step` calls to t = 0.60.
Outputs `T120.json` (SHA-256 `cb56ecf0…c096`), `T240.json` (`55a6ed5b…183f`), `S.json` (`d3ed468c…1cdd`).
M at t = 0 is round-off (≤ 7e-18) in every run.

## Total M at cell 0 (cells 0–2 max in brackets), t = 0.60

| | dt = 0.02 | dt = 0.01 | dt = 0.005 | dt = 0.0025 |
|---|---|---|---|---|
| N = 60 | | | −1.07e-6 (1.07e-6) | |
| N = 120 | −4.29e-8 (5.02e-8) | −9.80e-8 (9.80e-8) | −1.12e-7 (1.12e-7) | |
| N = 240 | | +5.44e-8 (5.44e-8) | +1.70e-8 (1.71e-8) | +7.84e-9 (1.52e-8) |

Annulus 5 ≤ r ≤ 10 at dt = 0.005: 1.25e-7 → 7.65e-9 → 4.66e-10 (N = 60 → 240), fourth order.

## Reading

- There is no spatial floor. At fixed dt the centre M falls 1.07e-6 → 1.12e-7 → 1.71e-8, and the
  dt → 0 extrapolations at cell 0 (Richardson, second order) are −1.17e-7 at N = 120 and +4.7e-9 at
  N = 240.
- The centre carries a time-stepping error that is second order in dt (cell-0 successive differences at
  N = 240: 3.74e-8, 0.92e-8, ratio 4.06; at N = 120: 5.51e-8, 1.40e-8, ratio 3.9) and has the opposite sign
  to the spatial error at N = 120. With dt = 0.06·dr the two cancel partly at N = 120 and the temporal part
  dominates at N = 240, which is the apparent floor.
- Its coefficient grows with resolution at the centre: error/dt² at dt = 0.06·dr is 1.9e-4 (N = 120) and
  5.0e-4 (N = 240). In the annulus the temporal part is not visible at these dt. So with dt ∝ dr the centre
  M converges more slowly than the interior.

Not settled: why the centre time error is only second order and grows with resolution (two resolutions
only; an order-reduction of the time integrator near r = 0, or a per-stage effect of the CMC lapse solve,
are the first candidates). An N = 480 time-step series would measure the growth.
