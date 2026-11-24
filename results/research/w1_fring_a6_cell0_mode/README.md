# Cell-0 fast mode in the A6 data: readings by the fixed methods

Methods: `PRERUN.md` (committed before any reading). Script: `research/w1_fring_a6/cell0_mode.py`. Numbers: `cell0_mode.json`. Stored A6 data only.

- **Vendor functions:** `enforce_algebraic_regularity` and `regularized_aa_prime` are defined in vendor/bb-palatini-unified-r0/vacuum_bssn.py and have no callers.
- **(1) Period, 0.0843 +/- 5% at both dt:** FAIL on M.
  - Dominant detrended cell-0 M period: 0.1649 at lambda .03 (+96% vs 0.0843) and 0.1901 at lambda .015 (+126%).
  - Against 3.685 x dr = 0.184 they sit at -10% and +3%.
  - The two dt differ by 15%.
  - Next peaks: lambda .03 0.1949 (0.59 of top), 0.0862 (0.52), 0.2366 (0.45); lambda .015 0.1534 (0.30), 0.2501 (0.23), 0.1377 (0.13).
  - Aa and a/b were not stored per step, so (1) can't be read on them.
- **(2) Amplitude largest at cell 0, lower by cells 1-2:** PASS.
  - At each arm's cell-0 period, amplitude relative to cell 0 for cells 1-5 is 0.38, 0.21, 0.04, 0.04, 0.03 at lambda .03 and 0.49, 0.13, 0.08, 0.05, 0.03 at lambda .015.
- **(3) Cell-0 Aa/r^2 departs from the cells 1-2 trend:** FAIL as fixed (every save required).
  - It departs at t 0.8235 and 1.6455 in both arms, with misfit 2-5x the cells 1-2 spread.
  - It doesn't depart at t 0.015 or 3.291, where the misfit is 0.03-0.9x the spread.
