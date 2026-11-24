# f_ring A6 runs: readings as written

Arms A1/A2 (ring) and B1/B2 (control), dr 0.05, lambda 0.03 and 0.015, to t 3.291, run Nov 24 2026 07:21-09:03 UTC on the pinned venv (Py 3.12.13 / NumPy 2.3.5), engine and vendor at #584's source. Declarations: `PRERUN.md` (committed before the full arms). Readings: `readings.json`; per-arm records and the A7 state saves: `arms/`.

## First check
Direct nine-point H, A1, cell 0, t 0.015 = 8.50297e-11 against 8.5030e-11: gap 3.7e-6. PASS. Every reading below follows it.

## A6 3(i), cells 0-5, every save
- H: PASS in all four arms at all four saves (largest ratio to bound: ring 0.41 at lambda .03 and 0.51 at lambda .015, both at t 3.291).
- M, ring: PASS at t 0.015 only; FAIL at 0.8235, 1.6455, 3.291 at both lambdas (worst |M| / bound: 155 at lambda .03, 135 at lambda .015, t 3.291).
- M, control: FAIL at every save; control at roundoff (|M| equals the sum of its term magnitudes, ~1e-14 at t 0.015).
- Largest ring lambda meeting (i) at every save: none (M). On H alone, both lambdas meet it.

## A6 3(ii) with A7 2(b,c), cell 0 (nine-point; seven-point agrees in every case)
| pair | t 0.015 | t 0.8235 | t 1.6455 | t 3.291 |
|---|---|---|---|---|
| ring A1-A2 | PASS, factor 3.94, q 2.98 | PASS, 3.27, q 2.71 | FAIL, 1.62, q 1.69 | FAIL, 0.84, q 0.75 |
| control B1-B2 | unresolved (control at roundoff) | FAIL, factor -1.85 | FAIL, -3.37 | FAIL, -1.97 |

Rate parts sit far above their floors at every resolved reading (floors <= 1e-14; ring rate parts 6.5e-11 to 5.3e-8). At t 3.291 the finer ring run has the larger H (production H cell 0: 2.50e-7 at lambda .03, 3.10e-7 at .015). Control negative factors mean the control's H changes sign between the two dt. Control H grows from ~3e-13 to ~1e-10.

## A6 3(iii), per arm over t 0.015 to 3.291 (max |H|, cells 0-5)
- No arm grows exponentially: the power-law fit has the smaller rms in all four (ring power slopes 1.34 at lambda .03 and 2.02 at .015; control 1.22 and 1.19). Growth over the span: ring x1.7e3 and x1.4e4, control x1.4e3 and x4.4e2.
- STOP by A6 3(iii) in every arm: no fall with dt refinement (resolved 3(ii) FAIL at t 3.291) while growing.
