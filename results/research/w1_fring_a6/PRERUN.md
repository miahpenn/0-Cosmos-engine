# f_ring A6 runs: declarations before any A6 step

Go: Jeremiah, "Go on A6", Nov 24 2026 07:11 UTC (relayed). Design: f_ring A6 (1Pc5Md4LCSX7_8ZHY3fAQDbsxl2soth5mPDMsZ69sOZ4) with A7 (1dr_BKBvWuk5nh_wVMMFOPSiPtAAAfpstJg9NHypnR9U); the docs rule. Engine at main as it stands; `origin_trace.run` refuses to start unless `engine/` and `vendor/` match #584's source. Team log not used. Script `research/w1_fring_a6/run.py`, gates `gates.py`, discrimination cases `synthetic.py`.

## Arms (A6 §4)
| arm | field | dr | lambda | dt | steps to t 3.291 |
|---|---|---|---|---|---|
| A1 | ring (corrected W1 slice) | 0.05 | 0.03 | 1.5e-3 | 2194 |
| A2 | ring | 0.05 | 0.015 | 7.5e-4 | 4388 |
| B1 | control (background only) | 0.05 | 0.03 | 1.5e-3 | 2194 |
| B2 | control | 0.05 | 0.015 | 7.5e-4 | 4388 |

No lambda 0.0075 arm: neither A6 nor A7 names one. Open edge, eight-node start, KO off, as the smoke. State saves (A7) at t 0.015, 0.8235, 1.6455, 3.291 in every arm: full geometry, scalars, matter, clock, plus the ledger totals.

## Recorded
- Every step: production H and M at cells 0-5.
- Each save, cells 0-5: production H; direct seven- and nine-point H; the connection term X^2 D_r C; adjusted direct H and production minus adjusted (the residual); B - 1 = X^2/a - 1 and U - 1 = sqrt(b)/X - 1; H terms (production Ricci, 16 pi rho, kinematic); M and its terms (D_r A_a + 3A_a/r, -(2/3) D_r K, other geometry, 8 pi j).
- A7 2(a) in-run ledger, every step, every arm: X (tangent, increment arithmetic, write; X is advanced linearly, so it has no integration piece), D = ln a + 2 ln b and T = (ln b - ln a)/3 (tangent, integration, increment arithmetic, write), with the conformal_ledger and x_update_ledger arithmetic (#587/#588). Each save checks the ledger telescopes to the stored endpoint.

## Input checks first
`research/w1_static_slice/checks.py` re-run on the corrected slice and its five faults before any step: the corrected slice passes on both grids; every fault is refused (`inputs.json`). Each A1/B1 run also checks its first 10 accepted steps row-for-row against the pinned origin trace (`results/research/w1_fring_origin_trace/{A,B}-0.05.json`).

## Gates and readings (nothing read before t 0.015)
1. **First check:** A1 direct nine-point H at cell 0, t 0.015, must equal 8.5030e-11 (areal trace) within 1e-3 relative. If it fails, nothing else is read.
2. **A6 3(i):** each arm, each save, cells 0-5: |H| <= 1e-3 (|Ricci| + |16 pi rho| + |kinematic|); |M| <= 1e-3 x the sum of its four terms' magnitudes. "Largest lambda meeting (i)" = the largest ring lambda passing at every save.
3. **A6 3(ii) under A7 2(b,c),** per dt pair (A1-A2, B1-B2) at each save, cell 0, nine-point (seven-point alongside):
   - Rate part = (coarse minus fine direct H) minus the write and increment-arithmetic curvature responses (X, D, T) of that pair, each on the frozen common backgrounds of the ledgers: mean a, b with the fine X for X; fine a, b with mean X for D and T.
   - Closure magnitude = |X closure| + |D,T closure| + |interaction| (full difference minus X-only minus D,T total). Floor = 10 x max(closure magnitude, one binary64 spacing of the rate part).
   - Below the floor: "unresolved". Otherwise the factor = (fine H + rate part) / fine H. Pass if >= 2; q = 1 + log2(factor) reported.
   - Raw production-H differences are reported, labelled as including writes (A7 2(d)).
4. **A6 3(iii),** per arm: stop if any value is non-finite, OR the growth is exponential, OR (ii) fails, resolved, at the last save while the arm grows.
   - Growth uses max |H| over cells 0-5 at every step from t 0.015 to 3.291. It fits ln h against t (exponential) and against ln t (power law).
   - "Exponential" = the exponential fit has the smaller rms AND its slope x span >= 1. "Growing" = h at the end > h at t 0.015.

## Short test before launch (t to 0.03 only, scratch, not results)
All four arms ran to t 0.03 and the combine step computed every reading; the A1 first check read 8.50297e-11 (gap 3.7e-6). No gate, threshold or definition above was changed after the test. The test showed the control's M and its terms, and the control's dt-pair H differences, sit at binary64 roundoff (~1e-14). So control 3(i) for M and control 3(ii) readings are reported as written, labelled "control at roundoff", and nothing is chosen from them; alongside each 3(ii) reading the report gives eps x sum|H terms| at cell 0 for scale. The lambda pick uses the ring only, as declared.

The 17 synthetic cases in `synthetic.json` show each gate passes and fails where it should.
