# B7 W2 cell 0-1 junction test, run r0 (RESEARCH ONLY): STOPPED at gate V

Follow-up to the origin check (PR #561). Go: Jeremiah, 2026-11-23 00:54:35 UTC ("Go on the joint test"). Design and
criteria fixed before the run: coordination-log entry 2026-11-23 ~01:00 UTC. Code: `even6_integral` and
`main_junction` in `research/b7_w2_origin_check_v1.py` (run with a third argument `junction`). Engine/vendor read-only.

**Result: gate V FAILED, so by the predeclared rule the swap into the research clone was not run.** No reading on the
junction.

| Gate V check (all grids dr 0.2/0.1/0.05, r_max 60) | Result |
|---|---|
| exact for 1, r^2, r^4 in every cell (rel <= 1e-12) | PASS (max rel 2.5e-14) |
| cell 1 / cell 0 r^6 error on dr 0.05 <= 10 | PASS (3.0) |
| r^6 max error over all cells falls with order >= 5.5 | **FAIL**: orders 5.57, -1.0 |

Disclosures:
- The failed check is a flaw in the gate, not evidence about the quadrature: the r^6 maximum over the whole domain
  sits at r of about 52-60, where F = r^7/7 is about 4e11 and the error is at the float64 floor (relative 9e-16 on
  dr 0.05). An absolute maximum on [0, 60] cannot show order there. The criterion was not changed after the run.
- Post-hoc, not gated: over r <= 5 the r^6 error falls 7.1e-4, 1.1e-5, 1.8e-7 (order 6.0, 6.0).
- Post-hoc, not gated: for pure r^6 the engine quadrature has the same cell 0 error as Q6 and no large cell 0 -> 1
  jump (dr 0.05: 8.9e-10 -> 2.0e-9). The r0 jump (8.9e-13 -> 2.1e-10) came from the r^2/r^4 part of that test
  polynomial: the engine's corrected trapezoids are not exact for r^4 beyond cell 0. Q6 is.

Nothing here changes W2 or the origin-check results. G4c stays FAIL. 3b (u = 1.0144) waits on W1; the -P_S S_r
production change stays on hold.
