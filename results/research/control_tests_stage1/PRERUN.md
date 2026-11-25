# Stage 1 PRERUN (written before launch)

**Go:** Jeremiah, 07:27:41 UTC (relayed by the coordinator): "If the readback finds nothing shared, go on stage 1 test". The readback is PR #605. Its finding is that the W1 wiggle comes from the areal build (`engine/v55_initial.py`, lines 373–430). Stage 1 doesn't use that path (`ctl.py`, lines 31–38; `t0_checks.py`, lines 17–19).

**Code:** PR #603 (`research/control_tests/test1.py`, `observer.py`, `mp_readback.py`), as merged. No engine or vendor change.

**Runtime overrides:** `scalar_system.V0 = 0` and `SD_QUARTIC_LAMBDA = 0`, set after `initialize()` (Jeremiah, 04:46: "Set V0 and lambda to 0").

**Arms (A13):**

| arm | dr | dt | r_max |
|---|---|---|---|
| 1a | 0.05 | 3.75e-4 | 8 |
| 1a | 0.025 | 3.75e-4 | 8 |
| 1a | 0.05 | 1.875e-4 | 8 |
| 1a | 0.05 | 3.75e-4 | 12 |
| 1b (K0 −0.1, c0 0.5) | 0.05 | 3.75e-4 | 8 |

- Every arm runs to t_end 2.1, with saves at t = 0.3k and an R-stage recorded step at each save.
- Input checks run at zero steps: 1b with K0 −0.1 / −0.017 × c0 0.5 / 0.9 × dr .05 / .025 / .0125, plus the three 1a grids.

**Gates as predeclared (none changed):**
- input checks A13 §1 and B1, B3 (i)–(iii), B4, B5, B6;
- G1b.1 identity at dr .05 across the saves, with R-stage readings;
- the 1a floors (for GK.1 in stage 2).
- A15 withdrew the 1b edge tolerance and every cross-grid 1b comparison. The G1b.1 order part reads UNDETERMINED by design.

**Discriminability:** checked on paper through A9–A15. B3 (iii) reads Kdot = 0 bitwise and can fail by design at edge cells, so it is scored on interior cells only. No gate is read against data before this run.

**Self-check:** on each arm's first save, the recording step must equal the production step bitwise, or the arm stops.

**Cost:** about 1.2 h serial, on 4 cores.
