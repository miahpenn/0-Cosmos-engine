# B7 W2 cell 0-1 junction test, rerun r1 (RESEARCH ONLY): NOT THE JUNCTION

Go: Jeremiah, 2026-11-23 00:58:37 UTC ("Go on the rerun"). Pre-run entry in the coordination log, 2026-11-23 ~01:05 UTC.
Only one thing changed from r0 (`../b7_w2_junction_r0/`): gate V's r^6 order check now uses cells with r <= 5, because
the r0 check over [0, 60] sat at the float64 floor. The change was logged before any swap result existed. The post-hoc
r <= 5 numbers from r0 had been seen; that is disclosed in the entry. Everything else is the same as r0: research clone
only, dr 0.2/0.1/0.05 at r_max 60, uniform dust, Hubble K_t, PS terms on, n_iter 6, the original G4c threshold against
the engine's PC2. The r0 gate reproduces with the code as changed. Run: `python research/b7_w2_origin_check_v1.py RING.npz OUT junction_r1`.

**Gate V: PASS** (exact through r^4 in every cell; r^6 order 5.99, 5.99 on r <= 5; cell1/cell0 = 3.0).

**Reading (predeclared): NOT THE JUNCTION.** Ring H01 with Q6 / engine = 0.957, 0.825, 0.714 (all >= 0.5).

| | dr 0.2 | dr 0.1 | dr 0.05 |
|---|---|---|---|
| ring/PC2, engine -> Q6 | 160x -> 153x | 42x -> 35x | 46x -> 33x |
| G4c with Q6 | FAIL | FAIL | FAIL |
| manufactured H01: engine / Q6 / exact B | 1.2e-7 / 1.1e-7 / 1.4e-8 | 1.0e-8 / 7.1e-9 / 9.0e-10 | 2.8e-9 / 4.4e-10 / 5.9e-11 |
| manufactured Q6 / exact B (needed <= 2 for "implicated") | 7.9 | 7.9 | 7.5 |
| B(Q6) - B(engine), cell 0 | -2e-16 | -2e-16 | -1e-16 |

Observations (not gated):
- Q6 and the engine give the same B in cell 0 to roundoff. Both reduce to an even quartic over [0, r_0]. The defect
  stays at cell 0 with either. The exact-B arm still removes it. So what is left sits in the first piece [0, r_0], or
  in what the curvature stencil does with it, and not in the junction or the trapezoids.
- PC2 rebuilt with Q6 drops 3-66x (1.0e-9 -> 3.5e-10; 6.9e-11 -> 1.0e-12). That is context only; G4c keeps the
  engine's PC2.

G4 overall stays UNTESTED. 3b (u = 1.0144) waits on W1; the -P_S S_r production change stays on hold.
