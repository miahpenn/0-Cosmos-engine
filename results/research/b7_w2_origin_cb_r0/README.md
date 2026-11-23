# B7 W2 coherent-B reference on the W1 slice, run r0 (RESEARCH ONLY)

Go: Jeremiah, 2026-11-23 02:18:35 UTC ("Go on coherent B"). Pre-run entry: coordination log, 2026-11-23 ~02:25 UTC
(design, gates and readings fixed before the run). Order per Aster's origin-gate review. Input and settings exactly as
3b (W1 A1_hubon_dr1, np.interp onto dr 0.2/0.1/0.05, r_max 60; Hubble K_t, PS terms, uniform dust, n_iter 6, same
evaluator). Run: `python research/b7_w2_origin_check_v1.py RING.npz OUT w1_cb W1_SLICE.npz`. Not covered: any engine
change; the -P_S S_r change stays on hold.

Reference B: (rB)' = src(r, B) by DOP853 (rtol 1e-13) over an even quintic spline through the same cell-centre field
values the evaluator sees. Uncertainty: cubic spline and rtol 1e-11. Footprint of H_0, H_1: B at cells 0-5 (measured).
Arms use B_override and change only the named cells. f = 1 - |h_X - h_F| / |h_E - h_F|, h = signed (H_0, H_1).

**Gates: all PASS** (I identity and 3b replay bitwise; V reference vs analytic B on the manufactured slice, 1e-5 or
better of the engine error; U budget; C cancellation, margin ~1e8; FP footprint exact).

**Reading (predeclared): NOT ASYMPTOTIC.** The class is COUPLED on all three grids (and with the cubic reference),
but f for cells 1-5 moves 0.61 -> 0.73 between dr 0.1 and dr 0.05, past the 0.1 agreement limit.

| | dr 0.2 | dr 0.1 | dr 0.05 |
|---|---|---|---|
| f, B_0 alone | 0.44 | 0.08 | 0.04 |
| f, cells 1-5 alone | 0.17 | 0.61 | 0.73 |
| f, cells 0-5 | 1.00 | 1.00 | 1.00 |
| H01 coherent B / engine B | 0.11 | 0.056 | 0.025 |
| coherent B argmax cell | 9 | 34 | 37 (r ~ 1.9, off centre) |
| G4c, engine B (x PC2) | 48 FAIL | 24 FAIL | 25 FAIL |
| G4c, coherent B (x PC2) | 16 FAIL | 6.0 PASS | 0.96 PASS |
| B_E - B_ref, cell 0 | -1.0e-10 | -1.7e-12 | -6.6e-14 |
| B_E - B_ref, cells 1-5 | -3.4e-10 to -5.5e-10 | -2.7e-11 to -3.4e-11 | -1.6e-12 to -2.2e-12 |

Observations (not gated):
- Coherent B removes the centre defect: 89-97.5% of H01, and the ring passes G4c on dr 0.1 and dr 0.05.
- The engine's B_0 is nearly right. From cell 1 on, the engine B sits at a near-constant offset below the reference,
  so the defect is a step between cells 0 and 1. Replacing B_0 alone does little; cells 1-5 carry most of it.

Disclosures:
- Gate U at dr 0.05 passed with little margin: u_h / |h_E - h_F| = 0.098 against the 0.1 limit. Per cell, u_B at
  cells 3-5 on dr 0.05 is about the size of B_E - B_ref, so per-cell attribution inside 3-5 is not resolved there.
- dr 0.2 is pre-asymptotic (as declared) and reads differently (B_0 alone 0.44).
- Swap arms' H_full is not meaningful (step at the footprint edge) and was not gated.
