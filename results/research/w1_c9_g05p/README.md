# C9 W1 input-check rerun (AMENDMENT 2): results

Run Nov 24 2026 ~02:28 UTC on Jeremiah's go (02:25 UTC). Script: research/w1_c9_g05p/run.py. Full record: report.json.
Declarations and Aster's clarification: PRERUN.md. Zero-step arrays only, no evolution.

**Outcome: input checks PASS on both grids.** native.npz sha256 16151c44...a98ea and the engine blobs are unchanged
from #576, so under Amendment 2, 5.2 the #576 readings stand and were not recomputed (G4c Run B 0.68x with K0-K5
passing, G1, E_ring, G4.2, G4.4, gate V).

| | dr 0.025 | dr 0.05 |
|---|---|---|
| G0.1-G0.4 | pass | pass |
| Entry (edge_tolerance) | accepted, 1.943e-7 | accepted, 1.949e-7 |
| last-5 max: D-D_s / PD | 2.86e-8 / 8.52e-8 | 3.20e-8 / 9.25e-8 |
| last-5 max: S / PS | 8.7e-10 / 2.3e-9 | 8.3e-10 / 3.3e-9 |
| G0.5'a P_tail / 4.207e-10 - 1 (<= 0.03) | +0.0161 | +0.0036 |
| G0.5'b m_D, m_S (<= 0.01) | 3.2e-4, 7.5e-4 | 1.3e-3, 3.0e-3 |
| G0.5'c last-5 dev D, S | 6.4e-4, 6.2e-4 | 1.2e-3, 3.1e-3 |
| G0.5'd D, PD vs bound | 2.86e-8 <= 3.85e-8; 8.52e-8 <= 1.01e-7 | 3.20e-8 <= 3.83e-8; 9.25e-8 <= 1.00e-7 |

Recorded only (no band): P_S,tail 5.44e-13 / 5.28e-13; f_snap(0) at r_cut 25: 0.6621 / 0.6618, at r_cut 60:
0.6712 / 0.6693 (a single-slice snapshot, not the period mean f_ring; archived 0.7415 is Read, Amendment 2, 3.3);
max |PD| 3.25e-5 / 3.24e-5 at cell 0; old G0.5 ratios as in report.json.

The #576 G0.5 FAIL stays on record. The fault-test baseline (the Amendment 2 discrimination check) read the same
unmodified arrays before this run; the edge_tolerance factor was fixed after it and before this run.
