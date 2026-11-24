# Amendment 4 predictions scored on stored data

Lumen²'s Amendment 4 (Drive 1tfEU7VNNJvFOx5acfRLcmA0e9D5J6_8zp-9x8K20t94) is scored as written, using the methods and thresholds in PRERUN.md, which were fixed before any score was read.
Data: A6 arms, the origin trace, and the #598 replay samples. No simulation and no engine change.
Script: `research/w1_fring_a6/amendment4_tests.py`. Full numbers: `amendment4.json`.

Input checks all pass:
- per_step equals the saves bitwise;
- the origin trace equals A1 at t .015 bitwise;
- the replay samples reproduce the #598 band M bitwise;
- synthetic PASS and FAIL cases are discriminated.

| | Score | Key numbers |
|---|---|---|
| P1 (M ≈ -(r/4) d_tH, \|Q\| ≤ .5) | **PARTIAL** | t .8235 and 1.6455: signs agree 16/16, \|Q\| ≤ .43 in 13 of 16 (cell r .125: 1.2-7.4). t 3.291: signs opposite at r .125 and .175 in both arms, \|Q\| 0.85-14.5. Overall 9/24 have \|Q\| > 1, 20/24 signs agree. All steps t ≥ .8235: median \|Q\| .86 (A1) / .73 (A2). |
| P2 (S_H A2/A1: ~.25 early, .8-1.25 late) | **PASS** | band max: .33 at t .015 and .99 at t 3.291. Per cell: early .08-2.5 (noisy), late .99-1.24. |
| P3 (cell-tied S3 vs space-tied S2, dr .05 vs .025, dt .00075, t .015) | **NEITHER AS WRITTEN** | dr .05 band H is a flat ~-2e-11 floor across cells 1-5 (peak at cell 3, r .175). dr .025 peaks at cell 1 (r .0375). Cell-1 growth .87, peak growth .79. |
| P4 (centre step) | skipped | needs a new run |
| P5 (band M grows with band d_tH, growths within 2x) | **PARTIAL** | Direction agrees in 4/4 cases (both grow; never M up with d_tH down). Magnitude is within 2x in 1/4: dt .000375 centred 1.98 vs 1.91 (in), backward 2.0 vs 5.3 (out); dt .0001875: 3.7 vs 1.7 and 3.6 vs 1.4 (ratios 2.1, 2.6, just out). |

Supplementary (not scored):
- **P3, H(t) - H(0).** Removing each grid's own t = 0 residual leaves a peak at cell 1 on both grids, with cell-1 growth 7.5x. That is S3's signature. The as-written comparison is masked by the initial-data floor.
- **P1/P2, late times.** At t 3.291, S_H ≈ d_tH: cell 1 S_H is -1.2e-7 against 4M/r of 1.3e-8. Late band H growth is mostly sourced, not M-driven, which is why the P1 link weakens late.
