# N1 fix: toy-S finite-k scan (research only)

`python research/pumping_finite_k_toyS_v1.py` (base main 29eac01, ~27 s). Edit: k_res^2 = 2 - 0.375 x0 at both scan sites; docstring corrected; output-only u* / edge flags added. Design and match rule in the coordination log first (2026-11-22, N1 FIX: DESIGN).
**Result: MATCH** with #538's archive (toy-s-rerun-2026-11-19/summary.txt): Gamma0 0.237 / 0.283 / 0.353 (archive 0.2372 / 0.2832 / 0.3527); the six compared rows (x0 2.0, 2.67; H 0.003, 0.0057, 0.01) agree to the printed digits (G 1.98e5, 678, 40.1; 3.05e9, 1.23e5, 863; archive 1.979e5, 678, 40.13; 3.052e9, 1.232e5, 863.4). The archived u* (1.125-1.3) lies inside the script's 0.7-1.4.
Scan flags: gamma0's u* 0.98 / 0.96 / 0.94, interior at all three x0. Main scan interior in 20 of 21 rows; EDGE (u* = 1.4, scan truncated) at x0 2.67, H 0.11, which is not an archived row. The script's own E-tolerance says NO for x0 1.6 at H <= 0.02 (not in the archive).
Per Aster's note this is a diagnostic: the match shows the corrected script reproduces the tables, not who ran what.
