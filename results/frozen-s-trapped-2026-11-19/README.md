# Frozen-S control and trapped-mode check (Nov 19 2026) — RESEARCH ONLY
Script: research/frozen_s_trapped_v1.py (toy S, planar lattice, seed 1). Raw output: summary.txt. Runtime 7 s.
Design and pass/fail were fixed in the coordination log before running.

## Frozen-S control (G_k: frozen / full; linear d0=1e-8 in brackets)
- H 0.0057 x0 2.0 (lin 679): d0 1e-5 679/679; 1e-3 680/646 linear; 1e-2 417/162 BACKREACTION; 1e-1 56/18.9 DETUNING
- H 0.0057 x0 2.67 (lin 1.2e5): 1e-5 1.35e5/8.2e4 linear; 1e-3 7858/1569 DETUNING; 1e-2 614/253 DETUNING; 1e-1 72/36 DETUNING
- H 0.01 x0 2.0 (lin 40.2): 1e-5..1e-2 linear; 1e-1 29.4/15.4 BACKREACTION
Reading: D's own detuning (nonlinear mass shift) cuts gain at large d0 in 5 of 12 rows, hardest at x0=2.67; back-reaction alone dominates in 2 rows (d0=1e-2 at x0=2.0, d0=0.1 at H=0.01). Back-reaction roughly halves the gain again where both act. Frozen-S exceeds linear at x0=2.67, d0=1e-5 (1.35e5 vs 1.2e5), unexplained.

## Trapped-mode check (peak-gain mode d_k*, zero-crossing omega_D; pass: omega_D <= 1.4642 in every linear window)
- x0 2.67, d0 1e-8 (k*=1.325): omega_D 1.303 / 1.243 / 1.191 in the three windows, all below sqrt2 and tracking omega_S (1.316/1.263/1.226). PASS.
- x0 2.0, d0 1e-8 (k*=1.276): window 1 reads 1.253 (below gap); windows 2 and 3 return NaN (<3 zero crossings found). INCONCLUSIVE for those windows; not a fail, not a pass.
- d0 1e-2 rows: x0=2.0 window 1 1.251; x0=2.67 picked k*=0.589 and gives NaN / 3.54 / NaN, noise-dominated, not interpretable.
- omega_S zero-crossing also NaN in several windows; the sqrt(1+0.375x) column is the usable comparison.
Caveat: NaN cause not investigated (windows, mode extraction). Fields were not saved; a re-read would need the series saved.
