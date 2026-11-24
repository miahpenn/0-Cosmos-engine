# Corrected W1 slice (static time) and its input checks

W1 derivation AMENDMENT 3 (Drive 1JU2EhMUs7dSio23VVu0WIPMt33iAY6u9gbD9EKE64aY), with the build path from AMENDMENT 4
(Drive 1J_dGcUVzAy9jceOwxByyvuQ1pyix9hHiC9owyoRNItM). Jeremiah's go: "Go", Nov 24 03:26 UTC (Lumen² channel).
Zero-step only: no evolution, no f_ring run.

## Build

`research/w1_static_slice/build.py` passes each grid's stored harmonic-balance vector (admitted native.npz, sha256
16151c44…a98ea) through the producer's own `mapped()` from Aster's native.py, unchanged. The solver's snapshot phase
is set per node to theta = omega_u H0 r^2/2 on every solver node (2383 / 4766 nodes, out to r about 119).
No Newton re-solve.

- **Positive control (A4 1):** theta = 0 reproduces native.npz **bitwise** (byte-equal S, PS, D, PD) on both grids. PASS.
- **Corrected slice:** `native_static.npz`, sha256 339652a2a7e1379b19812893e43bafcf9ab0211422f2b8685c951107b863a295.
  The per-grid array hashes and the provenance are in `build.json`.

## Input checks

These are run with the C9 gate script `research/w1_c9_g05p/run.py --static build.json`, which now adds G0.5'f.
'f is the static-frame flux ratio, in band [0.99, 1.01], as in A3 4.3. 'd uses the A3 4.5 PD bound. Without
`--static`, the script is unchanged for the admitted slice and only records 'f.

| corrected slice | G0.1-4, entry | 'a P_tail rel | 'b m_D / m_S | 'c | 'd | 'e ratio D / S | 'f D / S |
|---|---|---|---|---|---|---|---|
| dr 0.025 | PASS | +0.0161 | 3.3e-4 / 7.5e-4 | PASS | PASS | 0.756 / 0.717 | 0.9989 / 1.0058 PASS |
| dr 0.05 | PASS | +0.0036 | 1.3e-3 / 3.0e-3 | PASS | PASS | 0.757 / 0.722 | 0.9996 / **1.0106 FAIL** |

On dr 0.025 every check passes. On dr 0.05 every check passes except 'f for S, which reads 1.0106 against a 1.01
ceiling.

A3 section 5 says dr 0.05 is "record only for admission, as before". The script's combined outcome line still
requires both grids, so it reads REFUSED. Which reading governs admission is left to Jeremiah and Lumen².

Cause of the dr 0.05 miss (diagnostic, not a gate): the windowed flux mean has an end-of-window term. Moving the
window's upper edge from 58.5 to 59.8 shifts F_st by about 1.2% for D and S on both grids. That term is larger than
A3 4.3 allowed for, because the static-frame wavenumber is lower on the corrected slice. The engine's derivative is
4th order, so stencil error is negligible. The paper estimate of 0.2% / 0.4% stencil error used a 2nd-order size.

## Fault cases

Every fault fails 'f (`checks.json`):

| case | fails | 'f D / S (dr 0.025) |
|---|---|---|
| uncorrected (admitted) slice | 'f only | 1.245 / 1.286 |
| v^2 piece missing | 'a, 'c, 'f | 0.982 / 0.981 |
| wrong-sign phase | 'f only | 1.489 / 1.574 |
| D tail reversed | 'd, 'e, 'f | -0.514 / 1.006 |
| S tail reversed | 'e, 'f | 0.999 / -0.429 |

The v^2-missing case reads 0.98, not the paper's 0.94. It still fails 'f, and also fails 'a and 'c.

## Reproduce

```
python3 research/w1_static_slice/build.py V10_DIR native.py native.npz BUILD_DIR
python3 research/w1_static_slice/checks.py BUILD_DIR native.npz members.json OUTDIR
```

V10_DIR is the `v10/` folder of ringkit_v10_chk-w1-slice.tgz (Drive 16LGCwbEof7j0E84GSr21J9HsKDMpzDKn). native.py,
native.npz and members.json come from Aster's W1 bundle (Drive 1OUObtOwD7B51pDMiTacKJEjbavVBqv6D).

Environment: Python 3.13.18, numpy 2.5.3, scipy 1.18.1, vendor d6052d6.

## Next (each on Jeremiah's go)

1. Decide admission (dr 0.025 passes everything; dr 0.05 'f S is 0.06% over).
2. Rerun the zero-step C9 set (G1-G4) on the new bytes.
3. Then the f_ring run.
