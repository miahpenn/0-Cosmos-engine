# W1 initial-data entry (supplied_slice): checks

Go: Jeremiah, "Go on W1 entry" (projects channel, 2026-11-24 00:26 UTC). Design: Lumen², DESIGN-W1-initial-data-entry-Nov23-2026 (Drive 16eY2t9Mdz4NxUpY_SHcucH7lCyQOxa9bMVtkOhKcars). Gates from its section 6, posted in the thread before running, unchanged. Zero-step only; no evolution, no rerun, no eight-node integral, no uniform-phi guard (design 5b, open).

| Check | Result |
|---|---|
| C1 entry vs main 3213444 | **bit-identical**: every state array, scalars, matter and scalar_projection, 4 structures x dr 0.05/0.025 (dump sha256 80ac3a56… both) |
| C1 vs pre-#572 1a49290, max rel <= 1e-13 | **FAIL**, from #572's phi_r floor, not this patch: Aa (itself round-off, max 1.3e-15) moves by 3.4e-17 (2.9e-2 relative); hand_floor_reference dr 0.025 Lambda/B 4.6e-12 relative. Scalars bit-identical. |
| C2 scalar_projection | **PASS**, bit-identical vs main and vs 1a49290 |
| C3 round trip (rescaled_pocket as SuppliedSlice) | **PASS**, bit-identical at both dr |
| C4 momentum residual <= 1e-12 (W1) | **PASS**: 2.05e-16 (dr 0.05), 1.32e-16 (dr 0.025) |
| C4 j_S, j_D opposite sign at max abs(j_S) | **FAIL**: both negative (r 2.475: j_S -1.553e-9, j_D -4.88e-10; r 2.4625: -1.549e-9, -5.18e-10) |
| C5 refusals (8 cases) | **PASS**, each ValueError |
| C7 max abs(B6 - B5) <= 1e-13 (W1) | **PASS**: 0.0 at both dr (steps 5.8e-5, 3.0e-10, 1.4e-15, then 0) |

Not run: C6 (needs the eight-node integral in the engine), C8 (no archive comparison made), C9 (reruns; new go).

Discrimination probes (each must be flagged, all were): C1 uniform-phi-guard mutant vs main (bitwise flags) and H0 x(1+1e-12) (4.7 relative); C2 a 1e-10 change in one cell of S; C3 S x(1+1e-12); C4 residual without the P_S S_r term 3.9e-8; C7 second areal step 3.0e-10; C5 a valid slice and a matching H0_chart are accepted.

Inputs: W1 members `native.npz` (Aster, B7-W1-grid-aligned-comparison.zip, Drive 1OUObtOwD7B51pDMiTacKJEjbavVBqv6D; same file on branch aster/w1-profile cfae43c), sha256 16151c445c9dcfcc34cf61e28e268e991ab2bfea629beaa3aa52fadc4c5a98ea, r_max 60, phi = phi0 and Pi = pi0. The original run passed `H0_chart=0.0`. The original export script `native.py` in the linked archive applies the Hubble chart term to both momenta, using `H0=0.0056689681106181795` from the V10 `report.json` settings. The reproduction script now passes that value as provenance. On both grids the corrected entry accepts the original arrays and produces bit-identical geometry to the historical zero-metadata run. The input momenta already contain the chart term. edge_tolerance 1e-7 (members reach 9.3e-8 in the last five cells).

Reproduce (repo root, OPENBLAS_NUM_THREADS=1): `python3 research/w1_supplied_slice_entry/dump.py X.npz` in checkouts of 1a49290 and 3213444, then `python3 research/w1_supplied_slice_entry/checks.py before.npz main.npz native.npz outdir`. Python 3.13, NumPy 2.5.3, SciPy 1.18.1; vendor d6052d6.


## Review amendments

Amendment 1 uses current main for C1. The pre-#572 comparison and
the original C4 sign failure are retained in the original record. Its
replacement `max |j_S| > 1e-8` also failed. Amendment 2's 1000-fold
omission test passes post hoc, with ratios about 1.9e8 and 2.9e8.
This is evidence about the entry and its sensitivity to the signed
source. W1 admission and production reruns remain pending.
