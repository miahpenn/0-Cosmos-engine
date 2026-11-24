# B7 W1 eight-node areal integral: independent reproduction (research only)

Go: Jeremiah, "Go on eight-node" (projects channel, 2026-11-24 00:12 UTC). Zero-step only. No engine, initializer or gate change; the eight-node integral is substituted in memory by Aster's scripts, exactly as she ran them.

## Gates (written in the thread before running)

1. **G4c, predeclared and unchanged** (`research/b7_w2_origin_check_v1.py`): ring H_full max <= 10 x rescaled-pocket (PC2) H_full max, same grid, same integral, official post-lapse operator, at dr 0.05 and dr 0.025.
2. **Reproduction check:** on Aster's pinned engine commit, outputs match her archived files to 6 significant figures.

Observations only: direct areal curvature ratios (7- and 9-node), the coarse peak cell, the cells 0-5 quadrature share.

## Results

| | dr 0.05 | dr 0.025 |
|---|---:|---:|
| Official G4c ring/PC2 | 62.50x **FAIL** | 0.970x **PASS** |
| Direct areal, 9-node (obs.) | 2.41x | 0.797x |
| Direct areal, 7-node (obs.) | 1.91x | 0.671x |
| Official ring peak (obs.) | cell 38, r = 1.925 | cell 0 |
| Six-node frozen-source share of B discrepancy, cells 0-5 (obs.) | 98.85% | 99.85% |

- **Reproduction check: PASS.** On engine 229de81 / vendor d6052d6, all 196 numeric values in `full_summary.json`, `orders.json` and `probe_results.json` are bit-identical to Aster's archived files.
- **On current main 893c464 (includes #572, signed momentum):** every G4c ratio, peak cell and share above is unchanged. 72 values moved, all at round-off size (largest |H| change ~1e-12 in the homogeneous background, momentum maxima ~1e-16), consistent with #572's phi_r stencil round-off. Details in `comparison.txt`.

Not decided here: whether production runs at dr <= 0.025 or G4c switches to direct curvature. The archived linearly sampled W1 input (gate U stop at dr 0.025) was not used; inputs are Aster's fresh grid-aligned members.

## Inputs and reproduction

- `native.npz` from Aster's B7-W1-grid-aligned-comparison.zip (Drive 1OUObtOwD7B51pDMiTacKJEjbavVBqv6D), sha256 16151c445c9dcfcc34cf61e28e268e991ab2bfea629beaa3aa52fadc4c5a98ea. Not committed (600 KB).
- Archived references: B7-full-grid-paired-check.zip (Drive 12CvffjxI3QEyofCShsroq_QfwYjvzcrZ) summary.json sha256 47c34b40...; B7-origin-integral-repair-study.zip (Drive 1hPjLDNkWtp5_yDuefW7uAcE3gZDZu3aM) results.json 63b914c8..., orders.json 5428d3b5..., intervals.npz 5fa1348f....
- Handoff: "B7 W1 origin repair — evidence and production handoff" (Drive 1Fuoku6fJEZQcUiR4BM43Ore9PtYka8vab_XD31dUUdk). Patch B7-eight-node-areal.patch (Drive 1YHyR8_ckKVqo4gAlzPl73EXKu99Srfzk) was read, not applied.
- Scripts in `research/b7_w1_eight_node_repro/` are Aster's `probe.py`, `order_check.py`, `full.py`, `interval.py`, unchanged; `compare.py` is the comparison. Run from the repo root, each as `python3 <script> native.npz <outdir>`, in that order, with `OPENBLAS_NUM_THREADS=1`. Python 3.13.18, NumPy 2.5.3, SciPy 1.18.1.
