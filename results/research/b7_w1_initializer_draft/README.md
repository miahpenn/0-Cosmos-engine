# W1 through the production initializer — opt-in draft

This extends Aster's eight-node draft from base `229de812833849a8e57f597fa12b33fa20bd3df8`.
The production initializer now accepts a `W1ArealProfile` only under
`initial_structure="w1_research", outer_edge="open_edge"`. The four S/D fields
are already mapped to the engine chart; the input cell centres must exactly
match the requested grid. Nonfinite, complex, wrong-sized, shifted-grid,
missing and conflicting inputs are rejected. Data are copied, not resampled.
The archive COSMOS background and rest fluids are retained. Default starts,
wall initialization, evolution equations, vendor source and admission flags
are unchanged. Eight-node integration remains opt-in.

The saved signed-momentum draft is included because W1 has nonzero PS. One
implementation correction was necessary: phi is assigned exactly constant
in this initializer, so its radial derivative is supplied as exact zero.
Numerically differentiating that constant had introduced a ~2.46e-17 stencil
residual and changed Aa by ~8.67e-19 in the background control. The first
combined full-suite run had 3 failures, 1175 passes, 2 skips and 2 expected
failures; the unmodified base passed the affected controls. No tolerance or
historical research clone was changed to resolve this regression.

## W1 evidence

`native.npz` and `members.json` are recovered unchanged from
[B7-W1-grid-aligned-comparison.zip](https://drive.google.com/file/d/1OUObtOwD7B51pDMiTacKJEjbavVBqv6D/view).
The NPZ SHA-256 is `16151c445c9dcfcc34cf61e28e268e991ab2bfea629beaa3aa52fadc4c5a98ea`.
These are the previously solved aligned members, phase zero with the Hubble
chart term on; hatted outer radius 120.0064394643, not the archived 120.
No new harmonic-balance solve is performed.

Run `OPENBLAS_NUM_THREADS=1 python -m research.b7_w1_initializer_check`.
It checks the input and vendor hashes, calls the production kernel initializer,
and compares with the frozen research clone using the same eight-node rule.
All four fields and all nine checked geometry/gauge arrays, including the
CMC lapse, are bitwise identical on both grids. Time remains zero.

| dr | Official full-grid H maximum | Peak cell | Momentum maximum | Official G4c / PC2 |
|---|---:|---:|---:|---:|
| 0.05 | 6.483244797679359e-11 | 38 | 1.4447629651485532e-16 | 62.4984508876 — FAIL |
| 0.025 | 9.542197598483987e-12 | 0 | 1.9639900534604091e-16 | 0.9700073319 — PASS |

The PC2 control uses the same candidate initializer and all cells. This is an
integration check, not a new G4c acceptance criterion. The coarse failure
stands. The previous direct-areal curvature evidence is in the
[origin-repair handoff](https://docs.google.com/document/d/1Fuoku6fJEZQcUiR4BM43Ore9PtYka8vab_XD31dUUdk/edit);
no curvature operator is changed here. No W1 evolution or scientific
admission is claimed.

The checked runtime is Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, OpenBLAS
0.3.30 SkylakeX, one thread. This differs from the original report's runtime;
the stored W1 official H values nevertheless reproduce exactly. The report
records the actual source SHA and libraries; do not substitute the original
runtime label. Vendor remains `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

CI runs the complete repository suite, then this two-grid zero-step check,
and uploads the JSON result. Existing control and production-step smokes
remain in the workflow. W1 itself is never stepped by this check.

Local validation of executable source `bf000b3`: **1190 passed, 2 skipped,
2 xfailed** in 63.33 seconds. The warning is the installed anyio plugin's
pytest deprecation. `engine.run_controls` completed, and the workflow's
production smoke advanced the default N=40 state by one step with finite
positive geometry and finite diagnostics. The targeted input/clone checks
passed (27 tests). Logs retain the initial failed attempt and corrected run.
Hosted exact-head CI is separate; local success is not a hosted CI pass.
