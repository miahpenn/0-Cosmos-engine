# Eight-node areal integral in the engine: gate results

Gates predeclared in `GATES.md` (commit 1276bf6) before any run. Engine head: v55_initial.py blob d613ab44efa6919c6de208a45ac0c2ee8559d25a, production_kernel.py blob 849014d34c532b1f142def2b9b30ae02809e81d7. Python 3.13, NumPy 2.5.3, SciPy 1.18.1 (research scripts only), OPENBLAS_NUM_THREADS=1.

| Gate | Result | Reading |
|---|---|---|
| (a) suites | **PASS** | `pytest tests`: 1144 passed, 2 skipped, 2 xfailed (main 1e7bffc: 1133, same skips/xfails; +11 new). Root `pytest` with vendor: 1201 passed (main 1190). Mutation check: shifting the stencil nodes by 1e-4 fails 8 of the 11 new tests. |
| (b) default-off baselines | **PASS** | 192 arrays (geometry, scalars, matter; four structures; n 1200 and 2400) bit-identical to main with the flag off. Flag on moves 80 arrays: largest Lambda 2.04e-8 (n 1200), 2.54e-9 (n 2400); a 5.2e-10; X 1.3e-10 (n 1200). `gate_b.json`. |
| (c) G4c reproduction | **PASS** | Official ring/PC2 at dr 0.025: **0.970007x** (ring max 9.54e-12 at cell 0; control 9.84e-12). Engine integral vs Aster's `integrate(r,f,h,8)`: max abs difference 0.0 over 36 calls (bit-identical). Negative control, engine default integral: 26.27x at dr 0.025, 25.34x at dr 0.05. `gate_c_*.json`. |
| (d) dr 0.05 | **FAIL, on record** | 62.498x, peak cell 38. |
| (e) direct curvature (obs.) | recorded | Nine-node 2.408x / 0.797x; seven-node 1.911x / 0.671x (dr 0.05 / 0.025). |
| (f) order | **PASS** | Errors 5.90e-7 / 2.66e-9 / 1.07e-11 at dr 0.2 / 0.1 / 0.05; order 7.79 then 7.96. Cumulative rule: order 4.43 then 4.06 (fails, as required). `gates_fg.json`. |
| (g) edge noise gain | **PASS** | Max 4.666 at the last row (predicted 4.67); interior 1.32. |

## Observation, not a gate here (bears on C9 G3.2)

The background ("none") Hamiltonian max on the same operator is 2.91e-12 at dr 0.025 and 8.27e-13 at dr 0.05. The rescaled-pocket control is 9.84e-12 and 1.04e-12, i.e. 3.4x and 1.25x that floor. C9's G3.2 requires the control to sit at least 100x above it for the G4c ratio to count. As written, G3.2 would void the G4c reading at both spacings; the C9 reruns are not started until that is settled.

Scripts: `research/eight_node_engine/` (`dump_baselines.py`, `compare_baselines.py`, `gate_c.py`, `gates_fg.py`). Input native.npz sha256 16151c44...a98ea (Drive 1OUObtOwD7B51pDMiTacKJEjbavVBqv6D), not committed.
