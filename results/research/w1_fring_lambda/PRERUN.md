# dt-versus-dr origin trace: pre-run declarations

**Go:** Jeremiah's go, "Got it, thank you. Go", Nov 24 04:43 UTC in the Lumen² channel, answering "Go for both runs?". The coordinator verified it and relayed it.

**Input:** the corrected slice from #581, `results/research/w1_static_slice/native_static.npz`, with the same label and edge tolerance as the f_ring smoke. Engine main is unchanged since 6b5728d.

**Script:** `research/w1_fring_lambda/run.py`. It reuses the smoke's `build` and `residuals` unchanged.

**Runs:**

| run | arms | dr | dt | λ = dt/dr | steps | to t |
|---|---|---|---|---|---|---|
| 1 | A and B | 0.025 | 3.75e-4 | 0.015 | 40 | 0.015 |
| 2 | A and B | 0.05 | 7.5e-4 | 0.015 | 20 | 0.015 |

Both arms use the open edge and the eight-node start. Arm A is the ring; arm B is the "none" control.

**Saved every step:**
- pointwise H and M at cells 0 to 5;
- max|H| and max|M|, with their cell index and R.

Nothing beyond this. There is no gating.

**Comparison already on record:** the f_ring smoke ran at λ = 0.03. Its stored arm A max|H| at t 0.015 sits at cell 0 on both grids:

| dr | dt | max\|H\| at t 0.015 |
|---|---|---|
| 0.05 | 1.5e-3 | 1.431e-10 |
| 0.025 | 7.5e-4 | 1.410e-10 |

**Reading** (predeclared by Lumen² at 04:40 UTC, recorded here before any step):
- **If λ = dt/dr governs:** the two ring-arm cell-0 H values agree within 20% and sit well below ~1.4e-10.
- **If dr governs:** the dr 0.05 run stays near 1.4e-10.

**Also reported:**
- ring minus control at cell 0;
- a comparison against Aster's unverified lead of 2.706e-11 at (dr 0.025, dt 3.75e-4). No source for that number exists in the repo.
