# dt-versus-dr origin trace: results

PRERUN commit: 2307beb, pushed before any step. Script: `research/w1_fring_lambda/run.py`. Data: `arm_{A025,B025,A05,B05}.json` (every step, cells 0 to 5, max and location). There is no gating. All four arms stayed finite to t 0.015.

## Headline

**Halving dt cut the ring arm's cell-0 H by 5 to 8 times on both grids.** In the smoke, at λ 0.03, it sat at 1.4e-10. Here, at λ 0.015, it is 2.7e-11 at dr 0.025 and 1.8e-11 at dr 0.05. So the origin residual follows dt, not dr.

| ring arm, t 0.015 | dr 0.05 | dr 0.025 |
|---|---|---|
| λ 0.03 (smoke) | dt 1.5e-3: 1.431e-10 (cell 0) | dt 7.5e-4: 1.410e-10 (cell 0) |
| λ 0.015 (this run), H at cell 0 | dt 7.5e-4: **1.813e-11** | dt 3.75e-4: **2.706e-11** |
| λ 0.015, max\|H\| (cell) | 6.53e-11 (cell 37, the ring) | 2.706e-11 (cell 0) |
| ring minus control, H at cell 0 | 1.790e-11 | 2.735e-11 |
| control (arm B), H at cell 0 | 2.2e-13 | -2.8e-13 |

## Against the predeclared reading

- **"If dr governs, the dr 0.05 run stays near 1.4e-10":** it doesn't. It is 1.8e-11, 7.9 times lower.
- **"If λ governs, the two ring-arm cell-0 values agree within 20% and sit well below ~1.4e-10":** both sit well below, but they differ by 33% (dr 0.05 / dr 0.025 = 0.67). So the 20% condition is not met as written.

Neither branch holds exactly as declared. The 20% miss may not mean much, because at dr 0.05 cell-0 H changes sign during the run:

| t | 0 | 0.00375 | 0.0075 | 0.01125 | 0.015 |
|---|---|---|---|---|---|
| dr 0.05, cell-0 H | -2.10e-11 | -1.65e-11 | -7.2e-12 | +4.7e-12 | +1.81e-11 |
| dr 0.025, cell-0 H | 7.1e-12 | 1.0e-11 | 4.6e-12 | 1.36e-11 | 2.71e-11 |

So a single sample at t 0.015 compares two values both at the level of the slice's initial residual (about 2e-11 at dr 0.05 on cells 0 to 5).

## Aster's lead

Aster's 2.706e-11 at (dr 0.025, dt 3.75e-4) matches this run's cell-0 H, 2.7064e-11, to four figures. Before this run, no source for it existed in the repo.

## M

M stays at or below 5.5e-12 at cells 0 to 5 on both ring arms (max 3.6e-12 at cell 1 on dr 0.025; 5.5e-12 at cell 0 on dr 0.05). The controls stay at or below 2.3e-13.

Step wall time: 2.4 s at dr 0.025 and 1.2 s at dr 0.05, both arms.
