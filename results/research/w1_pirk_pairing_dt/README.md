# W1 PIRK pairing, timestep refinement series (code only, not run)

Design DESIGN-W1-PIRK-pairing-dt-refinement-Nov25-2026 (Drive 1MZrugyHkzJWhH3OeTl_VoW9LJjwaZEv7FAS3iRsw8GI), Amendment 1
(1u-9OTudwOuI74yEXIs8rjxIcYZrVYCwCA4F-nxjYVAo), Amendment 2 (1sgCtsP0--_qMa-fsjUDf8xF1_9ylNyAO1uIh0JVk2Y8).
No data from this series exists. The run needs Jeremiah's own go.

Driver and scorer: `research/w1_fring_smoke/pairing_dt_series.py`. Pins: engine 05c2a09a, vendor d6052d60 (checked by
`stage_momentum.provenance()` at the start of every run).

Runs: 2 grids (dr .025, .0125) x 2 arms (P production, Q primary-L2 endpoint for Aa and K) x 4 levels
(2dt0 20 steps, dt0 40, dt0/2 80, dt0/4 160, T = 0.015) x (plain, ±1-ulp floor on a, b, X, seed 20261125) = 32 runs,
1200 steps per grid.

Cost (estimate, not measured on this series): about 2.4 s/step at dr .025 and about 5 s/step at dr .0125, so about
8,900 CPU-seconds; about 40 to 45 minutes wall on 4 cores, fine-grid runs first.

When the go is given (engine venv active, Python 3.12/3.13 with the pinned packages):

    OUT=results/research/w1_pirk_pairing_dt/out
    results/research/w1_pirk_pairing_dt/run_dt_series.sh $OUT       # 4 workers, longest first; then the scorer writes $OUT/score.json

Input checks: C0 (P at dt0, step 1 reproduces the archived stages and end digest), C1 (step-1 primary predictor of P
and Q bitwise, every field), C2 (log2 of the band gap ratio max|X_pred - X_expl| per level pair >= 1.7; 2dt0 ratios count
only if the 2dt0/dt0 pair passes, A2.3), C3 (Q: recomputed primary L2 equals the substituted value bitwise on every step;
one failure stops the run). Gates G0, G1, G2 as in the design; a failed or non-finite 2dt0 run reads UNDETERMINED (A2.4).
