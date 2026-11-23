# B7 W1 outer-cell C = Λ − Λ_metric trace (research only)

Go: Jeremiah, 2026-11-23 22:47 UTC ("Go on the trace"), on the gate set fixed in Lumen²'s
DERIVATION-W1-outer-cell-C-predictions-Nov23-2026 (Drive 1SN10MQ6l0K4DaWt0eFRx7QJsRuUboaPJvsrGQTZbFYA),
Amendment 1 (1ju8hVnEDtIIrzAxH0lQu1QvsFGQJ1gLMEaMMz7dDU8M) and Amendment 2 (1tJGspPs8L4RF6sv5LWws-4sx6QvViRALAc00jSf_5PY).
No engine or production change. Gates were not changed after the run.

## How it was run

The driver `research/b7_w1_outer_cell_trace.py` needs Aster's draft initializer (`eight_node_areal`,
`W1ArealProfile`, fixture `results/research/b7_w1_initializer_draft/native.npz`), which is on branch
`aster/w1-profile` at cfae43c and not on main. It was run from a checkout of cfae43c:

    OPENBLAS_NUM_THREADS=1 python -m research.b7_w1_outer_cell_trace results/research/b7_w1_outer_cell_trace

numpy 2.5.3, scipy 1.17.0, mpmath 1.4.1. Eight runs: W1 and PC2 (rescaled pocket x=0.5), 1200 and 2400 cells
(dr 0.05, 0.025, r_max 60), dt 0.00075 and 0.000375, to t = 0.003, open edge. The production step is untouched;
the driver only wraps calls to record operands. Floor F = 1.5·(32/3)·eps/dr: 7.11e-14 (1200), 1.42e-13 (2400).

References, 60-digit mpmath on the last 8 cells:
- departure side-copy: a−1, b−1, X−1 from the initializer's own source departure, integrated with exact
  eight-node weights (matches stored a−1 to 2e-16 at setup), then advanced with the engine's binary64 rates.
  It never feeds back. Λ is the engine's stored Λ.
- arithmetic-only control: exact-weight arithmetic on the stored binary64 a, b.

## Results against the fixed gates

| Gate | Result |
|---|---|
| G0 | PASS. All eight peaks match Aster's to 3 digits: W1 1200 .522/.577, 2400 1.507/1.296; PC2 1200 .141/1.087, 2400 .986/1.977 (1e-12). |
| T1 | PASS. Worst stored-vs-rebuilt residual 0.50 of tolerance over every stage, outer three cells, a, b, X, Λ. The write changes only Aa[−1], bit for bit. |
| Setup | PASS. \|C64 − C_dep\| at the outer cell 0.14F–0.27F. |
| T2 | Mixed-stage piece is below F at both dt on both grids and both fixtures (largest 6.5e-18 against F 1.4e-13), so by the predeclared rule it does not matter at the edge. The K piece is exactly zero: K's L2 rate does not read Λ. The dt ratio is not applicable at this size. |
| T3 | FAIL. Size holds: worst per-stage \|ΔC64 − ΔC_dep\| is 0.30F (limit 2F). The 2400/1200 RMS ratio, required in [1.3, 3.0], is 1.80 (W1, dt .00075) and 1.40 (PC2, .00075) but 1.02 (W1, .000375) and 1.28 (PC2, .000375). |
| T4 | PASS. \|Aa_written,64 − Aa_written,dep\| ≤ 0.54F (limit 1.7F). Outer-cell momentum after each write ≤ 1.21e-12 (1200) and ≤ 2.74e-12 (2400). |

Report-only:
- Arithmetic control over departure difference (RMS): 0.24–0.73 at 1200, 0.63–1.00 at 2400. Lumen² predicted
  1/4–1/2. At 2400 the arithmetic part is about as large as the whole difference, which is consistent with the
  T3 ratio failing: the stage difference there is not dominated by the 1/dr storage term.
- Lag-1 sign correlation of the outer-cell stage difference: −0.49 to 0.35, no consistent pattern.
- The setup check is close to structural, because the initializer forms Λ with the same operator as Λ_metric
  (C64 at setup is 0 or 3e-23). The informative setup number is C64 − C_dep above.

Caveat on T3, not a gate change: each RMS uses 8 (dt .00075) or 16 (dt .000375) stage values, so the ratio
carries sampling noise of order tens of percent. The gate stays failed as written.

Full per-stage values are in `report.json` (`runs[].stages`, `runs[].steps`, `gates`).
