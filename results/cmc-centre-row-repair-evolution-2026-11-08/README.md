# CMC lapse centre-row repair: evolution comparison, old vs new row (2026-11-08)

Evidence for the repair PR that replaces the CMC lapse solver's centre row
(`engine/cmc_gauge.py`; diff in `centre_row_fix.diff`, SHA-256
`6ddf52815d24280d01cef845815d8fd3965bcd70520b4f783717b27c37c82bac`). The attribution is in
`results/issue4-gate-c-centre-lapse-defect-2026-11-08/`.

## What was run

`tools/cmc_centre_row_evolution_comparison.py` (SHA-256
`20ab4e39fbfb0df56848cdc1dbde5793f39c74b8ba3140d55c3f646ea6df6850`), from the repair worktree
(`main` at `afc4d96f` plus the fix; vendor `d6052d60`), Python 3.12.15, numpy 2.5.4, scipy 1.18.1.
For N = 60, 120, 240: `V55ProductionKernel().initialize(resolution=N, r_max=40)` (default initial data,
not the #420 lineage), then real `kernel.step` calls with dt = 0.06·dr (0.04, 0.02, 0.01) to t = 0.60,
once with the production solver from `main` (old row, loaded from `main`'s file) and once with the
repaired solver. Total constraints H − 16πρ and M − 8πj and the connection constraint C = Λ − Λ_metric
are recorded at cells 0–2 and in 5 ≤ r ≤ 10. Output: `compare_rows.json`, SHA-256
`188db0f43497074372c92ef1bddc208bdfacb990e4d4f9362ebe6131134513e0`.

## Results at t = 0.60

| Quantity, cells 0–2 | N = 60 old / new | N = 120 old / new | N = 240 old / new |
|---|---|---|---|
| connection constraint C | 7.1e-4 / 7.5e-4 | 8.6e-4 / 4.5e-4 | 9.1e-4 / **7.7e-5** |
| Hamiltonian H_total | 2.6e-4 / 5.0e-4 | 1.5e-2 / 9.4e-3 | 3.0e-2 / **5.9e-3** |
| momentum M_total | 1.0e-3 / 9.5e-4 | 6.9e-4 / 6.7e-4 | 7.4e-4 / 9.5e-4 |

- With the old row, the centre connection constraint does not converge (7.1e-4 → 9.1e-4) and the centre
  Hamiltonian constraint grows with resolution (2.6e-4 → 3.0e-2). With the new row the connection
  constraint converges (7.5e-4 → 7.7e-5) and the Hamiltonian is 5× smaller at N = 240.
- The centre momentum constraint is ~7e-4 to 1e-3 with either row and does not converge in either: it
  comes from a different, pre-existing centre source (consistent with the initializer/curvature mismatch
  recorded in PR #414), not from the lapse row.
- Away from the centre nothing changes: the 5–10 annulus values of M, H and C are identical to three
  digits for old and new. Final-state differences at r ≥ 5 are ≤ 5.3e-7 (N = 60), 6.0e-8, 7.9e-9 (N = 240)
  for every evolved field except the lapse, which shifts globally by 3.7e-5, 9.6e-6, 2.4e-6 (the old row's
  second-order error in the solve's normalization, now removed).
- New-row lapse at t = 0.60: α₁ − α₀ = 2.8e-5, 6.9e-6, 1.7e-6 (∝ h², a smooth even lapse); old row: 0.

## Reading

The repair does what the attribution predicted: it removes the lapse-row contribution to the centre
defects, the centre connection constraint now converges, and nothing changes away from the centre beyond
the expected shift of the lapse normalization. It does not fix everything at the centre: the
non-converging centre momentum constraint has another, pre-existing source and needs its own
investigation.
