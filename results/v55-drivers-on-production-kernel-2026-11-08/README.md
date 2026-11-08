# V5.5 drivers re-pointed to the production kernel: old vs new ledger (2026-11-08)

`engine/bidirectional_production.run` (the interface ledger) and `engine/coupled_engine` evolved their
own copy of the V5.5 system (`reference_pirk_unified` + `true_cmc_reference`). They are now thin drivers
over `V55ProductionKernel`. Why the old stack is not live: `docs/PROVENANCE_OLDER_V55_PATHS.md`.

`tools/v55_driver_ledger_snapshot.py` (SHA-256 `902edd32…5f63`) runs `bidirectional_production.run(40, 6.0, 0.06)`
from `main` at `41ee0cb4` (old drivers; `old_reference_stack.json`, SHA-256 `6e4d627f…6e71`) and from
this branch (`production_kernel.json`, SHA-256 `e2ccf534…b006`); vendor `d6052d60`, Python 3.12.15,
numpy 2.5.4, scipy 1.18.1. Default initial data (amplitude 0.01, width 7), r_max = 40, worldtube at r = 10.

| | old, t = 0 | new, t = 0 | old, t = 6 | new, t = 6 |
|---|---|---|---|---|
| phi_outer | 33.898 | 0.17906 | 34.584 | 0.19974 |
| Pi_outer | 0.17906 | 0.0034191 | 0.077283 | 0.0034745 |
| H_eff | 7.31e-2 | 5.67e-3 | 3.16e-2 | 5.54e-3 |
| M_MS (r = 10) | 2.33 | 4.71e-2 | 1.03 | 4.76e-2 |
| chi_sigma | 0.510 | 0.990 | 0.837 | 0.990 |
| rho_outer | 6.39e-4 | 3.84e-6 | 1.19e-4 | 3.67e-6 |
| work_pR | −0.503 | +1.99e-4 | −9.38e-2 | +1.92e-4 |
| flux_T | 4.8e-15 | −4.7e-19 | 2.33e-4 | −3.87e-4 |

Misner–Sharp balance residual max |dM/dt − (flux + work)| over 2 ≤ t ≤ 4: old 1.8e-5, new 1.1e-5.

Reading: the old initializer put the archive scale factor a = 33.8983 into the COSMOS scalar φ and the
archive φ into its momentum Π, so the cosmic potential energy, density and H₀ started far from the
archive operating point and the ledger quantities differ by one to two orders of magnitude (work changes
sign). On the production kernel the COSMOS scalar starts at the archive values (φ = 0.179055,
Π = 0.00341913). The two drivers now own no equations; their steps are bit-identical to
`V55ProductionKernel.step` (tested).
