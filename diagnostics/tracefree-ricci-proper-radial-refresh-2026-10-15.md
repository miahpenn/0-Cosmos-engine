# Trace-free Ricci / proper-radial reconciliation — 15 October 2026

**Disposition:** diagnostic-only; scientific admission remains **NOT ADMITTED**. No production-physics or vendor files were changed.

## Exact-head source snapshot for PR #274 face-metric validation

The previous workflow correctly rejected this PR because the production matter-flux file differs from its historical pin. Instead of relaxing that guard or treating the provenance failure as a numerical failure, the diagnostic is now pinned to the exact PR production snapshot below and rerun against it.

- Source snapshot commit (pre-audit-manifest update): `e6b501fb19c4f2e21ad9c655662840a3d2582d54`.
- Production-kernel blob: `9685e76bb5a6f26d6bf383e240de24f37523f265`.
- Vendor submodule gitlink: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

| Engine file | Required blob SHA at the new snapshot |
| --- | --- |
| `engine/production_kernel.py` | `9685e76bb5a6f26d6bf383e240de24f37523f265` |
| `engine/matter_system.py` | `e52a117690253ddad7b24d116145a65474a351d2` |
| `engine/matter_rhs.py` | `f41180efbcff2c72cf36477a765318ab539fac80` |
| `engine/scalar_system.py` | `ffc119b1a6f2d495d972323901f76a85783b0a1b` |
| `engine/v55_matter.py` | `0ee3c51fac7dc46c63feba2747a750542bcee3d0` |
| `engine/v55_pirk_adapter.py` | `eeee0492d80889d70e52b0946563eee5054c0f11` |
| `engine/cmc_gauge.py` | `7b49b5751abebd4faf0130208032cd2080b71f64` |
| `engine/v55_initial.py` | `7db2e1a10cf16fb40e4221b0df504fdf39aae611` |
| `engine/stress_energy.py` | `ac3e7518c62bb72e2218a61489fa486bb787809f` |
| `engine/valencia.py` | `03ed86398f424c03d9d6a3e646da9222a34fac5c` |
| `engine/worldtube.py` | `aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6` |

The workflow and script use this same source commit, kernel blob, and vendor gitlink. Both retain exact blob checks for the full engine manifest and reject any engine/vendor/dependency edits after the snapshot. This is a provenance refresh so the numerical diagnostic can actually execute on the face-metric implementation; it is not a waiver of a numerical check or a claim of scientific admission.

## Historical source pins and configuration (15 October run)

- Refreshed live-main / PR-base snapshot: `7860d11cd49ea3a3317d01ef1f4f4950c4db80b8`.
- Recorded production-source comparison baseline: `fcb32b6b4025d5674e0560b8fb24581933f7f23c`.
- Production-kernel blob: `9685e76bb5a6f26d6bf383e240de24f37523f265`.
- Vendor submodule gitlink: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

All six live-main engine blobs match their recorded hashes exactly:

| Engine file | Expected and observed blob SHA |
| --- | --- |
| `engine/production_kernel.py` | `9685e76bb5a6f26d6bf383e240de24f37523f265` |
| `engine/v55_initial.py` | `7db2e1a10cf16fb40e4221b0df504fdf39aae611` |
| `engine/v55_pirk_adapter.py` | `eeee0492d80889d70e52b0946563eee5054c0f11` |
| `engine/cmc_gauge.py` | `7b49b5751abebd4faf0130208032cd2080b71f64` |
| `engine/v55_matter.py` | `acc8f8ceae28d82511c9a5ad9d2ca670e09fb0f7` |
| `engine/stress_energy.py` | `ac3e7518c62bb72e2218a61489fa486bb787809f` |

The diagnostic uses resolutions N = 40, 60, 80, 120, 160; r_max = 80; CFL = 0.03; amplitude = 0.01; width = 7; D amplitude = 1e-10; radiation enabled; and one ordinary production step per resolution. This captures four actual geometry-stage L2 inputs: initial RHS, explicit predictor RHS, primary predictor RHS, and final implicit L2 RHS. The workflow checks out the exact PR-head SHA and fails closed on source, vendor, dependency, finiteness, and same-stage closure discrepancies.

## Historical measured numerical results (15 October source snapshot)

**Provenance note:** The measurements below belong to workflow run [38666458930](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38666458930) and its historical source pins; they are retained as prior evidence, not as results for the new face-metric snapshot above. Review the regenerated exact-head artifact for the new snapshot separately.


The first completed exact-head run was [diagnostic workflow 38666458930](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38666458930), on diagnostic implementation head `d93160349c64aea31b94fb3ef68add0151eff1fc`. It emitted report status `completed_diagnostic_only`, with an empty error list. The full repository regression suite passed **589 tests in 25.74 seconds**. The source/test/workflow code is unchanged by this report-only results addition; exact-head CI and the diagnostic artifact are being rerun for the final documentation-completed candidate.

The signed Ricci gap below is `vendor trace-free Ricci − proper-radial reconstruction`. Cell i=0 and i=1 values are reported separately. Annulus maxima are `max |gap|` on fixed **areal radius** 5 ≤ R_areal ≤ 8, not a coordinate-radius substitute. Values are from the initial RHS unless a column says final stage.

| N | i=0 signed gap | i=1 signed gap | Annulus max abs, initial (cells) | Annulus max abs, final L2 stage | Hessian max abs, initial / final |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 40 | +1.6562870e-3 | -2.0793835e-4 | 3.4972353e-5 (2) | 3.5054053e-5 | 1.4286555e-10 / 1.1042533e-10 |
| 60 | +1.6146273e-3 | -2.2420894e-4 | 1.8776969e-6 (2) | 1.9319766e-6 | 4.6860524e-11 / 4.7469005e-11 |
| 80 | +1.6023053e-3 | -2.2698875e-4 | 1.1088630e-6 (3) | 1.1389849e-6 | 1.7861021e-11 / 2.3474567e-11 |
| 120 | +1.5943981e-3 | -2.2742724e-4 | 3.4016485e-7 (5) | 3.5233915e-7 | 4.0025673e-12 / 1.2572294e-11 |
| 160 | +1.5918150e-3 | -2.2729395e-4 | 8.8532672e-8 (6) | 9.5874420e-8 | 1.3153626e-12 / 7.6923200e-12 |

At N=160, the centre-cell coordinate radius is 0.25, and the recorded areal radius is kept separately in the per-cell profiles. The cell-i=0 signed gap remains about +1.59e-3 across the sweep rather than tending to zero; by contrast the annulus maximum falls to 8.85e-8 on the initial RHS. The direct-D2 Hessian gap is small in absolute terms and declines with refinement on the initial RHS. These are distinct regional behaviors, not a single combined statistic.

Across all five resolutions and all four stages (20 stage records), every measured profile was finite, both reconstructed state digests were identical before/after each bounded step, and all stage checks passed. The maximum recorded same-stage A_a RHS reconstruction closure and K RHS reconstruction closure were both 0.0, compared with a roundoff bound of 1.4210855e-14. Each resolution advanced exactly one bounded step; this is not an extended campaign.

## Manufactured control

The regression suite includes a smooth constant-positive-curvature three-space with k = 1e-3. Since the radial and angular Ricci eigenvalues agree analytically, trace-free Ricci should vanish. The test checks both the vendor and proper-radial representations near the centre and requires their gap to be below 5e-6; it also verifies that the fixed R_areal annulus is populated. This test passed within the 589-test suite.

## Interpretation and limits

The diagnostic localizes a persistent centre Ricci-representation mismatch while the same representation gap becomes much smaller on the fixed physical-radius annulus under refinement. It does **not** establish which expression or operator is wrong, nor prove convergence of the full coupled solution. The independent geometric formula shares the pinned radial derivative family/grid; it is an independent representation, not an independent discretization. No fixed absolute Hessian cutoff is used to decide acceptance.

The workflow uploads strict-JSON `report.json`, signed `profiles.csv`, `profiles.npz`, source-pin text files, and `SHA256SUMS.txt` as artifact `tracefree-ricci-proper-radial-2026-10-15`. The exact ZIP digest and regenerated per-file digests are recorded with the final exact-head result in Issue #4. Keep this PR draft until its final exact-head workflow rerun and independent methodological review complete.

## Coordination

- Diagnostic implementation: [PR #220](https://github.com/miahpenn/0-Cosmos-engine/pull/220).
- Prior current-main initialized-slice observer: [PR #218](https://github.com/miahpenn/0-Cosmos-engine/pull/218).
- Earlier candidate refreshed into this current-main PR: [PR #219](https://github.com/miahpenn/0-Cosmos-engine/pull/219).
- Scientific admission remains **NOT ADMITTED**; the action-normalization/source provenance gate in Issue #210 is separate and remains open.
