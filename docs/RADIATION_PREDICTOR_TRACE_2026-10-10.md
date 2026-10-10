# 0* Cosmos Engine — Radiation Predictor Failure Trace

**Date:** 10 October 2026  
**Branch:** `physics/spatial-beta-covariant-source-closure`  
**Diagnostic workflow:** [Run #7](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938182)  
**CI on trigger commit:** [Run #576](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938171) — 127 tests passed in 59.58 s  
**Artifact:** [Download `0star-radiation-stage-trace.zip`](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938182/artifacts/11676486591)

## Provenance

| Role | Exact commit |
|---|---|
| Physics/numerics source baseline | `9eed1cf95c52f001c2b87afbbe6da37faa56d50e` |
| Instrumentation | `527b44931b059f4d118a86c898b91d79ca4d96f2` |
| Trace workflow | `ed0bc6f61103f2439556b27ad588116a149f4add` |
| Content-neutral trigger / run SHA | `c8a9967126ba07cd12e4da9c86b7066a703d3419` |
| Vendor submodule | `d6052d605673ce9d82cdc99b0df79855fa2ea215` |

Configuration: `N=80`, `r_max=80`, `dr=1`, `dt=CFL=0.03`, target `t=50`, `D=1e-10`, radiation ON, scalar profile amplitude `0.01`, width `7`. Runtime Python 3.12.15 / NumPy 2.5.3. A GitHub comparison from the source baseline to the trigger commit showed only three changed paths: trace code, trace workflow, and trace tests; production physics/numerics files were unchanged between those commits.

Source SHA-256 inventory (the full inventory is embedded in `trace.json`):

- `engine/production_kernel.py`: `2adfa35803078f72cf7f4139a0690c6ae5b8c4da9b9ca612de0bc2a255037f55`
- `engine/cmc_gauge.py`: `bfb760d64d2d5afda911d2d628bd939b84ba6ecc9da0b7f58b2a8bae7b5979c6`
- `engine/scalar_system.py`: `d9d4ab473bdf11669ca42be71a080866386bbbab88827d00d56e8a793b7a3a39`
- `engine/matter_system.py`: `0b47aabb986f6d7e14073d81fd345fcb40aa3b8354ff02a6cb4e16474a5afcd5`
- `engine/matter_rhs.py`: `1a9f50e325877c0a74756d8398154baa566ae198a3884f62b5b37bd822f6e4fa`
- `engine/valencia.py`: `dd01288447e99744f80ae6c8b7fd69a6ef318b799ee5cde127fd6295cf70fb56`
- `engine/v55_matter.py`: `e6deb37073a19ad841a7653cb37b82ef2d237ef310a62de94e86c9511dd12459`
- `engine/v55_pirk_adapter.py`: `4131b707813038f6a4f06efea8e9b69cc39f227739f2c43ba18f7543e84efc47`
- `engine/v55_initial.py`: `3687138660b93f4d62057e8775ac2ff1b4ce0176800061245cf07633fb823842`

## Captured event

The diagnostic workflow completed successfully because the failure was captured and the artifact uploaded; it did **not** reach its target time.

- Last accepted state: `t=46.98000000000102`, center proper time `tau=22.016115393070738`; 1,566 accepted steps.
- Failed attempted predictor: `t=47.01000000000102`, step `dt=0.03`.
- First failure stage: `explicit_predictor_first_CMC_lapse`.
- Outermost cell: `i=79`, `r=79.5`.
- Explicit admissibility error: `E>=|S|` violated with `E=3.97896173626637354e-8`, `|S|=3.97918851370748485e-8`, ratio `|S|/E=1.00005699412463422`.
- Predictor metric at this cell: `alpha=1`, `beta=0`, `gamma_rr_inv=1.4853571835116133`.
- Neighbor ratios on the attempted predictor: cell 77 `0.9301483`, cell 78 `0.9258852`, cell 79 `1.0000570`.
- Last accepted state at cell 79 remains admissible but nearly saturated: ratio `0.9983867169143378`; cone margin `C=+4.6971156703181716e-7`.

## Ordered cone-margin budget at cell 79

Here `C = U_E - sqrt(gamma_rr_inv) * |U_r|`. The budget evaluates the accepted state on the accepted metric, then a flux-only conservative update and a geometric-source update on that same accepted metric, and finally evaluates the resulting predictor state on its predictor metric.

| Ordered stage | Cone margin C | Increment |
|---|---:|---:|
| Accepted state / accepted metric | `+4.6971157e-7` | — |
| After flux update / accepted metric | `+1.4812240e-6` | `+1.0115125e-6` |
| After geometric sources / accepted metric | `+2.2171301e-6` | `+7.3590611e-7` |
| Predictor matter state / predictor metric | `−1.6900205e-8` | metric contribution `−2.2340303e-6` |
| Net accepted-to-predictor change | | `−4.8661177e-7` |

At the same cell, accepted-to-predictor metric values change as follows: `a: 0.60895534 → 0.60320350` (−0.9445%); `b: 1.28146656 → 1.28756177` (+0.4756%); `X: 0.94389565 → 0.94655832` (+0.2821%); `gamma_rr_inv=X²/a: 1.46306132 → 1.48535718` (+1.5239%). Outer lapse stays `alpha=1` and shift `beta=0`. The predictor conservative state has ratio `0.9925229671` when tested on the accepted metric but `1.0000569941` on the predictor metric.

The ledger is internally closed: across the retained 66 late-window predictor budgets, there were no missing source captures; maximum closure error and predictor `U_E` / `U_r` reconstruction errors were all zero. The run recorded 394 stage snapshots, 1,566 progress samples, and **zero instrumentation errors**. The full trace intentionally retains a two-unit late-time window, not a full trajectory checkpoint.

## Bounded conclusion

For this failing step, flux transport and geometric sources *increase* the cone margin when evaluated on the accepted metric. The final metric evaluation then reduces it by more than the margin remaining and is the decisive crossing in this ordered budget. The accepted-state margin had already been collapsing as the run approached the outer lightlike limit: accepted outer-cell ratios rose from about `0.99254` at `t=46.86` to `0.99839` at `t=46.98`.

This locates the crossing in the explicit predictor metric, but it does **not** prove the metric evolution is an incorrect physical law. The upstream source of the increasingly steep radial-metric change could involve the explicit geometry RHS, boundary/characteristic treatment, stage admissibility, or their coupling. No clipping, damping, floors, or boundary-law edits are justified by this result.

## Integrity

- Uploaded artifact ZIP SHA-256 (from GitHub metadata and independently verified): `bef03745ad4ea1ab62d8a42de9b957c57210571d933cf2728c96a6196b970c43`.
- `trace.json` SHA-256: `642a82b8b97ddc849b08f70b4379672b08e81f842d3ae8bdeaaa1d5374afd8a3`.
- `summary.json` SHA-256: `f36105276c2ef4fc7ee91ed8af5593792cf05ee27324149a3e17e4cca87fee22`.
- Both internal checksum entries matched the downloaded archive contents.

## Next diagnostic gate

Keep physics unchanged. On a new, separately pinned diagnostic, record accepted and first-predictor explicit RHS terms for `a`, `b`, and `X` at cells 75–79; reconstruct `Δgamma_rr_inv` from those terms; and audit the matching pinned light-characteristic boundary variables. The immediate question is whether the rise in `gamma_rr_inv` is supplied by the interior explicit geometry RHS, boundary reconstruction, or both. Only after that deterministic audit should a controlled resolution/CFL comparison be authorized. No stabilization or physical-law edit is authorized by this trace.