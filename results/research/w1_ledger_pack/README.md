# W1 ledger reproduction pack (Aster)

**What this adds:** Aster's eight research commits, 1377495f to b9c04d0b, built on PR #584 head 05c2a09. They are merged here unchanged, keeping the original Git objects, authorship and hashes.

**Source:** Jeremiah shared it on Nov 24 at 06:12 UTC: "See drive folder: W1 Ledger Reproduction Pack". The Drive folder is 1YWMzeAnLQ3wucrCH0vKYvetZfTNwjTyc, and it holds the git bundle `W1-ledger-research-b9c04d0b.bundle`.

**Design:** see f_ring A7 (Drive 1dr_BKBvWuk5nh_wVMMFOPSiPtAAAfpstJg9NHypnR9U).

Nothing has been rerun here. Reproduction waits for Jeremiah's go.

## Checks made on intake (no code run)

- **Bundle:** SHA-256 e2ff2f17…7ac2, matching the pack's SHA256SUMS. `git bundle verify` passes. It needs 05c2a09 as a prerequisite and has a single head, b9c04d0b.
- **Snapshot:** all 33 files under `snapshot/` in the pack match the bundle's b9c04d0b content by SHA-256.
- **Scope:** the pack touches only `research/w1_fring_smoke/*.py` (8 scripts) and `results/research/w1_fring_*` (7 result folders). It makes no change to engine, vendor, tests or workflows.
- **Pins:**
  - Engine source 05c2a09. Its `engine/` and `vendor/` are identical to main's.
  - Vendor d6052d6.
  - Input NPZ SHA-256 339652a2…a295, the corrected slice from #581.
  - `origin_trace.run` checks all three before any step and refuses to run otherwise.
- **Side effects:**
  - Subprocess is limited to read-only `git rev-parse` and `git diff --quiet`, in `origin_trace.run`.
  - There's no network access.
  - Files are written only by `x_endpoint` and `conformal_ledger`, as `np.savez` into a directory the caller names (`mktemp -d` in their replay commands). The other scripts print JSON to stdout.
  - Every `np.load` sets `allow_pickle=False`.

## Scripts

All are in `research/w1_fring_smoke/`.

| script | purpose | report |
|---|---|---|
| `origin_trace.py` | Signed origin trace through the first PIRK2 step. This is the base runner the others import. | `w1_fring_origin_trace` |
| `areal_trace.py` | Production Ricci compared with direct areal curvature (seven- and nine-point) on accepted slices. | `w1_fring_areal_trace` |
| `quarter_trace.py`, `quarter_stage_trace.py` | Nominal, half and quarter dt; predictor and corrector stages. | `w1_fring_quarter_trace` |
| `metric_balance.py` | Metric balance across predictor and corrector. | `w1_fring_metric_balance` |
| `x_update_ledger.py` | X write ledger at the first shared time. | `w1_fring_x_update` |
| `x_endpoint.py` | X writes and rates at t 0.015. | `w1_fring_x_endpoint` |
| `conformal_ledger.py` | Determinant and tracefree split of the a, b difference. | `w1_fring_conformal_ledger` |

## Before any reproduction

- **Environment:** Aster's runs used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and OpenBLAS 0.3.30 (SkylakeX, one thread). This container has Python 3.13.18, NumPy 2.5.3 and SciPy 1.18.1. Its longdouble has 63 mantissa bits, which meets the 64-bit-significand check the ledgers need. Bitwise readback comparisons may differ for that reason alone.
- **Source commit:** `origin_trace.run` diffs against 05c2a09. That commit stays reachable through `refs/pull/584/head`.
