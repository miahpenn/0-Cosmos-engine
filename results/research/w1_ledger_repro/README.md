# Reproduction of Aster's W1 ledgers (x_endpoint, conformal_ledger)

**Result:** both runs reproduced Aster's archived ledgers byte for byte. That includes all nine targets, on both Python stacks.

**Go:** Jeremiah, "Go", Nov 24 06:15 UTC, for reproducing Aster's numbers with the #587 pack. The installs were approved separately ("Yes, install them", 06:18 UTC). He agreed to 3.12.13, with an update to 3.12.14 only if the run disagreed (06:22 UTC). The pinned run matched, so no 3.12.14 rerun was made.

**What ran:** the pack's scripts unchanged, as its replay commands give them. For each script, `--divisor 1`, `2` and `4` (20, 40 and 80 steps at dr 0.025 to t 0.015), then the combine step. Engine 05c2a09, vendor d6052d6, NPZ 339652a2; `origin_trace.run` checks all three. Nothing else ran, and the A6 runs remain held.

## Environments

| | pinned (Aster's versions) | container |
|---|---|---|
| Python | **3.12.13** (Aster used 3.12.14, which this uv can't download) | 3.13.18 |
| NumPy / SciPy | 2.3.5 / 1.17.0 | 2.5.3 / 1.18.1 |
| OpenBLAS (wheel) | 0.3.30 | 0.3.34 |
| OpenBLAS core picked at runtime | SkylakeX | SkylakeX |
| `OPENBLAS_CORETYPE` | unset | unset |

- **Shared:** both runs had `OPENBLAS_NUM_THREADS=1` and longdouble with 63 mantissa bits.
- **CPU:** Intel Xeon @ 2.10GHz with AVX-512, on x86_64.
- **Aster's record:** Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, OpenBLAS 0.3.30 (SkylakeX, one thread).
- **Kernel:** the "Haswell" in the OpenBLAS build string is the build's minimum target. The core in use is chosen at runtime, and both runs chose SkylakeX.
- **Full details:** `pinned/env.json` and `container/env.json`.

## Targets (nine-point stencil, cell 0)

| target | archived | pinned: bitwise | container: rel. diff |
|---|---|---|---|
| H, nominal / half / quarter | 8.964956e-11 / 1.492382e-11 / -2.201987e-11 | PASS, PASS, PASS | 0, 0, 0 |
| X writes, nominal-half / half-quarter | 6.269787e-12 / 1.691661e-11 | PASS, PASS | 0, 0 |
| tracefree rate, nominal-half / half-quarter | 6.506347e-11 / 1.859661e-11 | PASS, PASS | 0, 0 |
| summed writes, nominal-half / half-quarter | 9.800459e-12 / 1.851571e-11 | PASS, PASS | 0, 0 |

## Whole ledgers

Every numeric leaf matches in both runs: x_endpoint 6,021 of 6,021 leaves and conformal_ledger 7,246 of 7,246. The two combine outputs are also byte-identical to the archived files. Their SHA-256 values are f99316a9…0a (x_endpoint) and 446fa04d…14 (conformal_ledger), the same as `results/research/w1_fring_x_endpoint/ledger.json` and `results/research/w1_fring_conformal_ledger/ledger.json`. The outputs are therefore not stored again; `*/output_sha256.txt` holds their hashes.

**Built-in checks:** all 12 replays exited 0. Each script compares every accepted state against Aster's archived readbacks and raises an error on any mismatch, and none was raised.

- **Wall time:** about 2 to 7 minutes per replay, run in parallel. See `*/replay_status.txt`.
- **Comparison:** `research/w1_ledger_repro/compare.py`, with results in `*/compare.json`.
