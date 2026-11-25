#!/usr/bin/env bash
# Runs the 32-run series with 4 workers (longest first), then scores it. Needs Jeremiah's go before use.
set -euo pipefail
OUT=${1:?output directory}
mkdir -p "$OUT"
export PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -m research.w1_fring_smoke.pairing_dt_series runs "$OUT" > "$OUT/runs.txt"
# A failed run is recorded in its JSON; keep the others going.
xargs -P 4 -L 1 sh -c 'python -m research.w1_fring_smoke.pairing_dt_series "$@" > "$5/log_$2_$3_$4${6:+_floor}.txt" 2>&1 || echo "failed: $*" >> "$5/failed.txt"' _ < "$OUT/runs.txt"
python -m research.w1_fring_smoke.pairing_dt_series score "$OUT"
