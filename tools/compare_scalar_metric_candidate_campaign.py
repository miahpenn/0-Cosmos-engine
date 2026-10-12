#!/usr/bin/env python3
"""Compare strict baseline and scalar-metric candidate t=12 campaign artifacts.

The candidate changes only the scalar Pi_t lapse-gradient contraction from
X^-2 alpha' f' to gamma^rr alpha' f', gamma^rr=X^2/a. This script compares
matched-resolution histories and saved checkpoint arrays. It does not choose
a physics interpretation or make changes to any input artifact.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import sys
import numpy as np

BASELINE_ROOT = Path(sys.argv[1] if len(sys.argv) > 1 else "baseline")
CANDIDATE_ROOT = Path(sys.argv[2] if len(sys.argv) > 2 else "runs/scalar-metric-candidate")
OUT = Path(sys.argv[3] if len(sys.argv) > 3 else "results/scalar-metric-candidate-comparison")
RESOLUTIONS = (40, 60, 80)
CHECKPOINT_TIMES = (2, 4, 6, 8, 10, 12)
DIAGNOSTICS = (
    "H_eff", "tau_rate", "hamiltonian_max", "hamiltonian_inner_absmax",
    "hamiltonian_full_absmax", "momentum_max", "momentum_inner_absmax",
    "momentum_full_absmax", "connection_max", "connection_inner_absmax",
    "connection_full_absmax", "determinant_constraint_max", "M_MS",
    "trapping_min", "lapse_min", "lapse_max", "cmc_residual_outer_max",
)

def finite_tree(x, label):
    if isinstance(x, dict):
        for k, v in x.items():
            finite_tree(v, f"{label}.{k}")
    elif isinstance(x, (list, tuple)):
        for i, v in enumerate(x):
            finite_tree(v, f"{label}[{i}]")
    elif isinstance(x, (int, float, np.number, np.ndarray)):
        if not bool(np.all(np.isfinite(np.asarray(x)))):
            raise ValueError(f"non-finite value at {label}")

def max_abs_diff(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch {a.shape} vs {b.shape}")
    delta = np.abs(a-b)
    return {
        "max_abs": float(np.max(delta)) if delta.size else 0.0,
        "rms": float(np.sqrt(np.mean(delta**2))) if delta.size else 0.0,
        "baseline_max_abs": float(np.max(np.abs(a))) if a.size else 0.0,
        "candidate_max_abs": float(np.max(np.abs(b))) if b.size else 0.0,
    }

def load_ledger(root, n):
    p = root / f"N{n}" / "ledgers.json"
    if not p.is_file():
        raise FileNotFoundError(p)
    obj = json.loads(p.read_text())
    finite_tree(obj, str(p))
    return obj

def compare_resolution(n):
    base = load_ledger(BASELINE_ROOT, n)
    cand = load_ledger(CANDIDATE_ROOT, n)
    br = base["result"]
    cr = cand["result"]
    if br["status"] != "completed" or cr["status"] != "completed":
        raise RuntimeError(f"N{n} non-completed result: baseline={br['status']} candidate={cr['status']}")
    if abs(float(br["t_final"])-12.0)>1e-9 or abs(float(cr["t_final"])-12.0)>1e-9:
        raise RuntimeError(f"N{n} did not reach t=12: {br['t_final']} / {cr['t_final']}")

    bh, ch = base["history"], cand["history"]
    btimes = np.asarray([float(x["t"]) for x in bh])
    ctimes = np.asarray([float(x["t"]) for x in ch])
    if btimes.shape != ctimes.shape:
        raise RuntimeError(f"N{n} history length differs: {len(btimes)} vs {len(ctimes)}")
    time_diff = max_abs_diff(btimes, ctimes)
    if time_diff["max_abs"] > 1e-12:
        raise RuntimeError(f"N{n} histories are not time-matched: {time_diff['max_abs']}")

    metric_deltas = {}
    for key in DIAGNOSTICS:
        pairs = [(row[key], row2[key]) for row, row2 in zip(bh, ch) if key in row and key in row2
                 and isinstance(row[key], (int, float)) and isinstance(row2[key], (int, float))
                 and math.isfinite(float(row[key])) and math.isfinite(float(row2[key]))]
        if pairs:
            metric_deltas[key] = max_abs_diff(
                [p[0] for p in pairs], [p[1] for p in pairs]
            )

    checkpoint_deltas = {}
    for t in CHECKPOINT_TIMES:
        name = f"checkpoint_t{t:012.6f}.npz"
        bp = BASELINE_ROOT / f"N{n}" / name
        cp = CANDIDATE_ROOT / f"N{n}" / name
        if not bp.is_file() or not cp.is_file():
            checkpoint_deltas[str(t)] = {"available": False, "baseline_exists": bp.is_file(), "candidate_exists": cp.is_file()}
            continue
        with np.load(bp, allow_pickle=False) as bz, np.load(cp, allow_pickle=False) as cz:
            bkeys, ckeys = set(bz.files), set(cz.files)
            if bkeys != ckeys:
                raise RuntimeError(f"N{n} t={t}: checkpoint keys differ: {sorted(bkeys^ckeys)}")
            arrays = {}
            for key in sorted(bkeys):
                x, y = np.asarray(bz[key]), np.asarray(cz[key])
                if x.dtype.kind in "fiu" and y.dtype.kind in "fiu":
                    arrays[key] = max_abs_diff(x, y)
            checkpoint_deltas[str(t)] = {"available": True, "fields": arrays}
    return {
        "resolution": n,
        "baseline_result": {k: br[k] for k in ("status","t_final","tau_final","steps","misner_sharp_current_max","misner_sharp_current_rms")},
        "candidate_result": {k: cr[k] for k in ("status","t_final","tau_final","steps","misner_sharp_current_max","misner_sharp_current_rms")},
        "history_rows": len(bh),
        "time_alignment_max_abs": time_diff["max_abs"],
        "diagnostic_time_series_deltas": metric_deltas,
        "checkpoint_field_deltas": checkpoint_deltas,
    }

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    result = {
        "schema": "scalar_inverse_metric_candidate_campaign_compare_v1",
        "status": "completed",
        "kind": "matched_t12_history_and_checkpoint_comparison",
        "baseline_run_id": 38208253687,
        "candidate_change": "replace X^-2 alpha' f' by (X^2/a) alpha' f' in scalar momentum RHS",
        "resolutions": list(RESOLUTIONS),
        "cases": [compare_resolution(n) for n in RESOLUTIONS],
        "interpretation_guardrail": "Comparison quantifies consequences of the single equation correction. A small difference does not invalidate the metric identity; it means this term is not the dominant cause of the existing constraint drift."
    }
    (OUT / "comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print("SCALAR_METRIC_CANDIDATE_COMPARISON " + json.dumps({
        "status": result["status"],
        "cases": [{
            "N": c["resolution"], "baseline_tau": c["baseline_result"]["tau_final"],
            "candidate_tau": c["candidate_result"]["tau_final"],
            "baseline_Ms_max": c["baseline_result"]["misner_sharp_current_max"],
            "candidate_Ms_max": c["candidate_result"]["misner_sharp_current_max"],
            "candidate_minus_baseline_diagnostic_deltas": c["diagnostic_time_series_deltas"],
            "checkpoint_t12": c["checkpoint_field_deltas"].get("12")
        } for c in result["cases"]
        ]}, sort_keys=True))

if __name__ == "__main__":
    main()
