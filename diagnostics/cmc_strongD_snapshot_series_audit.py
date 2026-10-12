"""Time-series audit of CMC lapse re-solve corrections on archived strong-D snapshots.

The 49 saved, full-state snapshots from the boundary-matched R80 trajectory are
re-solved one slice at a time with the recorded outer alpha held fixed. This is
not a re-evolution: the integral of the central lapse correction is only a
scale comparison against the archived R80-R160 proper-time discrepancy.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import hashlib

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import engine.cmc_gauge as cmc
from diagnostics.cmc_strongD_snapshot_jacobian import (
    load_state,
    profile_summary,
    residual_for_target,
    solve_map,
)

OUT = ROOT / "runs" / "cmc_strongD_snapshot_series"
REFERENCE_TAU_MISMATCH = 0.0106210937742
BASELINE_TAU_MISMATCH = 0.0117426496868
TOL = 1.0e-8


def read_rows(path):
    rows = []
    with Path(path).open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("state_snapshot") is None:
                continue
            if row.get("diagnostic_only") is not True:
                raise ValueError(f"row {line_number} is not marked diagnostic_only")
            rows.append(row)
    rows.sort(key=lambda item: float(item["t"]))
    # Keep only unique times in case an artifact carries duplicate labels.
    unique = []
    seen = set()
    for row in rows:
        t = float(row["t"])
        key = round(t, 10)
        if key not in seen:
            unique.append(row)
            seen.add(key)
    if len(unique) < 3:
        raise ValueError(f"expected a saved state series, found {len(unique)} unique snapshots")
    return unique


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: python diagnostics/cmc_strongD_snapshot_series_audit.py <profiles.jsonl>")
    profile_path = Path(sys.argv[1])
    rows = read_rows(profile_path)
    records = []
    finite = True
    same_target_closed = True

    for idx, row in enumerate(rows):
        state = load_state(row)
        grid = state.grid
        geom = state.geometry.copy()
        alpha_in = np.asarray(geom.alpha, dtype=float).copy()
        outer_alpha = float(alpha_in[-1])
        alpha_out, kdot_used = solve_map(state, geom, alpha_in, outer_alpha)
        geom_out = geom.copy()
        geom_out.alpha = alpha_out.copy()
        kdot_out = float(cmc.target_kdot(grid, geom_out, state.scalars, state.matter))
        kdot_in = float(cmc.target_kdot(grid, geom, state.scalars, state.matter))

        r_in, rel_in, _ = residual_for_target(
            state, geom, alpha_in, kdot_in, outer_alpha
        )
        r_same, rel_same, _ = residual_for_target(
            state, geom, alpha_out, kdot_in, outer_alpha
        )
        r_recomputed, rel_recomputed, _ = residual_for_target(
            state, geom_out, alpha_out, kdot_out, outer_alpha
        )
        delta = alpha_out - alpha_in

        rec = {
            "t": float(state.t),
            "tau_snapshot": float(state.tau),
            "outer_alpha_fixed": outer_alpha,
            "alpha_center_in": float(alpha_in[0]),
            "alpha_center_out": float(alpha_out[0]),
            "delta_alpha_center": float(delta[0]),
            "delta_alpha_rms": float(np.sqrt(np.mean(delta * delta))),
            "delta_alpha_max_abs": float(np.max(np.abs(delta))),
            "delta_alpha_max_radius": float((int(np.argmax(np.abs(delta))) + 0.5) * grid.dr),
            "input_target_kdot": kdot_in,
            "output_target_kdot": kdot_out,
            "delta_target_kdot": float(kdot_out - kdot_in),
            "input_cmc_residual_relative_max": rel_in,
            "mapped_same_target_residual_relative_max": rel_same,
            "mapped_recomputed_target_residual_relative_max": rel_recomputed,
            "input_cmc_residual_max_abs": float(np.max(np.abs(r_in))),
            "mapped_recomputed_residual_max_abs": float(np.max(np.abs(r_recomputed))),
            "all_finite": bool(
                np.all(np.isfinite(alpha_in)) and np.all(np.isfinite(alpha_out))
                and np.all(np.isfinite(delta)) and np.all(np.isfinite(r_in))
                and np.all(np.isfinite(r_same)) and np.all(np.isfinite(r_recomputed))
            ),
        }
        records.append(rec)
        finite = finite and rec["all_finite"]
        same_target_closed = same_target_closed and rel_same < 1.0e-10

    times = np.asarray([r["t"] for r in records], dtype=float)
    central_correction = np.asarray([r["delta_alpha_center"] for r in records], dtype=float)
    if np.any(np.diff(times) <= 0.0):
        raise ValueError("snapshot time series must be strictly increasing")
    central_correction_integral = float(np.trapz(central_correction, times))
    central_abs_integral = float(np.trapz(np.abs(central_correction), times))
    tau = np.asarray([r["tau_snapshot"] for r in records], dtype=float)
    tau_rates = np.diff(tau) / np.diff(times)

    payload = {
        "test": "strong-D saved-snapshot CMC lapse re-solve series",
        "provenance": {
            "source_run": "https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255",
            "source_artifact_id": 11710858125,
            "source_commit": "35072830f9e4d7029bb10243a0520b9a273a17dd",
            "profile_file": profile_path.name,
            "profile_file_sha256": hashlib.sha256(profile_path.read_bytes()).hexdigest(),
            "snapshots": len(records),
            "start_t": float(times[0]),
            "end_t": float(times[-1]),
            "snapshot_step_min": float(np.min(np.diff(times))),
            "snapshot_step_max": float(np.max(np.diff(times))),
            "snapshot_taus_are_archived_ledgers": True,
        },
        "production_equations_changed": False,
        "per_snapshot": records,
        "summary": {
            "central_delta_alpha_at_final_snapshot": float(central_correction[-1]),
            "central_delta_alpha_max_abs_over_series": float(np.max(np.abs(central_correction))),
            "central_delta_alpha_rms_over_series": float(np.sqrt(np.mean(central_correction**2))),
            "delta_alpha_rms_max_over_series": float(max(r["delta_alpha_rms"] for r in records)),
            "delta_alpha_max_abs_max_over_series": float(max(r["delta_alpha_max_abs"] for r in records)),
            "max_input_cmc_residual_relative": float(max(r["input_cmc_residual_relative_max"] for r in records)),
            "max_recomputed_target_residual_relative": float(max(r["mapped_recomputed_target_residual_relative_max"] for r in records)),
            "max_same_target_residual_relative": float(max(r["mapped_same_target_residual_relative_max"] for r in records)),
            "central_lapse_correction_integral_signed": central_correction_integral,
            "central_lapse_correction_integral_abs": central_abs_integral,
            "central_abs_integral_fraction_of_matched_tau_mismatch": central_abs_integral / REFERENCE_TAU_MISMATCH,
            "central_signed_integral_fraction_of_matched_tau_mismatch": central_correction_integral / REFERENCE_TAU_MISMATCH,
            "matched_R80_minus_R160_tau_at_t12": REFERENCE_TAU_MISMATCH,
            "baseline_R80_minus_R160_tau_at_t12": BASELINE_TAU_MISMATCH,
            "archived_tau_rate_min": float(np.min(tau_rates)),
            "archived_tau_rate_max": float(np.max(tau_rates)),
            "interpretation_limit": (
                "The correction integral holds each saved non-lapse state fixed and "
                "does not re-evolve the equations. It is a scale diagnostic, not a "
                "counterfactual trajectory or causal decomposition of accumulated tau."
            ),
        },
        "checks": {
            "all_snapshots_finite": bool(finite),
            "same_target_linear_equation_closes": bool(same_target_closed),
            "time_strictly_increasing": bool(np.all(np.diff(times) > 0.0)),
        },
    }
    payload["pass"] = all(payload["checks"].values())
    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / "cmc_strongD_snapshot_series.json"
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    print(f"RESULT_JSON={out_path}", flush=True)
    if not payload["pass"]:
        raise SystemExit("FAIL: snapshot series finite/closure/time gate failed")
    print("PASS: strong-D lapse re-solve series complete; saved-state fields were not evolved")


if __name__ == "__main__":
    main()
