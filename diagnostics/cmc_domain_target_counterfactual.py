"""Freeze-state test of R80 vs R160 global CMC target on matched snapshots.

At each shared saved time, keep the R80 geometry/scalars/matter fixed and the
outer alpha set to the R160 value at the common R=79.75 cell. Solve the lapse
once with R80's native proper-volume Kdot target and once with the archived
R160 full-domain Kdot target. No state evolution or production edit occurs.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import engine.cmc_gauge as cmc
from diagnostics.cmc_strongD_snapshot_jacobian import load_state, solve_map

OUT = ROOT / "runs" / "cmc_domain_target_counterfactual"
R80_N = 160
R160_N = 320
RMAX_R80 = 80.0
TAU_MISMATCH_T12 = 0.0106210937742


def read_rows(path):
    rows = []
    with Path(path).open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("state_snapshot") is not None:
                rows.append(row)
    return rows


def unique_t_map(rows, label):
    result = {}
    for row in rows:
        key = round(float(row["t"]), 8)
        if key in result:
            raise ValueError(f"{label}: duplicate snapshot time {key}")
        result[key] = row
    return result


def profile_error_metrics(prediction, reference, dr):
    d = np.asarray(prediction, dtype=float) - np.asarray(reference, dtype=float)
    i = int(np.argmax(np.abs(d)))
    return {
        "rms": float(np.sqrt(np.mean(d * d))),
        "max_abs": float(np.max(np.abs(d))),
        "max_abs_radius": float((i + 0.5) * dr),
        "center_signed": float(d[0]),
    }


def solve_with_fixed_target(state, geometry, alpha_input, outer_alpha, kdot):
    original = cmc.target_kdot

    def target_override(grid, geom, scalars, matter, outer_frac=0.20):
        return float(kdot)

    cmc.target_kdot = target_override
    try:
        return solve_map(state, geometry, alpha_input, outer_alpha)
    finally:
        cmc.target_kdot = original


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: python diagnostics/cmc_domain_target_counterfactual.py <r80-profiles.jsonl> <r160-profiles.jsonl>"
        )
    r80_path, r160_path = map(Path, sys.argv[1:])
    r80_rows = read_rows(r80_path)
    r160_rows = read_rows(r160_path)
    r80 = unique_t_map(r80_rows, "R80")
    r160 = unique_t_map(r160_rows, "R160")
    common = sorted(set(r80) & set(r160))
    if len(common) < 3:
        raise ValueError(f"insufficient matched snapshot times: {len(common)}")
    if len(common) != len(r80) or len(common) != len(r160):
        raise ValueError(
            f"R80/R160 time grids differ: R80={len(r80)}, R160={len(r160)}, common={len(common)}"
        )

    records = []
    all_finite = True
    all_native_positive = True
    forced_positive_count = 0
    for t_key in common:
        row80, row160 = r80[t_key], r160[t_key]
        state = load_state(row80)
        grid = state.grid
        geom = state.geometry.copy()
        alpha_in = np.asarray(geom.alpha, dtype=float).copy()
        alpha_ref_full = np.asarray(row160["state_snapshot"]["geometry"]["alpha"], dtype=float)
        if alpha_ref_full.shape != (R160_N,) or grid.n != R80_N:
            raise ValueError("unexpected resolution in archived full-state profiles")
        # R80 cell j=159 and R160 cell j=159 are the same coordinate r=79.75.
        outer_alpha = float(alpha_ref_full[R80_N - 1])
        target_ref = float(row160["target_Kdot_full_domain"])

        alpha_native, target_native = solve_map(
            state, geom, alpha_in, outer_alpha
        )
        native_positive = bool(np.all(np.isfinite(alpha_native)) and np.min(alpha_native) > 0.0)
        if not native_positive:
            all_native_positive = False

        forced = {
            "status": "not_run",
            "alpha_center": None,
            "alpha_min": None,
            "target_used": target_ref,
            "profile_error_vs_R160_core": None,
            "alpha_center_delta_vs_R160": None,
            "target_effect_center": None,
            "profile_error_reduction_fraction": None,
            "failure": None,
        }
        try:
            alpha_forced, target_used_forced = solve_with_fixed_target(
                state, geom, alpha_in, outer_alpha, target_ref
            )
            forced_positive_count += int(
                np.all(np.isfinite(alpha_forced)) and np.min(alpha_forced) > 0.0
            )
            forced_delta = float(alpha_forced[0] - alpha_ref_full[0])
            native_delta = float(alpha_native[0] - alpha_ref_full[0])
            ref_core = alpha_ref_full[:R80_N]
            native_error = profile_error_metrics(alpha_native, ref_core, grid.dr)
            forced_error = profile_error_metrics(alpha_forced, ref_core, grid.dr)
            denominator = native_error["rms"]
            reduction = (
                None if denominator == 0.0
                else float(1.0 - forced_error["rms"] / denominator)
            )
            forced.update({
                "status": "positive" if np.min(alpha_forced) > 0.0 else "non_positive",
                "alpha_center": float(alpha_forced[0]),
                "alpha_min": float(np.min(alpha_forced)),
                "target_used": float(target_used_forced),
                "profile_error_vs_R160_core": forced_error,
                "alpha_center_delta_vs_R160": forced_delta,
                "target_effect_center": float(alpha_forced[0] - alpha_native[0]),
                "profile_error_reduction_fraction": reduction,
                "native_alpha_center_delta_vs_R160": native_delta,
                "native_profile_error_vs_R160_core": native_error,
            })
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            forced["status"] = "failed"
            forced["failure"] = str(exc)[:1000]

        native_error = profile_error_metrics(
            alpha_native, alpha_ref_full[:R80_N], grid.dr
        )
        actual_error = profile_error_metrics(alpha_in, alpha_ref_full[:R80_N], grid.dr)
        rec = {
            "t": float(state.t),
            "tau_R80": float(row80["tau"]),
            "tau_R160": float(row160["tau"]),
            "tau_R80_minus_R160": float(row80["tau"] - row160["tau"]),
            "outer_alpha_from_R160_at_R79_75": outer_alpha,
            "outer_alpha_R80_stored": float(alpha_in[-1]),
            "outer_alpha_difference_stored_R80_minus_R160": float(alpha_in[-1] - outer_alpha),
            "Kdot_R80_native": float(target_native),
            "Kdot_R160_full_domain": target_ref,
            "Kdot_R80_minus_R160": float(target_native - target_ref),
            "alpha_center_R80_stored": float(alpha_in[0]),
            "alpha_center_R80_native_resolve": float(alpha_native[0]),
            "alpha_center_R160_stored": float(alpha_ref_full[0]),
            "actual_stored_alpha_core_error_vs_R160": actual_error,
            "native_solve_alpha_core_error_vs_R160": native_error,
            "native_solve_correction_center": float(alpha_native[0] - alpha_in[0]),
            "target_forced_solve": forced,
            "all_finite": bool(
                np.all(np.isfinite(alpha_in))
                and np.all(np.isfinite(alpha_native))
                and np.isfinite(target_native)
                and np.isfinite(target_ref)
            ),
        }
        records.append(rec)
        all_finite = all_finite and rec["all_finite"]

    times = np.asarray([r["t"] for r in records], dtype=float)
    actual_center_delta = np.asarray([
        r["alpha_center_R80_stored"] - r["alpha_center_R160_stored"]
        for r in records
    ], dtype=float)
    native_center_delta = np.asarray([
        r["alpha_center_R80_native_resolve"] - r["alpha_center_R160_stored"]
        for r in records
    ], dtype=float)
    forced_ok = forced_positive_count == len(records)
    forced_center_delta = None
    forced_center_integral = None
    forced_abs_integral = None
    if forced_ok:
        forced_center_delta = np.asarray([
            r["target_forced_solve"]["alpha_center_delta_vs_R160"]
            for r in records
        ], dtype=float)
        forced_center_integral = float(np.trapezoid(forced_center_delta, times))
        forced_abs_integral = float(np.trapezoid(np.abs(forced_center_delta), times))

    payload = {
        "test": "R80 frozen-state external global CMC target counterfactual",
        "provenance": {
            "R80_run": "https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38162318255",
            "R80_artifact_id": 11710858125,
            "R80_source_commit": "35072830f9e4d7029bb10243a0520b9a273a17dd",
            "R160_reference_artifact_id": 11694568292,
            "R80_profiles_sha256": hashlib.sha256(r80_path.read_bytes()).hexdigest(),
            "R160_profiles_sha256": hashlib.sha256(r160_path.read_bytes()).hexdigest(),
            "matched_snapshots": len(records),
            "start_t": float(times[0]),
            "end_t": float(times[-1]),
        },
        "production_equations_changed": False,
        "protocol": (
            "At each shared saved time, hold the R80 geometry/scalars/matter fixed; "
            "use R160 alpha at r=79.75 as the outer lapse for both solves; compare "
            "native R80 Kdot with R160 full-domain Kdot copied as a diagnostic target. "
            "No saved state is advanced."
        ),
        "per_snapshot": records,
        "summary": {
            "actual_tau_difference_at_final_time": float(records[-1]["tau_R80_minus_R160"]),
            "actual_central_alpha_difference_at_final_time": float(actual_center_delta[-1]),
            "max_abs_actual_central_alpha_difference": float(np.max(np.abs(actual_center_delta))),
            "max_abs_actual_central_alpha_difference_time": float(times[int(np.argmax(np.abs(actual_center_delta)))]),
            "integral_actual_stored_central_alpha_difference": float(np.trapezoid(actual_center_delta, times)),
            "integral_native_resolved_central_alpha_difference": float(np.trapezoid(native_center_delta, times)),
            "reference_tau_difference": TAU_MISMATCH_T12,
            "forced_target_positive_all_snapshots": bool(forced_ok),
            "forced_target_positive_snapshot_count": int(forced_positive_count),
            "forced_target_signed_central_alpha_difference_integral": forced_center_integral,
            "forced_target_absolute_central_alpha_difference_integral": forced_abs_integral,
            "forced_target_integral_fraction_of_tau_mismatch": (
                None if forced_center_integral is None
                else float(forced_center_integral / TAU_MISMATCH_T12)
            ),
            "max_abs_Kdot_R80_minus_R160": float(max(abs(r["Kdot_R80_minus_R160"]) for r in records)),
            "final_Kdot_R80_native": float(records[-1]["Kdot_R80_native"]),
            "final_Kdot_R160_full_domain": float(records[-1]["Kdot_R160_full_domain"]),
            "interpretation_limit": (
                "This fixed-state target swap isolates the gauge target in the lapse solve, "
                "not its feedback on later geometry/matter. It is a diagnostic counterfactual, "
                "not an evolved trajectory. A non-positive forced lapse means that target is "
                "not admissible on that frozen R80 slice with the matched outer boundary."
            ),
        },
        "checks": {
            "all_native_results_finite": bool(all_finite),
            "native_lapses_positive": bool(all_native_positive),
            "strictly_increasing_times": bool(np.all(np.diff(times) > 0.0)),
        },
    }
    payload["pass"] = all(payload["checks"].values())
    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / "cmc_domain_target_counterfactual.json"
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    print(f"RESULT_JSON={out_path}", flush=True)
    if not payload["pass"]:
        raise SystemExit("FAIL: native state or timing integrity gate")
    print("PASS: external-global-target frozen-state counterfactual completed")


if __name__ == "__main__":
    main()
