"""Diagnostic-only repair ON/OFF state-array comparison at registered times."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from engine.production_kernel import V55ProductionKernel
from engine.step1_hamiltonian_terms import NoRepairKernel

D_AMPLITUDE = 1.0e-4
RESOLUTION = 160
R_MAX = 80.0
CFL = 0.0075
FINAL_TIME = 22.5
SAMPLE_TIMES = (0.01, 4.0, 8.0, 12.0, 16.0, 20.0, 22.5)
TIME_TOLERANCE = 0.05


def _json_safe(value):
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)):
        if math.isnan(float(value)):
            return "NaN"
        if math.isinf(float(value)):
            return "Infinity" if value > 0 else "-Infinity"
        return float(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def _arrays_from_object(obj, prefix: str, out: dict[str, list[float]]) -> None:
    """Copy every numeric ndarray field recursively; never mutate the source."""
    if isinstance(obj, np.ndarray):
        if np.issubdtype(obj.dtype, np.number):
            out[prefix] = _json_safe(np.asarray(obj).tolist())
        return
    if isinstance(obj, dict):
        for key in sorted(obj):
            _arrays_from_object(obj[key], f"{prefix}.{key}", out)
        return
    if hasattr(obj, "__dict__"):
        for key, value in sorted(vars(obj).items()):
            _arrays_from_object(value, f"{prefix}.{key}", out)


def snapshot_state(state, sample_index: int, expected_target: float) -> dict:
    arrays: dict[str, list[float]] = {
        "grid.centers": np.asarray(state.grid.centers, dtype=float).tolist(),
        "grid.volumes": np.asarray(state.grid.volumes, dtype=float).tolist(),
    }
    for name in ("geometry", "scalars", "matter"):
        _arrays_from_object(getattr(state, name), name, arrays)
    return {
        "sample_index": int(sample_index),
        "expected_target_time": float(expected_target),
        "actual_t": float(state.t),
        "tau": float(state.tau),
        "e_folds": float(state.e_folds),
        "arrays": arrays,
    }


def run_case(*, repair: bool) -> dict:
    Kernel = V55ProductionKernel if repair else NoRepairKernel
    kernel = Kernel()
    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        D_amplitude=D_AMPLITUDE,
        include_radiation=True,
    )
    dt_nominal = CFL * state.grid.dr
    samples = [snapshot_state(state, 0, 0.0)]
    target_index = 0
    failure = None
    while state.t < FINAL_TIME - 1.0e-14:
        dt = min(dt_nominal, FINAL_TIME - state.t)
        try:
            state = kernel.step(state, dt)
            state.geometry.assert_finite_positive()
        except (FloatingPointError, ValueError) as exc:
            failure = {
                "t": float(state.t),
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
            break
        if (target_index < len(SAMPLE_TIMES)
                and state.t + 0.5 * dt >= SAMPLE_TIMES[target_index]):
            samples.append(snapshot_state(
                state, target_index + 1, SAMPLE_TIMES[target_index]
            ))
            target_index += 1
    return {
        "status": "numerical_failure" if failure else "completed",
        "failure": failure,
        "repair": bool(repair),
        "D_amplitude": D_AMPLITUDE,
        "resolution": RESOLUTION,
        "r_max": R_MAX,
        "cfl": CFL,
        "include_radiation": True,
        "requested_final_time": FINAL_TIME,
        "final_time": float(state.t),
        "sample_times_requested": [0.0, *SAMPLE_TIMES],
        "samples": samples,
    }


def admission_report(report: dict) -> dict:
    expected = [0.0, *SAMPLE_TIMES]
    samples = report.get("samples", [])
    schedule_ok = len(samples) == len(expected) and all(
        abs(float(row["actual_t"]) - target) <= TIME_TOLERANCE
        for row, target in zip(samples, expected)
    )
    config_ok = (
        report.get("D_amplitude") == D_AMPLITUDE
        and report.get("resolution") == RESOLUTION
        and report.get("r_max") == R_MAX
        and report.get("cfl") == CFL
        and report.get("include_radiation") is True
        and report.get("requested_final_time") == FINAL_TIME
        and abs(float(report.get("final_time", float("nan"))) - FINAL_TIME) <= 1e-10
    )
    arrays_ok = bool(samples)
    finite_ok = bool(samples)
    for row in samples:
        arrays = row.get("arrays", {})
        if not arrays or any(len(values) != RESOLUTION for values in arrays.values()):
            arrays_ok = False
        for values in arrays.values():
            try:
                if not np.all(np.isfinite(np.asarray(values, dtype=float))):
                    finite_ok = False
            except (TypeError, ValueError):
                finite_ok = False
    checks = {
        "completed_without_numerical_failure": (
            report.get("status") == "completed" and report.get("failure") is None
        ),
        "frozen_configuration_matches": config_ok,
        "all_eight_registered_samples_present_within_tolerance": schedule_ok,
        "all_state_arrays_cover_160_cells": arrays_ok,
        "all_recorded_array_values_finite": finite_ok,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "tolerance_notes": {
            "sample_time_absolute_tolerance": TIME_TOLERANCE,
            "no_threshold_on_state_differences": True,
            "no_physical_residual_threshold": True,
        },
    }


def compare_runs(on: dict, off: dict) -> dict:
    """Report full-array differences; deliberately apply no magnitude threshold."""
    rows = []
    for on_row, off_row in zip(on.get("samples", []), off.get("samples", [])):
        on_arrays = on_row["arrays"]
        off_arrays = off_row["arrays"]
        field_rows = {}
        for name in sorted(set(on_arrays) | set(off_arrays)):
            if name not in on_arrays or name not in off_arrays:
                field_rows[name] = {"comparable": False, "reason": "field missing on one side"}
                continue
            a = np.asarray(on_arrays[name], dtype=float)
            b = np.asarray(off_arrays[name], dtype=float)
            if a.shape != b.shape:
                field_rows[name] = {
                    "comparable": False, "on_shape": list(a.shape), "off_shape": list(b.shape)
                }
                continue
            delta = b - a
            flat = int(np.argmax(np.abs(delta))) if delta.size else 0
            field_rows[name] = {
                "comparable": True,
                "shape": list(a.shape),
                "bitwise_equal_values": bool(np.array_equal(a, b)),
                "max_abs_difference": float(np.max(np.abs(delta))) if delta.size else 0.0,
                "rms_difference": float(np.sqrt(np.mean(delta * delta))) if delta.size else 0.0,
                "max_difference_flat_index": flat,
                "signed_off_minus_on_at_max": float(delta.flat[flat]) if delta.size else 0.0,
                "on_value_at_max": float(a.flat[flat]) if delta.size else 0.0,
                "off_value_at_max": float(b.flat[flat]) if delta.size else 0.0,
            }
        rows.append({
            "sample_index": int(on_row["sample_index"]),
            "t_on": float(on_row["actual_t"]),
            "t_off": float(off_row["actual_t"]),
            "absolute_time_difference": abs(float(on_row["actual_t"]) - float(off_row["actual_t"])),
            "field_differences": field_rows,
        })
    return {
        "sample_count_compared": len(rows),
        "no_state_difference_threshold_applied": True,
        "samples": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="runs/state-array-ab/state-array-ab.json")
    args = parser.parse_args()

    # Sequential by design: do not run two expensive evolutions concurrently.
    on = run_case(repair=True)
    on["admission"] = admission_report(on)
    off = None
    if on["admission"]["passed"]:
        off = run_case(repair=False)
        off["admission"] = admission_report(off)

    result = {
        "diagnostic_only": True,
        "configuration": {
            "D_amplitude": D_AMPLITUDE,
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "cfl": CFL,
            "include_radiation": True,
            "final_time": FINAL_TIME,
            "target_times": [0.0, *SAMPLE_TIMES],
        },
        "repair_on": on,
        "repair_off": off,
        "comparison": compare_runs(on, off) if off is not None else None,
        "interpretation_limits": [
            "This is a state-array diagnostic, not a physics result.",
            "No state-difference or physical-residual threshold is applied.",
            "No production physics code or evolution equation is changed.",
            "If ON fails admission, OFF is not launched.",
        ],
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "repair_on_admission": on["admission"],
        "repair_off_admission": off["admission"] if off is not None else None,
        "samples_compared": result["comparison"]["sample_count_compared"] if off is not None else 0,
    }, indent=2))
    if not on["admission"]["passed"] or off is None or not off["admission"]["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
