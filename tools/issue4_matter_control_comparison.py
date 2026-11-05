#!/usr/bin/env python3
"""Create a pinned side-by-side comparison of both controls and accepted t=0.60."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def _load(path: Path) -> dict:
    payload = path.read_bytes()
    report = json.loads(payload)
    if report.get("status") != "diagnostic_evaluation_complete":
        raise RuntimeError(f"{path} is not a completed diagnostic evaluation")
    if report.get("production_steps_run") != 0 or report.get("state_advanced") is not False:
        raise RuntimeError(f"{path} violates zero-step contract")
    if report.get("zero_step_tripwire") != "armed_and_not_triggered":
        raise RuntimeError(f"{path} zero-step tripwire is not clean")
    for epsilon, channels in report["channels"].items():
        for channel in ("pointwise", "finite_volume"):
            values = channels[channel]
            for key in ("predictor", "stress_source", "residual_tangent_minus_predictor"):
                if not np.all(np.isfinite(np.asarray(values[key], dtype=np.float64))):
                    raise RuntimeError(f"{path} contains non-finite {epsilon}/{channel}/{key}")
    return {"report": report, "sha256": hashlib.sha256(payload).hexdigest()}


def _channel_summary(report: dict, channel: str, epsilon: str) -> dict:
    summary = report["summaries"][epsilon][channel]
    residual = summary["residual_tangent_minus_predictor"]
    predictor = summary["predictor"]
    return {
        "max_abs_residual_r_le_10": residual["max_abs_r_le_10"],
        "radius_at_max_abs_residual_r_le_10": residual["radius_at_max_r_le_10"],
        "max_abs_residual_annulus_5_to_10": residual["max_abs_annulus_5_to_10"],
        "rms_residual_annulus_5_to_10": residual["rms_annulus_5_to_10"],
        "max_abs_predictor_r_le_10": predictor["max_abs_r_le_10"],
        "radius_at_max_abs_predictor_r_le_10": predictor["radius_at_max_r_le_10"],
    }


def compare(arbitrary_path: Path, equation_path: Path, accepted_path: Path) -> dict:
    arbitrary = _load(arbitrary_path)
    equation = _load(equation_path)
    accepted = _load(accepted_path)
    reports = {
        "arbitrary_manufactured_family": arbitrary["report"],
        "equation_consistent_matter_direction": equation["report"],
        "accepted_t_0p60": accepted["report"],
    }
    epsilons = ("eps2e-05", "eps1e-05")
    rows = {}
    for label, report in reports.items():
        rows[label] = {
            "state_class": report["state_class"],
            "tangent_semantics": report["tangent_semantics"],
            "geometry_tangent_semantics": (
                "frozen_flat_background" if label != "accepted_t_0p60"
                else "accepted_full_metric_state_and_direction"
            ),
            "epsilon_summaries": {
                epsilon: {
                    channel: _channel_summary(report, channel, epsilon)
                    for channel in ("pointwise", "finite_volume")
                }
                for epsilon in epsilons
            },
            "production_steps_run": report["production_steps_run"],
            "state_advanced": report["state_advanced"],
            "zero_step_tripwire": report["zero_step_tripwire"],
        }

    deltas = {}
    for epsilon in epsilons:
        deltas[epsilon] = {}
        a = arbitrary["report"]["channels"][epsilon]
        e = equation["report"]["channels"][epsilon]
        for channel in ("pointwise", "finite_volume"):
            deltas[epsilon][channel] = {}
            for key in ("measured_tangent",):
                av = np.asarray(a[key], dtype=np.float64)
                ev = np.asarray(e[key], dtype=np.float64)
                deltas[epsilon][channel][key + "_max_abs_equation_minus_arbitrary"] = float(np.max(np.abs(ev-av)))
            for key in ("predictor", "residual_tangent_minus_predictor"):
                av = np.asarray(a[channel][key], dtype=np.float64)
                ev = np.asarray(e[channel][key], dtype=np.float64)
                deltas[epsilon][channel][key + "_max_abs_equation_minus_arbitrary"] = float(np.max(np.abs(ev-av)))

    return {
        "schema": "issue4-three-class-full-bianchi-comparison/v1",
        "status": "comparison_complete_diagnostic_only",
        "classes": rows,
        "manufactured_direction_deltas": deltas,
        "interpretation": {
            "arbitrary_family": "arbitrary manufactured control; not matter-equation consistent",
            "equation_consistent_family": (
                "matter direction satisfies pinned scalar spatial equations and fluid "
                "continuity/Valencia momentum equations, but geometry is frozen flat"
            ),
            "accepted_state": (
                "source-pinned accepted t=0.60 full-metric state; its report is evaluated "
                "through the same current full-Bianchi evaluator"
            ),
            "pointwise_invariance": (
                "compare pointwise residual arrays directly; if the residual is unchanged "
                "while tangent and source move together, this is a frozen-geometry "
                "directional cancellation, not a physical closure result"
            ),
            "admission": "all scientific admission and conservation-closure flags remain false",
        },
        "input_sha256": {
            "arbitrary_report": arbitrary["sha256"],
            "equation_consistent_report": equation["sha256"],
            "accepted_report": accepted["sha256"],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arbitrary", type=Path, required=True)
    parser.add_argument("--equation-consistent", type=Path, required=True)
    parser.add_argument("--accepted", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.arbitrary, args.equation_consistent, args.accepted)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({
        "status": result["status"],
        "classes": sorted(result["classes"]),
        "pointwise_residual_delta_eps1e-05": result["manufactured_direction_deltas"]["eps1e-05"]["pointwise"]["residual_tangent_minus_predictor_max_abs_equation_minus_arbitrary"],
    }, indent=2))


if __name__ == "__main__":
    main()
