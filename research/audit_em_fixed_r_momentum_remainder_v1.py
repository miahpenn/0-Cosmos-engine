#!/usr/bin/env python3
"""Report paired fixed-R momentum-remainder trends from saved ledger artifacts.

This tool performs no evolution and makes no physical-residual admission
decision. It validates artifact provenance and algebraic bookkeeping, then
reports neutral/charged differences across spatial and temporal sweeps.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
ROUND_OFF_ULPS = 128.0
TIME_ABS_TOL = 1.0e-12
RUN_GROUPS = (
    "runs",
    "fixed_N_temporal_refinement_runs",
    "domain_radius_sweep_extension_runs",
)
Q_NEUTRAL = 0.0
Q_CHARGED = 0.25
OBSERVATION_FIELDS = (
    "tube_fixed_R_mass_rate_residual",
    "tube_momentum_constraint_at_R_star",
    "tube_native_mass_gradient_identity_residual",
    "tube_hamiltonian_constraint_at_R_star",
    "tube_stress_flux_composite_interp",
    "tube_radial_pressure_work_composite_interp",
    "tube_fixed_R_mass_rate_finite_difference",
    "tube_predicted_fixed_R_mass_rate",
)


def _number(value: Any, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _optional_number(value: Any, label: str) -> float | None:
    if value is None:
        return None
    return _number(value, label)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read JSON artifact {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"artifact {path} must contain a JSON object")
    return value


def _validate_provenance(
    ledger: dict[str, Any], fixed_r: dict[str, Any]
) -> tuple[str, int]:
    provenance = ledger.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError("joined ledger is missing provenance")
    source_sha = provenance.get("source_head_sha")
    if not isinstance(source_sha, str) or len(source_sha) != 40:
        raise ValueError("joined ledger has no 40-character source_head_sha")
    try:
        int(source_sha, 16)
    except ValueError as exc:
        raise ValueError("source_head_sha is not hexadecimal") from exc

    pins = provenance.get("production_blob_pins")
    if not isinstance(pins, dict) or not pins:
        raise ValueError("joined ledger has no production_blob_pins")
    for name, pin in pins.items():
        if not isinstance(pin, dict):
            raise ValueError(f"invalid source pin record: {name}")
        if pin.get("match") is not True or pin.get("expected") != pin.get("observed"):
            raise ValueError(f"production source pin mismatch: {name}")

    if fixed_r.get("source_head_sha") != source_sha:
        raise ValueError("joined ledger and fixed-R report source heads differ")
    if fixed_r.get("production_source_pins_all_match") is not True:
        raise ValueError("fixed-R report did not verify production source pins")
    if fixed_r.get("vendor_gitlink_match") is not True:
        raise ValueError("fixed-R report did not verify the vendor gitlink")
    if fixed_r.get("physical_residual_threshold_applied") is not False:
        raise ValueError("fixed-R report must preserve unthresholded physical residuals")
    return source_sha, len(pins)


def _config_key(run: dict[str, Any]) -> tuple[Any, ...]:
    return (
        int(run.get("resolution", -1)),
        _optional_number(run.get("cfl"), "run.cfl"),
        _optional_number(run.get("nominal_dt"), "run.nominal_dt"),
        _optional_number(
            run.get("domain_radius_sweep_r_max"), "run.domain_radius_sweep_r_max"
        ),
        int(run.get("steps", -1)),
        _number(run.get("t_final"), "run.t_final"),
    )


def _config_dict(group: str, run: dict[str, Any]) -> dict[str, Any]:
    return {
        "run_group": group,
        "resolution": int(run["resolution"]),
        "cfl": _optional_number(run.get("cfl"), "run.cfl"),
        "nominal_dt": _optional_number(run.get("nominal_dt"), "run.nominal_dt"),
        "domain_radius_sweep_r_max": _optional_number(
            run.get("domain_radius_sweep_r_max"), "run.domain_radius_sweep_r_max"
        ),
        "steps": int(run["steps"]),
        "t_final": _number(run.get("t_final"), "run.t_final"),
    }


def _roundoff_check(
    label: str, observed: float, expected: float, operand_scale: float
) -> dict[str, Any]:
    residual = observed - expected
    scale = max(
        abs(operand_scale),
        abs(observed) + abs(expected),
        sys.float_info.min,
    )
    bound = ROUND_OFF_ULPS * sys.float_info.epsilon * scale
    ratio = 0.0 if residual == 0.0 else abs(residual) / bound
    if not math.isfinite(ratio) or ratio > 1.0:
        raise ArithmeticError(
            f"{label} failed operand-scaled roundoff check: residual={residual:.17e}, "
            f"bound={bound:.17e}, ratio={ratio:.17e}"
        )
    return {
        "passed": True,
        "residual": residual,
        "roundoff_bound": bound,
        "bound_ratio": ratio,
        "operand_scale": scale,
    }


def _paired_rows(
    group: str, run_neutral: dict[str, Any], run_charged: dict[str, Any]
) -> list[dict[str, Any]]:
    neutral_rows = run_neutral.get("accepted_slice_rows")
    charged_rows = run_charged.get("accepted_slice_rows")
    if not isinstance(neutral_rows, list) or not isinstance(charged_rows, list):
        raise ValueError(f"{group} run is missing accepted_slice_rows")
    if not neutral_rows or len(neutral_rows) != len(charged_rows):
        raise ValueError(f"{group} neutral/charged row counts differ or are empty")

    neutral_rows = sorted(neutral_rows, key=lambda row: _number(row.get("t"), "neutral.t"))
    charged_rows = sorted(charged_rows, key=lambda row: _number(row.get("t"), "charged.t"))
    config = _config_dict(group, run_neutral)
    output: list[dict[str, Any]] = []
    for index, (neutral, charged) in enumerate(zip(neutral_rows, charged_rows)):
        if not isinstance(neutral, dict) or not isinstance(charged, dict):
            raise ValueError(f"{group} accepted row {index} is not an object")
        neutral_t = _number(neutral.get("t"), f"{group}.neutral[{index}].t")
        charged_t = _number(charged.get("t"), f"{group}.charged[{index}].t")
        if abs(neutral_t - charged_t) > TIME_ABS_TOL:
            raise ValueError(
                f"{group} neutral/charged times do not match: "
                f"{neutral_t:.17e} versus {charged_t:.17e}"
            )
        row: dict[str, Any] = {
            **config,
            "time": 0.5 * (neutral_t + charged_t),
            "time_pair_abs_difference": abs(neutral_t - charged_t),
            "q_over_m_neutral": Q_NEUTRAL,
            "q_over_m_charged": Q_CHARGED,
        }
        for field in OBSERVATION_FIELDS:
            if field not in neutral or field not in charged:
                if field in ("tube_fixed_R_mass_rate_residual", "tube_momentum_constraint_at_R_star"):
                    raise ValueError(f"{group} accepted row is missing required field {field}")
                continue
            neutral_value = _number(neutral[field], f"neutral.{field}")
            charged_value = _number(charged[field], f"charged.{field}")
            row[f"neutral_{field}"] = neutral_value
            row[f"charged_{field}"] = charged_value
            row[f"delta_{field}"] = charged_value - neutral_value
        output.append(row)
    return output


def _collect_paired_observations(ledger: dict[str, Any]) -> list[dict[str, Any]]:
    all_rows: list[dict[str, Any]] = []
    settings = ledger.get("settings", {})
    if not isinstance(settings, dict):
        raise ValueError("joined ledger settings must be an object")
    charge_arms = settings.get("q_over_m_arms")
    if not isinstance(charge_arms, list) or not any(
        isinstance(value, (int, float)) and not isinstance(value, bool) and float(value) == Q_NEUTRAL
        for value in charge_arms
    ) or not any(
        isinstance(value, (int, float)) and not isinstance(value, bool) and float(value) == Q_CHARGED
        for value in charge_arms
    ):
        raise ValueError("joined ledger must declare q/m=0 and q/m=0.25 comparison arms")

    for group in RUN_GROUPS:
        runs = ledger.get(group)
        if not isinstance(runs, list):
            raise ValueError(f"joined ledger is missing run group {group}")
        configs: dict[tuple[Any, ...], dict[float, dict[str, Any]]] = {}
        for index, run in enumerate(runs):
            if not isinstance(run, dict):
                raise ValueError(f"{group}[{index}] is not an object")
            q_over_m = _number(run.get("q_over_m"), f"{group}[{index}].q_over_m")
            if q_over_m not in (Q_NEUTRAL, Q_CHARGED):
                continue
            key = _config_key(run)
            arms = configs.setdefault(key, {})
            if q_over_m in arms:
                raise ValueError(f"duplicate q/m={q_over_m} run for {group} config {key}")
            arms[q_over_m] = run
        for key, arms in configs.items():
            if set(arms) != {Q_NEUTRAL, Q_CHARGED}:
                raise ValueError(f"unpaired neutral/charged run for {group} config {key}")
            all_rows.extend(_paired_rows(group, arms[Q_NEUTRAL], arms[Q_CHARGED]))
    if not all_rows:
        raise ValueError("joined ledger contained no paired neutral/charged observations")
    return sorted(
        all_rows,
        key=lambda row: (
            row["run_group"],
            row["resolution"],
            row["cfl"] if row["cfl"] is not None else -1.0,
            row["domain_radius_sweep_r_max"]
            if row["domain_radius_sweep_r_max"] is not None
            else -1.0,
            row["time"],
        ),
    )


def _identity_reconciliation(fixed_r: dict[str, Any]) -> dict[str, Any]:
    comparison = fixed_r.get("charged_minus_neutral")
    arms = fixed_r.get("arms")
    if not isinstance(comparison, dict) or not isinstance(arms, dict):
        raise ValueError("fixed-R report is missing arms or charged_minus_neutral")
    try:
        neutral = arms["q_over_m_0"]["recomputed"]["interpolated"]
        charged = arms["q_over_m_0.25"]["recomputed"]["interpolated"]
    except (KeyError, TypeError) as exc:
        raise ValueError("fixed-R report is missing neutral/charged interpolants") from exc
    if not isinstance(neutral, dict) or not isinstance(charged, dict):
        raise ValueError("fixed-R interpolants must be objects")

    def v(source: dict[str, Any], key: str) -> float:
        return _number(source.get(key), f"fixed-R.{key}")

    coordinate = v(comparison, "coordinate_mass_rate")
    flux = v(comparison, "composite_flux")
    work = v(comparison, "composite_pressure_work")
    observed = v(comparison, "composite_fixed_R_residual_unthresholded")
    momentum_term = v(comparison, "momentum_constraint_mass_identity_term_composite_interp")
    radial_rr_term = v(comparison, "radial_rr_einstein_mass_identity_term_composite_interp")
    predicted = v(comparison, "geometric_identity_predicted_mass_balance_remainder_composite_interp")
    gap = v(comparison, "geometric_identity_minus_observed_mass_balance_remainder_composite_interp_unthresholded")

    arm_scale = sum(
        abs(v(source, key))
        for source in (neutral, charged)
        for key in (
            "coordinate_mass_rate",
            "composite_flux",
            "composite_pressure_work",
            "composite_fixed_R_residual_unthresholded",
            "momentum_constraint_mass_identity_term_composite_interp",
            "radial_rr_einstein_mass_identity_term_composite_interp",
            "geometric_identity_predicted_mass_balance_remainder_composite_interp",
            "geometric_identity_minus_observed_mass_balance_remainder_composite_interp_unthresholded",
        )
    )
    checks = {
        "paired_coordinate_mass_rate": _roundoff_check(
            "paired coordinate mass rate",
            coordinate,
            v(charged, "coordinate_mass_rate") - v(neutral, "coordinate_mass_rate"),
            arm_scale,
        ),
        "paired_flux": _roundoff_check(
            "paired composite flux",
            flux,
            v(charged, "composite_flux") - v(neutral, "composite_flux"),
            arm_scale,
        ),
        "paired_pressure_work": _roundoff_check(
            "paired composite pressure work",
            work,
            v(charged, "composite_pressure_work") - v(neutral, "composite_pressure_work"),
            arm_scale,
        ),
        "observed_fixed_R_formula": _roundoff_check(
            "observed fixed-R residual formula",
            observed,
            coordinate - flux - work,
            arm_scale,
        ),
        "paired_observed_fixed_R_residual": _roundoff_check(
            "paired fixed-R residual",
            observed,
            v(charged, "composite_fixed_R_residual_unthresholded")
            - v(neutral, "composite_fixed_R_residual_unthresholded"),
            arm_scale,
        ),
        "paired_momentum_identity_term": _roundoff_check(
            "paired momentum-constraint mass-identity term",
            momentum_term,
            v(charged, "momentum_constraint_mass_identity_term_composite_interp")
            - v(neutral, "momentum_constraint_mass_identity_term_composite_interp"),
            arm_scale,
        ),
        "paired_radial_RR_identity_term": _roundoff_check(
            "paired radial-RR mass-identity term",
            radial_rr_term,
            v(charged, "radial_rr_einstein_mass_identity_term_composite_interp")
            - v(neutral, "radial_rr_einstein_mass_identity_term_composite_interp"),
            arm_scale,
        ),
        "geometric_prediction_split": _roundoff_check(
            "geometric identity split",
            predicted,
            momentum_term + radial_rr_term,
            arm_scale,
        ),
        "paired_geometric_prediction": _roundoff_check(
            "paired geometric prediction",
            predicted,
            v(charged, "geometric_identity_predicted_mass_balance_remainder_composite_interp")
            - v(neutral, "geometric_identity_predicted_mass_balance_remainder_composite_interp"),
            arm_scale,
        ),
        "predicted_minus_observed_gap": _roundoff_check(
            "predicted-minus-observed remainder gap",
            gap,
            predicted - observed,
            arm_scale,
        ),
        "paired_gap": _roundoff_check(
            "paired geometric identity gap",
            gap,
            v(charged, "geometric_identity_minus_observed_mass_balance_remainder_composite_interp_unthresholded")
            - v(neutral, "geometric_identity_minus_observed_mass_balance_remainder_composite_interp_unthresholded"),
            arm_scale,
        ),
    }
    neutral_target = _number(
        arms["q_over_m_0"]["recomputed"].get("target_areal_radius"),
        "fixed-R neutral target_areal_radius",
    )
    charged_target = _number(
        arms["q_over_m_0.25"]["recomputed"].get("target_areal_radius"),
        "fixed-R charged target_areal_radius",
    )
    if abs(neutral_target - charged_target) > TIME_ABS_TOL:
        raise ValueError("neutral and charged fixed-R targets differ")
    return {
        "target_areal_radius": neutral_target,
        "neutral": {
            "composite_fixed_R_residual_unthresholded": v(
                neutral, "composite_fixed_R_residual_unthresholded"
            ),
            "momentum_constraint_mass_identity_term": v(
                neutral, "momentum_constraint_mass_identity_term_composite_interp"
            ),
            "radial_rr_einstein_mass_identity_term": v(
                neutral, "radial_rr_einstein_mass_identity_term_composite_interp"
            ),
        },
        "charged": {
            "composite_fixed_R_residual_unthresholded": v(
                charged, "composite_fixed_R_residual_unthresholded"
            ),
            "momentum_constraint_mass_identity_term": v(
                charged, "momentum_constraint_mass_identity_term_composite_interp"
            ),
            "radial_rr_einstein_mass_identity_term": v(
                charged, "radial_rr_einstein_mass_identity_term_composite_interp"
            ),
        },
        "charged_minus_neutral": {
            "coordinate_mass_rate": coordinate,
            "composite_flux": flux,
            "composite_pressure_work": work,
            "observed_fixed_R_residual_unthresholded": observed,
            "momentum_constraint_mass_identity_term": momentum_term,
            "radial_rr_einstein_mass_identity_term": radial_rr_term,
            "geometric_identity_predicted_remainder": predicted,
            "predicted_minus_observed_gap_unthresholded": gap,
        },
        "algebraic_roundoff_checks": checks,
        "all_algebraic_roundoff_checks_passed": all(
            check["passed"] for check in checks.values()
        ),
        "physical_residual_threshold_applied": False,
        "interpretation": (
            "This validates bookkeeping identities among the recorded operands only. "
            "The charged-minus-neutral fixed-R remainder and momentum-constraint "
            "contribution are unthresholded physical/numerical diagnostics; this "
            "does not establish global energy conservation."
        ),
    }


def build_report(
    ledger: dict[str, Any],
    fixed_r: dict[str, Any],
    ledger_sha256: str,
    fixed_r_sha256: str,
) -> dict[str, Any]:
    source_sha, pin_count = _validate_provenance(ledger, fixed_r)
    observations = _collect_paired_observations(ledger)
    identity = _identity_reconciliation(fixed_r)

    def final_rows(group: str) -> list[dict[str, Any]]:
        return [
            row for row in observations
            if row["run_group"] == group
            and abs(row["time"] - row["t_final"]) <= TIME_ABS_TOL
        ]

    spatial = final_rows("runs")
    temporal = final_rows("fixed_N_temporal_refinement_runs")
    radius = final_rows("domain_radius_sweep_extension_runs")
    return {
        "schema_version": SCHEMA_VERSION,
        "classification": "research-only paired fixed-R momentum-remainder trend audit",
        "provenance": {
            "source_head_sha": source_sha,
            "production_blob_pin_count": pin_count,
            "all_production_blob_pins_match": True,
            "vendor_gitlink_match": True,
            "joined_ledger_sha256": ledger_sha256,
            "fixed_r_recomputation_sha256": fixed_r_sha256,
        },
        "summary": {
            "paired_configuration_count": len(
                {
                    (
                        row["run_group"],
                        row["resolution"],
                        row["cfl"],
                        row["nominal_dt"],
                        row["domain_radius_sweep_r_max"],
                        row["steps"],
                        row["t_final"],
                    )
                    for row in observations
                }
            ),
            "paired_time_sample_count": len(observations),
            "spatial_refinement_final_sample_count": len(spatial),
            "temporal_refinement_final_sample_count": len(temporal),
            "domain_radius_final_sample_count": len(radius),
            "fixed_R_identity_algebraic_checks_passed": identity[
                "all_algebraic_roundoff_checks_passed"
            ],
            "physical_residual_threshold_applied": False,
            "interpretation": (
                "Paired observations are reported across matched saved times, spatial "
                "resolutions and timesteps. Algebraic checks fail closed on inconsistent "
                "artifact bookkeeping; physical residual magnitudes are never thresholded."
            ),
        },
        "spatial_refinement_final_time": spatial,
        "temporal_refinement_final_time": temporal,
        "domain_radius_final_time": radius,
        "paired_time_history": observations,
        "fixed_R_identity_reconciliation": identity,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="raw joined-ledger JSON")
    parser.add_argument(
        "--fixed-r-input", required=True, type=Path, help="independent fixed-R recomputation JSON"
    )
    parser.add_argument("--output", required=True, type=Path, help="output report JSON")
    args = parser.parse_args()
    try:
        ledger = _load_json(args.input)
        fixed_r = _load_json(args.fixed_r_input)
        report = build_report(
            ledger,
            fixed_r,
            _sha256(args.input),
            _sha256(args.fixed_r_input),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    except (ValueError, ArithmeticError, KeyError, TypeError) as exc:
        print(f"fixed-R momentum remainder audit failed closed: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report["summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
