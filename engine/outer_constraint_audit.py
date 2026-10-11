"""Full-grid Hamiltonian/momentum constraint audit. Diagnostic-only; no state repair added."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.step1_hamiltonian_terms import NoRepairKernel
from engine.stress_energy import assemble_total_stress_energy

# Frozen campaign configuration; keep aligned with the admitted Step 1 run.
D_AMPLITUDE = 1.0e-4
RESOLUTION = 160
R_MAX = 80.0
CFL = 0.0075
FINAL_TIME = 22.5
SAMPLE_TIMES = (0.01, 4.0, 8.0, 12.0, 16.0, 20.0, 22.5)
TIME_TOLERANCE = 0.05
DECOMPOSITION_TOLERANCE = 1.0e-12


def _weighted_l2(values: np.ndarray, volumes: np.ndarray, mask: np.ndarray) -> float:
    weights = np.asarray(volumes[mask], dtype=float)
    selected = np.asarray(values[mask], dtype=float)
    return float(np.sqrt(
        np.sum(weights * selected * selected)
        / max(float(np.sum(weights)), 1.0e-300)
    ))


def _summary(values: np.ndarray, radii: np.ndarray, volumes: np.ndarray,
             mask: np.ndarray) -> dict:
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        return {"cell_count": 0, "max_abs": None, "max_cell": None,
                "max_r": None, "volume_weighted_l2": None}
    local_index = int(np.argmax(np.abs(values[indices])))
    index = int(indices[local_index])
    return {
        "cell_count": int(indices.size),
        "max_abs": float(abs(values[index])),
        "max_cell": index,
        "max_r": float(radii[index]),
        "signed_value_at_max": float(values[index]),
        "volume_weighted_l2": _weighted_l2(values, volumes, mask),
    }


def record_sample(state, sample_index: int, expected_target: float) -> dict:
    """Record every cell and regional witnesses without modifying the state."""
    grid = state.grid
    geom = state.geometry
    _, vacuum, _ = adapter.vendor_modules()
    total = assemble_total_stress_energy(
        grid, geom, state.scalars, state.matter
    )
    raw = vacuum.constraints(grid, geom)

    radii = np.asarray(grid.centers, dtype=float)
    volumes = np.asarray(grid.volumes, dtype=float)
    raw_h = np.asarray(raw["hamiltonian"], dtype=float)
    raw_m = np.asarray(raw["momentum"], dtype=float)
    rho = np.asarray(total.rho, dtype=float)
    current = np.asarray(total.j, dtype=float)

    # Same physical residual definitions used by production_kernel diagnostics.
    H = raw_h - 16.0 * math.pi * rho
    M = raw_m - 8.0 * math.pi * current

    Ab = -0.5 * np.asarray(geom.Aa, dtype=float)
    geometry = vacuum.geometry_terms(grid, geom)
    curvature = np.asarray(geometry["R"], dtype=float)
    extrinsic_A = -(np.asarray(geom.Aa, dtype=float) ** 2 + 2.0 * Ab ** 2)
    extrinsic_K = (2.0 / 3.0) * np.asarray(geom.K, dtype=float) ** 2
    matter_source = -16.0 * math.pi * rho
    reconstructed_H = curvature + extrinsic_A + extrinsic_K + matter_source
    decomposition_error = reconstructed_H - H

    # Disjoint radial partitions; global maximum remains independently reported.
    masks = {
        "full_grid": np.ones(grid.n, dtype=bool),
        "inner_r_le_20": radii <= 20.0,
        "middle_20_lt_r_lt_64": (radii > 20.0) & (radii < 64.0),
        "outer_r_ge_64": radii >= 64.0,
    }
    center_cells = list(range(min(5, grid.n)))
    profile = []
    for i in range(grid.n):
        profile.append({
            "cell": int(i),
            "r": float(radii[i]),
            "raw_hamiltonian": float(raw_h[i]),
            "rho_total": float(rho[i]),
            "curvature_R": float(curvature[i]),
            "extrinsic_A": float(extrinsic_A[i]),
            "extrinsic_K": float(extrinsic_K[i]),
            "matter_source": float(matter_source[i]),
            "hamiltonian_residual_H": float(H[i]),
            "hamiltonian_decomposition_error": float(decomposition_error[i]),
            "raw_momentum": float(raw_m[i]),
            "momentum_density_j": float(current[i]),
            "momentum_residual_M": float(M[i]),
        })

    return {
        "sample_index": int(sample_index),
        "expected_target_time": float(expected_target),
        "actual_t": float(state.t),
        "tau": float(state.tau),
        "grid_cell_count": int(grid.n),
        "center_cells_0_to_4": [profile[i] for i in center_cells],
        "regional_summary": {
            name: {
                "H": _summary(H, radii, volumes, mask),
                "M": _summary(M, radii, volumes, mask),
            }
            for name, mask in masks.items()
        },
        "decomposition_error_max_abs_full_grid": float(np.max(np.abs(decomposition_error))),
        "profile": profile,
    }


def run_case(*, repair: bool = True) -> dict:
    Kernel = V55ProductionKernel if repair else NoRepairKernel
    kernel = Kernel()
    state = kernel.initialize(
        resolution=RESOLUTION,
        r_max=R_MAX,
        D_amplitude=D_AMPLITUDE,
        include_radiation=True,
    )
    dt_nominal = CFL * state.grid.dr
    samples = [record_sample(state, 0, 0.0)]
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
            samples.append(record_sample(
                state, target_index + 1, SAMPLE_TIMES[target_index]
            ))
            target_index += 1

    return {
        "status": "numerical_failure" if failure else "completed",
        "failure": failure,
        "diagnostic_only": True,
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
    samples = report.get("samples", [])
    expected_times = [0.0, *SAMPLE_TIMES]
    schedule_ok = len(samples) == len(expected_times)
    if schedule_ok:
        schedule_ok = all(
            abs(float(row["actual_t"]) - target) <= TIME_TOLERANCE
            for row, target in zip(samples, expected_times)
        )
    full_profile_ok = bool(samples) and all(
        row.get("grid_cell_count") == RESOLUTION
        and len(row.get("profile", [])) == RESOLUTION
        and [cell.get("cell") for cell in row.get("profile", [])] == list(range(RESOLUTION))
        for row in samples
    )
    finite_ok = True
    for row in samples:
        for cell in row.get("profile", []):
            for key, value in cell.items():
                if key != "cell" and not math.isfinite(float(value)):
                    finite_ok = False
                    break
            if not finite_ok:
                break
        if not finite_ok:
            break
    decomposition_ok = bool(samples) and all(
        math.isfinite(float(row["decomposition_error_max_abs_full_grid"]))
        and float(row["decomposition_error_max_abs_full_grid"]) <= DECOMPOSITION_TOLERANCE
        for row in samples
    )
    config_ok = (
        report.get("D_amplitude") == D_AMPLITUDE
        and report.get("resolution") == RESOLUTION
        and report.get("r_max") == R_MAX
        and report.get("cfl") == CFL
        and report.get("include_radiation") is True
        and report.get("requested_final_time") == FINAL_TIME
        and abs(float(report.get("final_time", float("nan"))) - FINAL_TIME) <= 1.0e-10
    )
    checks = {
        "completed_without_numerical_failure": (
            report.get("status") == "completed" and report.get("failure") is None
        ),
        "frozen_configuration_matches": config_ok,
        "all_registered_samples_present_within_tolerance": schedule_ok,
        "full_radial_profile_has_all_cells": full_profile_ok,
        "all_profile_values_finite": finite_ok,
        "full_grid_hamiltonian_decomposition_le_1e-12": decomposition_ok,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "tolerance_notes": {
            "sample_time_absolute_tolerance": TIME_TOLERANCE,
            "full_grid_decomposition_error_max_abs": DECOMPOSITION_TOLERANCE,
            "no_physical_residual_pass_threshold_is_applied": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repair", choices=("on", "off"), default="on")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    report = run_case(repair=args.repair == "on")
    report["admission"] = admission_report(report)
    default_output = f"runs/outer-constraint-audit/repair-{args.repair}.json"
    output = Path(args.output or default_output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["admission"], indent=2))
    if not report["admission"]["passed"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
