#!/usr/bin/env python3
"""Frozen-state audit of the CMC lapse target and production K-RHS consistency.

This diagnostic never advances the state. It repeatedly applies the current-main
CMC lapse solve to the same archived physical slice, while separating:
  * the second-order elliptic rows used by solve_cmc_lapse,
  * the vendor fourth-order lapse Laplacian in primary_l2_rhs,
  * the K-sector KO dissipation included by primary_l3_rhs,
  * the difference between the incoming-alpha target and the solved-alpha RHS.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from engine.cmc_gauge import solve_cmc_lapse
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import total_matter_projection
from engine import v55_pirk_adapter as adapter


EXPECTED_CAPTURE_HEAD = "a854c3d6eafde1a1bae8066b7cc212bac7511ffd"
EXPECTED_CAPTURE_RUN_ID = 38337821862
EXPECTED_VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
GEOMETRY_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")
SPECIES_NAMES = ("dark_matter", "baryons", "radiation")
PHYSICAL_GEOMETRY_FIELDS = tuple(x for x in GEOMETRY_FIELDS if x != "alpha")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if values.shape != weights.shape:
        raise ValueError("weighted_mean values and weights must have equal shapes")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(weights)):
        raise ValueError("weighted_mean inputs must be finite")
    total_weight = float(np.sum(weights))
    if not np.isfinite(total_weight) or total_weight <= 0.0:
        raise ValueError("weighted_mean requires a positive finite total weight")
    return float(np.sum(weights * values) / total_weight)


def profile_stats(values: np.ndarray, weights: np.ndarray, start: int = 3,
                  stop: int | None = None) -> dict:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if values.shape != weights.shape:
        raise ValueError("profile and weight arrays must have equal shapes")
    if stop is None:
        stop = values.size - 1
    if not (0 <= start < stop <= values.size):
        raise ValueError(f"invalid profile statistics slice [{start}:{stop}]")
    v = values[start:stop]
    w = weights[start:stop]
    if not np.all(np.isfinite(v)) or not np.all(np.isfinite(w)):
        raise ValueError("profile statistics require finite values and weights")
    return {
        "cell_start_inclusive": int(start),
        "cell_stop_exclusive": int(stop),
        "count": int(v.size),
        "max_abs": float(np.max(np.abs(v))),
        "max_signed": float(v[int(np.argmax(np.abs(v)))]),
        "rms": float(np.sqrt(np.mean(v * v))),
        "proper_volume_weighted_rms": float(np.sqrt(np.sum(w * v * v) / np.sum(w))),
    }


def find_artifact_dir(root: Path) -> Path:
    matches = list(root.rglob("stage_capture_profiles.npz"))
    if len(matches) != 1:
        raise RuntimeError(
            f"expected exactly one stage_capture_profiles.npz under {root}; "
            f"found {len(matches)}"
        )
    directory = matches[0].parent
    required = (
        "stage_capture_profiles.npz", "stage_capture.json", "report.json",
        "source_boundary_failure_accepted_state.npz",
        "source_boundary_failure_accepted_state.json", "SHA256SUMS.txt",
    )
    for name in required:
        if not (directory / name).is_file():
            raise RuntimeError(f"verified replay artifact lacks {name}")
    return directory


def verify_artifact(directory: Path) -> dict:
    checked = {}
    manifest = directory / "SHA256SUMS.txt"
    for line in manifest.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        expected, name = line.split(maxsplit=1)
        name = name.strip().lstrip("*")
        path = directory / name
        if not path.is_file():
            raise RuntimeError(f"manifest input missing: {name}")
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"manifest checksum mismatch for {name}: {actual}")
        checked[name] = actual
    if set(checked) != {
        "report.json", "source_boundary_failure_accepted_state.npz",
        "source_boundary_failure_accepted_state.json", "stage_capture.json",
        "stage_capture_profiles.npz",
    }:
        raise RuntimeError(f"unexpected replay manifest entries: {sorted(checked)}")
    report = json.loads((directory / "report.json").read_text(encoding="utf-8"))
    capture = json.loads((directory / "stage_capture.json").read_text(encoding="utf-8"))
    if report.get("provenance", {}).get("candidate_head") != EXPECTED_CAPTURE_HEAD:
        raise RuntimeError("frozen input is not from the expected exact observer head")
    if report.get("provenance", {}).get("input_artifact_run") != 38322337250:
        raise RuntimeError("frozen-state source run provenance changed")
    if report.get("provenance", {}).get("vendor_pin") != EXPECTED_VENDOR_COMMIT:
        raise RuntimeError("frozen replay's vendor provenance changed")
    if capture.get("accepted_state_count") != 27:
        raise RuntimeError("expected 27 accepted states in the frozen profile capture")
    if capture.get("geometry_source_term_record_count") != 703:
        raise RuntimeError("unexpected captured geometry-source stage count")
    return {"checksums": checked, "report": report, "capture": capture}


def load_terminal_state(directory: Path, capture: dict):
    """Rebuild containers and overwrite them with saved accepted state 026."""
    npz_path = directory / "stage_capture_profiles.npz"
    history = capture["accepted_state_history"]
    if len(history) != 27:
        raise RuntimeError("accepted-state history length differs from its header")
    terminal = history[-1]
    if abs(float(terminal["t"]) - 47.15300258518951) > 2.0e-14:
        raise RuntimeError(f"unexpected terminal input time {terminal['t']:.17g}")

    with np.load(npz_path, allow_pickle=False) as data:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=80, r_max=80.0, amplitude=0.01, width=7.0,
            D_amplitude=1.0e-10, include_radiation=True,
        )
        if not np.array_equal(
            np.asarray(state.grid.centers), np.asarray(data["accepted_026__r"])
        ):
            raise RuntimeError("terminal accepted-state grid differs from N=80/r_max=80")
        for name in GEOMETRY_FIELDS:
            key = f"accepted_026__geometry__{name}"
            if key not in data:
                raise RuntimeError(f"terminal geometry profile absent: {key}")
            setattr(state.geometry, name, np.asarray(data[key], dtype=float).copy())
        state.scalars = ScalarFields(*(
            np.asarray(data[f"accepted_026__scalars__{name}"], dtype=float).copy()
            for name in SCALAR_FIELDS
        ))
        for species_name in SPECIES_NAMES:
            species = getattr(state.matter, species_name)
            for field in SPECIES_FIELDS:
                key = f"accepted_026__matter__{species_name}__{field}"
                if key not in data:
                    raise RuntimeError(f"terminal matter profile absent: {key}")
                setattr(species, field, np.asarray(data[key], dtype=float).copy())
    state.t = float(terminal["t"])
    state.tau = float(terminal["tau"])
    return kernel, state


def physical_state_digest(state) -> str:
    """Digest non-lapse geometry and all matter/scalar profiles."""
    digest = hashlib.sha256()
    arrays = [("r", state.grid.centers)]
    arrays.extend((f"geometry.{name}", getattr(state.geometry, name))
                  for name in PHYSICAL_GEOMETRY_FIELDS)
    arrays.extend((f"scalars.{name}", getattr(state.scalars, name))
                  for name in SCALAR_FIELDS)
    for species_name in SPECIES_NAMES:
        species = getattr(state.matter, species_name)
        arrays.extend((f"matter.{species_name}.{name}", getattr(species, name))
                      for name in SPECIES_FIELDS)
    for name, values in arrays:
        array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(str(array.shape).encode("ascii") + b"\0")
        digest.update(array.tobytes())
    return digest.hexdigest()


def proper_volume_weights(grid, geometry) -> np.ndarray:
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    if np.any(a <= 0.0) or np.any(b <= 0.0) or np.any(X <= 0.0):
        raise RuntimeError("non-positive geometry in the CMC audit state")
    weights = r**2 * np.sqrt(np.maximum(a, 0.0)) * b / X**3
    weights *= float(grid.dr)
    return weights


def full_k_rhs(grid, geometry, scalars, matter, vacuum):
    l2 = vacuum.primary_l2_rhs(grid, geometry)
    l3_vacuum = vacuum.primary_l3_rhs(grid, geometry)
    l3_total = adapter.primary_l3_with_matter(grid, geometry, scalars, matter)
    return l2["K"] + l3_total["K"], l2["K"], l3_vacuum["K"], l3_total["K"] - l3_vacuum["K"]


def pde_residual_and_decomposition(grid, geometry, scalars, matter, vacuum,
                                   target_kdot: float) -> dict:
    """Split solved-lapse full K-RHS residual into stencil, KO and solve-row parts."""
    alpha = np.asarray(geometry.alpha, dtype=float)
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    inv = X * X / a
    ap = grid.cell_derivative_fourth(a, parity=1)
    bp = grid.cell_derivative_fourth(b, parity=1)
    Xp = grid.cell_derivative_fourth(X, parity=1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r
    lap2 = np.zeros_like(alpha)
    lap2[1:-1] = (
        (alpha[:-2] - 2.0 * alpha[1:-1] + alpha[2:]) / grid.dr**2
        + c[1:-1] * (alpha[2:] - alpha[:-2]) / (2.0 * grid.dr)
    )

    total = total_matter_projection(grid, geometry, scalars, matter)
    Q = (
        1.5 * geometry.Aa**2 + geometry.K**2 / 3.0
        + 4.0 * math.pi * (total["rho"] + total["pr"] + 2.0 * total["pt"])
    )
    # This is exactly the interior equation assembled by solve_cmc_lapse:
    # lap2(alpha) - Q*alpha/inv = -target_kdot/inv.
    second_order_residual = inv * lap2 - Q * alpha + target_kdot

    actual, l2_k, _l3_vacuum, _l3_matter = full_k_rhs(
        grid, geometry, scalars, matter, vacuum
    )
    ko_k = vacuum.ko_dissipation(grid, geometry.K, 1)
    # The equation above implies (-lap2 + Q*alpha) == target_kdot up to
    # second_order_residual; actual K RHS differs by the vendor L2 stencil and
    # the KO(K) term in primary_l3_rhs.
    lapse_stencil_delta = l2_k + lap2
    actual_minus_target = actual - target_kdot
    decomposition = lapse_stencil_delta + ko_k - second_order_residual
    closure = actual_minus_target - decomposition
    weights = proper_volume_weights(grid, geometry)
    mask = slice(3, grid.n - 1)
    terms = {
        "second_order_solve_row_residual": profile_stats(second_order_residual, weights, 1, grid.n - 1),
        "vendor_lapse_stencil_minus_solve_stencil": profile_stats(lapse_stencil_delta, weights, 1, grid.n - 1),
        "vendor_K_ko_dissipation": profile_stats(ko_k, weights, 1, grid.n - 1),
        "production_K_rhs_minus_target": profile_stats(actual_minus_target, weights, 1, grid.n - 1),
        "residual_decomposition_closure": profile_stats(closure, weights, 1, grid.n - 1),
        "outer_interior_K_rhs_minus_target": profile_stats(actual_minus_target, weights, mask.start, mask.stop),
        "weighted_average_production_K_rhs": weighted_mean(actual, weights),
        "weighted_average_target_difference": weighted_mean(actual - target_kdot, weights),
        "target_kdot": float(target_kdot),
        "max_abs_K_ko_over_outer_interior": float(np.max(np.abs(ko_k[mask]))),
    }
    scale = (
        np.abs(actual) + abs(target_kdot) + np.abs(lapse_stencil_delta)
        + np.abs(ko_k) + np.abs(second_order_residual)
    )
    bound = 64.0 * np.finfo(np.float64).eps * scale
    bad = np.abs(closure[1:-1]) > bound[1:-1]
    if np.any(bad):
        terms["decomposition_roundoff_failures"] = int(np.count_nonzero(bad))
        terms["decomposition_roundoff_max_ratio"] = float(
            np.max(np.abs(closure[1:-1][bad]) / bound[1:-1][bad])
        )
        raise RuntimeError("K-RHS residual decomposition failed its operand-scaled bound")
    terms["decomposition_roundoff_failures"] = 0
    terms["decomposition_roundoff_max_ratio"] = float(
        np.max(np.abs(closure[1:-1]) / np.maximum(bound[1:-1], np.finfo(float).tiny))
    )
    terms["vendor_KO_fraction_of_outer_interior_residual_max"] = (
        terms["vendor_K_ko_dissipation"]["max_abs"]
        / max(terms["outer_interior_K_rhs_minus_target"]["max_abs"], np.finfo(float).tiny)
    )
    return {
        "terms": terms,
        "arrays": {
            "second_order_residual": second_order_residual,
            "lapse_stencil_delta": lapse_stencil_delta,
            "ko_K": ko_k,
            "actual_minus_target": actual_minus_target,
            "decomposition_closure": closure,
            "actual_K_rhs": actual,
        },
    }


def run_audit(input_root: Path, output_dir: Path, iterations: int) -> dict:
    directory = find_artifact_dir(input_root)
    verified = verify_artifact(directory)
    kernel, state = load_terminal_state(directory, verified["capture"])
    if kernel is None:
        raise RuntimeError("failed to initialize current-main kernel containers")
    original_digest = physical_state_digest(state)
    _, vacuum, _ = adapter.vendor_modules()
    weights = proper_volume_weights(state.grid, state.geometry)
    records = []
    profile_arrays = {}
    converged = False
    convergence_limit = 128.0 * np.finfo(np.float64).eps

    for index in range(iterations):
        input_alpha = np.asarray(state.geometry.alpha, dtype=float).copy()
        input_rhs, input_l2, input_l3_vac, input_l3_matter = full_k_rhs(
            state.grid, state.geometry, state.scalars, state.matter, vacuum
        )
        input_target = weighted_mean(input_rhs, weights)
        solved_alpha, target = solve_cmc_lapse(
            state.grid, state.geometry, state.scalars, state.matter
        )
        if not np.isfinite(target) or not np.all(np.isfinite(solved_alpha)):
            raise RuntimeError(f"non-finite CMC solve output at iteration {index}")
        if np.any(solved_alpha <= 0.0):
            raise RuntimeError(f"non-positive CMC lapse at iteration {index}")
        target_delta = abs(float(target) - input_target)
        if target_delta > 128.0 * np.finfo(float).eps * max(
            abs(float(target)), abs(input_target), np.finfo(float).tiny
        ):
            raise RuntimeError("target_kdot does not equal the incoming RHS weighted mean")

        state.geometry.alpha = np.asarray(solved_alpha, dtype=float).copy()
        output_alpha = state.geometry.alpha.copy()
        output_rhs = pde_residual_and_decomposition(
            state.grid, state.geometry, state.scalars, state.matter, vacuum, float(target)
        )
        delta_alpha = output_alpha - input_alpha
        alpha_scale = max(float(np.max(np.abs(input_alpha))),
                          float(np.max(np.abs(output_alpha))), np.finfo(float).tiny)
        delta_max = float(np.max(np.abs(delta_alpha)))
        target_record = output_rhs["terms"]
        target_record = dict(target_record)
        target_record.update({
            "iteration": int(index),
            "input_alpha_min": float(np.min(input_alpha)),
            "input_alpha_max": float(np.max(input_alpha)),
            "output_alpha_min": float(np.min(output_alpha)),
            "output_alpha_max": float(np.max(output_alpha)),
            "max_abs_alpha_change": delta_max,
            "alpha_change_relative_to_profile_scale": delta_max / alpha_scale,
            "incoming_K_rhs_target": float(target),
            "independent_incoming_weighted_mean": input_target,
            "incoming_target_roundoff_delta": target_delta,
            "incoming_lapse_K_rhs": profile_stats(input_l2, weights),
            "incoming_vacuum_L3_K_rhs": profile_stats(input_l3_vac, weights),
            "incoming_matter_L3_K_difference": profile_stats(input_l3_matter, weights),
        })
        records.append(target_record)
        profile_arrays[f"alpha_input_{index:03d}"] = input_alpha
        profile_arrays[f"alpha_output_{index:03d}"] = output_alpha
        profile_arrays[f"incoming_K_rhs_{index:03d}"] = input_rhs
        profile_arrays[f"incoming_lapse_K_rhs_{index:03d}"] = input_l2
        profile_arrays[f"incoming_vacuum_L3_K_rhs_{index:03d}"] = input_l3_vac
        profile_arrays[f"incoming_matter_L3_K_difference_{index:03d}"] = input_l3_matter
        profile_arrays[f"output_K_rhs_{index:03d}"] = output_rhs["arrays"]["actual_K_rhs"]
        for name, values in output_rhs["arrays"].items():
            profile_arrays[f"{name}_{index:03d}"] = values
        if physical_state_digest(state) != original_digest:
            raise RuntimeError("frozen physical state changed during lapse-only repetitions")
        if delta_max <= convergence_limit * alpha_scale:
            converged = True
            break
    output_dir.mkdir(parents=True, exist_ok=True)
    profile_path = output_dir / "cmc_fixed_point_profiles.npz"
    np.savez_compressed(profile_path, **profile_arrays)
    report = {
        "schema": "cmc_lapse_frozen_state_target_operator_audit_v1",
        "diagnostic_only": True,
        "status": "completed_diagnostic_only",
        "input_artifact": {
            "run_id": EXPECTED_CAPTURE_RUN_ID,
            "expected_capture_head": EXPECTED_CAPTURE_HEAD,
            "artifact_directory": str(directory),
            "manifest_sha256": verified["checksums"],
        },
        "state": {
            "resolution": int(state.grid.n),
            "r_max": float(state.grid.centers[-1] + 0.5 * state.grid.dr),
            "t": float(state.t),
            "tau": float(state.tau),
            "accepted_state_count": int(verified["capture"]["accepted_state_count"]),
            "geometry_stage_count": int(verified["capture"]["geometry_source_term_record_count"]),
            "physical_state_digest_before": original_digest,
            "physical_state_digest_after": physical_state_digest(state),
        },
        "operator_description": {
            "target": "proper-volume average of vendor primary_l2_rhs[K] + primary_l3_rhs[K] + L3 matter correction, using incoming alpha",
            "lapse_solver": "second-order radial interior operator; alpha[0]=alpha[1] center row; alpha[-1]=1 outer row",
            "production_K_rhs": "vendor fourth-order -lap(alpha) plus curvature/matter L3 sources and KO(K)",
            "residual_decomposition": "production K RHS - target = vendor-lapse-stencil minus solve-stencil + KO(K) - second-order solve-row residual",
            "no_evolution": True,
        },
        "fixed_point": {
            "iterations_run": len(records),
            "iteration_limit": iterations,
            "converged_to_roundoff": converged,
            "roundoff_convergence_limit_factor": 128,
            "final_max_abs_alpha_change": records[-1]["max_abs_alpha_change"] if records else None,
            "iteration_records": records,
        },
        "interpretation_limits": [
            "The physical slice is frozen; only the lapse is replaced by repeated calls to the current-main CMC solver.",
            "This does not advance the trajectory or establish that any source term, boundary, gauge condition, or continuum equation is physically wrong.",
            "The decomposition separates production-vs-solve lapse stencils, K-sector KO dissipation, and solve-row residual; it is an operator audit, not a promotion or scientific-admission result.",
        ],
    }
    report_path = output_dir / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    files = [profile_path, report_path]
    (output_dir / "SHA256SUMS.txt").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in files),
        encoding="utf-8",
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=12)
    args = parser.parse_args()
    if args.iterations < 1 or args.iterations > 32:
        raise ValueError("iterations must be in [1, 32]")
    report = run_audit(args.input_dir, args.output_dir, args.iterations)
    print(json.dumps({
        "status": report["status"],
        "t": report["state"]["t"],
        "iterations_run": report["fixed_point"]["iterations_run"],
        "converged_to_roundoff": report["fixed_point"]["converged_to_roundoff"],
        "final_max_abs_alpha_change": report["fixed_point"]["final_max_abs_alpha_change"],
        "final_K_rhs_minus_target": report["fixed_point"]["iteration_records"][-1]["production_K_rhs_minus_target"],
        "artifact_output": str(args.output_dir),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
