#!/usr/bin/env python3
"""Capture exact source-resolved Lambda/M_r operands from a copied accepted state.

This diagnostic advances exactly one in-memory PIRK copy from a checksummed
t=0.60 endpoint profile. The caller's saved state is immutable and no trajectory
is written. Full radial state/term arrays are emitted for stage inspection.
No matter-source coefficient is inferred or fitted; scientific admission stays
false.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.matter_system import ConservedSpecies
from engine.production_kernel import LAMBDA_M, ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState, total_matter_projection

SOURCE_COMMIT = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_SCHEMA = "origin-cell-joint-shadow-evolved/v1"
EXPECTED_TIME = 0.60
R_MAX = 80.0
CFL = 0.03
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
MATTER_MAP = {
    "dark_matter": "dm",
    "baryons": "baryon",
    "radiation": "radiation",
}
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")
STAGE_NAMES = ("start_lapse_refreshed", "explicit_predictor", "primary_predictor")
L2_CALL_NAMES = (*STAGE_NAMES, "final_primary_L2")
EPS = float(np.finfo(np.float64).eps)


def _arr(value: Any) -> np.ndarray:
    return np.asarray(value, dtype=np.float64)


def _max_abs(value: Any) -> float:
    values = _arr(value)
    return float(np.max(np.abs(values))) if values.size else 0.0


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _state_digest(state: ProductionState) -> str:
    arrays: dict[str, np.ndarray] = {}
    for name in GEOMETRY_FIELDS:
        arrays["geometry." + name] = getattr(state.geometry, name)
    for name in SCALAR_FIELDS:
        arrays["scalars." + name] = getattr(state.scalars, name)
    for sector in (*MATTER_MAP.keys(), "charged_dust"):
        item = getattr(state.matter, sector, None)
        if item is None:
            continue
        for name in SPECIES_FIELDS:
            arrays["matter." + sector + "." + name] = getattr(item, name)
    digest = hashlib.sha256()
    for name in sorted(arrays):
        values = np.ascontiguousarray(np.asarray(arrays[name]))
        digest.update(name.encode("utf-8"))
        digest.update(str(values.dtype).encode("ascii"))
        digest.update(repr(values.shape).encode("ascii"))
        digest.update(values.tobytes())
    if state.electromagnetic is not None:
        for name in ("enclosed_charge_faces", "electric_field_faces"):
            if hasattr(state.electromagnetic, name):
                values = np.ascontiguousarray(np.asarray(getattr(state.electromagnetic, name)))
                digest.update(("electromagnetic." + name).encode("utf-8"))
                digest.update(str(values.dtype).encode("ascii"))
                digest.update(repr(values.shape).encode("ascii"))
                digest.update(values.tobytes())
    return digest.hexdigest()


def _copy_tree(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _copy_tree(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_copy_tree(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.copy()
    if isinstance(value, (float, int, np.floating, np.integer)):
        return float(value)
    return value


def _verify_source_pins(report: dict) -> dict:
    pins = report.get("source_pins", {}).get("source_pin_check", {})
    if pins.get("source_commit") != SOURCE_COMMIT:
        raise RuntimeError(f"accepted source anchor mismatch: {pins.get('source_commit')!r}")
    if pins.get("vendor_commit") != VENDOR_COMMIT or pins.get("vendor_checkout") != VENDOR_COMMIT:
        raise RuntimeError("artifact vendor pin is not the required exact gitlink")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("artifact does not attest a clean protected-path diff")
    blobs = pins.get("protected_source_blobs", {})
    if "engine/production_kernel.py" not in blobs:
        raise RuntimeError("artifact is missing the pinned production source map")
    for path, row in blobs.items():
        expected = row.get("source_blob")
        runtime = row.get("runtime_blob")
        if not expected or runtime != expected:
            raise RuntimeError(f"source/runtime blob mismatch in artifact for {path}")
        actual = subprocess.run(
            ["git", "rev-parse", f"HEAD:{path}"],
            cwd=ROOT, check=True, capture_output=True, text=True,
        ).stdout.strip()
        if actual != expected:
            raise RuntimeError(
                f"checked-out protected blob drift for {path}: expected {expected}, got {actual}"
            )
    vendor_actual = subprocess.run(
        ["git", "-C", str(ROOT / "vendor" / "bb-palatini-unified-r0"), "rev-parse", "HEAD"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout.strip()
    if vendor_actual != VENDOR_COMMIT:
        raise RuntimeError(f"vendor checkout drift: expected {VENDOR_COMMIT}, got {vendor_actual}")
    return {
        "source_commit": SOURCE_COMMIT,
        "vendor_commit": VENDOR_COMMIT,
        "protected_source_blobs_checked": len(blobs),
        "protected_blobs_match": True,
    }


def _verify_bundle(input_dir: Path, resolution: int) -> tuple[dict, dict[str, np.ndarray], dict]:
    input_dir = input_dir.resolve()
    required = ("report.json", "profiles.npz", "SHA256SUMS.txt")
    for name in required:
        if not (input_dir / name).is_file():
            raise FileNotFoundError(f"input artifact missing {name}: {input_dir}")
    expected_sums = {}
    for line in (input_dir / "SHA256SUMS.txt").read_text(encoding="utf-8").splitlines():
        parts = line.strip().split(maxsplit=1)
        if not parts:
            continue
        if len(parts) != 2:
            raise ValueError(f"malformed SHA256SUMS.txt line: {line!r}")
        expected_sums[parts[1].lstrip("*")] = parts[0]
    for name in ("report.json", "profiles.npz"):
        if name not in expected_sums:
            raise RuntimeError(f"checksum manifest does not cover {name}")
        actual = _file_sha256(input_dir / name)
        if actual != expected_sums[name]:
            raise RuntimeError(f"{name} SHA-256 mismatch: expected {expected_sums[name]}, got {actual}")

    report = json.loads((input_dir / "report.json").read_text(encoding="utf-8"))
    if report.get("schema") != EXPECTED_SCHEMA or report.get("status") != "completed_diagnostic_only":
        raise RuntimeError("input report schema/status does not match the accepted profile artifact")
    if report.get("scientific_admission") is not False:
        raise RuntimeError("input artifact has an unexpected scientific-admission state")
    provenance = _verify_source_pins(report)
    config = report.get("configuration", {})
    if int(config.get("resolution", -1)) != resolution:
        raise RuntimeError("report resolution does not match the requested artifact lane")
    if not math.isclose(float(config.get("r_max", math.nan)), R_MAX, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError("report radial extent mismatch")
    if not math.isclose(float(config.get("target_time", math.nan)), EXPECTED_TIME, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError("report target time is not 0.60")
    runs = report.get("runs", [])
    if len(runs) != 1 or int(runs[0].get("resolution", -1)) != resolution:
        raise RuntimeError("report must contain exactly the expected single-resolution run")
    if not math.isclose(float(runs[0].get("final_time", math.nan)), EXPECTED_TIME, rel_tol=0.0, abs_tol=1e-10):
        raise RuntimeError("saved endpoint time is not 0.60")
    if runs[0].get("all_final_arrays_finite") is not True:
        raise RuntimeError("input report does not attest finite final arrays")
    prefix = f"N{resolution}_full_continuous_metric_profile_final_"
    with np.load(input_dir / "profiles.npz", allow_pickle=False) as loaded:
        arrays = {k: np.asarray(loaded[k]).copy() for k in loaded.files}
    required_keys = (
        *(prefix + name for name in (*GEOMETRY_FIELDS, *SCALAR_FIELDS)),
        *(prefix + f"{species}_{field}" for species in MATTER_MAP.values() for field in SPECIES_FIELDS),
        prefix + "r", prefix + "t", prefix + "tau", prefix + "e_folds",
        prefix + "profile_C_lambda", prefix + "profile_Lambda_metric",
    )
    missing = [key for key in required_keys if key not in arrays]
    if missing:
        raise RuntimeError(f"profile NPZ missing required arrays: {missing}")
    for key in required_keys:
        if not np.all(np.isfinite(_arr(arrays[key]))):
            raise FloatingPointError(f"profile array contains non-finite values: {key}")
    return report, arrays, {
        "input_dir": str(input_dir),
        "sha256": expected_sums,
        "source_pin_check": provenance,
        "report_sha256": _file_sha256(input_dir / "report.json"),
        "profiles_sha256": _file_sha256(input_dir / "profiles.npz"),
    }


def _restore_profile_state(
    profile_arrays: dict[str, np.ndarray], resolution: int
) -> ProductionState:
    prefix = f"N{resolution}_full_continuous_metric_profile_final_"
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(resolution, R_MAX)
    geometry_values = {}
    for name in GEOMETRY_FIELDS:
        values = _arr(profile_arrays[prefix + name])
        if values.shape != grid.centers.shape or not np.all(np.isfinite(values)):
            raise ValueError(f"saved geometry field {name} has an unexpected shape or non-finite value")
        geometry_values[name] = values.copy()
    geometry = vacuum.VacuumState(**geometry_values)
    scalars = ScalarFields(*(
        _arr(profile_arrays[prefix + name]).copy() for name in SCALAR_FIELDS
    ))
    species = {}
    for sector, source_name in MATTER_MAP.items():
        species[sector] = ConservedSpecies(*(
            _arr(profile_arrays[prefix + f"{source_name}_{field}"]).copy()
            for field in SPECIES_FIELDS
        ))
    matter = V55MatterState(
        dark_matter=species["dark_matter"],
        baryons=species["baryons"],
        radiation=species["radiation"],
        charged_dust=None,
    )
    state = ProductionState(
        grid=grid,
        geometry=geometry,
        scalars=scalars,
        matter=matter,
        t=float(np.asarray(profile_arrays[prefix + "t"]).item()),
        tau=float(np.asarray(profile_arrays[prefix + "tau"]).item()),
        e_folds=float(np.asarray(profile_arrays[prefix + "e_folds"]).item()),
        electromagnetic=None,
    )
    saved_r = _arr(profile_arrays[prefix + "r"])
    if saved_r.shape != state.grid.centers.shape or not np.allclose(
        saved_r, _arr(state.grid.centers), rtol=0.0, atol=1e-12
    ):
        raise RuntimeError("saved profile radii do not match the pinned grid")
    if not math.isclose(state.t, EXPECTED_TIME, rel_tol=0.0, abs_tol=1e-10):
        raise RuntimeError("restored state time does not match t=0.60")
    return state


def _lambda_metric(grid: Any, geometry: Any) -> np.ndarray:
    r = _arr(grid.centers)
    a, b = _arr(geometry.a), _arr(geometry.b)
    d1 = lambda values, parity=1: _arr(grid.cell_derivative_fourth(values, parity=parity))
    return d1(a, 1) / (2.0 * a * a) - d1(b, 1) / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)


def _put_state_arrays(out: dict[str, np.ndarray], label: str, state: ProductionState) -> None:
    for name in GEOMETRY_FIELDS:
        out[f"{label}__geometry__{name}"] = _arr(getattr(state.geometry, name)).copy()
    for name in SCALAR_FIELDS:
        out[f"{label}__scalars__{name}"] = _arr(getattr(state.scalars, name)).copy()
    for sector, _ in MATTER_MAP.items():
        item = getattr(state.matter, sector)
        for name in SPECIES_FIELDS:
            out[f"{label}__matter__{sector}__{name}"] = _arr(getattr(item, name)).copy()


def _put_stage_arrays(
    out: dict[str, np.ndarray],
    label: str,
    geometry: Any,
    terms: dict,
    scalars: ScalarFields | None = None,
    matter: V55MatterState | None = None,
) -> None:
    for name in GEOMETRY_FIELDS:
        out[f"{label}__geometry__{name}"] = _arr(getattr(geometry, name)).copy()
    for section, values in terms.items():
        if isinstance(values, dict):
            for name, value in values.items():
                if isinstance(value, np.ndarray):
                    out[f"{label}__term__{section}__{name}"] = _arr(value).copy()
        elif isinstance(values, np.ndarray):
            out[f"{label}__term__{section}"] = _arr(values).copy()
    if scalars is not None:
        for name in SCALAR_FIELDS:
            out[f"{label}__scalars__{name}"] = _arr(getattr(scalars, name)).copy()
    if matter is not None:
        for sector in MATTER_MAP:
            item = getattr(matter, sector)
            for name in SPECIES_FIELDS:
                out[f"{label}__matter__{sector}__{name}"] = _arr(getattr(item, name)).copy()


def _region_summary(values: Any, r: np.ndarray) -> dict:
    values = _arr(values)
    mask10 = r <= 10.0 + 1e-12
    annulus = (r >= 5.0) & (r <= 10.0)
    return {
        "max_abs_full_grid": _max_abs(values),
        "max_abs_r_le_10": _max_abs(values[mask10]),
        "max_abs_fixed_annulus_r5_to10": _max_abs(values[annulus]) if np.any(annulus) else None,
        "signed_cells_0_1_2_6": {
            str(i): float(values[i]) for i in (0, 1, 2, 6) if i < values.size
        },
    }


def _raw_bound(*operands: np.ndarray, factor: float = 512.0) -> np.ndarray:
    scale = np.zeros_like(_arr(operands[0]))
    for values in operands:
        scale += np.abs(_arr(values))
    return factor * EPS * np.maximum(scale, np.finfo(np.float64).tiny)


def _capture_step(state: ProductionState, resolution: int, dt: float, out_dir: Path, report_bundle: dict) -> dict:
    if state.electromagnetic is not None or state.matter.charged_dust is not None:
        raise RuntimeError("this bounded capture is defined only for the archived EM-free six-sector state")
    kernel = V55ProductionKernel()
    _, vacuum, _ = adapter.vendor_modules()
    input_digest = _state_digest(state)
    stage_capture: dict[str, dict] = {}
    lambda2_capture: dict[str, dict] = {}
    stage_index = 0
    lambda2_index = 0

    original_geometry_stage = kernel._geometry_stage_terms
    original_lambda2 = vacuum.lambda_l2_rhs

    def record_geometry_stage(grid, geometry, scalars, matter, electromagnetic):
        nonlocal stage_index
        index = stage_index
        stage_index += 1
        if index >= len(STAGE_NAMES):
            raise RuntimeError("unexpected additional geometry stage in one-step capture")
        name = STAGE_NAMES[index]
        terms = original_geometry_stage(grid, geometry, scalars, matter, electromagnetic)
        stage_capture[name] = {
            "geometry": geometry.copy(),
            "scalars": scalars.copy(),
            "matter": matter.copy(),
            "terms": _copy_tree(terms),
        }
        return terms

    def record_lambda2(grid, geometry, *, lambda_m=LAMBDA_M):
        nonlocal lambda2_index
        index = lambda2_index
        lambda2_index += 1
        if index >= len(L2_CALL_NAMES):
            raise RuntimeError("unexpected additional Lambda L2 call")
        name = L2_CALL_NAMES[index]
        result = original_lambda2(grid, geometry, lambda_m=lambda_m)
        lambda2_capture[name] = {
            "geometry": geometry.copy(),
            "result": _arr(result).copy(),
            "lambda_m": float(lambda_m),
        }
        return result

    kernel._geometry_stage_terms = record_geometry_stage
    vacuum.lambda_l2_rhs = record_lambda2
    try:
        output = kernel.step(state, dt)
    finally:
        kernel._geometry_stage_terms = original_geometry_stage
        vacuum.lambda_l2_rhs = original_lambda2

    if stage_index != 3 or set(stage_capture) != set(STAGE_NAMES):
        raise RuntimeError(f"incomplete geometry stage capture: count={stage_index}, stages={sorted(stage_capture)}")
    if lambda2_index != 4 or set(lambda2_capture) != set(L2_CALL_NAMES):
        raise RuntimeError(f"incomplete Lambda L2 capture: count={lambda2_index}, calls={sorted(lambda2_capture)}")
    output_digest = _state_digest(state)
    if output_digest != input_digest:
        raise RuntimeError("input accepted-state arrays mutated during copied one-step capture")
    if output.t <= state.t:
        raise RuntimeError("copied-step output time did not advance as expected")

    r = _arr(state.grid.centers)
    raw: dict[str, np.ndarray] = {"r": r.copy()}
    _put_state_arrays(raw, "input", state)
    _put_state_arrays(raw, "copied_output", output)
    for name, entry in stage_capture.items():
        _put_stage_arrays(
            raw, "stage_" + name, entry["geometry"], entry["terms"],
            entry["scalars"], entry["matter"],
        )
    for name, entry in lambda2_capture.items():
        _put_stage_arrays(
            raw, "lambda2_" + name, entry["geometry"],
            {"lambda_l2": entry["result"]},
        )

    def stage_source(name: str) -> dict:
        entry = stage_capture[name]
        geometry, scalars, matter = entry["geometry"], entry["scalars"], entry["matter"]
        total = total_matter_projection(state.grid, geometry, scalars, matter)
        j = _arr(total["j"])
        m_geom = _arr(vacuum.momentum_constraint(state.grid, geometry))
        alpha_over_a = _arr(geometry.alpha) / _arr(geometry.a)
        l3_vacuum = _arr(vacuum.lambda_l3_rhs(state.grid, geometry))
        l3_matter = -8.0 * math.pi * LAMBDA_M * alpha_over_a * j
        l3_total = _arr(entry["terms"]["lambda_l3"])
        closure = l3_total - (l3_vacuum + l3_matter)
        bound = _raw_bound(l3_total, l3_vacuum, l3_matter)
        if not np.all(np.abs(closure) <= bound):
            raise FloatingPointError(f"{name} Lambda L3 (vacuum + matter) reconstruction exceeded the raw-operand bound")
        raw[f"source_{name}__j_total"] = j.copy()
        raw[f"source_{name}__M_r_geometry"] = m_geom.copy()
        raw[f"source_{name}__alpha_over_a"] = alpha_over_a.copy()
        raw[f"source_{name}__lambda_l3_vacuum"] = l3_vacuum.copy()
        raw[f"source_{name}__lambda_l3_matter"] = l3_matter.copy()
        raw[f"source_{name}__lambda_l3_total"] = l3_total.copy()
        raw[f"source_{name}__lambda_l3_reconstruction_residual"] = closure.copy()
        return {
            "j_total": j,
            "M_r_geometry": m_geom,
            "alpha_over_a": alpha_over_a,
            "lambda_l3_vacuum": l3_vacuum,
            "lambda_l3_matter": l3_matter,
            "lambda_l3_total": l3_total,
            "lambda_l3_reconstruction_residual": closure,
            "lambda_l3_raw_operand_bound": bound,
        }

    sources = {
        "start_lapse_refreshed": stage_source("start_lapse_refreshed"),
        "primary_predictor": stage_source("primary_predictor"),
    }
    final_entry = lambda2_capture["final_primary_L2"]
    geom_f = final_entry["geometry"]
    m_geom_f = _arr(vacuum.momentum_constraint(state.grid, geom_f))
    alpha_over_a_f = _arr(geom_f.alpha) / _arr(geom_f.a)
    l2_f = _arr(final_entry["result"])
    raw["source_final_primary_L2__M_r_geometry"] = m_geom_f.copy()
    raw["source_final_primary_L2__alpha_over_a"] = alpha_over_a_f.copy()
    raw["source_final_primary_L2__lambda_l2"] = l2_f.copy()

    stage0 = stage_capture["start_lapse_refreshed"]
    stagep = stage_capture["primary_predictor"]
    l2_0 = _arr(stage0["terms"]["lambda_l2"])
    l3_0 = _arr(stage0["terms"]["lambda_l3"])
    l3_p = _arr(stagep["terms"]["lambda_l3"])
    l2_calls_match = {}
    for name in STAGE_NAMES:
        stage_values = _arr(stage_capture[name]["terms"]["lambda_l2"])
        call_values = _arr(lambda2_capture[name]["result"])
        difference = stage_values - call_values
        l2_calls_match[name] = _max_abs(difference)
        if not np.array_equal(stage_values, call_values):
            if _max_abs(difference) > _max_abs(_raw_bound(stage_values, call_values)):
                raise FloatingPointError(f"recorded {name} Lambda L2 operand differs from its direct call")
    lambda_increment_predicted = 0.5 * dt * (l2_0 + l2_f + l3_0 + l3_p)
    lambda_increment_actual = _arr(output.geometry.Lambda) - _arr(state.geometry.Lambda)
    lambda_residual = lambda_increment_actual - lambda_increment_predicted
    lambda_bound = _raw_bound(
        _arr(output.geometry.Lambda), _arr(state.geometry.Lambda),
        0.5 * dt * _arr(np.abs(l2_0) + np.abs(l2_f) + np.abs(l3_0) + np.abs(l3_p)),
    )
    mask10 = r <= 10.0 + 1e-12
    if not np.all(np.abs(lambda_residual[mask10]) <= lambda_bound[mask10]):
        raise FloatingPointError("accepted copied-step Lambda update did not close within raw-operand-scaled roundoff for r<=10")

    # The native metric target is evaluated at the exact input and copied output;
    # no instantaneous chain-rate is substituted for this nonlinear finite difference.
    lambda_metric_input = _lambda_metric(state.grid, state.geometry)
    lambda_metric_output = _lambda_metric(state.grid, output.geometry)
    connection_input = _arr(state.geometry.Lambda) - lambda_metric_input
    connection_output = _arr(output.geometry.Lambda) - lambda_metric_output
    delta_connection_from_updates = (
        lambda_increment_actual - (lambda_metric_output - lambda_metric_input)
    )
    delta_connection_direct = connection_output - connection_input
    connection_residual = delta_connection_direct - delta_connection_from_updates
    connection_bound = _raw_bound(
        connection_input, connection_output, lambda_increment_actual,
        lambda_metric_output - lambda_metric_input,
    )
    if not np.all(np.abs(connection_residual) <= connection_bound):
        raise FloatingPointError("finite-step connection identity failed its arithmetic closure bound")

    source_m0 = LAMBDA_M * sources["start_lapse_refreshed"]["alpha_over_a"] * sources["start_lapse_refreshed"]["M_r_geometry"]
    source_mf = LAMBDA_M * alpha_over_a_f * m_geom_f
    source_j0 = sources["start_lapse_refreshed"]["lambda_l3_matter"]
    source_jp = sources["primary_predictor"]["lambda_l3_matter"]
    finite_step_source_combo = 0.5 * dt * (source_m0 + source_mf + source_j0 + source_jp)
    same_stage_start_pair = 0.5 * dt * (
        source_m0 + source_j0
        + LAMBDA_M * sources["primary_predictor"]["alpha_over_a"] * sources["primary_predictor"]["M_r_geometry"]
        + source_jp
    )
    raw["budget__lambda_increment_predicted"] = lambda_increment_predicted.copy()
    raw["budget__lambda_increment_actual"] = lambda_increment_actual.copy()
    raw["budget__lambda_increment_residual"] = lambda_residual.copy()
    raw["budget__lambda_increment_raw_bound"] = lambda_bound.copy()
    raw["budget__lambda_metric_input"] = lambda_metric_input.copy()
    raw["budget__lambda_metric_output"] = lambda_metric_output.copy()
    raw["budget__connection_input"] = connection_input.copy()
    raw["budget__connection_output"] = connection_output.copy()
    raw["budget__delta_connection_from_updates"] = delta_connection_from_updates.copy()
    raw["budget__delta_connection_direct"] = delta_connection_direct.copy()
    raw["budget__connection_residual"] = connection_residual.copy()
    raw["budget__connection_raw_bound"] = connection_bound.copy()
    raw["budget__source_combo_accepted_split"] = finite_step_source_combo.copy()
    raw["budget__source_combo_same_stage_counterfactual"] = same_stage_start_pair.copy()

    npz_path = out_dir / "stage_arrays.npz"
    np.savez_compressed(npz_path, **raw)
    source_report = {}
    for name, entry in sources.items():
        source_report[name] = {
            key: _region_summary(value, r)
            for key, value in entry.items()
            if isinstance(value, np.ndarray) and value.ndim == 1 and value.shape == r.shape
        }
    report = {
        "schema": "issue4-pirk-stage-resolved-momentum-source/v1",
        "status": "completed_diagnostic_only",
        "lane": "Issue #4 / full matter-coupled momentum-propagation identity",
        "scope": "one copied, in-memory PIRK step from the exact checksum-verified #420 t=0.60 endpoint profile; no accepted trajectory written",
        "input_artifact": {
            "profile_resolution": resolution,
            "input_time": float(state.t),
            "input_state_digest_recomputed": input_digest,
            "input_state_immutable": True,
            "source_report_sha256": report_bundle["report_sha256"],
            "profiles_npz_sha256": report_bundle["profiles_sha256"],
            "input_member_sha256": report_bundle["sha256"],
        },
        "pins": {
            **report_bundle["source_pin_check"],
            "source_pin_commit": SOURCE_COMMIT,
            "vendor_pin": VENDOR_COMMIT,
            "runtime_python": sys.version,
            "runtime_numpy": np.__version__,
            "runtime_platform": platform.platform(),
        },
        "step": {
            "dt": dt,
            "copied_output_time": float(output.t),
            "input_state_digest_after_copy": output_digest,
            "geometry_stage_names": list(STAGE_NAMES),
            "lambda_l2_call_names": list(L2_CALL_NAMES),
            "geometry_stage_count": len(stage_capture),
            "lambda_l2_call_count": len(lambda2_capture),
            "lambda_l2_stage_operand_max_differences": l2_calls_match,
        },
        "source_decomposition": {
            "lambda_m": LAMBDA_M,
            "matter_l3_definition": "-8*pi*lambda_m*(alpha/a)*j_total; diagnostic source term only",
            "M_r_geometry_definition": "pinned vendor momentum_constraint(grid, geometry), as consumed by vendor Lambda L2",
            "stages_with_independent_j_reconstruction": ["start_lapse_refreshed", "primary_predictor"],
            "stages_with_geometry_M_r": ["start_lapse_refreshed", "final_primary_L2"],
            "stage_arrays": source_report,
            "accepted_split_source_combo": _region_summary(finite_step_source_combo, r),
            "same_stage_pair_counterfactual": _region_summary(same_stage_start_pair, r),
            "final_primary_L2": {
                "M_r_geometry": _region_summary(m_geom_f, r),
                "alpha_over_a": _region_summary(alpha_over_a_f, r),
                "lambda_l2": _region_summary(l2_f, r),
            },
            "warning": "the accepted PIRK split samples final L2 at the final-primary mixed state and matter L3 at the primary predictor; their second-stage source operands are not interchangeable",
        },
        "closure": {
            "lambda_increment_update": _region_summary(lambda_residual, r),
            "lambda_increment_raw_operand_bound": _region_summary(lambda_bound, r),
            "lambda_increment_pass_r_le_10": True,
            "finite_step_connection_identity": _region_summary(connection_residual, r),
            "finite_step_connection_raw_operand_bound": _region_summary(connection_bound, r),
            "finite_step_connection_identity_pass_full_grid": True,
            "max_abs_lambda_l3_reconstruction_residual": max(
                _max_abs(sources[name]["lambda_l3_reconstruction_residual"])
                for name in sources
            ),
            "scientific_interpretation": "arithmetic/stage-accounting closure only; not the full M_r propagation identity or total stress-energy conservation",
        },
        "admission_flags": {
            "full_matter_coupled_propagation_identity_closed": False,
            "matter_source_sign_and_normalization_reviewed": False,
            "reference_metric_connection_coupling_convention_frozen": False,
            "total_stress_energy_sector_conservation_closed": False,
            "scientific_admission": False,
        },
        "limitations": [
            "No prior accepted interval ending at t=0.60 was replayed; this is one copied diagnostic step starting from the accepted t=0.60 state.",
            "No production equations, initializer, dependencies, vendor source, source pins, bounds, or scientific-admission flags are changed.",
            "The candidate source S_Q=-8*pi*alpha*Q_r remains provisional.",
            "Pointwise, finite-volume, and expanded stress divergences remain distinct.",
            "The same-stage source pairing is algebraic; differing PIRK stage states prevent promoting it to a finite-step conservation theorem.",
        ],
    }
    report["artifacts"] = {
        "stage_arrays.npz": _file_sha256(npz_path),
    }
    (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    sums = {
        "report.json": _file_sha256(out_dir / "report.json"),
        "stage_arrays.npz": _file_sha256(npz_path),
    }
    (out_dir / "SHA256SUMS.txt").write_text(
        "".join(f"{digest}  {name}\n" for name, digest in sums.items()),
        encoding="utf-8",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, help="directory with report.json, profiles.npz, and SHA256SUMS.txt from #420")
    parser.add_argument("--resolution", choices=(60, 120, 240), type=int, required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report, arrays, report_bundle = _verify_bundle(input_dir, args.resolution)
    state = _restore_profile_state(arrays, args.resolution)
    dt = CFL * R_MAX / args.resolution
    original_dt = float(report["runs"][0]["dt"])
    if not math.isclose(dt, original_dt, rel_tol=0.0, abs_tol=1e-12):
        # The copied step follows the same pinned CFL/dr rule, but records the
        # new dt explicitly rather than confusing it with the last historical step.
        if not math.isclose(original_dt, CFL * R_MAX / args.resolution, rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError("input trajectory timestep does not match the pinned CFL/dr rule")
    result = _capture_step(state, args.resolution, dt, out_dir, report_bundle)
    print(json.dumps({
        "resolution": args.resolution,
        "input_time": result["input_artifact"]["input_time"],
        "copied_output_time": result["step"]["copied_output_time"],
        "status": result["status"],
        "input_state_immutable": result["input_artifact"]["input_state_immutable"],
        "lambda_increment_pass_r_le_10": result["closure"]["lambda_increment_pass_r_le_10"],
        "connection_identity_pass_full_grid": result["closure"]["finite_step_connection_identity_pass_full_grid"],
        "scientific_admission": result["admission_flags"]["scientific_admission"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
