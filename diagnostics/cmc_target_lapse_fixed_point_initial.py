"""Initial-slice self-consistent CMC target/lapse map audit.

Snapshot-only diagnostic on the archived t=0 R80 and R160 states. It iterates
alpha <- solve_CMC(target_kdot(alpha)) with the spatial geometry, scalars and
matter frozen. There is no damping, clipping, fitted parameter or production
source change. Non-convergence is recorded as a valid diagnostic outcome.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import engine.cmc_gauge as cmc
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState

OUT = ROOT / "runs" / "cmc_target_lapse_fixed_point_initial"
R80_N, R80_MAX = 160, 80.0
R160_N, R160_MAX = 320, 160.0
AMPLITUDE, WIDTH, D_AMPLITUDE = 0.01, 7.0, 1.0e-4
INCLUDE_RADIATION = True
MAX_ITERATIONS = 20
STOP_RMS = 1.0e-13


def rms(values):
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_initial_row(path: Path, expected_n: int, expected_rmax: float):
    found = None
    with path.open("r", encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip():
                continue
            row = json.loads(line)
            if abs(float(row.get("t", -1.0))) <= 1.0e-12:
                if row.get("state_snapshot") is not None:
                    found = row
                    break
    if found is None:
        raise ValueError(f"{path}: no full t=0 state_snapshot found")
    snap = found["state_snapshot"]
    if expected_n == R80_N:
        n = len(snap["geometry"]["a"])
        radius = np.asarray(found["radius"], dtype=float)
    else:
        n = len(snap["geometry"]["a"])
        radius = np.asarray(found["radius_full"], dtype=float)
    if n != expected_n:
        raise ValueError(f"{path}: snapshot N={n}, expected {expected_n}")
    domain = found.get("domain", {})
    if expected_n == R80_N:
        full_n = int(domain.get("N", domain.get("full_N", expected_n)))
        full_rmax = float(domain.get("r_max", domain.get("full_r_max", expected_rmax)))
    else:
        full_n = int(found.get("full_N", domain.get("N", expected_n)))
        full_rmax = float(found.get("full_r_max", domain.get("r_max", expected_rmax)))
    if full_n != expected_n or abs(full_rmax - expected_rmax) > 1.0e-12:
        raise ValueError(f"{path}: unexpected domain metadata N={full_n}, r_max={full_rmax}")
    if radius.shape != (expected_n,) or not np.all(np.isfinite(radius)):
        raise ValueError(f"{path}: invalid radius array")
    return found


def restore_state(row, resolution: int, r_max: float):
    snap = row["state_snapshot"]
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        amplitude=AMPLITUDE,
        width=WIDTH,
        D_amplitude=D_AMPLITUDE,
        include_radiation=INCLUDE_RADIATION,
    )
    if len(state.grid.centers) != resolution:
        raise ValueError("reconstructed production grid has unexpected resolution")
    radius_key = "radius" if resolution == R80_N else "radius_full"
    radius_saved = np.asarray(row[radius_key], dtype=float)
    if not np.allclose(radius_saved, state.grid.centers, rtol=0.0, atol=1.0e-12):
        raise ValueError("archived t=0 radius values do not match production grid")
    for name, values in snap["geometry"].items():
        target = getattr(state.geometry, name)
        arr = np.asarray(values, dtype=float)
        if arr.shape != (resolution,) or not np.all(np.isfinite(arr)):
            raise ValueError(f"invalid archived geometry field {name}")
        target[:] = arr
    for name, values in snap["scalars"].items():
        target = getattr(state.scalars, name)
        arr = np.asarray(values, dtype=float)
        if arr.shape != (resolution,) or not np.all(np.isfinite(arr)):
            raise ValueError(f"invalid archived scalar field {name}")
        target[:] = arr
    for species, fields in snap["matter"].items():
        target_species = getattr(state.matter, species)
        for name, values in fields.items():
            arr = np.asarray(values, dtype=float)
            if arr.shape != (resolution,) or not np.all(np.isfinite(arr)):
                raise ValueError(f"invalid archived matter field {species}.{name}")
            getattr(target_species, name)[:] = arr
    state.t = float(row["t"])
    state.tau = float(row.get("tau", 0.0))
    state.e_folds = float(row.get("e_folds", 0.0))
    return state


def snapshot_target(state, alpha):
    geometry = state.geometry.copy()
    geometry.alpha = np.asarray(alpha, dtype=float).copy()
    return float(cmc.target_kdot(state.grid, geometry, state.scalars, state.matter))


def iterate_self_consistent_map(state, outer_alpha: float, domain_name: str):
    alpha = np.asarray(state.geometry.alpha, dtype=float).copy()
    previous_delta = None
    records = []
    status = "max_iterations_reached"
    failure = None
    for iteration in range(MAX_ITERATIONS):
        geom = state.geometry.copy()
        geom.alpha = alpha.copy()
        target_input = float(cmc.target_kdot(state.grid, geom, state.scalars, state.matter))
        try:
            alpha_out, target_used = cmc.solve_cmc_lapse(
                state.grid, geom, state.scalars, state.matter,
                outer_alpha=float(outer_alpha),
            )
        except Exception as exc:
            failure = {
                "iteration": iteration,
                "type": type(exc).__name__,
                "error": str(exc),
            }
            status = "map_step_failed"
            records.append({
                "iteration": iteration,
                "target_input": target_input,
                "failure": failure,
            })
            break
        alpha_out = np.asarray(alpha_out, dtype=float)
        if not np.all(np.isfinite(alpha_out)):
            status = "nonfinite_output"
            failure = {"iteration": iteration, "error": "nonfinite alpha output"}
            break
        delta = alpha_out - alpha
        delta_rms = rms(delta)
        delta_max = float(np.max(np.abs(delta)))
        target_output = snapshot_target(state, alpha_out)
        drift = float(target_output - target_input)
        row = {
            "iteration": iteration,
            "target_input": target_input,
            "target_used_by_solver": float(target_used),
            "target_output_alpha": target_output,
            "target_output_minus_input": drift,
            "alpha_delta_rms": delta_rms,
            "alpha_delta_max_abs": delta_max,
            "alpha_input_rms": rms(alpha),
            "alpha_output_rms": rms(alpha_out),
            "alpha_output_min": float(np.min(alpha_out)),
            "alpha_output_max": float(np.max(alpha_out)),
            "delta_contraction_ratio": (
                None if previous_delta is None or rms(previous_delta) == 0.0
                else float(delta_rms / rms(previous_delta))
            ),
            "all_finite": bool(np.all(np.isfinite(alpha_out))),
            "positive": bool(np.min(alpha_out) > 0.0),
        }
        records.append(row)
        if not row["positive"]:
            status = "nonpositive_output"
            failure = {"iteration": iteration, "error": "non-positive alpha output"}
            break
        if delta_rms <= STOP_RMS:
            status = "converged"
            alpha = alpha_out
            break
        previous_delta = delta.copy()
        alpha = alpha_out.copy()

    return {
        "domain": domain_name,
        "status": status,
        "converged": status == "converged",
        "iteration_count": len(records),
        "outer_alpha_boundary": float(outer_alpha),
        "initial_alpha_rms": rms(np.asarray(state.geometry.alpha)),
        "final_alpha_rms": rms(alpha),
        "final_alpha_delta_rms": records[-1].get("alpha_delta_rms") if records else None,
        "final_target": records[-1].get("target_output_alpha") if records else None,
        "failure": failure,
        "iterations": records,
    }


def fixed_target_control(state, outer_alpha: float, fixed_target: float):
    original = cmc.target_kdot
    rows = []
    alpha_input = np.asarray(state.geometry.alpha, dtype=float).copy()
    try:
        cmc.target_kdot = lambda *args, **kwargs: float(fixed_target)
        alpha_outs = []
        for k in range(3):
            geom = state.geometry.copy()
            geom.alpha = alpha_input.copy()
            alpha_out, used = cmc.solve_cmc_lapse(
                state.grid, geom, state.scalars, state.matter,
                outer_alpha=float(outer_alpha),
            )
            alpha_out = np.asarray(alpha_out, dtype=float)
            rows.append({
                "repeat": k,
                "target_used": float(used),
                "alpha_output_rms": rms(alpha_out),
                "alpha_min": float(np.min(alpha_out)),
                "positive": bool(np.min(alpha_out) > 0.0),
            })
            alpha_outs.append(alpha_out.copy())
            alpha_input = alpha_out.copy()
    finally:
        cmc.target_kdot = original
    return {
        "fixed_target": float(fixed_target),
        "repeats": rows,
        "repeat_1_vs_2_rms": rms(alpha_outs[1] - alpha_outs[0]),
        "repeat_2_vs_3_rms": rms(alpha_outs[2] - alpha_outs[1]),
        "repeat_outputs_identical_to_1e-12": bool(
            rms(alpha_outs[1] - alpha_outs[0]) < 1.0e-12
            and rms(alpha_outs[2] - alpha_outs[1]) < 1.0e-12
        ),
    }


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: python diagnostics/cmc_target_lapse_fixed_point_initial.py "
            "<r80-boundary-matched-profiles.jsonl> <r160-reference-profiles.jsonl>"
        )
    r80_path, r160_path = map(Path, sys.argv[1:])
    row80 = read_initial_row(r80_path, R80_N, R80_MAX)
    row160 = read_initial_row(r160_path, R160_N, R160_MAX)
    state80 = restore_state(row80, R80_N, R80_MAX)
    state160 = restore_state(row160, R160_N, R160_MAX)

    alpha160 = np.asarray(state160.geometry.alpha, dtype=float)
    # R80's outer cell is at r=79.75, exactly the corresponding R160 cell.
    outer_alpha80 = float(alpha160[R80_N - 1])
    outer_alpha160 = float(alpha160[-1])
    target80_initial = snapshot_target(state80, state80.geometry.alpha)
    target160_initial = snapshot_target(state160, state160.geometry.alpha)
    target160_archived = float(row160["target_Kdot_full_domain"])

    self80 = iterate_self_consistent_map(state80, outer_alpha80, "R80_boundary_matched")
    self160 = iterate_self_consistent_map(state160, outer_alpha160, "R160_native_boundary")

    fixed80 = fixed_target_control(state80, outer_alpha80, target160_archived)
    fixed160 = fixed_target_control(state160, outer_alpha160, target160_archived)

    payload = {
        "schema": "cmc_target_lapse_fixed_point_initial_v1",
        "status": "completed",
        "diagnostic_only": True,
        "production_physics_changed": False,
        "production_kernel_source_changed": False,
        "iterations_are_undamped": True,
        "max_iterations": MAX_ITERATIONS,
        "stop_alpha_delta_rms": STOP_RMS,
        "settings": {
            "R80": {"N": R80_N, "r_max": R80_MAX, "dr": float(state80.grid.dr), "outer_alpha": outer_alpha80},
            "R160": {"N": R160_N, "r_max": R160_MAX, "dr": float(state160.grid.dr), "outer_alpha": outer_alpha160},
            "final_time": 0.0,
        },
        "initial_targets": {
            "R80_native_target": target80_initial,
            "R160_recomputed_target": target160_initial,
            "R160_archived_target": target160_archived,
            "R80_minus_R160_target_gap": float(target80_initial - target160_initial),
            "R160_target_reconstruction_error": float(target160_initial - target160_archived),
        },
        "self_consistent_map": {
            "R80": self80,
            "R160": self160,
        },
        "fixed_target_control": {
            "target_source": "archived R160 full-domain target on the initial slice",
            "R80": fixed80,
            "R160": fixed160,
        },
        "interpretation_limit": (
            "This diagnostic iterates the actual alpha-dependent CMC target map on the archived t=0 slice "
            "with all non-lapse fields frozen and no damping. Non-convergence is a recorded outcome, not a "
            "reason to tune the iteration. The fixed-target control checks repeatability of the elliptic "
            "solve; alpha iterates are map probes, not full spacetime evolution."
        ),
        "inputs": {
            "R80_path": str(r80_path),
            "R80_sha256": sha256(r80_path),
            "R160_path": str(r160_path),
            "R160_sha256": sha256(r160_path),
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "report.json"
    out.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print("[CMC_TARGET_LAPSE_FIXED_POINT] FINAL " + json.dumps({
        "status": payload["status"],
        "R80_map_status": self80["status"],
        "R160_map_status": self160["status"],
        "R80_iterations": self80["iteration_count"],
        "R160_iterations": self160["iteration_count"],
        "R80_minus_R160_target_gap": payload["initial_targets"]["R80_minus_R160_target_gap"],
        "R80_fixed_target_repeatability": fixed80["repeat_outputs_identical_to_1e-12"],
        "R160_fixed_target_repeatability": fixed160["repeat_outputs_identical_to_1e-12"],
        "report": str(out),
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "failure.json").write_text(json.dumps({
            "status": "failed",
            "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        raise
