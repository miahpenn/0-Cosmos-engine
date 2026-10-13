#!/usr/bin/env python3
"""Read-only split of the production BSSN/PIRK L2 Aa RHS at real step stages.

The pinned vendor L2 Aa equation is split into its signed trace-free lapse
Hessian and trace-free Ricci pieces. Original production functions are wrapped
only to observe their actual inputs and outputs; returned objects and arrays
are never replaced or edited. Each resolution takes one ordinary production
step. This is a numerical attribution diagnostic, not a physics/admission gate.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SOURCE_COMMIT = "cb22285df4c578b1cf76931fbba7661b501d86e6"
SOURCE_KERNEL_BLOB = "2954e7a79656849b5f088053ea0f4ac97d4bb443"
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RESOLUTIONS = (40, 60, 80, 120, 160)
R_MAX = 80.0
CFL = 0.03
FIXED_RADII = (2.0, 3.0, 5.0, 8.0)
CENTRE_CELLS = tuple(range(6))
STAGE_NAMES = (
    "initial_rhs",
    "explicit_predictor_rhs",
    "primary_predictor_rhs",
    "final_implicit_l2_rhs",
)
OUT = Path("results/pirk-l2-centre-stage-decomposition")


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def _array(value) -> np.ndarray:
    return np.asarray(value, dtype=float).copy()


def _max_abs(values) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.max(np.abs(arr))) if arr.size else 0.0


def _finite(values) -> bool:
    return bool(np.all(np.isfinite(np.asarray(values, dtype=float))))


def _state_digest(state) -> str:
    h = hashlib.sha256()
    named = {}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        named["geometry." + name] = getattr(state.geometry, name)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        named["scalars." + name] = getattr(state.scalars, name)
    for sector in ("dark_matter", "baryons", "radiation"):
        bundle = getattr(state.matter, sector)
        for name in ("rest", "energy_t", "momentum_r"):
            named["matter." + sector + "." + name] = getattr(bundle, name)
    for name in sorted(named):
        arr = np.ascontiguousarray(np.asarray(named[name]))
        h.update(name.encode("utf-8"))
        h.update(str(arr.dtype).encode("ascii"))
        h.update(repr(arr.shape).encode("ascii"))
        h.update(arr.tobytes())
    return h.hexdigest()


def _interp_at_radius(radii, values, target: float) -> dict:
    radii = np.asarray(radii, dtype=float)
    values = np.asarray(values, dtype=float)
    if target < radii[0] or target > radii[-1]:
        return {"target_r": float(target), "in_domain": False}
    right = int(np.searchsorted(radii, target, side="left"))
    if right < radii.size and np.isclose(radii[right], target, rtol=0.0, atol=1e-14):
        left = right
        weight = 0.0
    else:
        left = right - 1
        if left < 0 or right >= radii.size:
            return {"target_r": float(target), "in_domain": False}
        weight = float((target - radii[left]) / (radii[right] - radii[left]))
    value = float((1.0 - weight) * values[left] + weight * values[right])
    return {
        "target_r": float(target),
        "sample_r": float(target),
        "value": value,
        "left_i": left,
        "right_i": right,
        "weight_right": weight,
    }


def _component_summary(capture: dict) -> dict:
    radii = capture["r"]
    lapse = capture["lapse_hessian_Aa"]
    ricci = capture["ricci_tracefree_Aa"]
    rhs_aa = capture["rhs_Aa"]
    sum_aa = lapse + ricci
    closure_aa = sum_aa - rhs_aa
    closure_k = -capture["lap_alpha"] - capture["rhs_K"]
    centre = []
    for i in CENTRE_CELLS:
        if i < radii.size:
            centre.append({
                "i": int(i),
                "r": float(radii[i]),
                "lapse_hessian_Aa": float(lapse[i]),
                "ricci_tracefree_Aa": float(ricci[i]),
                "sum_Aa": float(sum_aa[i]),
                "vendor_l2_Aa": float(rhs_aa[i]),
                "Aa_closure": float(closure_aa[i]),
                "K_lapse_laplacian": float(-capture["lap_alpha"][i]),
                "vendor_l2_K": float(capture["rhs_K"][i]),
                "K_closure": float(closure_k[i]),
            })
    fixed = []
    for target in FIXED_RADII:
        fixed.append({
            "target_r": float(target),
            "lapse_hessian_Aa": _interp_at_radius(radii, lapse, target),
            "ricci_tracefree_Aa": _interp_at_radius(radii, ricci, target),
            "sum_Aa": _interp_at_radius(radii, sum_aa, target),
            "vendor_l2_Aa": _interp_at_radius(radii, rhs_aa, target),
            "Aa_closure": _interp_at_radius(radii, closure_aa, target),
            "K_lapse_laplacian": _interp_at_radius(radii, -capture["lap_alpha"], target),
            "vendor_l2_K": _interp_at_radius(radii, capture["rhs_K"], target),
            "K_closure": _interp_at_radius(radii, closure_k, target),
        })
    aa_scale = max(1.0, _max_abs(lapse), _max_abs(ricci), _max_abs(rhs_aa))
    k_scale = max(1.0, _max_abs(capture["lap_alpha"]), _max_abs(capture["rhs_K"]))
    roundoff_aa = 64.0 * np.finfo(float).eps * aa_scale
    roundoff_k = 64.0 * np.finfo(float).eps * k_scale
    return {
        "stage": capture["stage"],
        "finite": all(_finite(capture[k]) for k in (
            "r", "alpha", "R", "Rr_r", "DrDr_alpha", "lap_alpha",
            "lapse_hessian_Aa", "ricci_tracefree_Aa", "rhs_Aa", "rhs_K"
        )),
        "lapse_hessian_Aa_max_abs": _max_abs(lapse),
        "ricci_tracefree_Aa_max_abs": _max_abs(ricci),
        "vendor_l2_Aa_max_abs": _max_abs(rhs_aa),
        "Aa_reconstruction_closure_max_abs": _max_abs(closure_aa),
        "Aa_roundoff_bound": roundoff_aa,
        "Aa_closure_within_roundoff_bound": bool(_max_abs(closure_aa) <= roundoff_aa),
        "K_lapse_laplacian_max_abs": _max_abs(capture["lap_alpha"]),
        "vendor_l2_K_max_abs": _max_abs(capture["rhs_K"]),
        "K_reconstruction_closure_max_abs": _max_abs(closure_k),
        "K_roundoff_bound": roundoff_k,
        "K_closure_within_roundoff_bound": bool(_max_abs(closure_k) <= roundoff_k),
        "Aa_closure_centre_i0_i5_max_abs": _max_abs(closure_aa[:6]),
        "Aa_closure_excluding_outer_two_max_abs": _max_abs(closure_aa[:-2]),
        "K_closure_centre_i0_i5_max_abs": _max_abs(closure_k[:6]),
        "centre_cells": centre,
        "fixed_coordinate_radius_samples": fixed,
    }


def _json_safe(value):
    """Normalize NumPy scalar/container types; strict JSON still rejects NaN/Inf."""
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def _write_outputs(report: dict, profiles: dict[str, np.ndarray]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(_json_safe(report), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **profiles)
    csv_path = OUT / "profiles.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("resolution", "stage", "i", "r", "lapse_hessian_Aa",
                         "ricci_tracefree_Aa", "l2_Aa_sum", "vendor_l2_Aa",
                         "Aa_closure", "minus_lap_alpha", "vendor_l2_K", "K_closure"))
        for key in sorted(profiles):
            if not key.endswith("_r"):
                continue
            prefix = key[:-2]
            r = profiles[key]
            required = ("lapse_hessian_Aa", "ricci_tracefree_Aa", "rhs_Aa",
                        "lap_alpha", "rhs_K")
            if not all(prefix + "_" + suffix in profiles for suffix in required):
                continue
            vals = {k: profiles[prefix + "_" + k] for k in required}
            n_tag = prefix.split("_")[0]
            stage = "_".join(prefix.split("_")[1:])
            for i, radius in enumerate(r):
                lapse = vals["lapse_hessian_Aa"][i]
                ricci = vals["ricci_tracefree_Aa"][i]
                rhs_aa = vals["rhs_Aa"][i]
                rhs_k = vals["rhs_K"][i]
                lap = vals["lap_alpha"][i]
                writer.writerow((n_tag, stage, i, float(radius), float(lapse),
                                 float(ricci), float(lapse + ricci), float(rhs_aa),
                                 float(lapse + ricci - rhs_aa), float(-lap),
                                 float(rhs_k), float(-lap - rhs_k)))
    checksum_paths = (report_path, npz_path, csv_path)
    lines = [f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}"
             for path in checksum_paths]
    (OUT / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUT / "SOURCE_COMMIT.txt").write_text(SOURCE_COMMIT + "\n", encoding="utf-8")
    (OUT / "SOURCE_KERNEL_BLOB.txt").write_text(SOURCE_KERNEL_BLOB + "\n", encoding="utf-8")
    (OUT / "VENDOR_COMMIT.txt").write_text(VENDOR_COMMIT + "\n", encoding="utf-8")


def main() -> int:
    head = _git("rev-parse", "HEAD")
    source_blob = _git("rev-parse", f"{SOURCE_COMMIT}:engine/production_kernel.py")
    head_blob = _git("rev-parse", "HEAD:engine/production_kernel.py")
    vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    if source_blob != SOURCE_KERNEL_BLOB or head_blob != SOURCE_KERNEL_BLOB:
        raise RuntimeError(
            f"production kernel provenance mismatch: expected {SOURCE_KERNEL_BLOB}; "
            f"source={source_blob}; diagnostic branch={head_blob}"
        )
    if vendor != VENDOR_COMMIT:
        raise RuntimeError(f"vendor pin mismatch: expected {VENDOR_COMMIT}, found {vendor}")

    from engine import v55_pirk_adapter as adapter
    from engine.production_kernel import V55ProductionKernel
    _, vacuum, _ = adapter.vendor_modules()

    reports = []
    profiles: dict[str, np.ndarray] = {}
    errors = []
    for resolution in RESOLUTIONS:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=resolution, r_max=R_MAX, amplitude=0.01,
            width=7.0, D_amplitude=1.0e-10, include_radiation=True,
        )
        dt = CFL * float(state.grid.dr)
        input_digest = _state_digest(state)
        captures = []
        active_stage = [None]
        capture_enabled = [False]
        original_stage_terms = adapter.geometry_stage_terms
        original_geometry_terms = vacuum.geometry_terms
        original_l2 = vacuum.primary_l2_rhs
        stage_names = STAGE_NAMES[:3]
        stage_counter = [0]

        def wrapped_geometry_terms(grid, geometry):
            terms = original_geometry_terms(grid, geometry)
            caller = sys._getframe(1)
            # Observe only primary L2 calls explicitly selected by the actual
            # stage-split adapter or the direct final implicit update. Other
            # constraint/diagnostic uses of the same helper remain untouched.
            if capture_enabled[0] and caller.f_code.co_name == "primary_l2_rhs":
                radii = _array(grid.centers)
                alpha = _array(geometry.alpha)
                captured = {
                    "stage": active_stage[0] or "final_implicit_l2_rhs",
                    "r": radii,
                    "alpha": alpha,
                    "R": _array(terms["R"]),
                    "Rr_r": _array(terms["Rr_r"]),
                    "DrDr_alpha": _array(terms["DrDr_alpha"]),
                    "lap_alpha": _array(terms["lap_alpha"]),
                }
                captured["lapse_hessian_Aa"] = -(
                    captured["DrDr_alpha"] - captured["lap_alpha"] / 3.0
                )
                captured["ricci_tracefree_Aa"] = captured["alpha"] * (
                    captured["Rr_r"] - captured["R"] / 3.0
                )
                captures.append(captured)
            return terms

        def wrapped_l2(grid, geometry):
            caller = sys._getframe(1)
            caller_name = caller.f_code.co_name
            should_capture = caller_name in ("geometry_stage_terms", "step")
            before = len(captures)
            previous_capture_state = capture_enabled[0]
            capture_enabled[0] = should_capture
            try:
                result = original_l2(grid, geometry)
            finally:
                capture_enabled[0] = previous_capture_state
            if not should_capture:
                return result
            if len(captures) != before + 1:
                raise RuntimeError(
                    f"{caller_name}: expected one geometry_terms call for selected L2 stage; "
                    f"captured {len(captures) - before}"
                )
            captured = captures[-1]
            captured["rhs_Aa"] = _array(result["Aa"])
            captured["rhs_K"] = _array(result["K"])
            captured["stage"] = active_stage[0] or "final_implicit_l2_rhs"
            return result

        def wrapped_stage_terms(grid, geometry, fields, matter, lambda_m=2.0):
            index = stage_counter[0]
            if index >= len(stage_names):
                raise RuntimeError("unexpected additional geometry stage during production step")
            active_stage[0] = stage_names[index]
            before = len(captures)
            try:
                result = original_stage_terms(
                    grid, geometry, fields, matter, lambda_m=lambda_m
                )
            finally:
                active_stage[0] = None
            if len(captures) != before + 1:
                raise RuntimeError(
                    f"{stage_names[index]} expected exactly one L2 split call, "
                    f"observed {len(captures) - before}"
                )
            stage_counter[0] += 1
            return result

        adapter.geometry_stage_terms = wrapped_stage_terms
        vacuum.geometry_terms = wrapped_geometry_terms
        vacuum.primary_l2_rhs = wrapped_l2
        try:
            accepted = kernel.step(state, dt)
        finally:
            adapter.geometry_stage_terms = original_stage_terms
            vacuum.geometry_terms = original_geometry_terms
            vacuum.primary_l2_rhs = original_l2

        output_digest = _state_digest(state)
        if input_digest != output_digest:
            raise AssertionError(
                f"N={resolution}: production step mutated its input state in place"
            )
        if len(captures) != 4:
            raise AssertionError(
                f"N={resolution}: expected 4 real L2 Aa/K evaluations; captured {len(captures)}"
            )
        expected_stages = list(STAGE_NAMES)
        observed_stages = [item["stage"] for item in captures]
        if observed_stages != expected_stages:
            raise AssertionError(
                f"N={resolution}: stage order mismatch: {observed_stages} != {expected_stages}"
            )
        accepted.geometry.assert_finite_positive()
        if not (np.isfinite(accepted.t) and accepted.t > state.t):
            raise AssertionError(f"N={resolution}: production step did not advance time")

        stage_reports = []
        all_closures_pass = True
        for capture in captures:
            summary = _component_summary(capture)
            stage_reports.append(summary)
            all_closures_pass = bool(
                all_closures_pass
                and summary["finite"]
                and summary["Aa_closure_within_roundoff_bound"]
                and summary["K_closure_within_roundoff_bound"]
            )
            prefix = f"N{resolution}_{capture['stage']}"
            for key in (
                "r", "alpha", "R", "Rr_r", "DrDr_alpha", "lap_alpha",
                "lapse_hessian_Aa", "ricci_tracefree_Aa", "rhs_Aa", "rhs_K"
            ):
                profiles[prefix + "_" + key] = _array(
                    capture["lap_alpha"] if key == "lap_alpha" else capture[key]
                )
        reports.append({
            "resolution": resolution,
            "r_max": R_MAX,
            "dr": float(state.grid.dr),
            "dt": float(dt),
            "source_state_sha256_before": input_digest,
            "source_state_sha256_after": output_digest,
            "source_state_unchanged": input_digest == output_digest,
            "accepted_time": float(accepted.t),
            "captured_actual_l2_calls": observed_stages,
            "stage_summaries": stage_reports,
            "all_component_reconstructions_within_roundoff": bool(all_closures_pass),
            "step_status": "completed",
        })
        if not all_closures_pass:
            errors.append(
                f"N={resolution}: L2 component reconstruction exceeded its roundoff bound"
            )

    status = "completed_diagnostic_only" if not errors else "diagnostic_closure_failure"
    report = {
        "schema": "v55_pirk_l2_centre_stage_decomposition_v1",
        "status": status,
        "source_commit": SOURCE_COMMIT,
        "source_production_kernel_blob": source_blob,
        "diagnostic_head_commit": head,
        "vendor_commit": vendor,
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "cfl": CFL,
        "requested_stage_order": list(STAGE_NAMES),
        "centre_cells": list(CENTRE_CELLS),
        "fixed_coordinate_radii": list(FIXED_RADII),
        "ordinary_production_steps_per_resolution": 1,
        "production_kernel_and_vendor_edited": False,
        "component_definitions": {
            "lapse_hessian_tracefree_Aa": "-(DrDr_alpha - lap_alpha/3)",
            "ricci_tracefree_Aa": "alpha*(Rr_r - R/3)",
            "vendor_l2_Aa": "lapse_hessian_tracefree_Aa + ricci_tracefree_Aa",
            "vendor_l2_K": "-lap_alpha",
        },
        "errors": errors,
        "resolution_reports": reports,
        "guardrails": [
            "The wrapper records vendor geometry_terms outputs but returns the original dictionary unchanged.",
            "The wrapper records primary_l2_rhs outputs but returns the original dictionary unchanged.",
            "No geometry/source array, RHS value, equation, coefficient, limiter, or boundary condition is edited.",
            "Stage decomposition is attribution of the coded discrete operator, not proof a term is physically incorrect.",
            "This run is one ordinary step per resolution, not a long campaign or admission test."
        ],
    }
    _write_outputs(report, profiles)
    print("PIRK_L2_CENTRE_STAGE_DECOMPOSITION_" + ("PASS" if not errors else "FAIL"))
    print(f"source={SOURCE_COMMIT} kernel_blob={source_blob} vendor={vendor} diagnostic_head={head}")
    print(f"resolutions={list(RESOLUTIONS)} steps_per_resolution=1")
    print("reconstruction_maxima=" + json.dumps([
        {
            "N": x["resolution"],
            "stages": [
                {
                    "stage": y["stage"],
                    "Aa": y["Aa_reconstruction_closure_max_abs"],
                    "Aa_bound": y["Aa_roundoff_bound"],
                    "K": y["K_reconstruction_closure_max_abs"],
                    "K_bound": y["K_roundoff_bound"],
                }
                for y in x["stage_summaries"]
            ],
        }
        for x in reports
    ], sort_keys=True))
    print(f"report_sha256={hashlib.sha256((OUT / 'report.json').read_bytes()).hexdigest()}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
