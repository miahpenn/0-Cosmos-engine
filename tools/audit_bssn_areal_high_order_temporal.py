#!/usr/bin/env python3
"""Research-only centred temporal discriminator with high-order fixed-R sampling.

Reuse the SHA-verified accepted t=0.60 state and take two forward copied kernel
steps for h=dt and h=dt/2. Compare forward and centred fixed-R mass rates with
endpoint/midpoint and time-averaged currents using both linear controls and
cubic interpolation. For each interpolation order, compare preassembled chi/F
with operand-wise interpolation followed by nonlinear reassembly.

No backward stepping, accepted-state mutation, source fitting, or production
changes are permitted. A passing diagnostic is not conservation admission.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))

from tools import audit_bssn_areal_misner_sharp_current as areal
from tools import audit_bssn_areal_operand_interpolation as operand_audit
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_HIGH_ORDER_TEMPORAL_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_HIGH_ORDER_TEMPORAL_OUT", f"results/bssn-areal-high-order-temporal-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
FIXED_RADII = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0)
METHODS = (
    ("composite_linear", "composite", "linear"),
    ("operand_linear", "operand", "linear"),
    ("composite_cubic", "composite", "cubic"),
    ("operand_cubic", "operand", "cubic"),
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(f"{name} contains non-finite values")
    return arr


def scalar_finite(name: str, value: Any) -> float:
    out = float(value)
    if not np.isfinite(out):
        raise RuntimeError(f"{name} is non-finite")
    return out


def query(method_name: str, profile: dict[str, Any], profile_arrays: dict[str, np.ndarray], R: float) -> dict[str, float]:
    _, source, order = next(spec for spec in METHODS if spec[0] == method_name)
    if source == "composite":
        return operand_audit.query_composite(profile, R, order)
    return operand_audit.query_operand(profile_arrays, R, order)


def error_summary(residual: np.ndarray, radii: np.ndarray, current: np.ndarray) -> dict[str, float]:
    residual = finite("residual", residual)
    current = finite("current scale", current)
    inner = radii <= 3.0
    annulus = (radii >= 2.0) & (radii <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("common radius samples do not cover required inner and annular regions")
    rms = float(np.sqrt(np.mean(residual[annulus] ** 2)))
    scale = float(np.sqrt(np.mean(current[annulus] ** 2)))
    return {
        "max_abs_R_le_3": float(np.max(np.abs(residual[inner]))),
        "rms_R_2_to_10": rms,
        "max_abs_R_2_to_10": float(np.max(np.abs(residual[annulus]))),
        "rms_relative_to_current_R_2_to_10": rms / max(scale, 1e-300),
    }


def profile_at_time(kernel, state):
    profile = areal.geometry_rate_and_mass(kernel, state)
    arrays = operand_audit.operand_arrays(profile)
    native_check = operand_audit.verify_native_recombination(profile, arrays)
    return profile, arrays, native_check


def temporal_arm(kernel, accepted, h: float) -> dict[str, Any]:
    accepted_digest = saved.tangent.state_digest(accepted)
    profile0, arrays0, check0 = profile_at_time(kernel, accepted)

    state1 = kernel.step(areal.clone_state(accepted), h)
    digest1 = saved.tangent.state_digest(state1)
    profile1, arrays1, check1 = profile_at_time(kernel, state1)

    state2 = kernel.step(areal.clone_state(state1), h)
    digest2 = saved.tangent.state_digest(state2)
    profile2, arrays2, check2 = profile_at_time(kernel, state2)

    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("saved accepted input was mutated during temporal probe")
    if saved.tangent.state_digest(state1) != digest1:
        raise RuntimeError("first shadow state mutated during second step or observation")
    if saved.tangent.state_digest(state2) != digest2:
        raise RuntimeError("second shadow state mutated during observation")
    if abs(float(state1.t) - (float(accepted.t) + h)) > 1e-10:
        raise RuntimeError("first shadow step did not advance exactly by h")
    if abs(float(state2.t) - (float(accepted.t) + 2.0 * h)) > 1e-10:
        raise RuntimeError("second shadow step did not advance exactly by h")
    if max(profile0["max_abs_beta"], profile1["max_abs_beta"], profile2["max_abs_beta"]) > 1e-14:
        raise RuntimeError("zero-shift precondition failed")

    lo = max(float(profile0["R"][0]), float(profile1["R"][0]), float(profile2["R"][0]))
    hi = min(float(profile0["R"][-1]), float(profile1["R"][-1]), float(profile2["R"][-1]), 10.0)
    if not hi > lo:
        raise RuntimeError("no common areal-radius interval across three forward states")
    radii = np.linspace(lo, hi, max(300, int(accepted.grid.n) * 2))

    samples: dict[str, dict[str, np.ndarray]] = {}
    profiles = (profile0, profile1, profile2)
    all_arrays = (arrays0, arrays1, arrays2)
    for method_name, _, _ in METHODS:
        chi = [np.empty_like(radii) for _ in range(3)]
        current = [np.empty_like(radii) for _ in range(3)]
        for i, radius in enumerate(radii):
            for j in range(3):
                q = query(method_name, profiles[j], all_arrays[j], float(radius))
                chi[j][i] = q["chi"]
                current[j][i] = q["current"]
        forward_rate = finite(f"{method_name}.forward_rate", -0.5 * radii * (chi[1] - chi[0]) / h)
        centered_rate = finite(f"{method_name}.centered_rate", -0.5 * radii * (chi[2] - chi[0]) / (2.0 * h))
        forward_average = finite(f"{method_name}.forward_current_average", 0.5 * (current[0] + current[1]))
        two_step_average = finite(f"{method_name}.two_step_current_average", (current[0] + 2.0 * current[1] + current[2]) / 4.0)
        forward_start = finite(f"{method_name}.forward_start_residual", forward_rate - current[0])
        forward_end = finite(f"{method_name}.forward_end_residual", forward_rate - current[1])
        centered_mid = finite(f"{method_name}.centered_midpoint_residual", centered_rate - current[1])
        centered_average = finite(f"{method_name}.centered_average_residual", centered_rate - two_step_average)
        rate_gap = finite(f"{method_name}.centered_minus_forward_rate", centered_rate - forward_rate)
        samples[method_name] = {
            "chi": chi,
            "current": current,
            "forward_rate": forward_rate,
            "centered_rate": centered_rate,
            "forward_current_average": forward_average,
            "two_step_current_average": two_step_average,
            "forward_residual_vs_endpoint_average": finite(f"{method_name}.forward_residual_avg", forward_rate - forward_average),
            "forward_residual_vs_start_current": forward_start,
            "forward_residual_vs_end_current": forward_end,
            "centered_residual_vs_midpoint_current": centered_mid,
            "centered_residual_vs_two_step_average": centered_average,
            "centered_minus_forward_rate": rate_gap,
        }

    result = {}
    for method_name, _, _ in METHODS:
        v = samples[method_name]
        result[method_name] = {
            "forward": {
                "residual_vs_endpoint_time_average": error_summary(v["forward_residual_vs_endpoint_average"], radii, v["forward_current_average"]),
                "residual_vs_start_instantaneous_current": error_summary(v["forward_residual_vs_start_current"], radii, v["current"][0]),
                "residual_vs_end_instantaneous_current": error_summary(v["forward_residual_vs_end_current"], radii, v["current"][1]),
            },
            "centered": {
                "residual_vs_midpoint_instantaneous_current": error_summary(v["centered_residual_vs_midpoint_current"], radii, v["current"][1]),
                "residual_vs_two_step_time_average": error_summary(v["centered_residual_vs_two_step_average"], radii, v["two_step_current_average"]),
            },
            "temporal_sampling": {
                "centered_minus_forward_rate": error_summary(v["centered_minus_forward_rate"], radii, v["centered_rate"]),
            },
            "fixed_R": {},
        }
        _, source, order = next(spec for spec in METHODS if spec[0] == method_name)
        for R in FIXED_RADII:
            if not lo <= R <= hi:
                continue
            queries = []
            for j in range(3):
                queries.append(query(method_name, profiles[j], all_arrays[j], float(R)))
            fwd_rate = -0.5 * R * (queries[1]["chi"] - queries[0]["chi"]) / h
            cen_rate = -0.5 * R * (queries[2]["chi"] - queries[0]["chi"]) / (2.0 * h)
            fwd_avg = 0.5 * (queries[0]["current"] + queries[1]["current"])
            cen_avg = (queries[0]["current"] + 2.0 * queries[1]["current"] + queries[2]["current"]) / 4.0
            result[method_name]["fixed_R"][str(R)] = {
                "forward_mass_rate": scalar_finite(f"{method_name}.{R}.forward_rate", fwd_rate),
                "forward_endpoint_average_current": scalar_finite(f"{method_name}.{R}.forward_avg", fwd_avg),
                "forward_residual_vs_endpoint_average": scalar_finite(f"{method_name}.{R}.forward_residual", fwd_rate - fwd_avg),
                "centered_mass_rate": scalar_finite(f"{method_name}.{R}.centered_rate", cen_rate),
                "centered_midpoint_current": scalar_finite(f"{method_name}.{R}.mid_current", queries[1]["current"]),
                "centered_two_step_average_current": scalar_finite(f"{method_name}.{R}.centered_avg", cen_avg),
                "centered_residual_vs_midpoint_current": scalar_finite(f"{method_name}.{R}.centered_mid_residual", cen_rate - queries[1]["current"]),
                "centered_residual_vs_two_step_average": scalar_finite(f"{method_name}.{R}.centered_avg_residual", cen_rate - cen_avg),
            }

    comparison = {}
    for order in ("linear", "cubic"):
        composite = samples[f"composite_{order}"]
        operand = samples[f"operand_{order}"]
        for key in ("forward_residual_vs_endpoint_average", "centered_residual_vs_midpoint_current", "centered_residual_vs_two_step_average"):
            delta = finite(f"{order}.{key}.operand_minus_composite", operand[key] - composite[key])
            annulus = (radii >= 2.0) & (radii <= 10.0)
            comparison[f"{order}_{key}_operand_minus_composite"] = {
                "rms_R_2_to_10": float(np.sqrt(np.mean(delta[annulus] ** 2))),
                "max_abs_R_2_to_10": float(np.max(np.abs(delta[annulus]))),
            }

    return {
        "step_h": float(h),
        "base_time": float(accepted.t),
        "first_shadow_time": float(state1.t),
        "second_shadow_time": float(state2.t),
        "centered_derivative_time": float(state1.t),
        "steps_advanced_from_saved_input": 2,
        "saved_accepted_input_immutable": True,
        "first_shadow_state_immutable_after_second_step": True,
        "second_shadow_state_immutable_during_observation": True,
        "common_areal_radius_interval": [float(lo), float(hi)],
        "common_samples": int(radii.size),
        "native_grid_reassembly": {"t0": check0, "t1": check1, "t2": check2},
        "methods": result,
        "operand_minus_composite": comparison,
    }


def main() -> int:
    if N not in (60, 120, 240):
        raise ValueError("resolution must be 60, 120 or 240")
    pins = origin_shadow.verify_source_pins()
    if pins.get("source_commit") != SOURCE or pins.get("vendor_commit") != VENDOR:
        raise RuntimeError("source/vendor pin mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected source/vendor/dependency diff is not empty")
    manifest = saved.check_artifact(FULL_DIR, N, full_profile=True)
    profile_arrays = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, profile_arrays, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not at t=0.60")
    digest = saved.tangent.state_digest(accepted)
    dt = CFL * float(accepted.grid.dr)
    arms = [{"arm": label, **temporal_arm(kernel, accepted, h)} for label, h in (("dt", dt), ("dt_half", dt / 2.0))]
    if saved.tangent.state_digest(accepted) != digest:
        raise RuntimeError("accepted input digest changed across temporal arms")

    report = {
        "schema": "bssn-areal-high-order-temporal-discriminator/v1",
        "status": "completed_diagnostic_only",
        "source_pins": {
            "source_commit": SOURCE,
            "vendor_commit": VENDOR,
            "live_pin_check": pins,
            "accepted_artifact_manifest": manifest,
            "accepted_state_digest_before_after_match": True,
        },
        "runtime": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "configuration": {
            "resolution": N,
            "accepted_time": TARGET_TIME,
            "dt_rule": "dt=0.03*dr; dt and dt/2",
            "fixed_R_mass_rate_forward": "-R*(chi(t+h)-chi(t))/(2*h)",
            "fixed_R_mass_rate_centered": "-R*(chi(t+2h)-chi(t))/(4*h), centred at t+h",
            "forward_current_average": "(F(t)+F(t+h))/2",
            "centered_midpoint_current": "F(t+h)",
            "centered_two_step_time_average": "(F(t)+2F(t+h)+F(t+2h))/4",
            "methods": "linear/cubic; preassembled chi/F versus operand-wise interpolation and reassembly",
            "state handling": "two copied forward steps only per timestep arm; accepted state never stepped backward or mutated",
        },
        "arms": arms,
        "interpretation_guard": [
            "The centred and forward rates are located at different times; their difference is reported as a temporal-sampling discriminator, not a formal order estimate by itself.",
            "A high-order residual is not a conservation proof. Inspect the dt and dt/2 arms and the spatial trend separately.",
            "The operand route interpolates the pinned observer's R_r and R_t fields; it does not independently derive R_r from interpolated primitive metric fields.",
            "No production equation, source normalization, initializer, vendor, or dependency was changed.",
        ],
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_input_immutable": True,
            "two_forward_shadow_steps_per_arm": True,
            "first_and_second_shadow_states_immutable": True,
            "native_grid_operand_reassembly_checked": True,
            "linear_low_order_control_included": True,
            "high_order_centered_temporal_sampling_completed": True,
            "integrated_current_identity_closed": False,
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_conservation_closed": False,
            "scientific_admission": False,
        },
        "scientific_admission": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SHA256SUMS.txt").write_text(f"{sha256(report_path)}  report.json\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "source_vendor_pins_match": True,
        "accepted_input_immutable": True,
        "scientific_admission": False,
        "arms": [{
            "arm": arm["arm"],
            "methods": {
                name: {
                    "forward_vs_endpoint_average": result["forward"]["residual_vs_endpoint_time_average"],
                    "centered_vs_midpoint": result["centered"]["residual_vs_midpoint_instantaneous_current"],
                    "centered_vs_two_step_average": result["centered"]["residual_vs_two_step_time_average"],
                    "centered_minus_forward_rate": result["temporal_sampling"]["centered_minus_forward_rate"],
                } for name, result in arm["methods"].items()
            },
            "operand_minus_composite": arm["operand_minus_composite"],
            "fixed_R_operand_cubic": arm["methods"]["operand_cubic"]["fixed_R"],
        } for arm in arms],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
