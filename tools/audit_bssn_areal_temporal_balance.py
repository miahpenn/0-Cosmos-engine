#!/usr/bin/env python3
"""Research-only temporal-sampling discriminator for the fixed-R MS balance.

Starting from the same SHA-verified accepted state at t=0.60, this probe takes
one and then a second FORWARD copied production-kernel step for each h in
{dt, dt/2}. It compares the existing one-step forward mass-rate estimate
with a two-step centred estimate, and checks endpoint instantaneous current
against the trapezoidal one-step and composite two-step time averages.

The archived accepted input is never stepped backward or mutated. No
production equation, initializer, vendor, or dependency is changed. Passing
this diagnostic establishes numerical bookkeeping only, not conservation.
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
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_TEMPORAL_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_TEMPORAL_OUT", f"results/bssn-areal-temporal-balance-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
FIXED_RADII = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(name: str, values: Any) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise RuntimeError(f"{name} contains non-finite values")
    return array


def interpolate(profile: dict[str, Any], field: str, radii: np.ndarray) -> np.ndarray:
    return finite(
        f"interpolated {field}",
        np.interp(radii, finite("profile R", profile["R"]), finite(field, profile[field])),
    )


def error_summary(error: np.ndarray, radii: np.ndarray, current: np.ndarray | None = None) -> dict[str, float]:
    err = finite("residual profile", error)
    inner = radii <= 3.0
    annulus = (radii >= 2.0) & (radii <= 10.0)
    if not np.any(inner) or not np.any(annulus):
        raise RuntimeError("common fixed-R samples do not cover required inner/annular regions")
    summary = {
        "max_abs_R_le_3": float(np.max(np.abs(err[inner]))),
        "rms_R_2_to_10": float(np.sqrt(np.mean(err[annulus] ** 2))),
        "max_abs_R_2_to_10": float(np.max(np.abs(err[annulus]))),
    }
    if current is not None:
        cur = finite("current profile for relative norm", current)
        scale = float(np.sqrt(np.mean(cur[annulus] ** 2)))
        summary["rms_relative_to_current_R_2_to_10"] = float(
            np.sqrt(np.mean(err[annulus] ** 2)) / max(scale, 1e-300)
        )
    return summary


def fixed_radius_sample(profile: dict[str, Any], radius: float) -> dict[str, float]:
    chi = float(np.interp(radius, profile["R"], profile["chi"]))
    current = float(np.interp(radius, profile["R"], profile["rhs"]))
    return {"chi": chi, "mass": 0.5 * radius * (1.0 - chi), "current": current}


def probe_temporal_window(kernel, accepted, h: float) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    accepted_digest_before = saved.tangent.state_digest(accepted)

    # Two forward steps only. Each call receives a copied state, and each
    # intermediate state is retained as an immutable observation input.
    state1 = kernel.step(areal.clone_state(accepted), h)
    state1_digest_before = saved.tangent.state_digest(state1)
    state2 = kernel.step(areal.clone_state(state1), h)
    state2_digest_before = saved.tangent.state_digest(state2)

    if abs(float(state1.t) - (float(accepted.t) + h)) > 1e-10:
        raise RuntimeError("first shadow step did not advance by exactly h")
    if abs(float(state2.t) - (float(accepted.t) + 2.0 * h)) > 1e-10:
        raise RuntimeError("second shadow step did not advance by exactly h")

    p0 = areal.geometry_rate_and_mass(kernel, accepted)
    p1 = areal.geometry_rate_and_mass(kernel, state1)
    p2 = areal.geometry_rate_and_mass(kernel, state2)

    if saved.tangent.state_digest(accepted) != accepted_digest_before:
        raise RuntimeError("accepted t=0.60 input was mutated by temporal probe")
    if saved.tangent.state_digest(state1) != state1_digest_before:
        raise RuntimeError("first shadow state was mutated during second-step/observer use")
    if saved.tangent.state_digest(state2) != state2_digest_before:
        raise RuntimeError("second shadow state was mutated by the observer")

    lo = max(float(p0["R"][0]), float(p1["R"][0]), float(p2["R"][0]))
    hi = min(float(p0["R"][-1]), float(p1["R"][-1]), float(p2["R"][-1]), 10.0)
    if not hi > lo:
        raise RuntimeError("no common fixed-areal-radius interval across the three states")
    radii = np.linspace(lo, hi, max(300, int(accepted.grid.n) * 2))

    chi0 = interpolate(p0, "chi", radii)
    chi1 = interpolate(p1, "chi", radii)
    chi2 = interpolate(p2, "chi", radii)
    f0 = interpolate(p0, "rhs", radii)
    f1 = interpolate(p1, "rhs", radii)
    f2 = interpolate(p2, "rhs", radii)

    # Primary rates are derived from chi at common R, preserving m=R(1-chi)/2.
    forward_rate = finite("forward accepted-step mass rate", -0.5 * radii * (chi1 - chi0) / h)
    centered_rate = finite("two-forward-step centred mass rate", -0.5 * radii * (chi2 - chi0) / (2.0 * h))

    # One-step endpoint average, plus the corresponding endpoint-instantaneous controls.
    forward_current_average = finite("one-step trapezoidal current", 0.5 * (f0 + f1))
    forward_residual_average = finite("forward minus endpoint-average current", forward_rate - forward_current_average)
    forward_residual_start = finite("forward minus start current", forward_rate - f0)
    forward_residual_end = finite("forward minus end current", forward_rate - f1)

    # The centred derivative is located at t+h. Compare it to the instantaneous
    # midpoint current and the composite-trapezoid average over [t,t+2h].
    two_step_current_average = finite("two-step time-averaged current", (f0 + 2.0 * f1 + f2) / 4.0)
    centered_residual_midpoint = finite("centred derivative minus midpoint current", centered_rate - f1)
    centered_residual_average = finite("centred derivative minus two-step average current", centered_rate - two_step_current_average)
    centered_residual_start = finite("centred derivative minus first endpoint current", centered_rate - f0)
    centered_residual_end = finite("centred derivative minus second endpoint current", centered_rate - f2)
    rate_sampling_gap = finite("centred minus forward rate", centered_rate - forward_rate)
    current_step_01 = finite("current change step 0 to 1", f1 - f0)
    current_step_12 = finite("current change step 1 to 2", f2 - f1)

    fixed_r: dict[str, Any] = {}
    for radius in FIXED_RADII:
        if not (lo <= radius <= hi):
            continue
        s0 = fixed_radius_sample(p0, radius)
        s1 = fixed_radius_sample(p1, radius)
        s2 = fixed_radius_sample(p2, radius)
        fwd = -0.5 * radius * (s1["chi"] - s0["chi"]) / h
        cen = -0.5 * radius * (s2["chi"] - s0["chi"]) / (2.0 * h)
        fwd_avg = 0.5 * (s0["current"] + s1["current"])
        cen_avg = (s0["current"] + 2.0 * s1["current"] + s2["current"]) / 4.0
        fixed_r[str(radius)] = {
            "forward_mass_rate": float(fwd),
            "forward_start_instantaneous_current": s0["current"],
            "forward_end_instantaneous_current": s1["current"],
            "forward_endpoint_time_average_current": float(fwd_avg),
            "forward_residual_vs_start_current": float(fwd - s0["current"]),
            "forward_residual_vs_end_current": float(fwd - s1["current"]),
            "forward_residual_vs_time_average": float(fwd - fwd_avg),
            "centered_mass_rate": float(cen),
            "centered_start_instantaneous_current": s0["current"],
            "centered_midpoint_instantaneous_current": s1["current"],
            "centered_end_instantaneous_current": s2["current"],
            "centered_two_step_time_average_current": float(cen_avg),
            "centered_residual_vs_midpoint_current": float(cen - s1["current"]),
            "centered_residual_vs_two_step_time_average": float(cen - cen_avg),
            "centered_minus_forward_rate": float(cen - fwd),
            "accepted_mass": s0["mass"],
            "midpoint_mass": s1["mass"],
            "second_shadow_mass": s2["mass"],
        }

    arrays = {
        "R": radii.copy(),
        "chi_t0": chi0.copy(), "chi_t1": chi1.copy(), "chi_t2": chi2.copy(),
        "current_t0": f0.copy(), "current_t1": f1.copy(), "current_t2": f2.copy(),
        "forward_mass_rate": forward_rate.copy(), "centered_mass_rate": centered_rate.copy(),
        "forward_current_time_average": forward_current_average.copy(),
        "two_step_current_time_average": two_step_current_average.copy(),
        "forward_residual_vs_time_average": forward_residual_average.copy(),
        "forward_residual_vs_start_current": forward_residual_start.copy(),
        "forward_residual_vs_end_current": forward_residual_end.copy(),
        "centered_residual_vs_midpoint_current": centered_residual_midpoint.copy(),
        "centered_residual_vs_two_step_time_average": centered_residual_average.copy(),
        "centered_residual_vs_start_current": centered_residual_start.copy(),
        "centered_residual_vs_end_current": centered_residual_end.copy(),
        "centered_minus_forward_rate": rate_sampling_gap.copy(),
        "current_change_01": current_step_01.copy(),
        "current_change_12": current_step_12.copy(),
    }

    report = {
        "step_h": float(h),
        "base_time": float(accepted.t),
        "first_shadow_time": float(state1.t),
        "second_shadow_time": float(state2.t),
        "forward_interval": [float(accepted.t), float(state1.t)],
        "centred_interval": [float(accepted.t), float(state2.t)],
        "centred_derivative_time": float(state1.t),
        "steps_advanced_from_saved_input": 2,
        "saved_accepted_input_immutable": True,
        "first_shadow_state_immutable_after_second_step": True,
        "second_shadow_state_immutable_during_observation": True,
        "common_areal_radius_interval": [float(lo), float(hi)],
        "common_samples": int(radii.size),
        "forward": {
            "residual_vs_endpoint_time_average": error_summary(forward_residual_average, radii, forward_current_average),
            "residual_vs_start_instantaneous_current": error_summary(forward_residual_start, radii, f0),
            "residual_vs_end_instantaneous_current": error_summary(forward_residual_end, radii, f1),
        },
        "centred": {
            "residual_vs_midpoint_instantaneous_current": error_summary(centered_residual_midpoint, radii, f1),
            "residual_vs_two_step_time_average": error_summary(centered_residual_average, radii, two_step_current_average),
            "residual_vs_start_instantaneous_current": error_summary(centered_residual_start, radii, f0),
            "residual_vs_end_instantaneous_current": error_summary(centered_residual_end, radii, f2),
        },
        "temporal_sampling": {
            "centered_minus_forward_rate": error_summary(rate_sampling_gap, radii, centered_rate),
            "current_change_0_to_1": error_summary(current_step_01, radii, f0),
            "current_change_1_to_2": error_summary(current_step_12, radii, f1),
        },
        "fixed_areal_radius": fixed_r,
        "max_abs_beta_all_three_states": float(max(p0["max_abs_beta"], p1["max_abs_beta"], p2["max_abs_beta"])),
    }
    return report, arrays


def main() -> int:
    if N not in (60, 120, 240):
        raise ValueError("resolution must be 60, 120 or 240")

    pin = origin_shadow.verify_source_pins()
    if pin.get("source_commit") != SOURCE or pin.get("vendor_commit") != VENDOR:
        raise RuntimeError("source/vendor pin mismatch")
    if pin.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected source diff not empty")

    manifest = saved.check_artifact(FULL_DIR, N, full_profile=True)
    profile_arrays = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, profile_arrays, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not at t=0.60")

    accepted_digest = saved.tangent.state_digest(accepted)
    dt = CFL * float(accepted.grid.dr)
    reports = []
    saved_profiles: dict[str, np.ndarray] = {}
    for label, h in (("dt", dt), ("dt_half", dt / 2.0)):
        arm_report, arm_profiles = probe_temporal_window(kernel, accepted, h)
        arm_report["arm"] = label
        reports.append(arm_report)
        for name, values in arm_profiles.items():
            saved_profiles[f"{label}_{name}"] = finite(f"{label}.{name}", values).copy()

    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted input digest changed across temporal arms")

    checks = {
        "source_pins_match": True,
        "accepted_artifact_sha_verified": True,
        "accepted_input_immutable_across_all_arms": True,
        "two_forward_shadow_steps_used_per_arm": all(x["steps_advanced_from_saved_input"] == 2 for x in reports),
        "first_shadow_states_immutable": all(x["first_shadow_state_immutable_after_second_step"] for x in reports),
        "second_shadow_states_immutable": all(x["second_shadow_state_immutable_during_observation"] for x in reports),
        "zero_shift_precondition": all(x["max_abs_beta_all_three_states"] <= 1e-14 for x in reports),
        "all_saved_temporal_profiles_finite": all(np.all(np.isfinite(x)) for x in saved_profiles.values()),
        "temporal_sampling_discriminator_completed": True,
        "integrated_current_identity_closed": False,
        "full_matter_coupled_propagation_identity_closed": False,
        "total_stress_energy_conservation_closed": False,
        "production_source_or_vendor_changed": False,
        "scientific_admission": False,
    }
    report = {
        "schema": "bssn-areal-temporal-sampling-discriminator/v1",
        "status": "completed_diagnostic_only",
        "source_pins": {
            "source_commit": SOURCE,
            "vendor_commit": VENDOR,
            "live_pin_check": pin,
            "full_profile_accepted_artifact": manifest,
        },
        "runtime": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "configuration": {
            "resolution": N,
            "accepted_time": TARGET_TIME,
            "dt_rule": "dt=0.03*dr; both dt and dt/2 arms",
            "accepted_input_reused_unchanged": True,
            "shadow_steps": "two forward copied production-kernel steps; no backward stepping",
            "fixed_R_mass_rate": "m_dot|R = -R*Delta(chi)|R/(2*Delta_t), interpolating chi to common R first",
            "forward_rate": "(m(t+h)-m(t))/h, evaluated from chi at common fixed R",
            "centred_rate": "(m(t+2h)-m(t))/(2h), centred at t+h and built from two forward steps",
            "one_step_current_average": "(F(t)+F(t+h))/2",
            "two_step_current_average": "(F(t)+2F(t+h)+F(t+2h))/4, composite trapezoid",
            "current": "F=4*pi*R^2*T_t^R from the existing pinned tensor transformation",
        },
        "probes": reports,
        "acceptance": checks,
        "interpretation_guard": [
            "Workflow success validates execution, pinned inputs and diagnostic bookkeeping only.",
            "The centred derivative and forward accepted-step difference are centred at different times; their difference is reported as a temporal-sampling discriminator, not a formal convergence order by itself.",
            "Instantaneous and time-averaged current residuals are separate measures; no source term is fitted or changed.",
            "Integrated current identity, full matter-coupled propagation and total stress-energy conservation remain open.",
        ],
        "scientific_admission": False,
    }

    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    profiles_path = OUT / "profiles.npz"
    np.savez_compressed(profiles_path, **saved_profiles)
    checksums_path = OUT / "SHA256SUMS.txt"
    checksums_path.write_text(
        f"{sha256(report_path)}  report.json\n{sha256(profiles_path)}  profiles.npz\n",
        encoding="utf-8",
    )

    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "source_pins_match": checks["source_pins_match"],
        "accepted_input_immutable": checks["accepted_input_immutable_across_all_arms"],
        "scientific_admission": False,
        "arms": {
            item["arm"]: {
                "forward_vs_endpoint_average": item["forward"]["residual_vs_endpoint_time_average"],
                "centered_vs_midpoint_current": item["centred"]["residual_vs_midpoint_instantaneous_current"],
                "centered_vs_two_step_average": item["centred"]["residual_vs_two_step_time_average"],
                "centered_minus_forward": item["temporal_sampling"]["centered_minus_forward_rate"],
                "fixed_R": item["fixed_areal_radius"],
            } for item in reports
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
