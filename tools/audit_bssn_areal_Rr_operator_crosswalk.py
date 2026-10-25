#!/usr/bin/env python3
"""Research-only crosswalk of three fixed-R areal-map radial derivatives.

Compare, at the same fixed R and with a four-point stencil anchored to the
accepted state throughout t, t+h and t+2h:
  1. derivative of the interpolated primitive b/X map (current raw-metric path);
  2. derivative of a cubic interpolant of the composed odd map R(r)=r*g(r^2);
  3. the existing fourth-order coordinate-grid derivative of native R(r),
     sampled at fixed R.

For all three options, hold the other primitive metric and stress operands,
including primitive-derived R_t, identical when rebuilding chi and the current.
This identifies whether the residual sensitivity comes from differentiating
b/X before composition, differentiating the composed R interpolant, or the
existing cell derivative. No option is designated physical authority.
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
from tools import audit_bssn_areal_primitive_metric_derivative as primitive
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_RR_CROSSWALK_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_RR_CROSSWALK_OUT", f"results/bssn-areal-Rr-crosswalk-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
FOUR_PI = 4.0 * np.pi
METHODS = (
    "primitive_bX_derivative",
    "composed_R_odd_cubic_derivative",
    "native_grid_fourth_derivative",
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


def scalar(name: str, value: Any) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise RuntimeError(f"{name} is non-finite")
    return result


def error_summary(values: np.ndarray, radii: np.ndarray) -> dict[str, float]:
    v = finite("summary values", values)
    inner = radii <= 3.0
    annulus = (radii >= 2.0) & (radii <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("common R samples do not cover inner and annular regions")
    return {
        "max_abs_R_le_3": float(np.max(np.abs(v[inner]))),
        "rms_R_2_to_10": float(np.sqrt(np.mean(v[annulus] ** 2))),
        "max_abs_R_2_to_10": float(np.max(np.abs(v[annulus]))),
    }


def composed_R_and_Rr(profile: dict[str, Any], r: float, stencil_start: int) -> tuple[float, float]:
    centers = finite("native coordinate centers", profile["r"])
    R_native = finite("native areal radius", profile["observer_R"])
    s = centers * centers
    g = R_native / centers
    nodes = s[stencil_start:stencil_start + 4]
    vals = g[stencil_start:stencil_start + 4]
    g_value, dg_ds = primitive.lagrange_value_and_derivative(nodes, vals, float(r) * float(r))
    Rvalue = float(r) * g_value
    Rr = g_value + 2.0 * float(r) * float(r) * dg_ds
    return scalar("composed R(r)", Rvalue), scalar("composed-interpolant Rr", Rr)


def recombine(q: dict[str, float], R_target: float, rr: float) -> dict[str, float]:
    a, X, alpha = q["a"], q["X"], q["alpha"]
    Rt, E, j, p_r = q["R_t_primitive_rate"], q["E"], q["j"], q["p_r"]
    if min(a, X, alpha, rr) <= 0:
        raise RuntimeError(f"nonpositive operand at R={R_target}: a={a}, X={X}, alpha={alpha}, Rr={rr}")
    g_rr = a / (X * X)
    chi = (X * X / a) * rr * rr - (Rt / alpha) ** 2
    j_up = j / g_rr
    T_t_R = -E * Rt - alpha * j_up * rr - (Rt / rr) * (j / alpha * Rt + p_r * rr)
    current = FOUR_PI * R_target * R_target * T_t_R
    return {
        "R_r": scalar("selected Rr", rr),
        "R_t": scalar("held primitive Rt", Rt),
        "chi": scalar("recomputed chi", chi),
        "current": scalar("recomputed transformed current", current),
    }


def temporal_arm(kernel, accepted, h: float) -> dict[str, Any]:
    accepted_digest = saved.tangent.state_digest(accepted)
    states = [accepted]
    states.append(kernel.step(areal.clone_state(states[0]), h))
    digest1 = saved.tangent.state_digest(states[1])
    states.append(kernel.step(areal.clone_state(states[1]), h))
    digest2 = saved.tangent.state_digest(states[2])
    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted input was mutated")
    if saved.tangent.state_digest(states[1]) != digest1:
        raise RuntimeError("first shadow state was mutated")
    if saved.tangent.state_digest(states[2]) != digest2:
        raise RuntimeError("second shadow state was mutated")
    if abs(float(states[1].t) - (float(accepted.t) + h)) > 1e-10:
        raise RuntimeError("first shadow time advance is incorrect")
    if abs(float(states[2].t) - (float(accepted.t) + 2.0 * h)) > 1e-10:
        raise RuntimeError("second shadow time advance is incorrect")

    profiles = [areal.geometry_rate_and_mass(kernel, state) for state in states]
    prims = [primitive.build_primitive_profile(kernel, p) for p in profiles]
    if max(p["max_abs_beta"] for p in profiles) > 1e-14:
        raise RuntimeError("zero-shift precondition failed")
    lo = max(float(p["R"][0]) for p in profiles)
    hi = min(float(p["R"][-1]) for p in profiles)
    hi = min(hi, 10.0)
    if not hi > lo:
        raise RuntimeError("no common fixed-R interval")
    radii = np.linspace(lo, hi, max(300, int(accepted.grid.n) * 2))

    samples: dict[str, list[dict[str, np.ndarray]]] = {
        method: [] for method in METHODS
    }
    derivative_gaps = {
        "primitive_minus_native": [],
        "composed_R_odd_minus_native": [],
        "primitive_minus_composed_R_odd": [],
        "composed_R_map_minus_fixed_R": [],
    }
    for j, (profile, prim) in enumerate(zip(profiles, prims)):
        by_method = {
            method: {key: np.empty_like(radii) for key in ("R_r", "chi", "current")}
            for method in METHODS
        }
        gaps = {key: np.empty_like(radii) for key in derivative_gaps}
        for i, R in enumerate(radii):
            anchor_r = primitive.invert_areal_radius(prims[0], float(R), "cubic")
            stencil = primitive.stencil_start_for_coordinate(prims[0], anchor_r)
            q = primitive.primitive_metric_query(
                prim,
                float(R),
                "cubic",
                stencil_start=stencil,
                seed=float(np.interp(R, prim["observer_R"], prim["r"])),
            )
            coordinate_r = q["coordinate_r"]
            composed_R, rr_composed = composed_R_and_Rr(prim, coordinate_r, stencil)
            rr_native = areal.interp4(float(R), profile["R"], profile["Rr"])
            rr_primitive = q["R_r_primitive_derivative"]
            rr_native = scalar("native grid derivative interpolated at fixed R", rr_native)
            rr_values = {
                "primitive_bX_derivative": rr_primitive,
                "composed_R_odd_cubic_derivative": rr_composed,
                "native_grid_fourth_derivative": rr_native,
            }
            for method, rr in rr_values.items():
                result = recombine(q, float(R), rr)
                by_method[method]["R_r"][i] = result["R_r"]
                by_method[method]["chi"][i] = result["chi"]
                by_method[method]["current"][i] = result["current"]
            gaps["primitive_minus_native"][i] = rr_primitive - rr_native
            gaps["composed_R_odd_minus_native"][i] = rr_composed - rr_native
            gaps["primitive_minus_composed_R_odd"][i] = rr_primitive - rr_composed
            gaps["composed_R_map_minus_fixed_R"][i] = composed_R - float(R)
        for method in METHODS:
            samples[method].append(by_method[method])
        for key in derivative_gaps:
            derivative_gaps[key].append(gaps[key])

    variant_report: dict[str, Any] = {}
    fixed_R_report: dict[str, Any] = {}
    for method, series in samples.items():
        chi0, chi1, chi2 = (x["chi"] for x in series)
        f0, f1, f2 = (x["current"] for x in series)
        forward_rate = finite(f"{method}.forward rate", -0.5 * radii * (chi1 - chi0) / h)
        centered_rate = finite(f"{method}.centered rate", -0.5 * radii * (chi2 - chi0) / (2.0 * h))
        forward_avg = finite(f"{method}.forward current avg", 0.5 * (f0 + f1))
        centered_avg = finite(f"{method}.centered current avg", (f0 + 2.0 * f1 + f2) / 4.0)
        variant_report[method] = {
            "forward_vs_endpoint_average": error_summary(forward_rate - forward_avg, radii),
            "centered_vs_midpoint_current": error_summary(centered_rate - f1, radii),
            "centered_vs_two_step_average": error_summary(centered_rate - centered_avg, radii),
            "centered_minus_forward_rate": error_summary(centered_rate - forward_rate, radii),
        }
        fixed_R_report[method] = {}
        for radius in (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0):
            if not lo <= radius <= hi:
                continue
            q0 = min(range(len(radii)), key=lambda idx: abs(float(radii[idx]) - radius))
            fixed_R_report[method][str(radius)] = {
                "sample_R": float(radii[q0]),
                "forward_residual": float((forward_rate - forward_avg)[q0]),
                "centered_midpoint_residual": float((centered_rate - f1)[q0]),
                "centered_average_residual": float((centered_rate - centered_avg)[q0]),
            }

    gap_report = {}
    for key, sequence in derivative_gaps.items():
        values = sequence[1]
        gap_report[key] = error_summary(values, radii)

    return {
        "step_h": float(h),
        "base_time": float(accepted.t),
        "first_shadow_time": float(states[1].t),
        "second_shadow_time": float(states[2].t),
        "steps_advanced_from_accepted_input": 2,
        "accepted_input_immutable": True,
        "first_shadow_state_immutable": True,
        "second_shadow_state_immutable": True,
        "common_areal_radius_interval": [float(lo), float(hi)],
        "common_samples": int(radii.size),
        "native_areal_map_max_abs_gap": max(p["native_areal_map_max_abs_gap"] for p in prims),
        "derivative_gaps_at_midpoint": gap_report,
        "variants": variant_report,
        "fixed_R": fixed_R_report,
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
    stored = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, stored, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not at t=0.60")
    digest = saved.tangent.state_digest(accepted)
    dt = CFL * float(accepted.grid.dr)
    arms = [{"arm": name, **temporal_arm(kernel, accepted, h)} for name, h in (("dt", dt), ("dt_half", dt / 2.0))]
    if saved.tangent.state_digest(accepted) != digest:
        raise RuntimeError("accepted input changed across temporal arms")
    report = {
        "schema": "bssn-areal-Rr-operator-crosswalk/v1",
        "status": "completed_diagnostic_only",
        "source_pins": {
            "source_commit": SOURCE,
            "vendor_commit": VENDOR,
            "live_pin_check": pins,
            "accepted_artifact_manifest": manifest,
            "accepted_state_digest_before_after_match": True,
        },
        "arms": arms,
        "interpretation_guard": [
            "This compares discrete derivative/reconstruction operators; no variant is declared physical authority.",
            "Smaller residual is not a method-admission criterion; the full resolution and timestep trend must be reviewed.",
            "The component crosswalk holds primitive-derived R_t and all other sampled operands fixed to isolate R_r.",
            "No production equation, source normalization, initializer, vendor/dependency, or admission flag changed.",
        ],
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_input_immutable": True,
            "forward_shadow_states_immutable": True,
            "all_three_Rr_operators_compared": True,
            "R_t_and_other_operands_held_fixed": True,
            "integrated_current_identity_closed": False,
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_conservation_closed": False,
            "scientific_admission": False,
        },
        "scientific_admission": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "report.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SHA256SUMS.txt").write_text(f"{sha256(path)}  report.json\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "scientific_admission": False,
        "arms": [{
            "arm": arm["arm"],
            "derivative_gaps_at_midpoint": arm["derivative_gaps_at_midpoint"],
            "variants": arm["variants"],
        } for arm in arms],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
