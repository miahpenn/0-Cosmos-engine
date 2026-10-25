#!/usr/bin/env python3
"""Exact fixed-R mass/current control on the target kernel's actual areal maps.

The accepted t=0.60 profile artifact is immutable. We take copied forward
kernel steps at dt and dt/2, capture only their areal-coordinate maps R(t,r),
and evaluate the exact analytic identity
    m(t,R)=(M0+M_RATE*t) f(R)
    chi(t,R)=1-2*m(t,R)/R
    F(t,R)=M_RATE*f(R)
on those same maps. No target stress tensor, physical source, or target mass
is used. The test isolates interpolation/query error using the target's real
map motion and temporal sample spacing, not its physical conservation law.
"""
from __future__ import annotations

import argparse
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
N = int(os.environ.get("AREAL_TARGET_MAP_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_TARGET_MAP_OUT", f"results/manufactured-areal-target-map-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
RMAX = 80.0
M0 = 0.30
M_RATE = 0.12
ANNULUS = (2.0, 10.0)
SHAPES = ("rational", "gaussian", "broad_rational")
METHODS = ("linear", "cubic")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(name: str, value: Any) -> np.ndarray:
    arr = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(f"{name} contains non-finite values")
    return arr


def shape_function(radius: np.ndarray, shape: str) -> np.ndarray:
    r = finite("shape radii", radius)
    if shape == "rational":
        out = r**3 / (7.0**3 + r**3)
    elif shape == "gaussian":
        out = np.exp(-((r / 6.0) ** 2))
    elif shape == "broad_rational":
        out = r**2 / (12.0**2 + r**2)
    else:
        raise ValueError(f"unsupported analytic shape: {shape}")
    return finite(f"{shape} shape function", out)


def sample(values: np.ndarray, coordinates: np.ndarray, query: np.ndarray, method: str) -> np.ndarray:
    coords = finite("map coordinates", coordinates)
    vals = finite("field values", values)
    q = finite("fixed-R queries", query)
    if coords.ndim != 1 or vals.shape != coords.shape or np.any(np.diff(coords) <= 0):
        raise RuntimeError("target areal map must be a monotone one-dimensional profile")
    if q[0] < coords[0] or q[-1] > coords[-1]:
        raise RuntimeError("fixed-R query falls outside a target map")
    if method == "linear":
        return finite("linear sample", np.interp(q, coords, vals))
    if method == "cubic":
        return finite("cubic sample", np.array([areal.interp4(float(x), coords, vals) for x in q]))
    raise ValueError(f"unsupported interpolation method: {method}")


def rms(values: np.ndarray) -> float:
    arr = finite("residual array", values)
    return float(np.sqrt(np.mean(arr * arr)))


def norm_stats(residual: np.ndarray, exact: np.ndarray) -> dict[str, float]:
    res = finite("balance residual", residual)
    truth = finite("exact current", exact)
    return {
        "rms": rms(res),
        "max_abs": float(np.max(np.abs(res))),
        "rms_relative_to_exact_current": rms(res) / max(rms(truth), 1e-300),
    }


def manufactured_profiles(R: np.ndarray, time: float, shape: str) -> tuple[np.ndarray, np.ndarray]:
    f = shape_function(R, shape)
    mass = finite("manufactured mass", (M0 + M_RATE * time) * f)
    chi = finite("manufactured chi", 1.0 - 2.0 * mass / R)
    current = finite("manufactured current", M_RATE * f)
    return chi, current


def observed_order(coarse: float, fine: float) -> float:
    if coarse <= 0 or fine <= 0:
        raise RuntimeError("observed order requires positive norms")
    return float(np.log2(coarse / fine))


def capture_maps(kernel, accepted, h: float) -> dict[str, Any]:
    accepted_before = saved.tangent.state_digest(accepted)
    start_profile = areal.geometry_rate_and_mass(kernel, accepted)
    start_R = finite("accepted areal map", start_profile["R"]).copy()
    state1 = kernel.step(areal.clone_state(accepted), h)
    state1_before = saved.tangent.state_digest(state1)
    p1 = areal.geometry_rate_and_mass(kernel, state1)
    R1 = finite("first shadow areal map", p1["R"]).copy()
    state2 = kernel.step(areal.clone_state(state1), h)
    state2_before = saved.tangent.state_digest(state2)
    p2 = areal.geometry_rate_and_mass(kernel, state2)
    R2 = finite("second shadow areal map", p2["R"]).copy()

    if abs(float(state1.t) - (float(accepted.t) + h)) > 1e-10:
        raise RuntimeError("first shadow map is not at t0+dt")
    if abs(float(state2.t) - (float(accepted.t) + 2.0 * h)) > 1e-10:
        raise RuntimeError("second shadow map is not at t0+2dt")
    if saved.tangent.state_digest(accepted) != accepted_before:
        raise RuntimeError("accepted t=0.60 input state mutated")
    if saved.tangent.state_digest(state1) != state1_before:
        raise RuntimeError("first shadow state mutated during observation")
    if saved.tangent.state_digest(state2) != state2_before:
        raise RuntimeError("second shadow state mutated during observation")
    maps = (start_R, R1, R2)
    if any(np.any(np.diff(x) <= 0.0) for x in maps):
        raise RuntimeError("target areal map is not monotone at all temporal levels")
    if max(float(np.max(np.abs(x.geometry.beta))) for x in (accepted, state1, state2)) > 1e-14:
        raise RuntimeError("target map control requires the pinned zero-shift branch")
    times = (float(accepted.t), float(state1.t), float(state2.t))
    return {
        "times": list(times),
        "R_maps": maps,
        "map_min": [float(x[0]) for x in maps],
        "map_max": [float(x[-1]) for x in maps],
        "map_rate_rms": float(rms((R1 - start_R) / h)),
        "map_rate_max_abs": float(np.max(np.abs((R1 - start_R) / h))),
        "input_state_immutable": True,
        "shadow_states_immutable": True,
    }


def compare_map_control(maps: dict[str, Any], query: np.ndarray, shape: str, method: str) -> dict[str, Any]:
    R0, R1, R2 = maps["R_maps"]
    t0, t1, t2 = maps["times"]
    h1, h2 = t1 - t0, t2 - t1
    if abs(h1 - h2) > 1e-10:
        raise RuntimeError("manufactured temporal steps are unequal")
    c0, F0 = manufactured_profiles(R0, t0, shape)
    c1, F1 = manufactured_profiles(R1, t1, shape)
    c2, F2 = manufactured_profiles(R2, t2, shape)
    chi0 = sample(c0, R0, query, method)
    chi1 = sample(c1, R1, query, method)
    chi2 = sample(c2, R2, query, method)
    f0 = sample(F0, R0, query, method)
    f1 = sample(F1, R1, query, method)
    f2 = sample(F2, R2, query, method)
    exact = M_RATE * shape_function(query, shape)
    forward_rate = finite("forward manufactured fixed-R rate", -0.5 * query * (chi1 - chi0) / h1)
    centered_rate = finite("centered manufactured fixed-R rate", -0.5 * query * (chi2 - chi0) / (2.0 * h1))
    forward_balance = finite("forward balance", forward_rate - 0.5 * (f0 + f1))
    centered_balance = finite("centered balance", centered_rate - f1)
    return {
        "forward_minus_endpoint_average_current": norm_stats(forward_balance, exact),
        "centered_minus_midpoint_current": norm_stats(centered_balance, exact),
        "forward_mass_rate_minus_exact_current": norm_stats(forward_rate - exact, exact),
        "centered_mass_rate_minus_exact_current": norm_stats(centered_rate - exact, exact),
        "interpolated_midpoint_current_minus_exact": norm_stats(f1 - exact, exact),
    }


def run_resolution() -> dict[str, Any]:
    if N not in (60, 120, 240):
        raise ValueError("resolution must be 60, 120, or 240")
    pins = origin_shadow.verify_source_pins()
    if pins.get("source_commit") != SOURCE or pins.get("vendor_commit") != VENDOR:
        raise RuntimeError("source/vendor pin mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected source diff not empty")
    manifest = saved.check_artifact(FULL_DIR, N, full_profile=True)
    arrays = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, arrays, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("wrong SHA-verified accepted input time")

    dt = CFL * float(accepted.grid.dr)
    n_query = max(300, int(accepted.grid.n) * 2)
    query = np.linspace(ANNULUS[0], ANNULUS[1], n_query)
    arms = []
    for label, factor in (("dt", 1.0), ("dt_half", 0.5)):
        h = factor * dt
        maps = capture_maps(kernel, accepted, h)
        if maps["map_min"][0] > ANNULUS[0] or maps["map_max"][0] < ANNULUS[1]:
            raise RuntimeError("accepted target map does not cover the requested annulus")
        results = {
            shape: {
                method: compare_map_control(maps, query, shape, method)
                for method in METHODS
            } for shape in SHAPES
        }
        arms.append({
            "step_arm": label,
            "dt": h,
            "map_motion": {k: v for k, v in maps.items() if k != "R_maps"},
            "results": results,
        })
    return {
        "resolution": N,
        "accepted_time": TARGET_TIME,
        "dr": float(accepted.grid.dr),
        "dt": dt,
        "accepted_state_digest": saved.tangent.state_digest(accepted),
        "accepted_artifact_manifest": manifest,
        "arms": arms,
        "input_state_immutable": all(x["map_motion"]["input_state_immutable"] for x in arms),
        "shadow_states_immutable": all(x["map_motion"]["shadow_states_immutable"] for x in arms),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(f"results/manufactured-areal-target-map-N{N}-2026-10-25/report.json"))
    args = parser.parse_args()
    result = run_resolution()
    report = {
        "schema": "manufactured-exact-fixed-R-target-map-control/v1",
        "status": "completed_diagnostic_only",
        "source_pin": SOURCE,
        "vendor_pin": VENDOR,
        "runtime": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform()},
        "analytic_identity": {
            "mass": "(M0 + M_RATE*t) * f(R)",
            "chi": "1 - 2*mass/R",
            "current": "M_RATE*f(R)",
            "identity": "-R/2 * partial_t chi|R == current",
            "M0": M0,
            "M_RATE": M_RATE,
            "shapes": {
                "rational": "R^3/(7^3+R^3)",
                "gaussian": "exp(-(R/6)^2)",
                "broad_rational": "R^2/(12^2+R^2)",
            },
        },
        "configuration": {
            "resolution": N,
            "accepted_input_time": TARGET_TIME,
            "cfl": CFL,
            "fixed_R_annulus": list(ANNULUS),
            "interpolation_methods": list(METHODS),
            "analytic_profile_on_target_map": True,
            "target_state_and_shadow_steps_only_supply_R_maps": True,
            "target_stress_tensor_or_physical_source_used": False,
        },
        "result": result,
        "production_changed": False,
        "full_matter_coupled_propagation_identity_closed": False,
        "total_stress_energy_sector_conservation_closed": False,
        "scientific_admission": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (args.output.parent / "SHA256SUMS.txt").write_text(f"{sha256(args.output)}  {args.output.name}\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "resolution": N,
        "accepted_state_immutable": result["input_state_immutable"],
        "shadow_states_immutable": result["shadow_states_immutable"],
        "arms": [{
            "step_arm": arm["step_arm"],
            "dt": arm["dt"],
            "map_motion": arm["map_motion"],
            "results": arm["results"],
        } for arm in result["arms"]],
        "scientific_admission": False,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
