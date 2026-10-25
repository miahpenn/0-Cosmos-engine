#!/usr/bin/env python3
"""Research-only R_r/R_t component ablation for the fixed-R mass/current balance.

Using the same pinned accepted state, reconstruct raw metric/stress operands with
a cubic stencil anchored to the accepted state and reused at t, t+h and t+2h.
Then independently choose the areal-map derivative R_r and metric-rate R_t
from either (a) that primitive-metric reconstruction or (b) the pinned observer's
R-indexed fields. This provides a four-arm component ablation:

    primitive R_r / primitive R_t
    observer  R_r / primitive R_t
    primitive R_r / observer  R_t
    observer  R_r / observer  R_t

Only those two geometric operands are swapped; the other primitive metric and
stress operands are held on the common raw-metric interpolation path. Neither
discrete derivative is declared the physical authority. No production equations
or source normalizations are modified.
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
N = int(os.environ.get("AREAL_GEOMETRY_ABLATION_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_GEOMETRY_ABLATION_OUT", f"results/bssn-areal-geometry-ablation-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
FOUR_PI = 4.0 * np.pi
FIXED_RADII = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0)
VARIANTS = {
    "primitive_Rr_primitive_Rt": ("primitive", "primitive"),
    "observer_Rr_primitive_Rt": ("observer", "primitive"),
    "primitive_Rr_observer_Rt": ("primitive", "observer"),
    "observer_Rr_observer_Rt": ("observer", "observer"),
}


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


def observer_at_fixed_R(profile: dict[str, Any], R: float, field: str) -> float:
    # Use the observer's own fourth-order interpolation at the query R; do not
    # silently downgrade this control to linear interpolation.
    result = areal.interp4(R, profile["R"], profile[field])
    return scalar(f"observer {field}@R={R}", result)


def recombine(query: dict[str, float], R: float, rr_source: str, rt_source: str) -> dict[str, float]:
    rr = query["R_r_primitive_derivative"] if rr_source == "primitive" else query["R_r_observer_cubic"]
    rt = query["R_t_primitive_rate"] if rt_source == "primitive" else query["R_t_observer_cubic"]
    a = query["a"]
    X = query["X"]
    alpha = query["alpha"]
    E = query["E"]
    j = query["j"]
    p_r = query["p_r"]
    if a <= 0.0 or X <= 0.0 or alpha <= 0.0 or rr <= 0.0:
        raise RuntimeError(f"invalid metric/stress operand at R={R}: a={a}, X={X}, alpha={alpha}, Rr={rr}")
    g_rr = a / (X * X)
    chi = (X * X / a) * rr * rr - (rt / alpha) ** 2
    j_up = j / g_rr
    T_t_R = -E * rt - alpha * j_up * rr - (rt / rr) * (j / alpha * rt + p_r * rr)
    current = FOUR_PI * R * R * T_t_R
    return {
        "R_r": scalar("selected Rr", rr),
        "R_t": scalar("selected Rt", rt),
        "chi": scalar("reassembled chi", chi),
        "mass": scalar("Misner-Sharp mass", 0.5 * R * (1.0 - chi)),
        "current": scalar("transformed current", current),
    }


def error_summary(values: np.ndarray, radii: np.ndarray) -> dict[str, float]:
    values = finite("summary profile", values)
    inner = radii <= 3.0
    annulus = (radii >= 2.0) & (radii <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("fixed-R sampling does not cover the inner and annular regions")
    return {
        "max_abs_R_le_3": float(np.max(np.abs(values[inner]))),
        "rms_R_2_to_10": float(np.sqrt(np.mean(values[annulus] ** 2))),
        "max_abs_R_2_to_10": float(np.max(np.abs(values[annulus]))),
    }


def compare_arm(kernel, accepted, h: float) -> dict[str, Any]:
    accepted_digest = saved.tangent.state_digest(accepted)
    states = [accepted]
    states.append(kernel.step(areal.clone_state(states[0]), h))
    digest1 = saved.tangent.state_digest(states[1])
    states.append(kernel.step(areal.clone_state(states[1]), h))
    digest2 = saved.tangent.state_digest(states[2])

    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted input mutated during copied forward stepping")
    if saved.tangent.state_digest(states[1]) != digest1:
        raise RuntimeError("first shadow state mutated during second step")
    if saved.tangent.state_digest(states[2]) != digest2:
        raise RuntimeError("second shadow state mutated during observation")
    if abs(float(states[1].t) - (float(accepted.t) + h)) > 1e-10:
        raise RuntimeError("first shadow step did not advance by h")
    if abs(float(states[2].t) - (float(accepted.t) + 2.0 * h)) > 1e-10:
        raise RuntimeError("second shadow step did not advance by h")

    profiles = [areal.geometry_rate_and_mass(kernel, state) for state in states]
    prims = [primitive.build_primitive_profile(kernel, p) for p in profiles]
    if max(p["max_abs_beta"] for p in profiles) > 1e-14:
        raise RuntimeError("zero-shift precondition failed")
    lo = max(float(p["R"][0]) for p in profiles)
    hi = min(float(p["R"][-1]) for p in profiles)
    hi = min(hi, 10.0)
    if not hi > lo:
        raise RuntimeError("no common fixed-areal-radius interval")
    radii = np.linspace(lo, hi, max(300, int(accepted.grid.n) * 2))

    sampled: dict[str, list[dict[str, np.ndarray]]] = {}
    for name, (rr_source, rt_source) in VARIANTS.items():
        by_time = []
        for time_index, (profile, prim) in enumerate(zip(profiles, prims)):
            result = {
                "R_r": np.empty_like(radii),
                "R_t": np.empty_like(radii),
                "chi": np.empty_like(radii),
                "mass": np.empty_like(radii),
                "current": np.empty_like(radii),
            }
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
                q["R_r_observer_cubic"] = observer_at_fixed_R(profile, float(R), "Rr")
                q["R_t_observer_cubic"] = observer_at_fixed_R(profile, float(R), "Rt")
                values = recombine(q, float(R), rr_source, rt_source)
                for field in result:
                    result[field][i] = values[field]
            by_time.append(result)
        sampled[name] = by_time

    method_summary: dict[str, Any] = {}
    fixed_results: dict[str, Any] = {}
    for name, series in sampled.items():
        chi0, chi1, chi2 = (s["chi"] for s in series)
        cur0, cur1, cur2 = (s["current"] for s in series)
        forward_rate = finite(f"{name}.forward rate", -0.5 * radii * (chi1 - chi0) / h)
        centered_rate = finite(f"{name}.centered rate", -0.5 * radii * (chi2 - chi0) / (2.0 * h))
        forward_avg = finite(f"{name}.forward current average", 0.5 * (cur0 + cur1))
        centered_avg = finite(f"{name}.centered current average", (cur0 + 2.0 * cur1 + cur2) / 4.0)
        forward_res = forward_rate - forward_avg
        centered_mid_res = centered_rate - cur1
        centered_avg_res = centered_rate - centered_avg
        method_summary[name] = {
            "forward_vs_endpoint_average": error_summary(forward_res, radii),
            "centered_vs_midpoint_current": error_summary(centered_mid_res, radii),
            "centered_vs_two_step_average": error_summary(centered_avg_res, radii),
            "centered_minus_forward_rate": error_summary(centered_rate - forward_rate, radii),
        }
        fixed_results[name] = {}
        for R in FIXED_RADII:
            if not lo <= R <= hi:
                continue
            idx = int(np.argmin(np.abs(radii - R)))
            fixed_results[name][str(R)] = {
                "sample_R": float(radii[idx]),
                "forward_residual": float(forward_res[idx]),
                "centered_midpoint_residual": float(centered_mid_res[idx]),
                "centered_average_residual": float(centered_avg_res[idx]),
                "R_r_midpoint": float(series[1]["R_r"][idx]),
                "R_t_midpoint": float(series[1]["R_t"][idx]),
            }

    ablation: dict[str, Any] = {}
    baseline = sampled["primitive_Rr_primitive_Rt"]
    for name, series in sampled.items():
        if name == "primitive_Rr_primitive_Rt":
            continue
        delta_chi_mid = finite(f"{name}.delta chi", series[1]["chi"] - baseline[1]["chi"])
        delta_current_mid = finite(f"{name}.delta current", series[1]["current"] - baseline[1]["current"])
        ablation[name] = {
            "delta_chi_at_midpoint": error_summary(delta_chi_mid, radii),
            "delta_current_at_midpoint": error_summary(delta_current_mid, radii),
            "delta_R_r_at_midpoint": error_summary(series[1]["R_r"] - baseline[1]["R_r"], radii),
            "delta_R_t_at_midpoint": error_summary(series[1]["R_t"] - baseline[1]["R_t"], radii),
        }
    # Include per-operand deltas of each one-operator replacement against the
    # all-primitive arm, to distinguish direct Rr effects from direct Rt effects.
    rr_only = sampled["observer_Rr_primitive_Rt"]
    rt_only = sampled["primitive_Rr_observer_Rt"]
    both = sampled["observer_Rr_observer_Rt"]
    coupling = {}
    for field in ("chi", "current"):
        drr = rr_only[1][field] - baseline[1][field]
        drt = rt_only[1][field] - baseline[1][field]
        dboth = both[1][field] - baseline[1][field]
        coupling[field] = {
            "Rr_only_delta": error_summary(drr, radii),
            "Rt_only_delta": error_summary(drt, radii),
            "both_delta": error_summary(dboth, radii),
            "nonadditive_coupling_delta": error_summary(dboth - drr - drt, radii),
        }
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
        "variants": method_summary,
        "fixed_R": fixed_results,
        "midpoint_component_ablation": ablation,
        "component_nonlinearity": coupling,
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
    artifact_arrays = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, artifact_arrays, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not at t=0.60")
    digest = saved.tangent.state_digest(accepted)
    dt = CFL * float(accepted.grid.dr)
    arms = [{"arm": arm, **compare_arm(kernel, accepted, h)} for arm, h in (("dt", dt), ("dt_half", dt / 2.0))]
    if saved.tangent.state_digest(accepted) != digest:
        raise RuntimeError("accepted state digest changed across both arms")

    report = {
        "schema": "bssn-areal-fixed-R-Rr-Rt-ablation/v1",
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
            "primitive_path": "fixed four-point cubic stencil anchored to accepted state for all times; derive R_r from b(r),X(r) interpolation; derive R_t from b_dot,X_dot; recompute chi and T_t^R from sampled stress operands",
            "observer_control": "observer R_r and R_t interpolated independently at the same areal radius using the observer's cubic interpolation",
            "component_arms": list(VARIANTS.keys()),
            "scope": "only R_r and R_t source is swapped; all other sampled metric and matter operands remain on the primitive interpolation path",
        },
        "arms": arms,
        "interpretation_guard": [
            "This is a discrete-operator ablation, not a source formula selection or conservation proof.",
            "The observer and primitive-derived R_r/R_t are competing numerical reconstructions; neither is declared physical authority here.",
            "The two-step centred derivative is located at t+h; the forward derivative and its endpoint-average current are located on [t,t+h].",
            "The all-primitive versus one-operator swap differences include nonlinear coupling in chi/current; the coupling ledger is reported explicitly.",
            "No production equation, source normalization, initializer, vendor, dependency, or scientific admission flag is changed.",
        ],
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_input_immutable": True,
            "forward_shadow_states_immutable": True,
            "primitive_metric_derivatives_recomputed": True,
            "observer_Rr_Rt_used_as_separate_controls": True,
            "component_ablation_complete": True,
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
        "scientific_admission": False,
        "arms": [{
            "arm": arm["arm"],
            "variants": arm["variants"],
            "component_nonlinearity": arm["component_nonlinearity"],
        } for arm in arms],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
