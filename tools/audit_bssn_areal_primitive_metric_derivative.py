#!/usr/bin/env python3
"""Research-only fixed-R reconstruction from primitive metric derivatives.

The earlier operand-wise diagnostic interpolated the observer's R_r and R_t.
This probe goes one step further: reconstruct b(r), X(r) from regular even
profiles, differentiate that interpolation analytically to obtain R_r, form
R_t from interpolated b_dot and X_dot, and then recompute chi and T_t^R from
the stress operands. It compares this primitive-metric path with the existing
observer's composite and operand-wise cubic interpolation.

The test uses the same SHA-verified accepted t=0.60 state and copied forward
steps (dt and dt/2); it does not mutate or step the saved state backward.
The analytic derivative of the local interpolant is a distinct discrete
operator, not presumed to be the physical authority. No production code or
normalization changes are made.
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
N = int(os.environ.get("AREAL_PRIMITIVE_DERIVATIVE_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_PRIMITIVE_DERIVATIVE_OUT", f"results/bssn-areal-primitive-derivative-N{N}-2026-10-25"))
TARGET_TIME = 0.60
CFL = 0.03
FOUR_PI = 4.0 * np.pi
FIXED_RADII = (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0)
METHODS = (
    ("observer_composite_cubic", "composite", "cubic"),
    ("observer_operand_cubic", "observer_operand", "cubic"),
    ("raw_metric_cubic", "raw_metric", "cubic"),
    ("raw_metric_cubic_fixed_stencil", "raw_metric_fixed_stencil", "cubic"),
    ("raw_metric_linear_control", "raw_metric", "linear"),
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


def lagrange_value_and_derivative(nodes: np.ndarray, values: np.ndarray, x: float) -> tuple[float, float]:
    nodes = finite("Lagrange nodes", nodes)
    values = finite("Lagrange values", values)
    if nodes.ndim != 1 or values.shape != nodes.shape or nodes.size != 4 or np.any(np.diff(nodes) <= 0):
        raise RuntimeError("primitive cubic interpolation requires four strictly increasing nodes")
    weights = np.ones(4, dtype=np.float64)
    deriv_weights = np.zeros(4, dtype=np.float64)
    for k in range(4):
        for ell in range(4):
            if ell != k:
                weights[k] *= (x - nodes[ell]) / (nodes[k] - nodes[ell])
        for m in range(4):
            if m == k:
                continue
            term = 1.0 / (nodes[k] - nodes[m])
            for ell in range(4):
                if ell != k and ell != m:
                    term *= (x - nodes[ell]) / (nodes[k] - nodes[ell])
            deriv_weights[k] += term
    return float(weights @ values), float(deriv_weights @ values)


def primitive_even(profile: dict[str, Any], field: str, r: float, order: str, stencil_start: int | None = None) -> tuple[float, float]:
    """Interpolate an even field in s=r^2; return value and d(field)/dr."""
    rs = profile["r"]
    nodes = rs * rs
    values = finite(f"even field {field}", profile[field])
    target = float(r) ** 2
    if order == "cubic":
        idx = int(np.searchsorted(nodes, target))
        start = min(max(idx - 2, 0), nodes.size - 4) if stencil_start is None else min(max(int(stencil_start), 0), nodes.size - 4)
        value, deriv_s = lagrange_value_and_derivative(nodes[start:start+4], values[start:start+4], target)
    elif order == "linear":
        value = float(np.interp(target, nodes, values))
        idx = int(np.searchsorted(nodes, target, side="right") - 1)
        idx = min(max(idx, 0), nodes.size - 2)
        deriv_s = float((values[idx+1] - values[idx]) / (nodes[idx+1] - nodes[idx]))
    else:
        raise ValueError(f"unsupported interpolation order {order}")
    deriv_r = 2.0 * float(r) * deriv_s
    if not np.isfinite(value) or not np.isfinite(deriv_r):
        raise RuntimeError(f"non-finite interpolated {field} at r={r}")
    return value, deriv_r


def primitive_odd(profile: dict[str, Any], field: str, r: float, order: str, stencil_start: int | None = None) -> float:
    """Interpolate an odd field as r times a regular even profile."""
    rs = profile["r"]
    regular = finite(f"regularized odd field {field}", profile[field] / rs)
    nodes = rs * rs
    target = float(r) ** 2
    if order == "cubic":
        idx = int(np.searchsorted(nodes, target))
        start = min(max(idx - 2, 0), nodes.size - 4) if stencil_start is None else min(max(int(stencil_start), 0), nodes.size - 4)
        value, _ = lagrange_value_and_derivative(nodes[start:start+4], regular[start:start+4], target)
    elif order == "linear":
        value = float(np.interp(target, nodes, regular))
    else:
        raise ValueError(f"unsupported interpolation order {order}")
    return scalar_finite(f"interpolated odd field {field}", float(r) * value)


def build_primitive_profile(kernel, profile: dict[str, Any]) -> dict[str, Any]:
    q = profile["state"]
    grid = q.grid
    r = finite("coordinate radial centers", grid.centers)
    if r.size < 4 or np.any(r <= 0) or np.any(np.diff(r) <= 0):
        raise RuntimeError("primitive derivative observer requires at least four positive radial centers")
    g = q.geometry
    terms = kernel._geometry_stage_terms(grid, g, q.scalars, q.matter, q.electromagnetic)
    explicit = terms["explicit"]
    projected = operand_audit.operand_arrays(profile)
    primitive = {
        "r": r,
        "a": finite("a", g.a),
        "b": finite("b", g.b),
        "X": finite("X", g.X),
        "alpha": finite("alpha", g.alpha),
        "b_dot": finite("b RHS", explicit["b"]),
        "X_dot": finite("X RHS", explicit["X"]),
        "E": finite("Eulerian energy", projected["E"]),
        "j": finite("covariant radial current", projected["j"]),
        "p_r": finite("radial pressure", projected["pr"]),
        "observer_operands": projected,
        "observer_R": finite("observer R", profile["R"]),
        "observer_Rr": finite("observer Rr", profile["Rr"]),
        "observer_Rt": finite("observer Rt", profile["Rt"]),
        "observer_chi": finite("observer chi", profile["chi"]),
        "observer_F": finite("observer F", profile["rhs"]),
        "state_time": float(q.t),
        "max_abs_beta": float(profile["max_abs_beta"]),
    }
    for name in ("a", "b", "X", "alpha"):
        if np.any(primitive[name] <= 0):
            raise RuntimeError(f"primitive geometry field {name} is not positive")
    # The regular even interpolation must reproduce the nodewise areal map.
    R_nodes = primitive["r"] * np.sqrt(primitive["b"]) / primitive["X"]
    map_gap = finite("native areal-radius map gap", R_nodes - primitive["observer_R"])
    max_map_gap = float(np.max(np.abs(map_gap)))
    if max_map_gap > 2e-11:
        raise RuntimeError(f"primitive areal map fails native-grid reconstruction: {max_map_gap}")
    primitive["native_areal_map_max_abs_gap"] = max_map_gap
    return primitive


def evaluate_areal_at_coordinate(profile: dict[str, Any], r: float, order: str, stencil_start: int | None = None) -> tuple[float, float, dict[str, float]]:
    b, b_r = primitive_even(profile, "b", r, order, stencil_start)
    X, X_r = primitive_even(profile, "X", r, order, stencil_start)
    if b <= 0 or X <= 0:
        raise RuntimeError("interpolated b or X is not positive")
    R = float(r) * np.sqrt(b) / X
    R_r = R * (1.0 / max(float(r), 1e-300) + 0.5 * b_r / b - X_r / X)
    if r == 0.0:
        R_r = np.sqrt(b) / X
    if not np.isfinite(R) or not np.isfinite(R_r) or R_r <= 0:
        raise RuntimeError(f"primitive areal map is non-finite or not increasing at r={r}: Rr={R_r}")
    return R, R_r, {"b": b, "b_r": b_r, "X": X, "X_r": X_r}


def invert_areal_radius(
    profile: dict[str, Any],
    R_target: float,
    order: str,
    stencil_start: int | None = None,
    seed: float | None = None,
) -> float:
    """Invert R(r); a supplied stencil is held fixed throughout the inversion."""
    R_nodes = profile["observer_R"]
    centers = profile["r"]
    if R_target < R_nodes[0] - 1e-12 or R_target > R_nodes[-1] + 1e-12:
        raise RuntimeError(f"fixed R={R_target} is outside primitive center-map range")
    if stencil_start is None:
        right = int(np.searchsorted(R_nodes, R_target, side="left"))
        if right == 0:
            lo, hi = 0.0, float(centers[0])
        elif right >= len(centers):
            lo, hi = float(centers[-2]), float(centers[-1])
        else:
            lo, hi = float(centers[right-1]), float(centers[right])
        Rlo = evaluate_areal_at_coordinate(profile, lo, order)[0]
        Rhi = evaluate_areal_at_coordinate(profile, hi, order)[0]
        if Rlo > R_target + 1e-10 or Rhi < R_target - 1e-10:
            raise RuntimeError(f"primitive map does not bracket R={R_target}: [{Rlo},{Rhi}]")
    else:
        center = float(np.interp(R_target, R_nodes, centers)) if seed is None else float(seed)
        span = float(centers[1] - centers[0]) * 1.25
        lo = max(0.0, center - span)
        hi = min(float(centers[-1]), center + span)
        bracketed = False
        for _ in range(8):
            Rlo = evaluate_areal_at_coordinate(profile, lo, order, stencil_start)[0]
            Rhi = evaluate_areal_at_coordinate(profile, hi, order, stencil_start)[0]
            if Rlo <= R_target <= Rhi:
                bracketed = True
                break
            span *= 2.0
            lo = max(0.0, center - span)
            hi = min(float(centers[-1]), center + span)
        if not bracketed:
            raise RuntimeError(f"fixed-stencil primitive map did not bracket R={R_target}: [{Rlo},{Rhi}]")
    for _ in range(48):
        mid = 0.5 * (lo + hi)
        Rmid = evaluate_areal_at_coordinate(profile, mid, order, stencil_start)[0]
        if Rmid < R_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def primitive_metric_query(
    profile: dict[str, Any],
    R_target: float,
    order: str,
    stencil_start: int | None = None,
    seed: float | None = None,
) -> dict[str, float]:
    r = invert_areal_radius(profile, R_target, order, stencil_start, seed)
    R_actual, R_r, geo = evaluate_areal_at_coordinate(profile, r, order, stencil_start)
    a, _ = primitive_even(profile, "a", r, order, stencil_start)
    X, _ = primitive_even(profile, "X", r, order, stencil_start)
    alpha, _ = primitive_even(profile, "alpha", r, order, stencil_start)
    b, _ = primitive_even(profile, "b", r, order, stencil_start)
    b_dot, _ = primitive_even(profile, "b_dot", r, order, stencil_start)
    X_dot, _ = primitive_even(profile, "X_dot", r, order, stencil_start)
    E, _ = primitive_even(profile, "E", r, order, stencil_start)
    p_r, _ = primitive_even(profile, "p_r", r, order, stencil_start)
    j = primitive_odd(profile, "j", r, order, stencil_start)

    R_t = R_target * (0.5 * b_dot / b - X_dot / X)
    g_rr = a / (X * X)
    chi = (X * X / a) * R_r * R_r - (R_t / alpha) ** 2
    j_up = j / g_rr
    T_t_R = -E * R_t - alpha * j_up * R_r - (R_t / R_r) * (j / alpha * R_t + p_r * R_r)
    F = FOUR_PI * R_target * R_target * T_t_R
    return {
        "coordinate_r": scalar_finite("inverted coordinate radius", r),
        "R_reconstructed": scalar_finite("reconstructed areal radius", R_actual),
        "R_r_primitive_derivative": scalar_finite("primitive-derived R_r", R_r),
        "R_t_primitive_rate": scalar_finite("primitive-derived R_t", R_t),
        "chi": scalar_finite("primitive-derived chi", chi),
        "current": scalar_finite("primitive-derived fixed-R current", F),
        "a": scalar_finite("a sampled", a),
        "b": scalar_finite("b sampled", b),
        "X": scalar_finite("X sampled", X),
        "alpha": scalar_finite("alpha sampled", alpha),
        "E": scalar_finite("E sampled", E),
        "j": scalar_finite("j sampled", j),
        "p_r": scalar_finite("p_r sampled", p_r),
        "R_r_observer_interp": scalar_finite("observer Rr for discrepancy", np.interp(R_target, profile["observer_R"], profile["observer_Rr"])),
        "R_t_observer_interp": scalar_finite("observer Rt for discrepancy", np.interp(R_target, profile["observer_R"], profile["observer_Rt"])),
    }


def observer_query(profile: dict[str, Any], primitive: dict[str, Any], R_target: float, source: str, order: str) -> dict[str, float]:
    if source == "composite":
        if order == "cubic":
            chi = areal.interp4(R_target, profile["R"], profile["chi"])
            current = areal.interp4(R_target, profile["R"], profile["rhs"])
        else:
            chi = np.interp(R_target, profile["R"], profile["chi"])
            current = np.interp(R_target, profile["R"], profile["rhs"])
        return {"chi": scalar_finite("composite chi", chi), "current": scalar_finite("composite current", current)}
    if source == "observer_operand":
        # Existing operand audit interpolates the pinned observer's R_r and R_t,
        # rather than differentiating the primitive-metric interpolation.
        return operand_audit.query_operand(primitive["observer_operands"], R_target, order)
    raise ValueError(f"unknown source {source}")


def stencil_start_for_coordinate(profile: dict[str, Any], r: float) -> int:
    nodes = profile["r"] * profile["r"]
    idx = int(np.searchsorted(nodes, float(r) ** 2))
    return min(max(idx - 2, 0), nodes.size - 4)


def queries_for_method(
    profiles: tuple[dict[str, Any], dict[str, Any], dict[str, Any]],
    primitives: tuple[dict[str, Any], dict[str, Any], dict[str, Any]],
    R: float,
    method: str,
) -> list[dict[str, float]]:
    _, source, order = next(spec for spec in METHODS if spec[0] == method)
    if source == "raw_metric":
        return [primitive_metric_query(primitives[j], R, order) for j in range(3)]
    if source == "raw_metric_fixed_stencil":
        # Anchor the four-point interpolation stencil to the accepted state and
        # reuse it for all three time slices, preventing time-varying stencil hops.
        anchor_r = invert_areal_radius(primitives[0], R, "cubic")
        start = stencil_start_for_coordinate(primitives[0], anchor_r)
        return [
            primitive_metric_query(
                primitives[j],
                R,
                "cubic",
                stencil_start=start,
                seed=float(np.interp(R, primitives[j]["observer_R"], primitives[j]["r"])),
            )
            for j in range(3)
        ]
    return [observer_query(profiles[j], primitives[j], R, source, order) for j in range(3)]


def residual_summary(residual: np.ndarray, radii: np.ndarray, current: np.ndarray) -> dict[str, float]:
    residual = finite("mass-current residual", residual)
    current = finite("current scale", current)
    inner = radii <= 3.0
    annulus = (radii >= 2.0) & (radii <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("common samples do not cover both required regions")
    rms = float(np.sqrt(np.mean(residual[annulus] ** 2)))
    scale = float(np.sqrt(np.mean(current[annulus] ** 2)))
    return {
        "max_abs_R_le_3": float(np.max(np.abs(residual[inner]))),
        "rms_R_2_to_10": rms,
        "max_abs_R_2_to_10": float(np.max(np.abs(residual[annulus]))),
        "rms_relative_to_current_R_2_to_10": rms / max(scale, 1e-300),
    }


def temporal_arm(kernel, accepted, h: float) -> dict[str, Any]:
    accepted_digest = saved.tangent.state_digest(accepted)
    profile0 = areal.geometry_rate_and_mass(kernel, accepted)
    state1 = kernel.step(areal.clone_state(accepted), h)
    state1_digest = saved.tangent.state_digest(state1)
    profile1 = areal.geometry_rate_and_mass(kernel, state1)
    state2 = kernel.step(areal.clone_state(state1), h)
    state2_digest = saved.tangent.state_digest(state2)
    profile2 = areal.geometry_rate_and_mass(kernel, state2)
    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted state mutated during primitive derivative probe")
    if saved.tangent.state_digest(state1) != state1_digest:
        raise RuntimeError("first shadow state mutated during second step/observation")
    if saved.tangent.state_digest(state2) != state2_digest:
        raise RuntimeError("second shadow state mutated during observation")
    if abs(float(state1.t) - (float(accepted.t) + h)) > 1e-10:
        raise RuntimeError("first shadow step failed exact h advance")
    if abs(float(state2.t) - (float(accepted.t) + 2.0*h)) > 1e-10:
        raise RuntimeError("second shadow step failed exact h advance")

    p0 = build_primitive_profile(kernel, profile0)
    p1 = build_primitive_profile(kernel, profile1)
    p2 = build_primitive_profile(kernel, profile2)
    if max(p0["max_abs_beta"], p1["max_abs_beta"], p2["max_abs_beta"]) > 1e-14:
        raise RuntimeError("zero-shift precondition failed")

    lo = max(float(profile0["R"][0]), float(profile1["R"][0]), float(profile2["R"][0]))
    hi = min(float(profile0["R"][-1]), float(profile1["R"][-1]), float(profile2["R"][-1]), 10.0)
    if not hi > lo:
        raise RuntimeError("no common fixed-R interval")
    radii = np.linspace(lo, hi, max(300, int(accepted.grid.n) * 2))
    profiles = (profile0, profile1, profile2)
    primitives = (p0, p1, p2)
    sampled: dict[str, dict[str, list[np.ndarray]]] = {}
    for method, _, _ in METHODS:
        vals: dict[str, list[np.ndarray]] = {
            "chi": [np.empty_like(radii) for _ in range(3)],
            "current": [np.empty_like(radii) for _ in range(3)],
        }
        for i, R in enumerate(radii):
            queries = queries_for_method(profiles, primitives, float(R), method)
            for j, q in enumerate(queries):
                vals["chi"][j][i] = q["chi"]
                vals["current"][j][i] = q["current"]
        chi, current = vals["chi"], vals["current"]
        forward_rate = finite(f"{method}.forward rate", -0.5 * radii * (chi[1] - chi[0]) / h)
        centered_rate = finite(f"{method}.centered rate", -0.5 * radii * (chi[2] - chi[0]) / (2.0*h))
        forward_avg = finite(f"{method}.forward current average", 0.5 * (current[0] + current[1]))
        centered_avg = finite(f"{method}.centered current average", (current[0] + 2.0*current[1] + current[2]) / 4.0)
        sampled[method] = {
            "chi": chi,
            "current": current,
            "forward_rate": [forward_rate],
            "centered_rate": [centered_rate],
            "forward_avg": [forward_avg],
            "centered_avg": [centered_avg],
            "forward_residual_avg": [finite(f"{method}.forward residual", forward_rate - forward_avg)],
            "centered_residual_mid": [finite(f"{method}.centered midpoint residual", centered_rate - current[1])],
            "centered_residual_avg": [finite(f"{method}.centered avg residual", centered_rate - centered_avg)],
            "rate_gap": [finite(f"{method}.rate gap", centered_rate - forward_rate)],
        }

    method_summary: dict[str, Any] = {}
    for method, _, _ in METHODS:
        s = sampled[method]
        method_summary[method] = {
            "forward_vs_endpoint_average": residual_summary(s["forward_residual_avg"][0], radii, s["forward_avg"][0]),
            "centered_vs_midpoint_current": residual_summary(s["centered_residual_mid"][0], radii, s["current"][1]),
            "centered_vs_two_step_average": residual_summary(s["centered_residual_avg"][0], radii, s["centered_avg"][0]),
            "centered_minus_forward_rate": residual_summary(s["rate_gap"][0], radii, s["centered_rate"][0]),
            "fixed_R": {},
        }
        for R in FIXED_RADII:
            if not lo <= R <= hi:
                continue
            qs = queries_for_method(profiles, primitives, R, method)
            forward_rate = -0.5 * R * (qs[1]["chi"] - qs[0]["chi"]) / h
            centered_rate = -0.5 * R * (qs[2]["chi"] - qs[0]["chi"]) / (2.0*h)
            forward_avg = 0.5 * (qs[0]["current"] + qs[1]["current"])
            centered_avg = (qs[0]["current"] + 2.0*qs[1]["current"] + qs[2]["current"]) / 4.0
            row = {
                "forward_mass_rate": scalar_finite("fixed-R forward mass rate", forward_rate),
                "forward_endpoint_average_current": scalar_finite("fixed-R forward avg current", forward_avg),
                "forward_residual": scalar_finite("fixed-R forward residual", forward_rate - forward_avg),
                "centered_mass_rate": scalar_finite("fixed-R centred mass rate", centered_rate),
                "centered_midpoint_current": scalar_finite("fixed-R midpoint current", qs[1]["current"]),
                "centered_two_step_average_current": scalar_finite("fixed-R two-step current average", centered_avg),
                "centered_residual_vs_midpoint": scalar_finite("fixed-R centred midpoint residual", centered_rate - qs[1]["current"]),
                "centered_residual_vs_average": scalar_finite("fixed-R centred average residual", centered_rate - centered_avg),
            }
            if method.startswith("raw_metric"):
                row.update({
                    "R_r_primitive_at_midpoint": qs[1]["R_r_primitive_derivative"],
                    "R_r_observer_interpolated_at_midpoint": qs[1]["R_r_observer_interp"],
                    "R_r_gap_primitive_minus_observer": qs[1]["R_r_primitive_derivative"] - qs[1]["R_r_observer_interp"],
                    "R_t_primitive_at_midpoint": qs[1]["R_t_primitive_rate"],
                    "R_t_observer_interpolated_at_midpoint": qs[1]["R_t_observer_interp"],
                    "R_t_gap_primitive_minus_observer": qs[1]["R_t_primitive_rate"] - qs[1]["R_t_observer_interp"],
                    "coordinate_r_midpoint": qs[1]["coordinate_r"],
                })
            method_summary[method]["fixed_R"][str(R)] = row

    direct_comparisons = {}
    for raw_name in ("raw_metric_cubic", "raw_metric_cubic_fixed_stencil"):
        raw = sampled[raw_name]
        controls = ["observer_composite_cubic", "observer_operand_cubic"]
        controls.append("raw_metric_cubic" if raw_name.endswith("fixed_stencil") else "raw_metric_cubic_fixed_stencil")
        for control in controls:
            if control == raw_name:
                continue
            base = sampled[control]
            for key in ("chi", "current", "forward_rate", "centered_rate", "forward_residual_avg", "centered_residual_mid", "centered_residual_avg"):
                if key in ("chi", "current"):
                    a = raw[key][1]
                    b = base[key][1]
                else:
                    a = raw[key][0]
                    b = base[key][0]
                delta = finite(f"{raw_name} minus {control} {key}", a - b)
                annulus = (radii >= 2.0) & (radii <= 10.0)
                direct_comparisons[f"{raw_name}_minus_{control}_{key}"] = {
                    "rms_R_2_to_10": float(np.sqrt(np.mean(delta[annulus]**2))),
                    "max_abs_R_2_to_10": float(np.max(np.abs(delta[annulus]))),
                }

    return {
        "step_h": float(h),
        "base_time": float(accepted.t),
        "first_shadow_time": float(state1.t),
        "second_shadow_time": float(state2.t),
        "accepted_input_immutable": True,
        "first_shadow_state_immutable": True,
        "second_shadow_state_immutable": True,
        "common_areal_radius_interval": [float(lo), float(hi)],
        "common_samples": int(radii.size),
        "native_areal_map_max_abs_gap": max(p["native_areal_map_max_abs_gap"] for p in primitives),
        "methods": method_summary,
        "raw_metric_vs_observer_comparisons": direct_comparisons,
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
        raise RuntimeError("accepted input is not t=0.60")
    digest = saved.tangent.state_digest(accepted)
    dt = CFL * float(accepted.grid.dr)
    arms = [{"arm": label, **temporal_arm(kernel, accepted, h)} for label, h in (("dt", dt), ("dt_half", dt/2.0))]
    if saved.tangent.state_digest(accepted) != digest:
        raise RuntimeError("accepted input digest changed across both arms")
    report = {
        "schema": "bssn-areal-primitive-metric-derivative/v1",
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
            "raw_metric_path": "invert R=r*sqrt(b)/X by bracketing; interpolate even primitive fields in r^2; derive R_r by differentiating the same b/X interpolants; derive R_t from b_dot/X_dot; interpolate j as r times an even profile; then recompute chi and F from complete sampled stress operands",
            "controls": "existing composite cubic, observer operand cubic (existing pinned R_r/R_t), and primitive-metric linear interpolation control",
            "no_backward_steps": True,
            "accepted_state_never_mutated": True,
        },
        "arms": arms,
        "interpretation_guard": [
            "Primitive-derived R_r is a different discrete derivative operator from the pinned observer's fourth-order cell derivative; neither is declared the physical authority by this experiment.",
            "Reconstructing the exact same continuum expression with a different spatial operator is a discriminator, not a source change or a conservation proof.",
            "Linear primitive interpolation is a deliberately lower-order control; its derivative need not agree with cubic interpolation.",
            "Rates at t+h/2 and t+h are located at different times, so their difference is reported as a sampling measure, not as a formal truncation order.",
            "No production equation, initializer, vendor, dependency, source normalization, or scientific admission flag was changed.",
        ],
        "acceptance": {
            "source_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_input_immutable": True,
            "forward_shadow_states_immutable": True,
            "primitive_areal_map_reconstructs_native_nodes": True,
            "raw_R_r_derived_from_primitive_metric_interpolants": True,
            "stress_operands_reassembled_after_interpolation": True,
            "linear_low_order_control_included": True,
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
                    "forward_RMS": res["forward_vs_endpoint_average"]["rms_R_2_to_10"],
                    "centered_midpoint_RMS": res["centered_vs_midpoint_current"]["rms_R_2_to_10"],
                    "centered_average_RMS": res["centered_vs_two_step_average"]["rms_R_2_to_10"],
                    "rate_gap_RMS": res["centered_minus_forward_rate"]["rms_R_2_to_10"],
                } for name, res in arm["methods"].items()
            },
            "raw_metric_vs_observer": arm["raw_metric_vs_observer_comparisons"],
            "fixed_R_raw_metric_cubic": arm["methods"]["raw_metric_cubic"]["fixed_R"],
            "fixed_R_raw_metric_cubic_fixed_stencil": arm["methods"]["raw_metric_cubic_fixed_stencil"]["fixed_R"],
        } for arm in arms],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
