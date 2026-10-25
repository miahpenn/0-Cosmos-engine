#!/usr/bin/env python3
"""Research-only fixed-state compatibility audit for the areal-map radial derivative.

Compare the native-grid tangent D_r(R_t) of R_r=D_r(R) with the time
tangent of a fixed-stencil primitive b/X interpolant for R_r. The second
tangent is assembled analytically from the same interpolants of b, X, b_dot,
and X_dot. Then hold the transformed matter operands fixed and split the
change in the fixed-R Misner-Sharp mass tangent into:
  (i) the spatial R_r reconstruction change, including its induced chi/current;
  (ii) the R_r time-tangent change alone.

This is a diagnostic/control, not a selection of physical authority. It uses
only a SHA-verified accepted state at t=0.60; no production code is changed.
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

from engine.v55_matter import total_matter_projection
from tools import audit_bssn_areal_misner_sharp_current as areal
from tools import audit_bssn_areal_primitive_metric_derivative as primitive
from tools import audit_full_profile_propagation_candidate as saved
from tools import audit_origin_cell_joint_shadow_one_step as origin_shadow

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
N = int(os.environ.get("AREAL_RR_TANGENT_RESOLUTION", "120"))
FULL_DIR = Path(os.environ.get("FULL_PROFILE_ARTIFACT_DIR", f"reference/full/N{N}"))
OUT = Path(os.environ.get("AREAL_RR_TANGENT_OUT", f"results/bssn-areal-Rr-time-tangent-N{N}-2026-10-25"))
TARGET_TIME = 0.60
FOUR_PI = 4.0 * np.pi


def finite(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(f"{name} contains non-finite values")
    return arr


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def summary(values: np.ndarray, fixed_R: np.ndarray) -> dict[str, float]:
    values = finite("summary values", values)
    inner = fixed_R <= 3.0
    annulus = (fixed_R >= 2.0) & (fixed_R <= 10.0)
    if not inner.any() or not annulus.any():
        raise RuntimeError("common fixed-R samples miss required inner/annular region")
    return {
        "max_abs_R_le_3": float(np.max(np.abs(values[inner]))),
        "rms_R_2_to_10": float(np.sqrt(np.mean(values[annulus] ** 2))),
        "max_abs_R_2_to_10": float(np.max(np.abs(values[annulus]))),
    }


def fixed_R_sample(native_R: np.ndarray, values: np.ndarray, fixed_R: np.ndarray, name: str) -> np.ndarray:
    native_R = finite("native areal coordinates", native_R)
    values = finite(name, values)
    return finite(
        f"fixed-R samples of {name}",
        np.asarray([areal.interp4(float(target), native_R, values) for target in fixed_R], dtype=np.float64),
    )


def main() -> int:
    if N not in (60, 120, 240):
        raise ValueError("resolution must be 60, 120 or 240")
    pins = origin_shadow.verify_source_pins()
    if pins.get("source_commit") != SOURCE or pins.get("vendor_commit") != VENDOR:
        raise RuntimeError("source/vendor pin mismatch")
    if pins.get("protected_engine_vendor_dependency_diff_empty") is not True:
        raise RuntimeError("protected production/vendor/dependency diff is not empty")

    manifest = saved.check_artifact(FULL_DIR, N, full_profile=True)
    arrays = saved.load_profile_arrays(FULL_DIR)
    prefix = f"N{N}_full_continuous_metric_profile_final_"
    kernel, accepted = saved.reconstruct_state(N, arrays, prefix)
    origin_shadow._install_regular_row(kernel)
    if abs(float(accepted.t) - TARGET_TIME) > 1e-12:
        raise RuntimeError("accepted input is not the pinned t=0.60 state")
    accepted_digest = saved.tangent.state_digest(accepted)

    observed = areal.geometry_rate_and_mass(kernel, accepted)
    q = observed["state"]
    grid = q.grid
    g = q.geometry
    r = finite("coordinate centers", grid.centers)
    R = finite("areal radius", observed["R"])
    Rt = finite("native areal time rate", observed["Rt"])
    Rr_native = finite("native D4 areal radial derivative", observed["Rr"])
    if np.any(np.diff(R) <= 0.0) or np.any(Rr_native <= 0.0):
        raise RuntimeError("native areal map is not strictly monotone")

    # Build the primitive b/X profile and evaluate a *fixed coordinate-centred*
    # cubic stencil at each native cell. Stencils stay fixed while differentiating
    # the interpolants in time, so the rate comparison does not include stencil reselection.
    prim = primitive.build_primitive_profile(kernel, observed)
    n = r.size
    if n < 4:
        raise RuntimeError("at least four native centers are required")
    Rr_primitive = np.empty(n, dtype=np.float64)
    Rr_t_primitive = np.empty(n, dtype=np.float64)
    Rt_primitive = np.empty(n, dtype=np.float64)
    R_primitive = np.empty(n, dtype=np.float64)
    for i, ri in enumerate(r):
        start = min(max(i - 2, 0), n - 4)
        b, b_r = primitive.primitive_even(prim, "b", float(ri), "cubic", start)
        X, X_r = primitive.primitive_even(prim, "X", float(ri), "cubic", start)
        bdot, bdot_r = primitive.primitive_even(prim, "b_dot", float(ri), "cubic", start)
        Xdot, Xdot_r = primitive.primitive_even(prim, "X_dot", float(ri), "cubic", start)
        if min(b, X) <= 0.0:
            raise RuntimeError(f"nonpositive primitive metric factor at cell {i}")
        Ri = float(ri) * np.sqrt(b) / X
        Bi = 1.0 / float(ri) + 0.5 * b_r / b - X_r / X
        Rti = Ri * (0.5 * bdot / b - Xdot / X)
        Bti = (
            0.5 * (bdot_r / b - b_r * bdot / (b * b))
            - (Xdot_r / X - X_r * Xdot / (X * X))
        )
        R_primitive[i] = Ri
        Rt_primitive[i] = Rti
        Rr_primitive[i] = Ri * Bi
        Rr_t_primitive[i] = Rti * Bi + Ri * Bti

    R_map_gap = finite("primitive minus native R map", R_primitive - R)
    Rt_map_gap = finite("primitive minus native Rt", Rt_primitive - Rt)
    Rr_gap = finite("primitive Rr minus native-grid Rr", Rr_primitive - Rr_native)
    Rr_t_native = finite(
        "native-grid Rr time tangent D4(Rt)",
        grid.cell_derivative_fourth(Rt, parity=-1),
    )
    Rr_t_gap = finite(
        "primitive time tangent of Rr minus D4(Rt)",
        Rr_t_primitive - Rr_t_native,
    )

    terms = kernel._geometry_stage_terms(grid, g, q.scalars, q.matter, q.electromagnetic)
    explicit = terms["explicit"]
    l2, l3 = terms["primary_l2"], terms["primary_l3"]
    adot = finite("a RHS", explicit["a"])
    bdot_native = finite("b RHS", explicit["b"])
    Xdot_native = finite("X RHS", explicit["X"])
    Kdot = finite("K RHS", finite("K L2", l2["K"]) + finite("K L3", l3["K"]))
    Aadot = finite("Aa RHS", finite("Aa L2", l2["Aa"]) + finite("Aa L3", l3["Aa"]))
    a = finite("a", g.a)
    X = finite("X", g.X)
    alpha = finite("lapse", g.alpha)
    gamma_inv = finite("gamma^rr", X * X / a)
    gamma_inv_dot = finite("gamma^rr tangent", gamma_inv * (2.0 * Xdot_native / X - adot / a))
    Ktheta = finite("Ktheta", g.K / 3.0 - g.Aa / 2.0)
    Ktheta_dot = finite("Ktheta tangent", Kdot / 3.0 - Aadot / 2.0)

    projection = total_matter_projection(grid, g, q.scalars, q.matter)
    rho = finite("Eulerian energy", projection["rho"])
    j = finite("covariant radial current", projection["j"])
    pr = finite("radial pressure", projection["pr"])

    def variant(Rr: np.ndarray, Rr_t: np.ndarray, label: str) -> dict[str, np.ndarray]:
        chi = finite(f"{label} chi", gamma_inv * Rr * Rr - (R * Ktheta) ** 2)
        chi_dot = finite(
            f"{label} chi tangent",
            gamma_inv_dot * Rr * Rr
            + 2.0 * gamma_inv * Rr * Rr_t
            - 2.0 * R * Rt * Ktheta * Ktheta
            - 2.0 * R * R * Ktheta * Ktheta_dot,
        )
        # Use the same native-grid D4 operator for the spatial derivative of each
        # reconstructed chi field. This keeps the control bounded to Rr/Rr_t.
        chi_r = finite(f"{label} D4(chi)", grid.cell_derivative_fourth(chi, parity=1))
        mdot_fixed_R = finite(
            f"{label} instantaneous fixed-R mass tangent",
            -0.5 * R * (chi_dot - (Rt / Rr) * chi_r),
        )
        # Recompute the transformed current from the same matter operands.
        TtR = finite(
            f"{label} transformed current",
            -rho * Rt
            - alpha * gamma_inv * j * Rr
            - (Rt / Rr) * (j * Rt / alpha + pr * Rr),
        )
        current = finite(f"{label} current RHS", FOUR_PI * R * R * TtR)
        return {"Rr": Rr, "Rr_t": Rr_t, "chi": chi, "chi_dot": chi_dot,
                "chi_r": chi_r, "mass_rate": mdot_fixed_R, "current": current,
                "residual": finite(f"{label} residual", mdot_fixed_R - current)}

    baseline = variant(Rr_native, Rr_t_native, "native_D4_Rr")
    hybrid = variant(Rr_primitive, Rr_t_native, "primitive_Rr_native_tangent_control")
    primitive_matched = variant(Rr_primitive, Rr_t_primitive, "primitive_Rr_matched_tangent")

    # Isolate the Rr time-tangent contribution with Rr, chi, chi_r and current
    # held fixed between the hybrid and matched primitive variants.
    expected_rate_delta = finite(
        "analytic mass-rate delta from Rr tangent alone",
        -R * gamma_inv * Rr_primitive * (Rr_t_primitive - Rr_t_native),
    )
    measured_rate_delta = finite(
        "measured mass-rate delta from Rr tangent alone",
        primitive_matched["mass_rate"] - hybrid["mass_rate"],
    )
    tangent_component_closure = finite(
        "Rr tangent contribution closure", measured_rate_delta - expected_rate_delta
    )
    current_change = finite("primitive minus native current", primitive_matched["current"] - baseline["current"])
    mass_rate_change = finite("primitive minus native fixed-R mass rate", primitive_matched["mass_rate"] - baseline["mass_rate"])
    residual_change = finite("primitive minus native residual", primitive_matched["residual"] - baseline["residual"])
    delta_closure = finite("mass-rate/current delta closure", residual_change - mass_rate_change + current_change)

    low = max(float(R[0]), 0.5)
    high = min(float(R[-1]), 10.0)
    if not high > low:
        raise RuntimeError("accepted state has no common fixed-R analysis interval")
    fixed_R = np.linspace(low, high, max(300, 2 * N))
    metrics = {
        "primitive_R_map_minus_native": summary(fixed_R_sample(R, R_map_gap, fixed_R, "R_map_gap"), fixed_R),
        "primitive_Rt_minus_native": summary(fixed_R_sample(R, Rt_map_gap, fixed_R, "Rt_map_gap"), fixed_R),
        "primitive_Rr_minus_native": summary(fixed_R_sample(R, Rr_gap, fixed_R, "Rr_gap"), fixed_R),
        "primitive_Rr_tangent_minus_D4_Rt": summary(fixed_R_sample(R, Rr_t_gap, fixed_R, "Rr_t_gap"), fixed_R),
        "native_fixed_R_mass_rate_minus_current": summary(fixed_R_sample(R, baseline["residual"], fixed_R, "native residual"), fixed_R),
        "primitive_fixed_R_mass_rate_minus_current": summary(fixed_R_sample(R, primitive_matched["residual"], fixed_R, "primitive residual"), fixed_R),
        "delta_mass_rate_primitive_minus_native": summary(fixed_R_sample(R, mass_rate_change, fixed_R, "mass-rate change"), fixed_R),
        "delta_current_primitive_minus_native": summary(fixed_R_sample(R, current_change, fixed_R, "current change"), fixed_R),
        "delta_residual_primitive_minus_native": summary(fixed_R_sample(R, residual_change, fixed_R, "residual change"), fixed_R),
        "Rr_tangent_only_mass_rate_change": summary(fixed_R_sample(R, measured_rate_delta, fixed_R, "Rr tangent-only mass-rate change"), fixed_R),
        "Rr_tangent_only_closure_error": summary(fixed_R_sample(R, tangent_component_closure, fixed_R, "Rr tangent closure"), fixed_R),
        "delta_budget_closure_error": summary(fixed_R_sample(R, delta_closure, fixed_R, "delta closure"), fixed_R),
    }

    fixed_rows: dict[str, Any] = {}
    for radius in (1.0, 2.0, 3.0, 5.0, 8.0, 10.0):
        if low <= radius <= high:
            fixed_rows[str(radius)] = {}
            for label, data in (("native", baseline), ("primitive_Rr_native_tangent", hybrid),
                                ("primitive_Rr_matched_tangent", primitive_matched)):
                fixed_rows[str(radius)][label] = {
                    name: float(areal.interp4(radius, R, values))
                    for name, values in (
                        ("Rr", data["Rr"]), ("Rr_t", data["Rr_t"]),
                        ("mass_rate", data["mass_rate"]), ("current", data["current"]),
                        ("residual", data["residual"]),
                    )
                }

    # Repeat the tangent-only comparison using the exact fixed-R-anchored
    # stencil policy from audit_bssn_areal_Rr_operator_crosswalk.py. This avoids
    # choosing a different cubic stencil merely because a target falls near a cell.
    fixed_R_anchored_rows: dict[str, Any] = {}
    chi_native_grid = finite(
        "native chi from Rr and Ktheta",
        gamma_inv * Rr_native * Rr_native - (R * Ktheta) ** 2,
    )
    chi_r_native_grid = finite(
        "native D4(chi) for fixed-R transport",
        grid.cell_derivative_fourth(chi_native_grid, parity=1),
    )
    tangent_closure_abs: list[float] = []
    for target in (1.0, 2.0, 3.0, 5.0, 8.0, 10.0):
        if not low <= target <= high:
            continue
        anchor_r = primitive.invert_areal_radius(prim, float(target), "cubic")
        stencil_start = primitive.stencil_start_for_coordinate(prim, anchor_r)
        q = primitive.primitive_metric_query(
            prim,
            float(target),
            "cubic",
            stencil_start=stencil_start,
            seed=float(np.interp(float(target), prim["observer_R"], prim["r"])),
        )
        rq = float(q["coordinate_r"])
        bq, bq_r = primitive.primitive_even(prim, "b", rq, "cubic", stencil_start)
        Xq, Xq_r = primitive.primitive_even(prim, "X", rq, "cubic", stencil_start)
        bdq, bdq_r = primitive.primitive_even(prim, "b_dot", rq, "cubic", stencil_start)
        Xdq, Xdq_r = primitive.primitive_even(prim, "X_dot", rq, "cubic", stencil_start)
        if min(bq, Xq, rq) <= 0.0:
            raise RuntimeError(f"invalid anchored primitive metric at R={target}: r={rq}, b={bq}, X={Xq}")
        Rq = rq * np.sqrt(bq) / Xq
        Bq = 1.0 / rq + 0.5 * bq_r / bq - Xq_r / Xq
        Bdotq = (
            0.5 * (bdq_r / bq - bq_r * bdq / (bq * bq))
            - (Xdq_r / Xq - Xq_r * Xdq / (Xq * Xq))
        )
        Rtq = Rq * (0.5 * bdq / bq - Xdq / Xq)
        Rr_anchored = Rq * Bq
        Rr_dot_anchored = Rtq * Bq + Rq * Bdotq

        Rr_native_at_R = float(areal.interp4(float(target), R, Rr_native))
        Rr_dot_native_at_R = float(areal.interp4(float(target), R, Rr_t_native))
        Rt_native_at_R = float(areal.interp4(float(target), R, Rt))
        gamma_at_R = float(areal.interp4(float(target), R, gamma_inv))
        gamma_dot_at_R = float(areal.interp4(float(target), R, gamma_inv_dot))
        Ktheta_at_R = float(areal.interp4(float(target), R, Ktheta))
        Ktheta_dot_at_R = float(areal.interp4(float(target), R, Ktheta_dot))
        chi_r_at_R = float(areal.interp4(float(target), R, chi_r_native_grid))
        current_native_at_R = float(areal.interp4(float(target), R, observed["rhs"]))
        current_anchored = float(q["current"])

        # The transport derivative is held fixed during the tangent-only A/B
        # comparison. Its displayed primitive-Rr variant changes only the
        # 1/Rr transport coefficient; D4(chi) itself is the same measured field.
        def anchored_mass_terms(Rr_value: float, Rr_dot_value: float) -> dict[str, float]:
            chi_dot_value = (
                gamma_dot_at_R * Rr_value * Rr_value
                + 2.0 * gamma_at_R * Rr_value * Rr_dot_value
                - 2.0 * float(target) * Rt_native_at_R * Ktheta_at_R * Ktheta_at_R
                - 2.0 * float(target) * float(target) * Ktheta_at_R * Ktheta_dot_at_R
            )
            coordinate_tangent_term = -0.5 * float(target) * chi_dot_value
            fixed_R_transport_term = (
                0.5 * float(target) * Rt_native_at_R / Rr_value * chi_r_at_R
            )
            return {
                "chi_dot": float(chi_dot_value),
                "coordinate_tangent_term": float(coordinate_tangent_term),
                "fixed_R_transport_term": float(fixed_R_transport_term),
                "fixed_R_mass_rate": float(coordinate_tangent_term + fixed_R_transport_term),
            }

        native_terms = anchored_mass_terms(Rr_native_at_R, Rr_dot_native_at_R)
        primitive_Rr_native_tangent_terms = anchored_mass_terms(Rr_anchored, Rr_dot_native_at_R)
        primitive_matched_terms = anchored_mass_terms(Rr_anchored, Rr_dot_anchored)
        tangent_delta_measured = (
            primitive_matched_terms["fixed_R_mass_rate"]
            - primitive_Rr_native_tangent_terms["fixed_R_mass_rate"]
        )
        tangent_delta_expected = (
            -float(target) * gamma_at_R * Rr_anchored
            * (Rr_dot_anchored - Rr_dot_native_at_R)
        )
        tangent_delta_closure = tangent_delta_measured - tangent_delta_expected
        tangent_closure_abs.append(abs(tangent_delta_closure))
        fixed_R_anchored_rows[str(target)] = {
            "coordinate_r_anchor": float(anchor_r),
            "coordinate_r_query": rq,
            "stencil_start": int(stencil_start),
            "reconstructed_R_minus_target": float(Rq - float(target)),
            "Rr_native_D4": Rr_native_at_R,
            "Rr_primitive_fixed_R_stencil": float(Rr_anchored),
            "delta_Rr_primitive_minus_native": float(Rr_anchored - Rr_native_at_R),
            "Rr_dot_native_D4_Rt": Rr_dot_native_at_R,
            "Rr_dot_primitive_fixed_R_stencil": float(Rr_dot_anchored),
            "delta_Rr_dot_primitive_minus_native": float(Rr_dot_anchored - Rr_dot_native_at_R),
            "native": {
                **native_terms,
                "current": current_native_at_R,
                "residual": float(native_terms["fixed_R_mass_rate"] - current_native_at_R),
            },
            "primitive_Rr_native_tangent_control": {
                **primitive_Rr_native_tangent_terms,
                "current": current_anchored,
                "residual": float(primitive_Rr_native_tangent_terms["fixed_R_mass_rate"] - current_anchored),
            },
            "primitive_Rr_matched_tangent": {
                **primitive_matched_terms,
                "current": current_anchored,
                "residual": float(primitive_matched_terms["fixed_R_mass_rate"] - current_anchored),
            },
            "tangent_only_mass_rate_delta_measured": float(tangent_delta_measured),
            "tangent_only_mass_rate_delta_expected": float(tangent_delta_expected),
            "tangent_only_delta_closure_error": float(tangent_delta_closure),
            "current_primitive_minus_native": float(current_anchored - current_native_at_R),
            "transport_control_note": "D4(chi) is held fixed between the tangent-only pair; the primitive_Rr transport term uses that same chi_r with its own Rr denominator.",
        }

    # Explicitly expose the fixed-R evaluation-order commutator.  The existing
    # fixed_R lane evaluates each native-grid mass rate first and then interpolates
    # the composite fields; the anchored lane interpolates native operands to the
    # fixed-R sample and then reconstructs the mass rate.  Since this expression
    # is nonlinear, these are distinct measurements even when the current matches.
    # Keep both; do not use the commutator to select a physical derivative.
    native_sampling_order_commutator: dict[str, Any] = {}
    sampling_rate_deltas: list[float] = []
    sampling_current_deltas: list[float] = []
    sampling_residual_deltas: list[float] = []
    sampling_closure_errors: list[float] = []
    for target, anchored_row in fixed_R_anchored_rows.items():
        composite_row = fixed_rows[target]["native"]
        reconstructed_row = anchored_row["native"]
        delta_rate = float(reconstructed_row["fixed_R_mass_rate"] - composite_row["mass_rate"])
        delta_current = float(reconstructed_row["current"] - composite_row["current"])
        delta_residual = float(reconstructed_row["residual"] - composite_row["residual"])
        closure_error = float(delta_residual - delta_rate + delta_current)
        sampling_rate_deltas.append(delta_rate)
        sampling_current_deltas.append(delta_current)
        sampling_residual_deltas.append(delta_residual)
        sampling_closure_errors.append(closure_error)
        native_sampling_order_commutator[target] = {
            "evaluate_native_grid_then_interpolate": {
                "mass_rate": float(composite_row["mass_rate"]),
                "current": float(composite_row["current"]),
                "residual": float(composite_row["residual"]),
            },
            "interpolate_native_operands_then_evaluate": {
                "mass_rate": float(reconstructed_row["fixed_R_mass_rate"]),
                "current": float(reconstructed_row["current"]),
                "residual": float(reconstructed_row["residual"]),
            },
            "interpolate_operands_minus_interpolate_composite": {
                "mass_rate": delta_rate,
                "current": delta_current,
                "residual": delta_residual,
            },
            "algebraic_budget_closure_error": closure_error,
            "interpretation_guard": (
                "Sampling-order commutator only; not a physical residual and not a "
                "selector for either derivative/interpolation convention."
            ),
        }
    native_sampling_order_summary = {
        "sample_count": len(native_sampling_order_commutator),
        "max_abs_mass_rate_delta": float(np.max(np.abs(sampling_rate_deltas))),
        "rms_mass_rate_delta": float(np.sqrt(np.mean(np.square(sampling_rate_deltas)))),
        "max_abs_current_delta": float(np.max(np.abs(sampling_current_deltas))),
        "max_abs_residual_delta": float(np.max(np.abs(sampling_residual_deltas))),
        "max_abs_algebraic_budget_closure_error": float(np.max(np.abs(sampling_closure_errors))),
    }

    if saved.tangent.state_digest(accepted) != accepted_digest:
        raise RuntimeError("accepted input mutated during Rr tangent audit")
    if float(np.max(np.abs(R_map_gap))) > 2e-11:
        raise RuntimeError("primitive areal map does not reproduce native nodes")
    if float(np.max(np.abs(Rt_map_gap))) > 2e-11:
        raise RuntimeError("primitive areal-time rate does not reproduce native RHS nodes")

    report = {
        "schema": "bssn-areal-Rr-time-tangent-compatibility/v1",
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
            "accepted_time": float(accepted.t),
            "fixed_R_interval": [float(low), float(high)],
            "fixed_R_sample_count": int(fixed_R.size),
            "primitive_interpolant": "four-point cubic in r^2 with coordinate-cell stencil held fixed",
            "native_tangent": "Rr=D4(R), Rr_dot=D4(Rt)",
            "primitive_tangent": "Rr=R*(1/r+b_r/(2b)-X_r/X); analytically differentiate this fixed cubic interpolant using b_dot and X_dot",
            "current_control": "same accepted-state rho,j,p_r,alpha,gamma^rr,R_t; report both native-Rr and primitive-Rr currents",
            "tangent_only_control": "hold primitive Rr and chi_r fixed; swap only D4(Rt) for the analytic tangent of primitive Rr",
            "fixed_R_anchored_crosswalk": "also samples six fixed areal radii using the accepted-state anchored four-point r^2 stencil from the existing primitive-Rr crosswalk; D4(chi) is held common to isolate the Rr time-tangent contribution",
            "native_fixed_R_sampling_order": "explicitly compares evaluate-native-grid-then-interpolate against interpolate-native-operands-then-evaluate; the commutator is bookkeeping and does not select a physical stencil",
        },
        "profile_summary_R_2_to_10": metrics,
        "fixed_R": fixed_rows,
        "fixed_R_anchored_stencil_crosswalk": fixed_R_anchored_rows,
        "native_fixed_R_sampling_order_commutator": native_sampling_order_commutator,
        "native_fixed_R_sampling_order_summary": native_sampling_order_summary,
        "acceptance": {
            "source_vendor_pins_match": True,
            "accepted_artifact_sha_verified": True,
            "accepted_state_immutable": True,
            "primitive_map_and_Rt_native_node_reconstruction_passed": True,
            "Rr_tangent_component_closure": float(np.max(np.abs(tangent_component_closure))) <= 2e-12,
            "mass_rate_current_delta_closure": float(np.max(np.abs(delta_closure))) <= 2e-12,
            "fixed_R_anchored_tangent_delta_closure_passed": bool(max(tangent_closure_abs, default=0.0) <= 2e-12),
            "areal_Rr_time_tangent_compatibility_closed": False,
            "integrated_current_identity_closed": False,
            "full_matter_coupled_propagation_identity_closed": False,
            "total_stress_energy_conservation_closed": False,
            "scientific_admission": False,
        },
        "scientific_admission": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "SHA256SUMS.txt").write_text(
        f"{sha256(report_path)}  report.json\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
