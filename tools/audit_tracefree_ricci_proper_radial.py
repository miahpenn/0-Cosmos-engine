#!/usr/bin/env python3
"""Independent proper-radial Ricci / direct-D2 Hessian cross-check at actual PIRK L2 inputs.

This diagnostic captures the exact a,b,X,Lambda,alpha arrays at the four real
production L2 stage calls while preserving all wrapped production return values.
It computes a second curvature route from proper radial distance and areal radius,
which does not use the evolved Lambda curvature expression. A direct-D2 lapse
Hessian is independently assembled and compared with the captured vendor term.
No production, vendor, physics, gauge, or campaign parameter is edited.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import math
import os

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SOURCE_COMMIT = os.environ.get("SOURCE_COMMIT", "e6b501fb19c4f2e21ad9c655662840a3d2582d54")
SOURCE_KERNEL_BLOB = os.environ.get("SOURCE_KERNEL_BLOB", "9685e76bb5a6f26d6bf383e240de24f37523f265")
VENDOR_COMMIT = os.environ.get("VENDOR_COMMIT", "d6052d605673ce9d82cdc99b0df79855fa2ea215")
SOURCE_BLOBS = {
    "engine/production_kernel.py": "9685e76bb5a6f26d6bf383e240de24f37523f265",
    "engine/matter_system.py": "e52a117690253ddad7b24d116145a65474a351d2",
    "engine/matter_rhs.py": "f41180efbcff2c72cf36477a765318ab539fac80",
    "engine/scalar_system.py": "ffc119b1a6f2d495d972323901f76a85783b0a1b",
    "engine/v55_matter.py": "0ee3c51fac7dc46c63feba2747a750542bcee3d0",
    "engine/v55_pirk_adapter.py": "eeee0492d80889d70e52b0946563eee5054c0f11",
    "engine/cmc_gauge.py": "7b49b5751abebd4faf0130208032cd2080b71f64",
    "engine/v55_initial.py": "7db2e1a10cf16fb40e4221b0df504fdf39aae611",
    "engine/stress_energy.py": "ac3e7518c62bb72e2218a61489fa486bb787809f",
    "engine/valencia.py": "03ed86398f424c03d9d6a3e646da9222a34fac5c",
    "engine/worldtube.py": "aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6",
}
NATIVE_CANDIDATE_PATHS = frozenset({
    "engine/native_geometry/__init__.py",
    "engine/native_geometry/constraints.py",
    "engine/native_geometry/geometry.py",
    "engine/native_geometry/operators.py",
    "engine/native_geometry/evolution.py",
})

# Owner-authorized, exact-blob exception for the first bounded CMC runtime seam.
# Historical SOURCE_BLOBS remain immutable and are checked against SOURCE_COMMIT.
# Only these exact head blobs may replace the historical head copy; any further
# production edit changes the blob and fails closed until separately authorized.
AUTHORIZED_RUNTIME_HEAD_BLOBS = {
    "engine/production_kernel.py": "cc4ae0f58dd5df4e05f596f6e677526586956f26",
    "engine/v55_pirk_adapter.py": "be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3",
}
AUTHORIZED_RUNTIME_PATHS = frozenset(AUTHORIZED_RUNTIME_HEAD_BLOBS)

RESOLUTIONS = (40, 60, 80, 120, 160)
R_MAX = 80.0
CFL = 0.03
FIXED_RADII = (2.0, 3.0, 5.0, 8.0)
STAGE_NAMES = (
    "initial_rhs",
    "explicit_predictor_rhs",
    "primary_predictor_rhs",
    "final_implicit_l2_rhs",
)
OUT = Path("results/tracefree-ricci-proper-radial-2026-10-15")


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()




def _unexpected_source_pinned_paths(
    changed_paths: list[str], native_paths_in_source: set[str]
) -> list[str]:
    """Allow additive native files or the two exact-blob authorized runtime paths."""
    unexpected = []
    for path in changed_paths:
        if path in AUTHORIZED_RUNTIME_PATHS:
            # Exact candidate blob verification happens independently in main().
            continue
        if path not in NATIVE_CANDIDATE_PATHS:
            unexpected.append(path)
        elif path in native_paths_in_source:
            unexpected.append(f"{path} (native candidate exists in pinned source)")
    return unexpected


def _source_blob_matches(path: str, expected: str, baseline_blob: str, head_blob: str) -> bool:
    """Require the immutable historical baseline and either its blob or one exact approved head blob."""
    if baseline_blob != expected:
        return False
    authorized_head = AUTHORIZED_RUNTIME_HEAD_BLOBS.get(path)
    return head_blob == expected or (
        authorized_head is not None and head_blob == authorized_head
    )


def _arr(x) -> np.ndarray:
    return np.asarray(x, dtype=float).copy()


def _finite(x) -> bool:
    return bool(np.all(np.isfinite(np.asarray(x, dtype=float))))


def _max_abs(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.max(np.abs(x))) if x.size else 0.0


def _json_finite_safe(value):
    """Make nested report values strict-JSON-safe without hiding failure flags.

    Non-finite numeric witnesses become JSON null; companion finite/check flags
    and the diagnostic errors remain intact so a failed run can still be audited.
    """
    if isinstance(value, np.ndarray):
        return _json_finite_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_finite_safe(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): _json_finite_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_finite_safe(item) for item in value]
    return value


def _write_strict_json_report(path: Path, payload) -> None:
    """Write a report even when failure diagnostics include NaN/Inf witnesses."""
    safe_payload = _json_finite_safe(payload)
    path.write_text(
        json.dumps(safe_payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


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
        a = np.ascontiguousarray(np.asarray(named[name]))
        h.update(name.encode("utf-8"))
        h.update(str(a.dtype).encode("ascii"))
        h.update(repr(a.shape).encode("ascii"))
        h.update(a.tobytes())
    return h.hexdigest()


def _d4(values, h: float, parity: int, order: int = 1) -> np.ndarray:
    """Independent port of the pinned fourth-order cell-centred derivative."""
    u = np.asarray(values, dtype=float)
    n = u.size
    if u.ndim != 1 or n < 5 or parity not in (-1, 1) or order not in (1, 2):
        raise ValueError("invalid independent derivative input")
    ext = np.empty(n + 4)
    ext[2:-2] = u
    ext[1] = parity * u[0]
    ext[0] = parity * u[1]
    out = np.empty_like(u)
    for i in range(n - 2):
        s = ext[i:i + 5]
        if order == 1:
            out[i] = (s[0] - 8*s[1] + 8*s[3] - s[4]) / (12*h)
        else:
            out[i] = (-s[0] + 16*s[1] - 30*s[2] + 16*s[3] - s[4]) / (12*h*h)
    for i in (n - 2, n - 1):
        idx = np.arange(n - 5, n)
        offsets = idx.astype(float) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[order] = float(math.factorial(order))
        out[i] = float(np.linalg.solve(matrix, target) @ u[idx]) / h**order
    return out


def _interp(r, values, target: float) -> dict:
    r = np.asarray(r, dtype=float)
    v = np.asarray(values, dtype=float)
    if target < r[0] or target > r[-1]:
        return {"target_r": float(target), "in_domain": False}
    right = int(np.searchsorted(r, target, side="left"))
    if right < r.size and np.isclose(r[right], target, rtol=0.0, atol=1e-14):
        left, weight = right, 0.0
    else:
        left = right - 1
        if left < 0 or right >= r.size:
            return {"target_r": float(target), "in_domain": False}
        weight = float((target - r[left]) / (r[right] - r[left]))
    return {
        "target_r": float(target),
        "value": float((1.0-weight)*v[left] + weight*v[right]),
        "left_i": int(left), "right_i": int(right), "weight_right": float(weight),
    }


def _peak(values, r, mask=None) -> dict:
    v = np.asarray(values, dtype=float)
    r = np.asarray(r, dtype=float)
    ids = np.arange(v.size) if mask is None else np.flatnonzero(mask)
    if ids.size == 0:
        return {"max_abs": None, "signed": None, "i": None, "r": None}
    i = int(ids[np.argmax(np.abs(v[ids]))])
    return {"max_abs": float(abs(v[i])), "signed": float(v[i]),
            "i": i, "r": float(r[i])}


def _independent_metric_crosscheck(capture: dict, h: float) -> tuple[dict, dict]:
    """Independent curvature via dl^2 + R(l)^2 dOmega^2; direct-D2 lapse Hessian."""
    r = np.asarray(capture["r"], dtype=float)
    a = np.asarray(capture["a"], dtype=float)
    b = np.asarray(capture["b"], dtype=float)
    X = np.asarray(capture["X"], dtype=float)
    alpha = np.asarray(capture["alpha"], dtype=float)
    Lambda = np.asarray(capture["Lambda"], dtype=float)
    q = X / np.sqrt(a)              # dr/dl
    areal = r * np.sqrt(b) / X

    # R_areal is odd; q and alpha are even.
    areal_r = _d4(areal, h, -1, 1)
    areal_rr = _d4(areal, h, -1, 2)
    q_r = _d4(q, h, +1, 1)
    areal_prime_l = q * areal_r
    areal_second_l = q*q*areal_rr + q*q_r*areal_r
    ricci_radial_metric = -2.0 * areal_second_l / areal
    ricci_scalar_metric = (
        -4.0 * areal_second_l / areal
        + 2.0 * (1.0 - areal_prime_l**2) / areal**2
    )
    ricci_tf_metric = ricci_radial_metric - ricci_scalar_metric / 3.0
    ricci_tf_vendor = np.asarray(capture["Rr_r"]) - np.asarray(capture["R"]) / 3.0
    ricci_gap = ricci_tf_vendor - ricci_tf_metric

    # Independent scalar Hessian: H^l_l=q^2 alpha''+q q' alpha',
    # H^theta_theta=(R_areal'/R_areal) q alpha'. Uses direct D2(alpha).
    alpha_r = _d4(alpha, h, +1, 1)
    alpha_rr = _d4(alpha, h, +1, 2)
    hess_radial = q*q*alpha_rr + q*q_r*alpha_r
    hess_angular = (areal_prime_l / areal) * q * alpha_r
    hess_lap = hess_radial + 2.0*hess_angular
    hess_tf_metric = -(hess_radial - hess_lap / 3.0)
    hess_tf_vendor = -(
        np.asarray(capture["DrDr_alpha"]) - np.asarray(capture["lap_alpha"]) / 3.0
    )
    hess_gap = hess_tf_metric - hess_tf_vendor

    a_r = _d4(a, h, +1, 1)
    b_r = _d4(b, h, +1, 1)
    metric_connection = a_r/(2.0*a*a) - b_r/(a*b) + 2.0/r*(1.0/b - 1.0/a)
    lambda_gap = Lambda - metric_connection

    arrays = {
        "metric_areal_radius": areal,
        "metric_areal_prime_proper": areal_prime_l,
        "metric_areal_second_proper": areal_second_l,
        "ricci_radial_metric": ricci_radial_metric,
        "ricci_scalar_metric": ricci_scalar_metric,
        "ricci_tracefree_metric": ricci_tf_metric,
        "ricci_tracefree_vendor": ricci_tf_vendor,
        "ricci_tracefree_gap": ricci_gap,
        "ricci_tracefree_Aa_metric": alpha*ricci_tf_metric,
        "ricci_tracefree_Aa_gap": alpha*ricci_gap,
        "hessian_tracefree_metric": hess_tf_metric,
        "hessian_tracefree_vendor": hess_tf_vendor,
        "hessian_tracefree_gap": hess_gap,
        "metric_connection": metric_connection,
        "lambda_metric_gap": lambda_gap,
    }
    finite = all(_finite(v) for v in arrays.values())
    centre_mask = np.arange(r.size) < min(6, r.size)
    centre_first_two_mask = np.arange(r.size) < min(2, r.size)
    fixed_annulus_mask = (areal >= 5.0) & (areal <= 8.0)
    annulus_count = int(np.count_nonzero(fixed_annulus_mask))
    exterior_mask = r >= 5.0
    centre = []
    for i in range(min(6, r.size)):
        centre.append({
            "i": int(i), "r_coordinate": float(r[i]), "r_areal": float(areal[i]),
            "vendor_tf_ricci": float(ricci_tf_vendor[i]),
            "metric_tf_ricci": float(ricci_tf_metric[i]),
            "vendor_minus_metric_tf_ricci": float(ricci_gap[i]),
            "alpha_weighted_ricci_gap": float(alpha[i]*ricci_gap[i]),
            "vendor_tf_hessian": float(hess_tf_vendor[i]),
            "direct_d2_tf_hessian": float(hess_tf_metric[i]),
            "direct_d2_minus_vendor_hessian": float(hess_gap[i]),
            "lambda_minus_metric_connection": float(lambda_gap[i]),
        })
    fixed = []
    for target in FIXED_RADII:
        fixed.append({
            "target_r_coordinate": float(target),
            "vendor_minus_metric_tf_ricci": _interp(r, ricci_gap, target),
            "alpha_weighted_ricci_gap": _interp(r, alpha*ricci_gap, target),
            "direct_d2_minus_vendor_hessian": _interp(r, hess_gap, target),
            "lambda_minus_metric_connection": _interp(r, lambda_gap, target),
        })
    hess_max = _max_abs(hess_gap)
    summary = {
        "finite": bool(finite),
        "ricci_gap_vendor_minus_metric_max_abs_full": _max_abs(ricci_gap),
        "ricci_gap_vendor_minus_metric_peak_full": _peak(ricci_gap, r),
        "ricci_gap_max_abs_centre_i0_i5": _max_abs(ricci_gap[centre_mask]),
        "ricci_gap_peak_centre_i0_i5": _peak(ricci_gap, r, centre_mask),
        "ricci_gap_max_abs_centre_i0_i1": _max_abs(ricci_gap[centre_first_two_mask]),
        "ricci_gap_peak_centre_i0_i1": _peak(ricci_gap, r, centre_first_two_mask),
        "ricci_gap_fixed_areal_annulus_R5_R8": {
            "count": annulus_count,
            "max_abs": _max_abs(ricci_gap[fixed_annulus_mask]) if annulus_count else None,
            "peak": _peak(ricci_gap, r, fixed_annulus_mask),
            "vendor_tf_max_abs": _max_abs(ricci_tf_vendor[fixed_annulus_mask]) if annulus_count else None,
            "proper_radial_tf_max_abs": _max_abs(ricci_tf_metric[fixed_annulus_mask]) if annulus_count else None,
        },
        "ricci_gap_max_abs_r_ge_5": _max_abs(ricci_gap[exterior_mask]) if np.any(exterior_mask) else None,
        "ricci_gap_peak_r_ge_5": _peak(ricci_gap, r, exterior_mask),
        "alpha_weighted_ricci_gap_max_abs_full": _max_abs(alpha*ricci_gap),
        "direct_d2_hessian_gap_max_abs_full": hess_max,
        "direct_d2_hessian_gap_peak_full": _peak(hess_gap, r),
        "direct_d2_hessian_gap_max_abs_centre_i0_i5": _max_abs(hess_gap[centre_mask]),
        "direct_d2_hessian_gap_max_abs_centre_i0_i1": _max_abs(hess_gap[centre_first_two_mask]),
        "direct_d2_hessian_gap_fixed_areal_annulus_R5_R8": _max_abs(hess_gap[fixed_annulus_mask]) if annulus_count else None,
        "direct_d2_hessian_gap_max_abs_r_ge_5": _max_abs(hess_gap[exterior_mask]) if np.any(exterior_mask) else None,
        "direct_d2_hessian_gap_relative_to_larger_operator_full": float(
            hess_max / max(_max_abs(hess_tf_metric), _max_abs(hess_tf_vendor), np.finfo(float).tiny)
        ),
        "lambda_metric_connection_gap_max_abs_full": _max_abs(lambda_gap),
        "lambda_metric_connection_gap_peak_full": _peak(lambda_gap, r),
        "centre_cells": centre,
        "centre_cells_first_two": centre[:2],
        "fixed_areal_annulus_definition": "5 <= R_areal <= 8",
        "fixed_coordinate_radius_samples": fixed,
    }
    return arrays, summary


def main() -> int:
    head = _git("rev-parse", "HEAD")
    source_pins = {}
    for path, expected in SOURCE_BLOBS.items():
        baseline_blob = _git("rev-parse", f"{SOURCE_COMMIT}:{path}")
        head_file_blob = _git("rev-parse", f"HEAD:{path}")
        authorized_head = AUTHORIZED_RUNTIME_HEAD_BLOBS.get(path)
        source_pins[path] = {
            "historical_expected_blob": expected,
            "baseline_blob": baseline_blob,
            "diagnostic_head_blob": head_file_blob,
            "authorized_runtime_head_blob": authorized_head,
            "matches": _source_blob_matches(path, expected, baseline_blob, head_file_blob),
        }
        if not source_pins[path]["matches"]:
            raise RuntimeError(f"source provenance mismatch for {path}: {source_pins[path]}")
    changed = _git("diff", "--name-only", SOURCE_COMMIT, "HEAD", "--", "engine", "vendor", ".gitmodules", "requirements.txt")
    changed_paths = changed.splitlines() if changed else []
    native_paths_in_source = set()
    for path in NATIVE_CANDIDATE_PATHS:
        try:
            _git("cat-file", "-e", f"{SOURCE_COMMIT}:{path}")
        except subprocess.CalledProcessError:
            # The source pin is already reachable and all production pins above
            # passed; a missing candidate path means this new file is additive.
            continue
        native_paths_in_source.add(path)
    unexpected = _unexpected_source_pinned_paths(changed_paths, native_paths_in_source)
    if unexpected:
        raise RuntimeError(
            "production/vendor/dependency paths changed from current-main source pin:\n"
            + "\n".join(unexpected)
        )
    source_blob = _git("rev-parse", f"{SOURCE_COMMIT}:engine/production_kernel.py")
    head_blob = _git("rev-parse", "HEAD:engine/production_kernel.py")
    vendor = _git("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
    authorized_kernel_head = AUTHORIZED_RUNTIME_HEAD_BLOBS.get("engine/production_kernel.py")
    kernel_head_matches = head_blob == SOURCE_KERNEL_BLOB or head_blob == authorized_kernel_head
    if source_blob != SOURCE_KERNEL_BLOB or not kernel_head_matches:
        raise RuntimeError(
            f"production kernel provenance mismatch: historical={SOURCE_KERNEL_BLOB}; "
            f"source={source_blob}; diagnostic branch={head_blob}; "
            f"authorized_head={authorized_kernel_head}"
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
            if capture_enabled[0] and caller.f_code.co_name == "primary_l2_rhs":
                captures.append({
                    "stage": active_stage[0] or "final_implicit_l2_rhs",
                    "r": _arr(grid.centers),
                    "alpha": _arr(geometry.alpha),
                    "a": _arr(geometry.a),
                    "b": _arr(geometry.b),
                    "X": _arr(geometry.X),
                    "Lambda": _arr(geometry.Lambda),
                    "R": _arr(terms["R"]),
                    "Rr_r": _arr(terms["Rr_r"]),
                    "DrDr_alpha": _arr(terms["DrDr_alpha"]),
                    "lap_alpha": _arr(terms["lap_alpha"]),
                })
            return terms

        def wrapped_l2(grid, geometry):
            caller_name = sys._getframe(1).f_code.co_name
            should_capture = caller_name in ("geometry_stage_terms", "step")
            before = len(captures)
            previous = capture_enabled[0]
            capture_enabled[0] = should_capture
            try:
                result = original_l2(grid, geometry)
            finally:
                capture_enabled[0] = previous
            if not should_capture:
                return result
            if len(captures) != before + 1:
                raise RuntimeError(
                    f"{caller_name}: expected one geometry_terms call for selected L2 stage; "
                    f"captured {len(captures)-before}"
                )
            captures[-1]["rhs_Aa"] = _arr(result["Aa"])
            captures[-1]["rhs_K"] = _arr(result["K"])
            captures[-1]["lapse_hessian_Aa"] = -(
                captures[-1]["DrDr_alpha"] - captures[-1]["lap_alpha"]/3.0
            )
            captures[-1]["ricci_tracefree_Aa"] = captures[-1]["alpha"]*(
                captures[-1]["Rr_r"] - captures[-1]["R"]/3.0
            )
            return result

        def wrapped_stage_terms(grid, geometry, fields, matter, lambda_m=2.0, **stage_kwargs):
            index = stage_counter[0]
            if index >= len(stage_names):
                raise RuntimeError("unexpected additional geometry stage during production step")
            active_stage[0] = stage_names[index]
            before = len(captures)
            try:
                result = original_stage_terms(
                    grid, geometry, fields, matter,
                    lambda_m=lambda_m, **stage_kwargs,
                )
            finally:
                active_stage[0] = None
            if len(captures) != before + 1:
                raise RuntimeError(
                    f"{stage_names[index]} expected one selected L2 call, "
                    f"observed {len(captures)-before}"
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
            raise AssertionError(f"N={resolution}: production step mutated its input state")
        if len(captures) != 4:
            raise AssertionError(f"N={resolution}: expected 4 actual L2 captures, got {len(captures)}")
        observed = [x["stage"] for x in captures]
        if observed != list(STAGE_NAMES):
            raise AssertionError(f"N={resolution}: stage order mismatch: {observed}")
        accepted.geometry.assert_finite_positive()
        if not (np.isfinite(accepted.t) and accepted.t > state.t):
            raise AssertionError(f"N={resolution}: bounded production step did not advance")

        stage_reports = []
        all_stage_checks_pass = True
        profile_keys = (
            "r", "alpha", "a", "b", "X", "Lambda", "R", "Rr_r",
            "DrDr_alpha", "lap_alpha", "lapse_hessian_Aa", "ricci_tracefree_Aa",
            "rhs_Aa", "rhs_K", "metric_areal_radius", "metric_areal_prime_proper",
            "metric_areal_second_proper", "ricci_radial_metric", "ricci_scalar_metric",
            "ricci_tracefree_metric", "ricci_tracefree_vendor", "ricci_tracefree_gap",
            "ricci_tracefree_Aa_metric", "ricci_tracefree_Aa_gap",
            "hessian_tracefree_metric", "hessian_tracefree_vendor",
            "hessian_tracefree_gap", "metric_connection", "lambda_metric_gap",
        )
        for capture in captures:
            h = float(capture["r"][1] - capture["r"][0])
            metric_arrays, metric_summary = _independent_metric_crosscheck(capture, h)
            capture.update(metric_arrays)
            l2_aa_sum = capture["lapse_hessian_Aa"] + capture["ricci_tracefree_Aa"]
            aa_closure = l2_aa_sum - capture["rhs_Aa"]
            k_closure = -capture["lap_alpha"] - capture["rhs_K"]
            finite = all(_finite(capture[k]) for k in profile_keys)
            aa_bound = 64.0*np.finfo(float).eps*max(
                1.0, _max_abs(capture["lapse_hessian_Aa"]),
                _max_abs(capture["ricci_tracefree_Aa"]), _max_abs(capture["rhs_Aa"])
            )
            k_bound = 64.0*np.finfo(float).eps*max(
                1.0, _max_abs(capture["lap_alpha"]), _max_abs(capture["rhs_K"])
            )
            stage_checks = {
                "stage": capture["stage"], "finite": bool(finite and metric_summary["finite"]),
                "Aa_reconstruction_closure_max_abs": _max_abs(aa_closure),
                "Aa_roundoff_bound": aa_bound,
                "Aa_closure_within_roundoff_bound": bool(_max_abs(aa_closure) <= aa_bound),
                "K_reconstruction_closure_max_abs": _max_abs(k_closure),
                "K_roundoff_bound": k_bound,
                "K_closure_within_roundoff_bound": bool(_max_abs(k_closure) <= k_bound),
                "independent_metric_crosscheck": metric_summary,
            }
            stage_reports.append(stage_checks)
            stage_pass = (
                stage_checks["finite"]
                and stage_checks["Aa_closure_within_roundoff_bound"]
                and stage_checks["K_closure_within_roundoff_bound"]
            )
            all_stage_checks_pass = bool(all_stage_checks_pass and stage_pass)
            for key in profile_keys:
                profiles[f"N{resolution}_{capture['stage']}_{key}"] = _arr(capture[key])

        if not all_stage_checks_pass:
            errors.append(f"N={resolution}: nonfinite or same-stage RHS reconstruction failed; Ricci/Hessian gaps are reported without a fixed physical-error cutoff")
        reports.append({
            "resolution": resolution, "r_max": R_MAX, "dr": float(state.grid.dr), "dt": float(dt),
            "source_state_sha256_before": input_digest, "source_state_sha256_after": output_digest,
            "source_state_unchanged": input_digest == output_digest,
            "accepted_time": float(accepted.t), "captured_actual_l2_calls": observed,
            "stage_summaries": stage_reports,
            "all_checks_pass": bool(all_stage_checks_pass), "step_status": "completed",
        })

    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "schema": "v55_tracefree_ricci_proper_radial_refresh_v1",
        "status": "completed_diagnostic_only" if not errors else "diagnostic_crosscheck_failure",
        "source_commit": SOURCE_COMMIT,
        "source_production_kernel_blob": SOURCE_KERNEL_BLOB,
        "source_blob_manifest": source_pins,
        "diagnostic_head_commit": _git("rev-parse", "HEAD"),
        "vendor_commit": VENDOR_COMMIT,
        "resolutions": list(RESOLUTIONS), "r_max": R_MAX, "cfl": CFL,
        "requested_stage_order": list(STAGE_NAMES),
        "ordinary_production_steps_per_resolution": 1,
        "production_kernel_and_vendor_edited": bool(head_blob != source_blob or vendor != VENDOR_COMMIT),
        "authorized_runtime_source_transition": bool(head_blob in AUTHORIZED_RUNTIME_HEAD_BLOBS.values()),
        "metric_method": {
            "areal_radius": "R_areal = r*sqrt(b)/X",
            "proper_radial_factor": "q = X/sqrt(a), d/dl = q*d/dr",
            "radial_ricci": "R^l_l = -2*R_areal''(l)/R_areal",
            "scalar_ricci": "R = -4*R_areal''(l)/R_areal + 2*(1 - R_areal'(l)^2)/R_areal^2",
            "trace_free_radial_ricci": "R^l_l - R/3",
            "connection_usage": "Lambda is captured and the metric-connection residual is reported; Lambda is not used in the independent metric-Ricci formula",
            "regional_split": "first two centre cells (i=0,1) separated from fixed areal-radius annulus 5<=R_areal<=8; legacy centre-six/r>=5/coordinate-radius samples retained",
        },
        "hessian_method": {
            "radial_hessian": "q^2*D2(alpha) + q*D1(q)*D1(alpha)",
            "angular_hessian": "(D_l R_areal / R_areal)*q*D1(alpha)",
            "trace_free_Aa_contribution": "-(H^l_l - (H^l_l+2*H^theta_theta)/3)",
            "acceptance_rule": "reported unthresholded; no fixed absolute Hessian discrepancy cut is applied",
        },
        "errors": errors,
        "resolution_reports": reports,
        "guardrails": [
            "Original production geometry_terms and primary_l2_rhs outputs are returned unchanged.",
            "The cross-check captures the exact a,b,X,Lambda,alpha arrays at the four real L2 input stages.",
            "The independent Ricci discrepancy is retained as evidence and is never forced to vanish.",
            "The independent metric-Ricci route is not substituted into BSSN production.",
            "Historical production/vendor pins were verified; the report records whether the exact owner-authorized CMC runtime source blobs are active.",
            "One ordinary production step per resolution is not a long campaign or scientific admission.",
        ],
    }
    report_path = OUT / "report.json"
    _write_strict_json_report(report_path, report)
    npz_path = OUT / "profiles.npz"
    np.savez_compressed(npz_path, **profiles)
    csv_path = OUT / "profiles.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow((
            "resolution", "stage", "i", "r_coordinate", "alpha", "a", "b", "X", "Lambda",
            "R_areal", "vendor_tf_ricci", "metric_tf_ricci", "vendor_minus_metric_tf_ricci",
            "alpha_weighted_ricci_gap", "vendor_tf_hessian", "direct_d2_tf_hessian",
            "direct_d2_minus_vendor_hessian", "lambda_minus_metric_connection",
            "vendor_l2_Aa", "Aa_reconstruction_closure", "vendor_l2_K", "K_reconstruction_closure",
        ))
        suffixes = (
            "alpha", "a", "b", "X", "Lambda", "metric_areal_radius",
            "ricci_tracefree_vendor", "ricci_tracefree_metric", "ricci_tracefree_gap",
            "ricci_tracefree_Aa_gap", "hessian_tracefree_vendor", "hessian_tracefree_metric",
            "hessian_tracefree_gap", "lambda_metric_gap", "rhs_Aa", "lapse_hessian_Aa",
            "ricci_tracefree_Aa", "rhs_K", "lap_alpha",
        )
        for key in sorted(profiles):
            if not key.endswith("_r"):
                continue
            prefix = key[:-2]
            if not all(prefix + "_" + x in profiles for x in suffixes):
                continue
            arrays = {name: profiles[prefix + "_" + name] for name in suffixes}
            resolution = int(prefix.split("_")[0][1:])
            stage = "_".join(prefix.split("_")[1:])
            r = profiles[key]
            for i in range(r.size):
                Aa_close = arrays["lapse_hessian_Aa"][i] + arrays["ricci_tracefree_Aa"][i] - arrays["rhs_Aa"][i]
                K_close = -arrays["lap_alpha"][i] - arrays["rhs_K"][i]
                writer.writerow((
                    resolution, stage, i, float(r[i]), float(arrays["alpha"][i]),
                    float(arrays["a"][i]), float(arrays["b"][i]), float(arrays["X"][i]),
                    float(arrays["Lambda"][i]), float(arrays["metric_areal_radius"][i]),
                    float(arrays["ricci_tracefree_vendor"][i]), float(arrays["ricci_tracefree_metric"][i]),
                    float(arrays["ricci_tracefree_gap"][i]), float(arrays["ricci_tracefree_Aa_gap"][i]),
                    float(arrays["hessian_tracefree_vendor"][i]), float(arrays["hessian_tracefree_metric"][i]),
                    float(arrays["hessian_tracefree_gap"][i]), float(arrays["lambda_metric_gap"][i]),
                    float(arrays["rhs_Aa"][i]), float(Aa_close), float(arrays["rhs_K"][i]), float(K_close),
                ))
    sums_paths = (report_path, npz_path, csv_path)
    (OUT / "SHA256SUMS.txt").write_text(
        "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in sums_paths),
        encoding="utf-8",
    )
    (OUT / "SOURCE_COMMIT.txt").write_text(SOURCE_COMMIT + "\n", encoding="utf-8")
    (OUT / "SOURCE_KERNEL_BLOB.txt").write_text(SOURCE_KERNEL_BLOB + "\n", encoding="utf-8")
    (OUT / "VENDOR_COMMIT.txt").write_text(VENDOR_COMMIT + "\n", encoding="utf-8")

    print("TRACEFREE_RICCI_PROPER_RADIAL_REFRESH_" + ("PASS" if not errors else "FAIL"))
    print(f"source={SOURCE_COMMIT} kernel_blob={SOURCE_KERNEL_BLOB} vendor={VENDOR_COMMIT} head={report['diagnostic_head_commit']}")
    print("resolution_summary=" + json.dumps(_json_finite_safe([
        {
            "N": item["resolution"],
            "stages": [
                {
                    "stage": stage["stage"],
                    "ricci_gap_i0": stage["independent_metric_crosscheck"]["centre_cells"][0]["vendor_minus_metric_tf_ricci"],
                    "ricci_gap_i1": stage["independent_metric_crosscheck"]["centre_cells"][1]["vendor_minus_metric_tf_ricci"],
                    "ricci_gap_first2_max": stage["independent_metric_crosscheck"]["ricci_gap_max_abs_centre_i0_i1"],
                    "ricci_gap_fixed_areal_R5_R8": stage["independent_metric_crosscheck"]["ricci_gap_fixed_areal_annulus_R5_R8"]["max_abs"],
                    "fixed_areal_annulus_count": stage["independent_metric_crosscheck"]["ricci_gap_fixed_areal_annulus_R5_R8"]["count"],
                    "ricci_gap_rge5": stage["independent_metric_crosscheck"]["ricci_gap_max_abs_r_ge_5"],
                    "hessian_gap_max": stage["independent_metric_crosscheck"]["direct_d2_hessian_gap_max_abs_full"],
                    "lambda_metric_gap_max": stage["independent_metric_crosscheck"]["lambda_metric_connection_gap_max_abs_full"],
                } for stage in item["stage_summaries"]
            ],
        } for item in reports
    ]), sort_keys=True, allow_nan=False))
    print(f"report_sha256={hashlib.sha256(report_path.read_bytes()).hexdigest()}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
