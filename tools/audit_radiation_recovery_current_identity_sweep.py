#!/usr/bin/env python3
"""Paired radiation-recovery/current-identity sweep; diagnostic-only."""
from __future__ import annotations

import json
import math
from pathlib import Path
from fractions import Fraction
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import audit_current_source_momentum_constraint_tangent as audit
from engine.matter_system import KAPPA, _metric_arrays, initialize_radiation
from engine.v55_matter import metric_slice_from_q
from engine.valencia import FluidConserved, recover_radiation

EXPECTED_INITIAL = "a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0"
EXPECTED_ACCEPTED = "c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a"
TOLS = (1.0e-8, 1.0e-10, 1.0e-12, 1.0e-13, 1.0e-14)
EPS = np.finfo(np.float64).eps


def traced_pressure(E: float, S2: float, tol: float) -> tuple[float, int, str, float, float]:
    """Mirror the pinned recovery bisection to expose iteration/exit class."""
    if E < 0.0:
        raise ValueError("negative energy")
    if S2 == 0.0:
        return E / 3.0, 0, "zero_velocity_analytic", E / 3.0, E / 3.0
    if math.sqrt(max(S2, 0.0)) >= E:
        raise ValueError("radiation state outside strict E>|S| cone")

    def residual(p: float) -> tuple[float, float]:
        z = E + p
        v2 = S2 / (z * z)
        if v2 >= 1.0:
            return math.inf, v2
        return (3.0 * p - E) + 4.0 * p * v2 / (1.0 - v2), v2

    lo, hi = 0.0, E / 3.0
    f_lo, _ = residual(lo)
    f_hi, _ = residual(hi)
    if f_hi < 0.0:
        hi = math.nextafter(hi, math.inf)
        f_hi, _ = residual(hi)
    if not math.isfinite(f_hi) or f_lo > 0.0 or f_hi < 0.0:
        raise ValueError("no physical pressure bracket")
    for iteration in range(1, 81):
        mid = 0.5 * (lo + hi)
        f_mid, v2 = residual(mid)
        if abs(f_mid) <= tol * E:
            return mid, iteration, "converged", lo, hi
        if f_mid >= 0.0:
            hi = mid
        else:
            lo = mid
    mid = 0.5 * (lo + hi)
    f_mid, v2 = residual(mid)
    if not math.isfinite(f_mid) or v2 >= 1.0:
        raise ValueError("iteration-limit invalid endpoint")
    return mid, 80, "iteration_limit_return", lo, hi



def _down(x: float) -> float:
    if not math.isfinite(x):
        raise ValueError("non-finite interval endpoint")
    y = math.nextafter(x, -math.inf)
    if not math.isfinite(y):
        raise ValueError("outward rounding overflow")
    return y


def _up(x: float) -> float:
    if not math.isfinite(x):
        raise ValueError("non-finite interval endpoint")
    y = math.nextafter(x, math.inf)
    if not math.isfinite(y):
        raise ValueError("outward rounding overflow")
    return y


def _pt(x: float) -> tuple[float, float]:
    if not math.isfinite(x):
        raise ValueError("non-finite point interval")
    return (x, x)


def _iadd(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    return (_down(a[0] + b[0]), _up(a[1] + b[1]))


def _isub(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    return (_down(a[0] - b[1]), _up(a[1] - b[0]))


def _imul(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    values = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    if not all(math.isfinite(x) for x in values):
        raise ValueError("non-finite interval product")
    return (_down(min(values)), _up(max(values)))


def _idiv(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    if b[0] <= 0.0 <= b[1]:
        raise ValueError("interval divisor contains zero")
    values = (a[0] / b[0], a[0] / b[1], a[1] / b[0], a[1] / b[1])
    if not all(math.isfinite(x) for x in values):
        raise ValueError("non-finite interval quotient")
    return (_down(min(values)), _up(max(values)))


def _isqrt(a: tuple[float, float]) -> tuple[float, float]:
    if a[0] < 0.0 or not math.isfinite(a[1]):
        raise ValueError("invalid square-root interval")
    return (_down(math.sqrt(a[0])) if a[0] > 0.0 else 0.0,
            _up(math.sqrt(a[1])))


def exact_residual_sign(E: float, S2: float, pressure: float) -> int:
    """Exact sign of f(p) for the exact rational values represented by binary64 inputs."""
    e = Fraction.from_float(E)
    s2 = Fraction.from_float(S2)
    p = Fraction.from_float(pressure)
    denominator = (e + p) * (e + p) - s2
    if denominator <= 0:
        raise ValueError("exact residual denominator is non-positive")
    value = 3 * p - e + (4 * p * s2) / denominator
    return (value > 0) - (value < 0)


def residual_interval(E: float, S2: float, pressure: float) -> tuple[float, float]:
    """Outward-rounded enclosure of f(p) for the radiation pressure equation."""
    e, s2, p = _pt(E), _pt(S2), _pt(pressure)
    z = _iadd(e, p)
    v2 = _idiv(s2, _imul(z, z))
    one_minus_v2 = _isub(_pt(1.0), v2)
    first = _isub(_imul(_pt(3.0), p), e)
    second = _idiv(_imul(_imul(_pt(4.0), p), v2), one_minus_v2)
    return _iadd(first, second)


def identity_residual_interval(E: float, pressure: float, v2: float) -> tuple[float, float]:
    """Enclose the production identity residual using q.pressure and q.v2."""
    e, p, vv = _pt(E), _pt(pressure), _pt(v2)
    first = _isub(_imul(_pt(3.0), p), e)
    second = _idiv(_imul(_imul(_pt(4.0), p), vv),
                   _isub(_pt(1.0), vv))
    return _iadd(first, second)


def identity_correction_interval(E: float, pressure: float, v2: float,
                                 gamma_rr: float, Sr: float) -> tuple[float, float]:
    """Enclose the signed correction using the production result's p and v2."""
    f_value = identity_residual_interval(E, pressure, v2)
    root_metric = _isqrt(_pt(gamma_rr))
    root_v2 = _isqrt(_pt(v2))
    sign = -1.0 if Sr < 0.0 else (1.0 if Sr > 0.0 else 0.0)
    return _idiv(_imul(_imul(_imul(f_value, root_metric), _pt(sign)),
                       root_v2), _pt(KAPPA))


def interval_certificate(E: float, S2: float, p_trace: float, trace_lo: float,
                         trace_hi: float, pressure: float, v2: float,
                         gamma_rr: float, Sr: float, observed: float,
                         roundoff_bound: float) -> dict:
    """Fail-closed residual/bracket/pressure-error certificate for one cell."""
    result = {"certified": False, "failure": None}
    try:
        if not all(math.isfinite(x) for x in
                   (E, S2, p_trace, trace_lo, trace_hi, pressure, v2, gamma_rr, Sr, observed, roundoff_bound)):
            raise ValueError("non-finite input")
        if E <= 0.0 or S2 < 0.0 or gamma_rr <= 0.0 or not (0.0 <= v2 < 1.0):
            raise ValueError("outside physical certificate domain")
        if not (0.0 <= pressure <= E / 3.0):
            raise ValueError("production pressure outside physical bracket")

        f_lo = residual_interval(E, S2, 0.0)
        hi = E / 3.0
        f_hi = residual_interval(E, S2, hi)
        upper_endpoint = hi
        upper_branch = "E_over_3"
        if f_hi[0] < 0.0:
            upper_endpoint = math.nextafter(hi, math.inf)
            f_hi = residual_interval(E, S2, upper_endpoint)
            upper_branch = "nextafter_E_over_3"
        bracket_certified = f_lo[1] <= 0.0 and f_hi[0] >= 0.0
        if not bracket_certified:
            raise ValueError("endpoint residual signs not certified")

        # Certify S2 <= 8 E^2 / 9 using outward-rounded arithmetic.
        threshold = _idiv(_imul(_pt(8.0), _imul(_pt(E), _pt(E))), _pt(9.0))
        derivative_condition = _pt(S2)[1] <= threshold[0]
        if not derivative_condition:
            raise ValueError("cannot certify f'(p) >= 3 over [0,E/3]")

        root_residual = residual_interval(E, S2, pressure)
        rho_f = max(abs(root_residual[0]), abs(root_residual[1]))
        pressure_error_bound = _up(rho_f / 3.0)
        trace_lo_residual = residual_interval(E, S2, trace_lo)
        trace_hi_residual = residual_interval(E, S2, trace_hi)
        trace_bracket_certified = trace_lo_residual[1] <= 0.0 and trace_hi_residual[0] >= 0.0
        trace_lo_exact_sign = exact_residual_sign(E, S2, trace_lo)
        trace_hi_exact_sign = exact_residual_sign(E, S2, trace_hi)
        trace_bracket_exact_signs_certified = trace_lo_exact_sign <= 0 and trace_hi_exact_sign >= 0
        trace_analytic_zero_velocity = S2 == 0.0 and trace_lo == E / 3.0 and trace_hi == E / 3.0
        trace_classification_certified = (
            trace_bracket_certified or trace_bracket_exact_signs_certified
            or trace_analytic_zero_velocity)
        trace_residual = residual_interval(E, S2, p_trace)
        identity_residual = identity_residual_interval(E, pressure, v2)
        correction = identity_correction_interval(
            E, pressure, v2, gamma_rr, Sr)
        correction_covers_observed = correction[0] <= observed <= correction[1]
        correction_covers_with_bound = (
            correction[0] - roundoff_bound <= observed <= correction[1] + roundoff_bound)
        result.update({
            "certified": trace_classification_certified,
            "pressure_root_certified": True,
            "trace_classification_certified": trace_classification_certified,
            "failure": None if trace_classification_certified else "traced solver endpoint signs not certified",
            "endpoint_residual_at_zero": list(f_lo),
            "endpoint_residual_at_upper": list(f_hi),
            "upper_endpoint_branch": upper_branch,
            "upper_endpoint": upper_endpoint,
            "bracket_certified": True,
            "derivative_lower_bound": 3.0,
            "derivative_condition_S2_le_8E2_over_9": True,
            "production_pressure_residual_interval": list(root_residual),
            "rho_f": rho_f,
            "pressure_error_bound_rho_over_3": pressure_error_bound,
            "trace_bracket_lower_residual_interval": list(trace_lo_residual),
            "trace_bracket_upper_residual_interval": list(trace_hi_residual),
            "trace_bracket_signs_certified": trace_classification_certified,
            "trace_bracket_interval_signs_certified": trace_bracket_certified,
            "trace_bracket_exact_signs_certified": trace_bracket_exact_signs_certified,
            "trace_lo_pressure": trace_lo,
            "trace_hi_pressure": trace_hi,
            "trace_lo_pressure_hex": trace_lo.hex(),
            "trace_hi_pressure_hex": trace_hi.hex(),
            "trace_lo_exact_residual_sign": trace_lo_exact_sign,
            "trace_hi_exact_residual_sign": trace_hi_exact_sign,
            "trace_analytic_zero_velocity_exception": trace_analytic_zero_velocity,
            "trace_pressure_residual_interval": list(trace_residual),
            "production_identity_residual_interval_using_q_v2": list(identity_residual),
            "trace_pressure_delta": pressure - p_trace,
            "identity_correction_interval": list(correction),
            "correction_interval_covers_observed_signed_difference": correction_covers_observed,
            "correction_interval_covers_observed_with_existing_roundoff_bound": correction_covers_with_bound,
            "existing_strict_roundoff_bound_used_for_coverage": roundoff_bound,
        })
    except (ValueError, OverflowError, ZeroDivisionError) as exc:
        result["failure"] = str(exc)
    return result


def inspect_state(label, metrics, state, tol: float) -> dict:
    rows = []
    for i, metric in enumerate(metrics):
        U = FluidConserved(float(state.rest[i]), float(state.energy_t[i]), float(state.momentum_r[i]))
        E = float(U.energy_t / metric.sqrt_gamma)
        Sr = float(U.momentum_r / metric.sqrt_gamma)
        S2 = float(metric.gamma_rr_inv * Sr * Sr)
        p_trace, iterations, exit_class, trace_lo, trace_hi = traced_pressure(E, S2, tol)
        q = recover_radiation(metric, U, tol=tol)
        pressure_delta = float(q.pressure - p_trace)
        v2 = float(q.v2)
        f_p = float((3.0 * q.pressure - E) + 4.0 * q.pressure * v2 / (1.0 - v2)) if v2 != 0.0 else float(3.0 * q.pressure - E)
        direct = Sr / KAPPA
        projected = float((q.rho + q.pressure) * q.lorentz()**2 * metric.gamma_rr * q.v_r / KAPPA)
        observed = projected - direct
        sign_v = 0.0 if v2 == 0.0 else math.copysign(math.sqrt(v2), Sr)
        predicted = f_p * math.sqrt(metric.gamma_rr) * sign_v / KAPPA
        remainder = observed - predicted
        wrong_sign_prediction = -predicted
        roundoff_bound = 128.0 * EPS * (abs(projected) + abs(direct)) + 1.0e-30
        certificate = interval_certificate(E, S2, p_trace, trace_lo, trace_hi, q.pressure, v2, metric.gamma_rr, Sr, observed, roundoff_bound)
        rows.append({
            "cell": i, "E": E, "S_r": Sr, "S2": S2,
            "pressure": float(q.pressure), "v2": v2, "f_p": f_p,
            "iterations": iterations, "exit_class": exit_class,
            "trace_pressure_delta": pressure_delta,
            "j_direct": direct, "j_projected": projected,
            "observed_signed_difference": observed,
            "predicted_signed_difference": predicted,
            "arithmetic_remainder": remainder,
            "wrong_sign_prediction": wrong_sign_prediction,
            "wrong_sign_remainder": observed - wrong_sign_prediction,
            "strict_roundoff_bound": roundoff_bound,
            "interval_certificate": certificate,
            "strict_identity_passed": bool(abs(observed) <= roundoff_bound),
        })
    finite = all(all(math.isfinite(v) for k, v in row.items() if isinstance(v, (int, float))) for row in rows)
    wrong_sign_detected = all(
        abs(row["wrong_sign_remainder"]) > abs(row["arithmetic_remainder"])
        for row in rows if row["predicted_signed_difference"] != 0.0
    )
    return {
        "state": label, "tol": tol, "cells": len(rows), "all_finite": finite,
        "max_abs_observed_difference": max(abs(x["observed_signed_difference"]) for x in rows),
        "max_abs_predicted_difference": max(abs(x["predicted_signed_difference"]) for x in rows),
        "max_abs_arithmetic_remainder": max(abs(x["arithmetic_remainder"]) for x in rows),
        "max_strict_roundoff_bound": max(x["strict_roundoff_bound"] for x in rows),
        "strict_identity_passed_all_cells": all(x["strict_identity_passed"] for x in rows),
        "max_trace_pressure_delta": max(abs(x["trace_pressure_delta"]) for x in rows),
        "exit_class_counts": {name: sum(x["exit_class"] == name for x in rows) for name in sorted({x["exit_class"] for x in rows})},
        "wrong_sign_negative_control_detected": wrong_sign_detected,
        "per_cell": rows,
    }


def main() -> int:
    kernel = audit.V55ProductionKernel()
    initial = kernel.initialize(resolution=60, r_max=80.0, include_radiation=True)
    initial_digest = audit.state_digest(initial)
    if initial_digest != EXPECTED_INITIAL:
        raise RuntimeError(f"initial state digest mismatch: {initial_digest}")
    accepted = kernel.step(initial, 0.015)
    accepted_digest = audit.state_digest(accepted)
    if accepted_digest != EXPECTED_ACCEPTED:
        raise RuntimeError(f"accepted state digest mismatch: {accepted_digest}")
    accepted_before = audit.state_digest(accepted)
    metric_slice = metric_slice_from_q(accepted.grid, accepted.geometry)
    metrics = _metric_arrays(metric_slice)
    r = np.asarray(accepted.grid.centers, dtype=np.float64)
    profile = 1.0 + 0.1 * np.exp(-((r / 7.0) ** 2))
    v = 0.012 * np.tanh(r / 6.0)
    manufactured = initialize_radiation(metrics, 8.0e-9 * profile, v)

    paired = []
    for state_name, state in (("accepted", accepted.matter.radiation), ("manufactured", manufactured)):
        for tol in TOLS:
            paired.append(inspect_state(state_name, metrics, state, tol))
    accepted_after = audit.state_digest(accepted)
    if accepted_before != accepted_after:
        raise RuntimeError("diagnostic mutated accepted state")
    trends = {}
    for state_name in ("accepted", "manufactured"):
        group = [x for x in paired if x["state"] == state_name]
        trends[state_name] = {
            "tolerance_order": [x["tol"] for x in group],
            "observed_maxima": [x["max_abs_observed_difference"] for x in group],
            "predicted_maxima": [x["max_abs_predicted_difference"] for x in group],
            "remainder_maxima": [x["max_abs_arithmetic_remainder"] for x in group],
            "strict_identity_pass": [x["strict_identity_passed_all_cells"] for x in group],
        }
    # Explicit same-state production-path A/B: compare the existing default
    # tolerance with a tighter diagnostic tolerance without changing the default.
    accepted_default = next(x for x in paired if x["state"] == "accepted" and x["tol"] == 1.0e-12)
    accepted_tight = next(x for x in paired if x["state"] == "accepted" and x["tol"] == 1.0e-14)
    tight_by_cell = {row["cell"]: row for row in accepted_tight["per_cell"]}
    ab_rows = []
    for default_row in accepted_default["per_cell"]:
        tight_row = tight_by_cell[default_row["cell"]]
        bound = default_row["strict_roundoff_bound"]
        ab_rows.append({
            "cell": default_row["cell"],
            "default_tol": 1.0e-12,
            "tight_tol": 1.0e-14,
            "default_pressure": default_row["pressure"],
            "tight_pressure": tight_row["pressure"],
            "pressure_change_tight_minus_default": tight_row["pressure"] - default_row["pressure"],
            "default_v2": default_row["v2"],
            "tight_v2": tight_row["v2"],
            "default_j_direct": default_row["j_direct"],
            "tight_j_direct": tight_row["j_direct"],
            "default_j_projected": default_row["j_projected"],
            "tight_j_projected": tight_row["j_projected"],
            "default_signed_difference": default_row["observed_signed_difference"],
            "tight_signed_difference": tight_row["observed_signed_difference"],
            "default_predicted_signed_difference": default_row["predicted_signed_difference"],
            "tight_predicted_signed_difference": tight_row["predicted_signed_difference"],
            "default_remainder": default_row["arithmetic_remainder"],
            "tight_remainder": tight_row["arithmetic_remainder"],
            "unchanged_strict_roundoff_bound": bound,
            "default_fails_strict_bound": not default_row["strict_identity_passed"],
            "tight_passes_strict_bound": tight_row["strict_identity_passed"],
            "default_pressure_residual": default_row["f_p"],
            "tight_pressure_residual": tight_row["f_p"],
            "default_wrong_sign_control_detected": abs(default_row["wrong_sign_remainder"]) > abs(default_row["arithmetic_remainder"]),
            "tight_wrong_sign_control_detected": abs(tight_row["wrong_sign_remainder"]) > abs(tight_row["arithmetic_remainder"]),
        })
    ab_all_default_fail = all(row["default_fails_strict_bound"] for row in ab_rows)
    ab_all_tight_pass = all(row["tight_passes_strict_bound"] for row in ab_rows)
    ab_controls_detected = all(row["default_wrong_sign_control_detected"] and row["tight_wrong_sign_control_detected"] for row in ab_rows)
    ab_pressure_correction_accounted = all(
        abs(row["default_remainder"]) <= row["unchanged_strict_roundoff_bound"]
        and abs(row["tight_remainder"]) <= row["unchanged_strict_roundoff_bound"]
        for row in ab_rows
    )
    production_path_ab = {
        "status": "diagnostic_only",
        "same_accepted_state": True,
        "arms": {"default": 1.0e-12, "tight_diagnostic": 1.0e-14},
        "cell_count": len(ab_rows),
        "per_cell": ab_rows,
        "summary": {
            "all_default_cells_fail_unchanged_strict_bound": ab_all_default_fail,
            "all_tight_cells_pass_unchanged_strict_bound": ab_all_tight_pass,
            "both_wrong_sign_controls_detected_in_both_arms": ab_controls_detected,
            "pressure_residual_correction_accounts_within_bound_in_both_arms": ab_pressure_correction_accounted,
            "tolerance_cause_supported_by_preregistered_rule": bool(
                ab_all_default_fail and ab_all_tight_pass and ab_controls_detected and ab_pressure_correction_accounted
            ),
        },
        "decision_rule": "Support tolerance as cause only if default fails, tight passes the unchanged bound, per-cell signed correction accounts for the observed discrepancy within that bound in both arms, and wrong-sign negative controls remain detected in both arms. Otherwise tolerance alone is not established.",
    }
    report = {
        "schema": "issue4_radiation_recovery_paired_current_identity_v1",
        "status": "diagnostic_only",
        "provenance": audit.provenance(),
        "configuration": {"resolution": 60, "r_max": 80.0, "dt": 0.015, "production_steps_advanced": 1, "tolerances": list(TOLS)},
        "state_identity": {"initial_sha256": initial_digest, "accepted_expected_sha256": EXPECTED_ACCEPTED, "accepted_before_sha256": accepted_before, "accepted_after_sha256": accepted_after, "accepted_state_unchanged": accepted_before == accepted_after},
        "identity": "j_projected - j_direct = f(p) * sqrt(gamma_rr) * sign(S_r) * sqrt(v2) / KAPPA",
        "bound_policy": "unchanged strict bound: 128*eps*(abs(j_projected)+abs(j_direct))+1e-30; no recovery-aware relaxation",
        "paired_results": paired,
        "production_path_ab": production_path_ab,
        "interval_certificate_summary": {
            "certified_cells": sum(row["interval_certificate"]["certified"] for group in paired for row in group["per_cell"]),
            "pressure_root_certified_cells": sum(row["interval_certificate"].get("pressure_root_certified", False) for group in paired for row in group["per_cell"]),
            "trace_bracket_certified_cells": sum(row["interval_certificate"].get("trace_classification_certified", False) for group in paired for row in group["per_cell"]),
            "total_cells": sum(group["cells"] for group in paired),
            "all_cells_certified": all(row["interval_certificate"]["certified"] for group in paired for row in group["per_cell"]),
            "all_pressure_roots_certified": all(row["interval_certificate"].get("pressure_root_certified", False) for group in paired for row in group["per_cell"]),
            "all_trace_brackets_certified": all(row["interval_certificate"].get("trace_classification_certified", False) for group in paired for row in group["per_cell"]),
            "interval_trace_brackets_certified_cells": sum(row["interval_certificate"].get("trace_bracket_interval_signs_certified", False) for group in paired for row in group["per_cell"]),
            "exact_trace_brackets_certified_cells": sum(row["interval_certificate"].get("trace_bracket_exact_signs_certified", False) for group in paired for row in group["per_cell"]),
            "all_raw_observed_differences_covered": all(row["interval_certificate"].get("correction_interval_covers_observed_signed_difference", False) for group in paired for row in group["per_cell"]),
            "all_differences_covered_with_existing_roundoff_bound": all(row["interval_certificate"].get("correction_interval_covers_observed_with_existing_roundoff_bound", False) for group in paired for row in group["per_cell"] if row["interval_certificate"].get("pressure_root_certified", False)),
        },
        "tolerance_trends": trends,
        "gates": {"recovery_error_source_isolated": False, "projection_step_size_plateau_established": False, "independent_full_propagation_identity_closed": False, "scientific_admission": False, "acceptance_oracle_satisfied": False},
        "interpretation_limit": "The pressure-residual identity is a diagnostic decomposition, not by itself proof that recovery causes the strict discrepancy or that full momentum propagation closes. Solver classification and negative control are retained; all scientific gates remain false.",
    }
    out = Path("results/accepted-state-raw-momentum-witness/radiation-recovery-current-identity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(_json_safe(report), indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "state_identity": report["state_identity"], "trends": trends, "gates": report["gates"]}, sort_keys=True), flush=True)
    return 0


def _json_safe(value):
    """Normalize nested NumPy scalars before strict JSON serialization."""
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    return value


if __name__ == "__main__":
    raise SystemExit(main())
