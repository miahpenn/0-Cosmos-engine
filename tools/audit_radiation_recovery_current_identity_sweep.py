#!/usr/bin/env python3
"""Paired radiation-recovery/current-identity sweep; diagnostic-only."""
from __future__ import annotations

import json
import math
from pathlib import Path
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


def traced_pressure(E: float, S2: float, tol: float) -> tuple[float, int, str]:
    """Mirror the pinned recovery bisection to expose iteration/exit class."""
    if E < 0.0:
        raise ValueError("negative energy")
    if S2 == 0.0:
        return E / 3.0, 0, "zero_velocity_analytic"
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
            return mid, iteration, "converged"
        if f_mid >= 0.0:
            hi = mid
        else:
            lo = mid
    mid = 0.5 * (lo + hi)
    f_mid, v2 = residual(mid)
    if not math.isfinite(f_mid) or v2 >= 1.0:
        raise ValueError("iteration-limit invalid endpoint")
    return mid, 80, "iteration_limit_return"


def inspect_state(label, metrics, state, tol: float) -> dict:
    rows = []
    for i, metric in enumerate(metrics):
        U = FluidConserved(float(state.rest[i]), float(state.energy_t[i]), float(state.momentum_r[i]))
        E = float(U.energy_t / metric.sqrt_gamma)
        Sr = float(U.momentum_r / metric.sqrt_gamma)
        S2 = float(metric.gamma_rr_inv * Sr * Sr)
        p_trace, iterations, exit_class = traced_pressure(E, S2, tol)
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
            "strict_identity_passed": abs(observed) <= roundoff_bound,
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
    report = {
        "schema": "issue4_radiation_recovery_paired_current_identity_v1",
        "status": "diagnostic_only",
        "provenance": audit.provenance(),
        "configuration": {"resolution": 60, "r_max": 80.0, "dt": 0.015, "production_steps_advanced": 1, "tolerances": list(TOLS)},
        "state_identity": {"initial_sha256": initial_digest, "accepted_expected_sha256": EXPECTED_ACCEPTED, "accepted_before_sha256": accepted_before, "accepted_after_sha256": accepted_after, "accepted_state_unchanged": accepted_before == accepted_after},
        "identity": "j_projected - j_direct = f(p) * sqrt(gamma_rr) * sign(S_r) * sqrt(v2) / KAPPA",
        "bound_policy": "unchanged strict bound: 128*eps*(abs(j_projected)+abs(j_direct))+1e-30; no recovery-aware relaxation",
        "paired_results": paired,
        "tolerance_trends": trends,
        "gates": {"recovery_error_source_isolated": False, "projection_step_size_plateau_established": False, "independent_full_propagation_identity_closed": False, "scientific_admission": False, "acceptance_oracle_satisfied": False},
        "interpretation_limit": "The pressure-residual identity is a diagnostic decomposition, not by itself proof that recovery causes the strict discrepancy or that full momentum propagation closes. Solver classification and negative control are retained; all scientific gates remain false.",
    }
    out = Path("results/accepted-state-raw-momentum-witness/radiation-recovery-current-identity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "state_identity": report["state_identity"], "trends": trends, "gates": report["gates"]}, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
