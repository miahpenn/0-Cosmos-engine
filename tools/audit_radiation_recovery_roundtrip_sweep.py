#!/usr/bin/env python3
"""Radiation primitive-recovery round-trip sweep; diagnostic-only, no production edits."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tools import audit_current_source_momentum_constraint_tangent as audit
from engine.matter_system import (
    KAPPA, _metric_arrays, initialize_radiation,
)
from engine.v55_matter import metric_slice_from_q
from engine.valencia import FluidConserved, recover_radiation

EXPECTED_INITIAL = "a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0"
EXPECTED_ACCEPTED = "c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a"
TOLS = (1.0e-8, 1.0e-10, 1.0e-12, 1.0e-13, 1.0e-14)


def roundtrip(metrics, state, tol: float) -> dict:
    n = len(metrics)
    rel_energy = np.full(n, np.nan)
    rel_momentum = np.full(n, np.nan)
    rel_current = np.full(n, np.nan)
    residual_norm = np.full(n, np.nan)
    failures = []
    for i, metric in enumerate(metrics):
        U = FluidConserved(
            rest=float(state.rest[i]),
            energy_t=float(state.energy_t[i]),
            momentum_r=float(state.momentum_r[i]),
        )
        E_in = U.energy_t / metric.sqrt_gamma
        S_in = U.momentum_r / metric.sqrt_gamma
        try:
            q = recover_radiation(metric, U, tol=tol)
            W = q.lorentz()
            h = q.rho + q.pressure
            E_out = h * W * W - q.pressure
            S_out = h * W * W * metric.gamma_rr * q.v_r
            denom_E = max(abs(E_in), 1.0e-300)
            denom_S = max(abs(S_in), abs(E_in), 1.0e-300)
            rel_energy[i] = abs(E_out - E_in) / denom_E
            rel_momentum[i] = abs(S_out - S_in) / denom_S
            rel_current[i] = abs((S_out - S_in) / KAPPA)
            residual_norm[i] = abs(E_out - E_in) / denom_E + abs(S_out - S_in) / denom_S
        except (ValueError, ArithmeticError, OverflowError) as exc:
            failures.append({"cell": i, "error": f"{type(exc).__name__}: {exc}"})
    finite = np.isfinite(rel_energy) & np.isfinite(rel_momentum) & np.isfinite(rel_current)
    return {
        "tol": tol,
        "cells": n,
        "failures": failures,
        "failure_count": len(failures),
        "all_finite": bool(np.all(finite)),
        "max_relative_energy_roundtrip_error": float(np.nanmax(rel_energy)) if np.any(finite) else None,
        "max_relative_momentum_roundtrip_error": float(np.nanmax(rel_momentum)) if np.any(finite) else None,
        "max_absolute_current_error": float(np.nanmax(rel_current)) if np.any(finite) else None,
        "max_combined_roundtrip_error": float(np.nanmax(residual_norm)) if np.any(finite) else None,
        "cell_0_relative_momentum_error": float(rel_momentum[0]) if finite[0] else None,
        "cell_1_relative_momentum_error": float(rel_momentum[1]) if finite[1] else None,
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

    metric_slice = metric_slice_from_q(accepted.grid, accepted.geometry)
    metrics = _metric_arrays(metric_slice)
    r = np.asarray(accepted.grid.centers, dtype=np.float64)
    profile = 1.0 + 0.1 * np.exp(-((r / 7.0) ** 2))
    v = 0.012 * np.tanh(r / 6.0)
    manufactured = initialize_radiation(metrics, 8.0e-9 * profile, v)

    results = {
        "accepted_state": {
            str(tol): roundtrip(metrics, accepted.matter.radiation, tol)
            for tol in TOLS
        },
        "manufactured_state": {
            str(tol): roundtrip(metrics, manufactured, tol)
            for tol in TOLS
        },
    }
    report = {
        "schema": "issue4_radiation_recovery_roundtrip_sweep_v1",
        "status": "diagnostic_only",
        "provenance": audit.provenance(),
        "configuration": {
            "resolution": 60, "r_max": 80.0, "dt": 0.015,
            "production_steps_advanced": 1, "recovery_tolerances": list(TOLS),
        },
        "state_identity": {
            "initial_sha256": initial_digest,
            "accepted_sha256": accepted_digest,
            "accepted_state_unchanged": audit.state_digest(accepted) == accepted_digest,
        },
        "results": results,
        "gates": {
            "recovery_error_source_isolated": False,
            "projection_step_size_plateau_established": False,
            "independent_full_propagation_identity_closed": False,
            "scientific_admission": False,
            "acceptance_oracle_satisfied": False,
        },
        "interpretation_limit": (
            "This measures recovery round-trip error only. It does not change the strict identity bound, "
            "prove the projection mismatch is solely recovery error, or establish propagation closure."
        ),
    }
    out = Path("results/accepted-state-raw-momentum-witness/radiation-recovery-roundtrip")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps(report, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
