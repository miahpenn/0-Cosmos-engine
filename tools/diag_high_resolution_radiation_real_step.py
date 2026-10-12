#!/usr/bin/env python3
"""Bounded N=240/320 radiation real-stage admissibility capture.

Each resolution receives exactly one normal production kernel step on the
unchanged PR #10 RHS/admissibility repair candidate. This checks whether the
older combined CMC A/B's N=240/320 "radiation RHS created inadmissible state"
survives after the artificial unit-step guard is removed. It performs no
retry, projection, floor, clamp, source edit, or time refinement.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import matter_system as ms
from engine import production_kernel as production_module
from engine.production_kernel import V55ProductionKernel
from engine.v55_matter import metric_slice_from_q

RESOLUTIONS = (240, 320)
R_MAX = 80.0
CFL = 0.03
VENDOR_COMMIT = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
OUT = Path("artifacts/high-resolution-radiation-real-step")


def finite_or_none(value):
    x = float(value)
    return x if math.isfinite(x) else None


def state_snapshot(metric, radiation, stage):
    metrics = ms._metric_arrays(metric)
    sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics], dtype=float)
    gamma_rr_inv = np.asarray([m.gamma_rr_inv for m in metrics], dtype=float)
    U_E = np.asarray(radiation.energy_t, dtype=float)
    U_r = np.asarray(radiation.momentum_r, dtype=float)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        E = U_E / sqrt_gamma
        abs_S = np.sqrt(np.maximum(gamma_rr_inv, 0.0)) * np.abs(U_r) / sqrt_gamma
        ratio = abs_S / np.abs(E)
    margin = U_E - np.sqrt(np.maximum(gamma_rr_inv, 0.0)) * np.abs(U_r)

    def finite_min(arr):
        valid = np.isfinite(arr)
        if not np.any(valid):
            return {"value": None, "i": None}
        idxs = np.flatnonzero(valid)
        i = int(idxs[np.argmin(arr[valid])])
        return {"value": finite_or_none(arr[i]), "i": i}

    def finite_max(arr):
        valid = np.isfinite(arr)
        if not np.any(valid):
            return {"value": None, "i": None}
        idxs = np.flatnonzero(valid)
        i = int(idxs[np.argmax(arr[valid])])
        return {"value": finite_or_none(arr[i]), "i": i}

    def cell(i):
        i = int(i)
        m = metrics[i]
        return {
            "i": i,
            "r": finite_or_none(metric.r[i]),
            "U_E": finite_or_none(U_E[i]),
            "U_r": finite_or_none(U_r[i]),
            "E": finite_or_none(E[i]),
            "abs_S": finite_or_none(abs_S[i]),
            "abs_S_over_E": finite_or_none(ratio[i]),
            "conserved_cone_margin": finite_or_none(margin[i]),
            "a": finite_or_none(metric.a[i]),
            "b": finite_or_none(metric.b[i]),
            "X": finite_or_none(metric.X[i]),
            "alpha": finite_or_none(metric.alpha[i]),
            "beta": finite_or_none(metric.beta[i]),
            "Aa": finite_or_none(metric.Aa[i]),
            "K": finite_or_none(metric.K[i]),
            "gamma_rr": finite_or_none(m.gamma_rr),
            "gamma_rr_inv": finite_or_none(m.gamma_rr_inv),
            "sqrt_gamma": finite_or_none(m.sqrt_gamma),
        }

    ratio_valid = np.isfinite(ratio)
    top_indices = np.flatnonzero(ratio_valid)
    top_indices = sorted(top_indices, key=lambda i: float(ratio[i]), reverse=True)[:5]
    valid_margin = np.isfinite(margin)
    metric_valid = (
        np.all(np.isfinite(sqrt_gamma))
        and np.all(sqrt_gamma > 0.0)
        and np.all(np.isfinite(gamma_rr_inv))
        and np.all(gamma_rr_inv > 0.0)
    )
    conservative_finite = (
        np.all(np.isfinite(radiation.rest))
        and np.all(np.isfinite(U_E))
        and np.all(np.isfinite(U_r))
    )
    return {
        "stage": stage,
        "cells": int(len(U_E)),
        "all_conservative_values_finite": bool(conservative_finite),
        "all_metric_cells_valid": bool(metric_valid),
        "full_domain": {
            "min_cone_margin": finite_min(margin),
            "max_absS_over_E": finite_max(ratio),
            "all_cells_satisfy_E_nonnegative_and_cone": bool(
                conservative_finite and metric_valid
                and np.all(E >= 0.0) and np.all(abs_S <= E)
            ),
            "nonfinite_ratio_cells": int(np.count_nonzero(~ratio_valid)),
        },
        "centre_first_three_cells": {
            "min_cone_margin": finite_min(margin[:3]),
            "max_absS_over_E": finite_max(ratio[:3]),
            "all_cells_satisfy_E_nonnegative_and_cone": bool(
                np.all(E[:3] >= 0.0) and np.all(abs_S[:3] <= E[:3])
            ),
        },
        "cell_0": cell(0),
        "cell_1": cell(1) if len(U_E) > 1 else None,
        "cell_2": cell(2) if len(U_E) > 2 else None,
        "top_five_absS_over_E": [cell(i) for i in top_indices],
    }


def one_step(n):
    kernel = V55ProductionKernel()
    initial = kernel.initialize(
        resolution=n, r_max=R_MAX, include_radiation=True
    )
    dt = CFL * float(initial.grid.dr)
    initial_metric = metric_slice_from_q(initial.grid, initial.geometry)
    initial_snapshot = state_snapshot(
        initial_metric, initial.matter.radiation, "initial_state"
    )
    captures = []
    original_validator = production_module.validate_species_admissibility

    def validator_spy(metric, radiation, species, *, stage):
        if species is ms.Species.RADIATION:
            snap = state_snapshot(metric, radiation, stage)
            try:
                result = original_validator(metric, radiation, species, stage=stage)
            except Exception as exc:
                entry = {
                    "validation": "rejected",
                    "snapshot": snap,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "validator_diagnostic": getattr(exc, "diagnostic", None),
                }
                captures.append(entry)
                raise
            else:
                captures.append({"validation": "passed", "snapshot": snap})
                return result
        return original_validator(metric, radiation, species, stage=stage)

    production_module.validate_species_admissibility = validator_spy
    final = None
    error = None
    try:
        final = kernel.step(initial, dt)
    except Exception as exc:
        error = {
            "type": type(exc).__name__,
            "message": str(exc),
            "validator_diagnostic": getattr(exc, "diagnostic", None),
        }
    finally:
        production_module.validate_species_admissibility = original_validator

    result = {
        "N": int(n),
        "r_max": R_MAX,
        "dr": float(initial.grid.dr),
        "dt": dt,
        "CFL": CFL,
        "vendor_commit": VENDOR_COMMIT,
        "initial": initial_snapshot,
        "actual_kernel_validator_calls": captures,
    }
    if final is None:
        diag = error.get("validator_diagnostic") if error else None
        if diag and diag.get("stage"):
            result["step_status"] = "actual_stage_admissibility_rejection"
        elif error and "radiation RHS created inadmissible state" in error["message"]:
            result["step_status"] = "unit_step_guard_still_active_or_wrong_source"
        else:
            result["step_status"] = "other_kernel_failure"
        result["step_error"] = error
    else:
        result["step_status"] = "completed"
        final_metric = metric_slice_from_q(final.grid, final.geometry)
        result["final"] = state_snapshot(
            final_metric, final.matter.radiation, "accepted_final_state"
        )
        result["geometry_max_abs_delta"] = {
            key: finite_or_none(np.max(np.abs(
                getattr(final.geometry, key) - getattr(initial.geometry, key)
            )))
            for key in ("a", "b", "X", "alpha", "Aa", "K", "Lambda")
        }
    return result


def main():
    cases = [one_step(n) for n in RESOLUTIONS]
    failures = [
        c for c in cases
        if c["step_status"] in ("other_kernel_failure", "unit_step_guard_still_active_or_wrong_source")
    ]
    result = {
        "schema": "high_resolution_radiation_real_step_v1",
        "status": "diagnostic_capture_complete" if not failures else "diagnostic_harness_or_unexpected_failure",
        "diagnostic_only": True,
        "source_branch": "diag/high-resolution-radiation-real-step-after-rhs-contract-2026-10-12",
        "source_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "candidate_base": "6f6fd9d81af739cacce16fcd428bdfeac996e343",
        "pinned_geometry_vendor": VENDOR_COMMIT,
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "dt_rule": "dt = 0.03 * dr; one production kernel step per resolution",
        "question": "Does the N=240/320 centre-cell rejection from the old combined CMC A/B persist at actual production admissibility stages after the dt=1 RHS extraction guard is bypassed?",
        "interpretation_limits": [
            "A captured actual-stage rejection is a real bounded integration failure, not proof of which equation is wrong.",
            "A successful single step does not prove high-resolution campaign stability or production readiness.",
            "This test does not change or tune the cone condition, matter source, geometry RHS, gauge, or boundary.",
            "The earlier CMC A/B artifact is pinned to a source where species_rhs calls evolve_species(dt=1.0) and the radiation output-state guard is unconditional; its N=240/320 exception is therefore not admitted as a real finite-step trajectory failure.",
        ],
        "cases": cases,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.json"
    report_path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    digest = hashlib.sha256(report_path.read_bytes()).hexdigest()
    (OUT / "report.sha256").write_text(digest + "  report.json\n")
    print(json.dumps({
        "status": result["status"],
        "report": str(report_path),
        "sha256": digest,
        "cases": [{
            "N": c["N"],
            "step_status": c["step_status"],
            "validator_stages": [
                {"stage": x["snapshot"]["stage"], "validation": x["validation"]}
                for x in c["actual_kernel_validator_calls"]
            ],
            "error": c.get("step_error"),
        } for c in cases],
    }, indent=2, allow_nan=False))
    if failures:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
