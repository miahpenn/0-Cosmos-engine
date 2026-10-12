#!/usr/bin/env python3
"""Audit discrete propagation of Lambda - Lambda_metric over one production step.

Reads t=10 checkpoints from the already-completed t=12 campaign, plus the exact
production initializer at t=0. Runs one ordinary production step on copies of
these four states and records the actual split-scheme Lambda increment versus
the change in the metric-derived radial connection from (a,b). No physics
coefficient or production equation is changed.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path("source").resolve()))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()
R_CUT = 5.0
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED_FIELDS = ("rest", "energy_t", "momentum_r")


def d1(grid, arr, parity):
    return np.asarray(grid.cell_derivative_fourth(arr, parity=parity), dtype=float)


def lambda_metric(grid, geometry):
    """Metric-derived contracted radial conformal connection from a,b."""
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    return (
        d1(grid, a, +1)/(2.0*a**2)
        - d1(grid, b, +1)/(a*b)
        + 2.0/r*(1.0/b - 1.0/a)
    )


def peak(values, r, mask):
    ids = np.flatnonzero(mask)
    if ids.size == 0:
        return None
    idx = int(ids[np.argmax(np.abs(values[ids]))])
    return {"i": idx, "r": float(r[idx]), "signed": float(values[idx]),
            "abs": float(abs(values[idx]))}


def restore(path: Path, kernel):
    with np.load(path, allow_pickle=False) as archive:
        data = {k: archive[k].copy() for k in archive.files}
    n = int(data["r"].size)
    dr = float(data["r"][1] - data["r"][0])
    state = kernel.initialize(
        resolution=n, r_max=n*dr, include_radiation=True
    )
    if not np.allclose(state.grid.centers, data["r"], rtol=0.0,
                       atol=1.e-12*max(1.0, n*dr)):
        raise ValueError(f"{path}: checkpoint grid mismatch")
    for name in GEOMETRY_FIELDS:
        setattr(state.geometry, name, np.asarray(data[name], dtype=float).copy())
    state.scalars = ScalarFields(*(
        np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS
    ))
    species = {}
    for name, prefix in SPECIES.items():
        species[name] = ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{field}"], dtype=float).copy()
            for field in CONSERVED_FIELDS
        ))
    state.matter = V55MatterState(
        dark_matter=species["dark_matter"],
        baryons=species["baryons"],
        radiation=species["radiation"],
    )
    state.t = float(data["t"])
    state.tau = float(data["tau"])
    state.e_folds = float(data["e_folds"])
    state.history = []
    state.handoffs = []
    return state


def _finite(label, arr):
    arr = np.asarray(arr, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise FloatingPointError(f"non-finite values in {label}")
    return arr


def audit_case(kernel, state, label):
    grid = state.grid
    r = np.asarray(grid.centers, dtype=float)
    n = len(r)
    dt = 0.03 * float(grid.dr)
    a0 = np.asarray(state.geometry.a, dtype=float).copy()
    b0 = np.asarray(state.geometry.b, dtype=float).copy()
    lam0 = np.asarray(state.geometry.Lambda, dtype=float).copy()
    lm0 = lambda_metric(grid, state.geometry)
    c0 = lam0 - lm0

    stage_terms = []
    outside_lambda_l2 = []
    in_geometry_terms = [False]
    original_terms = adapter.geometry_stage_terms
    original_lambda_l2 = vacuum.lambda_l2_rhs

    def traced_terms(g, geometry, fields, matter, *args, **kwargs):
        in_geometry_terms[0] = True
        try:
            result = original_terms(g, geometry, fields, matter, *args, **kwargs)
        finally:
            in_geometry_terms[0] = False
        stage_terms.append({
            "call_index": len(stage_terms),
            "geometry_t": float(state.t),
            "a": np.asarray(geometry.a, dtype=float).copy(),
            "b": np.asarray(geometry.b, dtype=float).copy(),
            "Lambda": np.asarray(geometry.Lambda, dtype=float).copy(),
            "explicit_a": _finite("explicit_a", result["explicit"]["a"]).copy(),
            "explicit_b": _finite("explicit_b", result["explicit"]["b"]).copy(),
            "lambda_l2": _finite("lambda_l2", result["lambda_l2"]).copy(),
            "lambda_l3": _finite("lambda_l3", result["lambda_l3"]).copy(),
        })
        return result

    def traced_lambda_l2(g, geometry, *args, **kwargs):
        result = original_lambda_l2(g, geometry, *args, **kwargs)
        if not in_geometry_terms[0]:
            outside_lambda_l2.append({
                "a": np.asarray(geometry.a, dtype=float).copy(),
                "b": np.asarray(geometry.b, dtype=float).copy(),
                "Lambda": np.asarray(geometry.Lambda, dtype=float).copy(),
                "rhs": _finite("final_lambda_l2_rhs", result).copy(),
            })
        return result

    adapter.geometry_stage_terms = traced_terms
    vacuum.lambda_l2_rhs = traced_lambda_l2
    try:
        final = kernel.step(state, dt)
    finally:
        adapter.geometry_stage_terms = original_terms
        vacuum.lambda_l2_rhs = original_lambda_l2

    if len(stage_terms) != 3:
        raise ValueError(f"{label}: expected 3 geometry stage evaluations, got {len(stage_terms)}")
    if len(outside_lambda_l2) < 1:
        raise ValueError(f"{label}: final direct lambda_l2 RHS was not observed")
    # The production split uses geometry_stage_terms at g0, the explicit
    # predictor and the primary predictor; the final direct call is lambda_l2
    # on final_u_old_v with updated a,b,X,Aa,K and the primary Lambda.
    l20 = stage_terms[0]["lambda_l2"]
    l30 = stage_terms[0]["lambda_l3"]
    l31 = stage_terms[2]["lambda_l3"]
    l2final = outside_lambda_l2[-1]["rhs"]
    lambda_l2_contribution_rate = 0.5*(l20 + l2final)
    lambda_l3_contribution_rate = 0.5*(l30 + l31)
    rhs_lambda_rate = lambda_l2_contribution_rate + lambda_l3_contribution_rate

    a1 = _finite("final_a", final.geometry.a)
    b1 = _finite("final_b", final.geometry.b)
    lam1 = _finite("final_Lambda", final.geometry.Lambda)
    lm1 = lambda_metric(grid, final.geometry)
    c1 = lam1 - lm1

    actual_lambda_rate = (lam1 - lam0)/dt
    metric_lambda_rate = (lm1 - lm0)/dt
    constraint_rate = (c1 - c0)/dt
    rhs_closure = actual_lambda_rate - rhs_lambda_rate
    connection_closure = constraint_rate - (actual_lambda_rate - metric_lambda_rate)

    # The predictor/final metric-connection changes can be inspected separately.
    a_predictor = a0 + dt*stage_terms[0]["explicit_a"]
    b_predictor = b0 + dt*stage_terms[0]["explicit_b"]
    predictor_geom = state.geometry.copy()
    predictor_geom.a = a_predictor
    predictor_geom.b = b_predictor
    lm_predictor = lambda_metric(grid, predictor_geom)
    predictor_delta_lm = lm_predictor - lm0

    masks = {
        "first_four_cells": np.arange(n) < min(4, n),
        "inner_r_lt_5": r < R_CUT,
        "exterior_r_ge_5": r >= R_CUT,
        "full_domain": np.ones(n, dtype=bool),
    }
    arrays = {
        "C_before": c0,
        "C_after": c1,
        "delta_C": c1-c0,
        "delta_C_rate": constraint_rate,
        "Lambda_rate_actual": actual_lambda_rate,
        "Lambda_rate_rhs_split": rhs_lambda_rate,
        "Lambda_rate_rhs_l2": lambda_l2_contribution_rate,
        "Lambda_rate_rhs_l3": lambda_l3_contribution_rate,
        "Lambda_metric_rate": metric_lambda_rate,
        "Lambda_metric_rate_predictor": predictor_delta_lm/dt,
        "closure_actual_vs_rhs": rhs_closure,
        "closure_C_identity": connection_closure,
    }
    for name, values in arrays.items():
        _finite(name, values)

    summary = {
        "case": label,
        "resolution": n,
        "base_t": float(state.t),
        "dt": float(dt),
        "returned_t": float(final.t),
        "stage_geometry_call_count": len(stage_terms),
        "observed_final_lambda_l2_call_count": len(outside_lambda_l2),
        "interpretation_limit": (
            "One production step from a copied state; rate decomposition is "
            "numerical bookkeeping and does not change production equations."
        ),
        "regions": {
            region: {
                name: peak(values, r, mask)
                for name, values in arrays.items()
            }
            for region, mask in masks.items()
        },
        "cells_0_to_4": {
            str(i): {
                "r": float(r[i]),
                **{name: float(values[i]) for name, values in arrays.items()},
            }
            for i in range(min(5, n))
        },
        "component_identity": {
            "max_abs_actual_lambda_rate_minus_split_rhs": float(np.max(np.abs(rhs_closure))),
            "max_abs_connection_rate_identity_error": float(np.max(np.abs(connection_closure))),
            "expected_roundoff_scale": 1.e-10,
            "pass": bool(
                np.max(np.abs(rhs_closure)) < 1.e-10
                and np.max(np.abs(connection_closure)) < 1.e-10
            ),
        },
    }
    profiles = {"r": r}
    profiles.update(arrays)
    return summary, f"N{n}_{label}", profiles


def main():
    kernel = V55ProductionKernel()
    cases = []
    all_profiles = {}
    for n in (60, 80):
        initial = kernel.initialize(resolution=n, r_max=80.0, include_radiation=True)
        summary, key, profiles = audit_case(kernel, initial, "INITIAL")
        cases.append(summary)
        all_profiles.update({f"{key}_{k}": v for k, v in profiles.items()})

        checkpoint = Path("runs/t12") / f"N{n}" / "checkpoint_t00010.000000.npz"
        if not checkpoint.exists():
            raise FileNotFoundError(checkpoint)
        saved = restore(checkpoint, kernel)
        summary, key, profiles = audit_case(kernel, saved, "SAVED_T10")
        cases.append(summary)
        all_profiles.update({f"{key}_{k}": v for k, v in profiles.items()})

    report = {
        "schema": "v55_lambda_connection_rate_budget_v1",
        "status": "completed",
        "diagnostic_only": True,
        "production_steps_run": 4,
        "source_commit": "6d8b3f2790644591326df57c55908863f994cb0c",
        "campaign_run_id": 38208253687,
        "cases": cases,
        "claim_limit": (
            "This tests discrete constraint propagation on one step from frozen states. "
            "It does not assert that any single RHS term is the continuum cause of a "
            "multi-step Hamiltonian constraint violation."
        ),
    }
    out = Path("results/lambda-connection-budget")
    out.mkdir(parents=True, exist_ok=True)
    (out/"lambda_connection_budget.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n"
    )
    np.savez_compressed(out/"lambda_connection_profiles.npz", **all_profiles)
    print("LAMBDA_CONNECTION_BUDGET "+json.dumps({
        "status":report["status"],
        "cases":[{"case":x["case"],"N":x["resolution"],"t":x["base_t"],
                  "dt":x["dt"],"closure":x["component_identity"]["max_abs_connection_rate_identity_error"],
                  "pass":x["component_identity"]["pass"]} for x in cases]
    },sort_keys=True))


if __name__ == "__main__":
    main()
