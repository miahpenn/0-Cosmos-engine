#!/usr/bin/env python3
"""Frozen-state audit of the production Lambda RHS against the discrete metric connection.

This is a diagnostic only. It evaluates the exact pinned production RHS at four
frozen states (N=60/80 at t=0 and archived t~10), then takes one ordinary
production step from each independent state to compare the measured change in
C = Lambda - Lambda_metric with the instantaneous RHS budget. No RHS is altered.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
CAMPAIGN_RUN_ID = 38208253687
LAMBDA_M = 2.0
R_CUT = 5.0

import sys
SOURCE = Path("source").resolve()
sys.path.insert(0, str(SOURCE))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine.stress_energy import assemble_total_stress_energy
from engine import v55_pirk_adapter as adapter

_, vacuum, _ = adapter.vendor_modules()

GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED_FIELDS = ("rest", "energy_t", "momentum_r")


def d1(grid, values, parity=1):
    return np.asarray(grid.cell_derivative_fourth(
        np.asarray(values, dtype=float), parity=parity
    ), dtype=float)


def metric_connection(grid, geometry):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    return d1(grid, a, +1) / (2.0 * a**2) - d1(grid, b, +1) / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)


def metric_connection_rhs(grid, geometry, explicit_rhs):
    """Exact time derivative of the discrete Lambda_metric functional under a_t,b_t."""
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    da = np.asarray(explicit_rhs["a"], dtype=float)
    db = np.asarray(explicit_rhs["b"], dtype=float)
    ap = d1(grid, a, +1)
    bp = d1(grid, b, +1)
    dap = d1(grid, da, +1)
    dbp = d1(grid, db, +1)

    # Differentiate Lambda_metric = D(a)/(2a^2) - D(b)/(ab)
    #                           + 2/r * (1/b - 1/a).
    # D is a fixed linear finite-difference operator, so D(a_t) and D(b_t)
    # are the exact derivatives of the discrete derivative terms.
    return (
        dap / (2.0 * a**2)
        - ap * da / a**3
        - dbp / (a * b)
        + bp * da / (a**2 * b)
        + bp * db / (a * b**2)
        + 2.0 / r * (-db / b**2 + da / a**2)
    )


def restore_checkpoint(path: Path, kernel):
    with np.load(path, allow_pickle=False) as archive:
        data = {key: archive[key].copy() for key in archive.files}
    n = int(data["r"].size)
    dr = float(data["r"][1] - data["r"][0])
    r_max = float(n * dr)
    state = kernel.initialize(resolution=n, r_max=r_max, include_radiation=True)
    if not np.allclose(state.grid.centers, data["r"], rtol=0.0,
                       atol=1.0e-12 * max(1.0, r_max)):
        raise ValueError(f"{path}: checkpoint grid differs from reconstructed grid")
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


def summarize(values, radius, mask):
    v = np.asarray(values, dtype=float)
    selected = np.asarray(mask, dtype=bool) & np.isfinite(v)
    indices = np.flatnonzero(selected)
    if not indices.size:
        return {"count": 0}
    local = v[indices]
    peak_local = int(np.argmax(np.abs(local)))
    peak_index = int(indices[peak_local])
    return {
        "count": int(indices.size),
        "max_abs": float(abs(v[peak_index])),
        "peak_signed": float(v[peak_index]),
        "peak_r": float(radius[peak_index]),
        "peak_i": peak_index,
        "rms": float(np.sqrt(np.mean(local**2))),
        "mean_abs": float(np.mean(np.abs(local))),
    }


def analyze_case(state, label, kernel):
    # Production step always solves the CMC lapse and zeros shift/driver state
    # on its start slice. Mirror that preparation for the frozen RHS evaluation.
    grid = state.grid
    geom = state.geometry.copy()
    geom.alpha = kernel._solve_lapse(grid, geom, state.scalars, state.matter)[0]
    geom.beta.fill(0.0)
    geom.B.fill(0.0)

    terms = adapter.geometry_stage_terms(
        grid, geom, state.scalars, state.matter, lambda_m=LAMBDA_M
    )
    radius = np.asarray(grid.centers, dtype=float)
    metric_lam = metric_connection(grid, geom)
    c_constraint = np.asarray(geom.Lambda, dtype=float) - metric_lam

    lam_rhs_l2 = np.asarray(terms["lambda_l2"], dtype=float)
    lam_rhs_l3 = np.asarray(terms["lambda_l3"], dtype=float)
    lam_rhs = lam_rhs_l2 + lam_rhs_l3
    lam_metric_rhs = metric_connection_rhs(grid, geom, terms["explicit"])

    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    raw = vacuum.constraints(grid, geom)
    m_geometry = np.asarray(raw["momentum"], dtype=float)
    m_source = 8.0 * math.pi * np.asarray(total.j, dtype=float)
    m_residual = m_geometry - m_source
    h_residual = np.asarray(raw["hamiltonian"], dtype=float) - 16.0 * math.pi * np.asarray(total.rho, dtype=float)

    # The vendor Lambda evolution contains lambda_m*alpha*M_geometry/a and
    # the local adapter adds -8*pi*lambda_m*alpha*j/a. Isolate that documented
    # momentum-constraint addition; do not remove it from any evolved run.
    momentum_addition = LAMBDA_M * np.asarray(geom.alpha, dtype=float) * m_residual / np.asarray(geom.a, dtype=float)
    lam_rhs_without_momentum_addition = lam_rhs - momentum_addition

    c_rhs_full = lam_rhs - lam_metric_rhs
    c_rhs_without_momentum_addition = lam_rhs_without_momentum_addition - lam_metric_rhs

    rhs_arrays = {
        "connection_constraint_C": c_constraint,
        "hamiltonian_residual_H": h_residual,
        "momentum_residual_M": m_residual,
        "lambda_rhs_l2": lam_rhs_l2,
        "lambda_rhs_l3_including_matter": lam_rhs_l3,
        "lambda_rhs_total": lam_rhs,
        "metric_connection_rhs_from_a_b": lam_metric_rhs,
        "momentum_constraint_addition_lambda_m_alpha_M_over_a": momentum_addition,
        "connection_constraint_rhs_full": c_rhs_full,
        "connection_constraint_rhs_without_momentum_addition": c_rhs_without_momentum_addition,
    }
    for name, arr in rhs_arrays.items():
        if arr.shape != radius.shape or not np.all(np.isfinite(arr)):
            raise FloatingPointError(f"{label}: nonfinite/malformed {name}")

    # One ordinary production step from this input only. The RHS above was
    # already recorded on a copy, so state inputs remain frozen and independent.
    dt = 0.03 * float(grid.dr)
    result = kernel.step(state, dt)
    c_after = np.asarray(result.geometry.Lambda, dtype=float) - metric_connection(grid, result.geometry)
    c_step_rate = (c_after - c_constraint) / dt
    rhs_arrays["one_step_measured_delta_C_over_dt"] = c_step_rate

    masks = {
        "full_domain": np.ones_like(radius, dtype=bool),
        "r_lt_5": radius < R_CUT,
        "r_ge_5": radius >= R_CUT,
        "formal_stencil_i_ge_2": np.arange(radius.size) >= 2,
    }
    summaries = {
        field: {region: summarize(values, radius, mask) for region, mask in masks.items()}
        for field, values in rhs_arrays.items()
    }
    first_cells = {}
    for i in range(min(5, radius.size)):
        first_cells[str(i)] = {
            "r": float(radius[i]),
            **{name: float(values[i]) for name, values in rhs_arrays.items()},
        }

    # Algebraic closure checks for the RHS split. These should be roundoff-level
    # and detect accidental sign/factor loss in the diagnostic itself.
    split_closure = float(np.max(np.abs(
        c_rhs_full - c_rhs_without_momentum_addition - momentum_addition
    )))
    total_lambda_closure = float(np.max(np.abs(lam_rhs - lam_rhs_l2 - lam_rhs_l3)))

    return {
        "case": label,
        "resolution": int(radius.size),
        "r_max": float(radius.size * grid.dr),
        "dr": float(grid.dr),
        "base_t": float(state.t),
        "dt_one_step": dt,
        "returned_t": float(result.t),
        "stage_semidiscrete_rhs": {
            "interpretation": (
                "Frozen-input instantaneous RHS comparison. connection_constraint_rhs_full "
                "is d(Lambda - Lambda_metric)/dt implied by the production Lambda RHS and "
                "the exact derivative of the discrete metric-connection functional under the "
                "production explicit a,b RHS. The without-momentum-addition column is diagnostic "
                "algebra only; it is not an evolved counterfactual."
            ),
            "summaries": summaries,
            "first_five_cells": first_cells,
            "closures": {
                "momentum_addition_split_max_abs": split_closure,
                "Lambda_RHS_split_max_abs": total_lambda_closure,
            },
        },
        "one_step": {
            "steps_taken": 1,
            "physics_or_RHS_modified": False,
            "measured_delta_C_over_dt_summaries": {
                region: summarize(c_step_rate, radius, mask)
                for region, mask in masks.items()
            },
            "measured_delta_C_over_dt_first_five_cells": {
                str(i): float(c_step_rate[i]) for i in range(min(5, radius.size))
            },
            "input_and_output_C_first_five_cells": {
                str(i): {
                    "input_C": float(c_constraint[i]),
                    "output_C": float(c_after[i]),
                    "delta_C": float(c_after[i] - c_constraint[i]),
                }
                for i in range(min(5, radius.size))
            },
        },
    }


def main():
    kernel = V55ProductionKernel()
    cases = []
    for n in (60, 80):
        initial = kernel.initialize(resolution=n, r_max=80.0, include_radiation=True)
        cases.append(analyze_case(initial, f"N{n}_INITIAL_T0", kernel))
        checkpoint = Path(f"runs/t12/N{n}/checkpoint_t00010.000000.npz")
        if not checkpoint.exists():
            raise FileNotFoundError(f"required checkpoint missing: {checkpoint}")
        saved = restore_checkpoint(checkpoint, kernel)
        cases.append(analyze_case(saved, f"N{n}_SAVED_T10", kernel))

    report = {
        "schema": "lambda_metric_connection_rhs_consistency_v1",
        "status": "completed",
        "kind": "diagnostic_only_frozen_rhs_plus_one_step_per_input_state",
        "source_commit": SOURCE_COMMIT,
        "campaign_run_id": CAMPAIGN_RUN_ID,
        "settings": {
            "resolutions": [60, 80],
            "r_max": 80.0,
            "cfl_for_one_step": 0.03,
            "lambda_m": LAMBDA_M,
            "checkpoint_inputs": "Archived checkpoints near t=10 from the admitted execution campaign",
            "equation_changes": False,
            "production_main_modified": False,
        },
        "cases": cases,
        "interpretation_guardrail": (
            "This audit identifies semidiscrete RHS consistency and momentum-constraint "
            "addition contributions. A nonzero result is not alone proof of an incorrect "
            "continuum equation because the connection RHS deliberately uses the momentum "
            "constraint in this BSSN formulation. Do not flip signs, alter lambda_m, replace "
            "Ricci, or change production equations from this audit alone."
        ),
    }
    out = Path("results/lambda-metric-connection-rhs")
    out.mkdir(parents=True, exist_ok=True)
    path = out / "lambda_metric_connection_rhs.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("LAMBDA_METRIC_CONNECTION_RHS_AUDIT " + json.dumps({
        "status": report["status"],
        "source_commit": SOURCE_COMMIT,
        "cases": [{
            "case": case["case"],
            "dt": case["dt_one_step"],
            "formal_stencil_Cdot": case["stage_semidiscrete_rhs"]["summaries"]["connection_constraint_rhs_full"]["formal_stencil_i_ge_2"],
            "formal_stencil_Cdot_without_M": case["stage_semidiscrete_rhs"]["summaries"]["connection_constraint_rhs_without_momentum_addition"]["formal_stencil_i_ge_2"],
            "formal_stencil_M": case["stage_semidiscrete_rhs"]["summaries"]["momentum_residual_M"]["formal_stencil_i_ge_2"],
            "step_delta_C_over_dt": case["one_step"]["measured_delta_C_over_dt_summaries"]["formal_stencil_i_ge_2"],
            "i2": case["stage_semidiscrete_rhs"]["first_five_cells"].get("2"),
            "closures": case["stage_semidiscrete_rhs"]["closures"],
        } for case in cases]
    }, sort_keys=True))


if __name__ == "__main__":
    main()
