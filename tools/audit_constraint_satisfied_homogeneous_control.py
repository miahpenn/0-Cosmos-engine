#!/usr/bin/env python3
"""Source-pinned homogeneous constraint-satisfying control for Issue #4.

This is a null control, not a nontrivial propagation oracle. It creates a
regular spatially flat, homogeneous FLRW slice with H=M=C_Lambda=det=0
under the repository's native discrete constraints, measures the production
semi-discrete momentum tangent by epsilon-halved directional perturbations,
and performs bounded dt/dt-half checks. No source physics is changed.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import initialize_from_archive
from engine import v55_pirk_adapter as adapter
from tools.audit_current_source_momentum_constraint_tangent import (
    SOURCE_BASELINE,
    EXPECTED_VENDOR,
    provenance,
    state_digest,
    observe_constraint,
    rhs_direction,
    tangent_for_epsilon,
    one_step_check,
    regional_summaries,
    fixed_radius_samples,
    finite,
)

SCHEMA = "issue4_homogeneous_constraint_satisfied_control_v1"
PHI0 = 0.179055
PI0 = 0.00341913
MACHINE_EPS = float(np.finfo(np.float64).eps)


def max_abs(values: Any) -> float:
    arr = np.asarray(values, dtype=np.float64)
    return float(np.max(np.abs(arr))) if arr.size else 0.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def make_homogeneous_state(
    kernel: V55ProductionKernel,
    resolution: int,
    r_max: float,
) -> tuple[ProductionState, dict[str, Any]]:
    """Construct a regular flat FLRW slice, normalizing K from total rho."""
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(resolution, r_max)
    geom = vacuum.flat_state(grid)

    # All conformal metric/connection fields are exactly flat and regular.
    geom.a[:] = 1.0
    geom.b[:] = 1.0
    geom.X[:] = 1.0
    geom.alpha[:] = 1.0
    geom.beta.fill(0.0)
    geom.Aa.fill(0.0)
    geom.K.fill(0.0)
    geom.Lambda.fill(0.0)
    geom.B.fill(0.0)

    scalars = ScalarFields(
        S=np.zeros(resolution, dtype=np.float64),
        PS=np.zeros(resolution, dtype=np.float64),
        D=np.zeros(resolution, dtype=np.float64),
        PD=np.zeros(resolution, dtype=np.float64),
        phi=np.full(resolution, PHI0, dtype=np.float64),
        Pi=np.full(resolution, PI0, dtype=np.float64),
    )
    matter = initialize_from_archive(grid, geom, include_radiation=True)
    rho_pre = np.asarray(
        assemble_total_stress_energy(grid, geom, scalars, matter).rho,
        dtype=np.float64,
    )
    if rho_pre.shape != (resolution,) or not np.all(np.isfinite(rho_pre)):
        raise FloatingPointError("total energy density is malformed")
    rho_scale = max(max_abs(rho_pre), np.finfo(float).tiny)
    rho_nonuniformity = max_abs(rho_pre - float(np.mean(rho_pre))) / rho_scale
    if rho_nonuniformity > 1.0e-12:
        raise AssertionError(
            f"homogeneous fixture energy density is not homogeneous: {rho_nonuniformity}"
        )
    rho0 = float(np.mean(rho_pre))
    if rho0 <= 0.0:
        raise AssertionError("homogeneous control requires strictly positive total density")

    # With R=0 and A_ij=0, H_geom = 2 K^2/3; choose K so
    # 2 K^2/3 = 16 pi rho_total. The negative branch matches contraction.
    K0 = -math.sqrt(24.0 * math.pi * rho0)
    geom.K[:] = K0
    state = ProductionState(
        grid=grid, geometry=geom, scalars=scalars, matter=matter
    )

    # Recompute the production CMC lapse after the homogeneous constraint solve.
    state.geometry.alpha = np.asarray(
        kernel._solve_lapse(grid, state.geometry, scalars, matter)[0],
        dtype=np.float64,
    )
    state.geometry.beta.fill(0.0)
    state.geometry.B.fill(0.0)
    state.geometry.assert_finite_positive()

    return state, {
        "density_mean": rho0,
        "density_relative_nonuniformity": float(rho_nonuniformity),
        "K_homogeneous": float(K0),
        "parity_class": (
            "a,b,X,alpha,K,phi,Pi and fluid scalars even; Aa=Lambda=beta=0"
        ),
        "metric_construction": "a=b=X=alpha=1, beta=Aa=Lambda=B=0 before CMC re-solve",
    }


def constraint_witness(state: ProductionState) -> dict[str, Any]:
    _, vacuum, _ = adapter.vendor_modules()
    grid, geom = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=np.float64)
    raw = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, state.scalars, state.matter)
    H = finite(
        "homogeneous Hamiltonian residual",
        np.asarray(raw["hamiltonian"], dtype=np.float64)
        - 16.0 * math.pi * np.asarray(total.rho, dtype=np.float64),
    )
    M = finite(
        "homogeneous momentum residual",
        np.asarray(raw["momentum"], dtype=np.float64)
        - 8.0 * math.pi * np.asarray(total.j, dtype=np.float64),
    )
    d1 = lambda values: grid.cell_derivative_fourth(values, parity=1)
    a = np.asarray(geom.a, dtype=np.float64)
    b = np.asarray(geom.b, dtype=np.float64)
    X = np.asarray(geom.X, dtype=np.float64)
    Lambda_metric = (
        d1(a) / (2.0 * a * a)
        - d1(b) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    C = finite(
        "homogeneous connection residual",
        np.asarray(geom.Lambda, dtype=np.float64) - Lambda_metric,
    )
    det = finite("homogeneous determinant residual", raw["determinant"])
    M_observer = observe_constraint(state)["M_total"]
    return {
        "max_abs": {
            "H": max_abs(H),
            "M_direct": max_abs(M),
            "M_observer": max_abs(M_observer),
            "C_Lambda": max_abs(C),
            "determinant": max_abs(det),
            "alpha_minus_one": max_abs(np.asarray(geom.alpha) - 1.0),
            "a_b2_minus_one": max_abs(a * b * b - 1.0),
            "b_minus_X2": max_abs(b - X * X),
        },
        "first_four_cells": [
            {
                "i": i,
                "r": float(r[i]),
                "H": float(H[i]),
                "M": float(M[i]),
                "C_Lambda": float(C[i]),
                "determinant": float(det[i]),
                "alpha": float(geom.alpha[i]),
            }
            for i in range(min(4, r.size))
        ],
        "fixed_radius_samples": {
            str(radius): {
                "H": float(np.interp(radius, r, H)),
                "M": float(np.interp(radius, r, M)),
                "C_Lambda": float(np.interp(radius, r, C)),
                "determinant": float(np.interp(radius, r, det)),
            }
            for radius in (2.0, 3.0, 5.0, 8.0)
            if r[0] <= radius <= r[-1]
        },
        "constraint_definitions": {
            "H": "vendor geometric Hamiltonian minus 16*pi times assembled total rho",
            "M": "vendor geometric momentum constraint minus 8*pi times assembled total j",
            "C_Lambda": "evolved Lambda minus production metric-derived Lambda identity",
            "determinant": "native vendor determinant constraint output",
        },
    }


def run_resolution(
    kernel: V55ProductionKernel,
    resolution: int,
    r_max: float,
    epsilons: tuple[float, float],
    dt: float,
) -> dict[str, Any]:
    state, construction = make_homogeneous_state(kernel, resolution, r_max)
    initial_digest = state_digest(state)
    r = np.asarray(state.grid.centers, dtype=np.float64)
    constraints = constraint_witness(state)
    rates = rhs_direction(kernel, state)
    tangents = [
        tangent_for_epsilon(kernel, state, rates, eps)
        for eps in epsilons
    ]
    coarse, fine = tangents
    delta = np.asarray(coarse["Mdot"]) - np.asarray(fine["Mdot"])
    m0 = observe_constraint(state)["M_total"]

    # Each accepted-step check gets an isolated copy: neither the full step nor
    # the half-step is allowed to alter the frozen witness state.
    full_step = one_step_check(
        kernel, copy.deepcopy(state), m0, np.asarray(fine["Mdot"]), dt
    )
    half_step = one_step_check(
        kernel, copy.deepcopy(state), m0, np.asarray(fine["Mdot"]), 0.5 * dt
    )
    final_digest = state_digest(state)
    if final_digest != initial_digest:
        raise RuntimeError("observer or bounded step checks mutated the frozen control state")

    return {
        "resolution": int(resolution),
        "r_max": float(r_max),
        "dr": float(state.grid.dr),
        "r0": float(r[0]),
        "r1": float(r[1]),
        "initial_state_sha256": initial_digest,
        "construction": construction,
        "constraints": constraints,
        "tangent_by_epsilon": [
            {
                "epsilon": item["epsilon"],
                "Mdot_max_abs": max_abs(item["Mdot"]),
                "Mdot_cells_0_3": [float(x) for x in np.asarray(item["Mdot"])[:4]],
                "Mdot_geometric_total_cells_0_3": [
                    float(x) for x in np.asarray(item["Mdot_geom"])[:4]
                ],
                "Mdot_reconstructed_cells_0_3": [
                    float(x) for x in np.asarray(item["Mdot_reconstructed"])[:4]
                ],
                "Mdot_geometry_term_contributions": {
                    name: {
                        "cells_0_3": [float(x) for x in np.asarray(values)[:4]],
                        "fixed_radius_samples": fixed_radius_samples(values, r),
                    }
                    for name, values in item["Mdot_geom_terms"].items()
                },
                "Mdot_current_sector_contributions": {
                    name: {
                        "cells_0_3": [float(x) for x in np.asarray(values)[:4]],
                        "fixed_radius_samples": fixed_radius_samples(values, r),
                    }
                    for name, values in item["Mdot_sector"].items()
                },
                "Mdot_fixed_radius_samples": fixed_radius_samples(item["Mdot"], r),
                "Mdot_regional_summaries": regional_summaries(item["Mdot"], r),
                "arithmetic_closure": {
                    "Mdot_split_ok": bool(item["closure_ok"]),
                    "geometric_split_ok": bool(item["Mdot_geom_tangent_closure_ok"]),
                    "current_projection_ok": bool(item["plus_minus_current_projection_closure_ok"]),
                    "max_abs": float(item["closure_max_abs"]),
                    "max_bound": float(item["closure_max_bound"]),
                },
            }
            for item in tangents
        ],
        "epsilon_halving": {
            "coarse_epsilon": float(epsilons[0]),
            "fine_epsilon": float(epsilons[1]),
            "Mdot_difference_max_abs": max_abs(delta),
            "Mdot_difference_cells_0_3": [float(x) for x in delta[:4]],
            "Mdot_difference_fixed_radius_samples": fixed_radius_samples(delta, r),
        },
        "bounded_step_checks": {
            "dt": full_step,
            "dt_half": half_step,
            "interpretation": (
                "null-control consistency only; a failed finite-step check is reported "
                "separately from algebraic constraint/tangent closure"
            ),
        },
        "state_immutability": {
            "initial_sha256": initial_digest,
            "final_sha256": final_digest,
            "unchanged": initial_digest == final_digest,
        },
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    pins = provenance()
    if pins.get("production_source_baseline_anchor") != SOURCE_BASELINE:
        raise RuntimeError("production source baseline anchor mismatch")
    if pins.get("vendor_commit") != EXPECTED_VENDOR:
        raise RuntimeError("vendor pin mismatch")
    epsilons = tuple(float(x) for x in args.epsilons)
    if len(epsilons) != 2 or any(not math.isfinite(x) or x <= 0 for x in epsilons):
        raise ValueError("two finite positive directional epsilons are required")
    if epsilons[0] <= epsilons[1]:
        raise ValueError("epsilon order must be coarse then fine")
    if not args.resolutions or any(n < 8 for n in args.resolutions):
        raise ValueError("resolutions must each be at least eight")
    if not math.isfinite(args.r_max) or args.r_max <= 10:
        raise ValueError("r_max must be finite and greater than 10")
    if not math.isfinite(args.dt) or args.dt <= 0:
        raise ValueError("dt must be finite and positive")

    kernel = V55ProductionKernel()
    rows = [
        run_resolution(kernel, n, args.r_max, epsilons, args.dt)
        for n in args.resolutions
    ]
    constraint_tolerance = 5.0e-11
    alpha_tolerance = 5.0e-10
    all_constraints_controlled = all(
        row["constraints"]["max_abs"]["H"] <= constraint_tolerance
        and row["constraints"]["max_abs"]["M_direct"] <= constraint_tolerance
        and row["constraints"]["max_abs"]["M_observer"] <= constraint_tolerance
        and row["constraints"]["max_abs"]["C_Lambda"] <= constraint_tolerance
        and row["constraints"]["max_abs"]["determinant"] <= constraint_tolerance
        and row["constraints"]["max_abs"]["alpha_minus_one"] <= alpha_tolerance
        and row["constraints"]["max_abs"]["a_b2_minus_one"] <= 1.0e-14
        and row["constraints"]["max_abs"]["b_minus_X2"] <= 1.0e-14
        for row in rows
    )
    all_tangent_closures = all(
        item["arithmetic_closure"]["Mdot_split_ok"]
        and item["arithmetic_closure"]["geometric_split_ok"]
        and item["arithmetic_closure"]["current_projection_ok"]
        for row in rows for item in row["tangent_by_epsilon"]
    )
    all_state_unchanged = all(row["state_immutability"]["unchanged"] for row in rows)
    all_steps_completed = all(
        row["bounded_step_checks"]["dt"].get("completed") is True
        and row["bounded_step_checks"]["dt_half"].get("completed") is True
        for row in rows
    )
    all_null_tangents_small = all(
        max(item["Mdot_max_abs"] for item in row["tangent_by_epsilon"]) <= 5.0e-9
        for row in rows
    )
    return {
        "schema": SCHEMA,
        "status": (
            "completed_diagnostic_only"
            if all_constraints_controlled and all_tangent_closures
            and all_state_unchanged and all_steps_completed and all_null_tangents_small
            else "diagnostic_gate_failed"
        ),
        "scientific_admission": False,
        "production_physics_changed": False,
        "production_steps": "one dt and one dt/2 check per resolution; no campaign",
        "source_pins": pins,
        "configuration": {
            "resolutions": [int(x) for x in args.resolutions],
            "r_max": float(args.r_max),
            "epsilons": list(epsilons),
            "dt": float(args.dt),
            "scalar_operating_point": {"phi": PHI0, "Pi": PI0},
            "radiation": True,
            "charged_dust": False,
            "radial_maxwell": False,
        },
        "acceptance_gates": {
            "constraint_tolerance": constraint_tolerance,
            "alpha_minus_one_tolerance": alpha_tolerance,
            "max_null_Mdot_tolerance": 5.0e-9,
            "all_constraints_controlled": bool(all_constraints_controlled),
            "all_tangent_closures_pass": bool(all_tangent_closures),
            "all_frozen_states_unchanged": bool(all_state_unchanged),
            "all_dt_and_dt_half_checks_completed": bool(all_steps_completed),
            "all_null_Mdot_below_tolerance": bool(all_null_tangents_small),
        },
        "resolutions": rows,
        "interpretation_limits": [
            "This is an independent homogeneous flat-FLRW null control; all spatial gradients of H and C_Lambda vanish.",
            "The measured M tangent is derived from the production coupled RHS by symmetric directional differences, not from the independent inhomogeneous BSSN propagation law.",
            "Passing verifies a baseline parity/normalization/immutability control only; it cannot close Issue #4 or validate the hybrid initializer.",
            "No production or vendor source is modified and no scientific admission is made.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=[60, 80, 120, 160])
    parser.add_argument("--r-max", type=float, default=80.0)
    parser.add_argument("--epsilons", nargs=2, type=float, default=[2.0e-4, 1.0e-4])
    parser.add_argument("--dt", type=float, default=0.015)
    parser.add_argument(
        "--output", type=Path,
        default=Path("artifacts/issue4-homogeneous-control/report.json"),
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        report = run(args)
        code = 0 if report["status"] == "completed_diagnostic_only" else 2
    except Exception as exc:
        report = {
            "schema": SCHEMA,
            "status": "failed_closed",
            "scientific_admission": False,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        code = 1
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
