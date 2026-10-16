#!/usr/bin/env python3
"""Diagnostic CMC-lapse centre-series and lapse-profile sensitivity audit.

Uses the exact pinned production source and vendor. The baseline arm is the
repository production initializer. A second counterfactual replaces only the
lapse profile by a locally flattened profile before evaluating the RHS; all
metric, scalar, and matter fields are otherwise held fixed. This is not a CMC
solution and is not a production change.

The user-supplied hybrid B/Newton state is not present as a pinned array artifact
in the repository, so this run must not be described as a reproduction of that
hybrid state. scientific_admission remains False.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse
from engine.matter_system import ConservedSpecies
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import V55MatterState, total_matter_projection
from tools import audit_current_source_momentum_constraint_tangent as tangent


SCHEMA = "cmc_lapse_centre_series_and_profile_sensitivity_v1"
SOURCE_BASELINE = "7f47669a1c0edfafef9fc5f27b2db80840d3f035"
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
DEFAULT_RESOLUTIONS = (60, 80, 120, 160)
DEFAULT_EPSILONS = (2.0e-5, 1.0e-5)
MACHINE_EPS = float(np.finfo(np.float64).eps)


def finite(name: str, values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if not bool(np.all(np.isfinite(arr))):
        raise FloatingPointError(f"non-finite values in {name}")
    return arr


def compact_profile(values: Any, r: np.ndarray) -> dict[str, Any]:
    arr = finite("compact profile", values)
    return {
        "cell_0": float(arr[0]),
        "cell_1": float(arr[1]) if arr.size > 1 else None,
        "radius_cell_0": float(r[0]),
        "radius_cell_1": float(r[1]) if r.size > 1 else None,
        "fixed_radius_samples": tangent.fixed_radius_samples(arr, r),
    }


def c_m_c_residual(state: ProductionState) -> dict[str, Any]:
    """Evaluate the same tridiagonal row equation used by solve_cmc_lapse."""
    grid, geom = state.grid, state.geometry
    r = finite("CMC radii", grid.centers)
    h = float(grid.dr)
    a = finite("CMC a", geom.a)
    b = finite("CMC b", geom.b)
    X = finite("CMC X", geom.X)
    alpha = finite("CMC alpha", geom.alpha)
    if h <= 0.0 or alpha.shape != r.shape or r.size < 4:
        raise ValueError("invalid grid/lapse shape for CMC residual")

    d1 = lambda f: grid.cell_derivative_fourth(f, parity=1)
    ap, bp, Xp = d1(a), d1(b), d1(X)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r
    inv = X * X / a
    total = total_matter_projection(grid, geom, state.scalars, state.matter)
    Q = (
        1.5 * np.asarray(geom.Aa, dtype=np.float64) ** 2
        + np.asarray(geom.K, dtype=np.float64) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"], dtype=np.float64)
            + np.asarray(total["pr"], dtype=np.float64)
            + 2.0 * np.asarray(total["pt"], dtype=np.float64)
        )
    )
    # Re-solve only to obtain the CMC target rate used by the production solver.
    alpha_check, used_kdot = solve_cmc_lapse(
        grid, geom.copy(), state.scalars, state.matter
    )
    m = Q / inv
    residual = np.zeros_like(alpha)
    # Solver centre row: alpha[0] - alpha[1] = 0; outer row: alpha[-1] = 1.
    residual[0] = alpha[0] - alpha[1]
    residual[-1] = alpha[-1] - 1.0
    for i in range(1, r.size - 1):
        d2 = (alpha[i - 1] - 2.0 * alpha[i] + alpha[i + 1]) / (h * h)
        d1_centered = (alpha[i + 1] - alpha[i - 1]) / (2.0 * h)
        lhs = d2 + c[i] * d1_centered - m[i] * alpha[i]
        residual[i] = lhs + float(used_kdot) / inv[i]
    alpha_difference = alpha - np.asarray(alpha_check, dtype=np.float64)
    return {
        "target_kdot": float(used_kdot),
        "max_abs_solver_row_residual_interior": float(np.max(np.abs(residual[1:-1]))),
        "centre_row_residual": float(residual[0]),
        "outer_normalization_residual": float(residual[-1]),
        "alpha_rerun_max_abs_difference": float(np.max(np.abs(alpha_difference))),
        "solver_row_residual_witness": compact_profile(residual, r),
    }


def lapse_series(state: ProductionState) -> dict[str, Any]:
    """Save CMC alpha values/derivatives using the exact vendor geometry terms."""
    grid, geom = state.grid, state.geometry
    r = finite("lapse-series radii", grid.centers)
    alpha = finite("lapse-series alpha", geom.alpha)
    if np.any(alpha <= 0.0):
        raise FloatingPointError("lapse-series alpha must be positive")
    _, vacuum, _ = adapter.vendor_modules()
    vendor_geometry = vacuum.geometry_terms(grid, geom)
    alpha_r = finite(
        "cell derivative alpha",
        grid.cell_derivative_fourth(alpha, parity=1),
    )
    alpha_hessian_rr = finite(
        "vendor DrDr_alpha", vendor_geometry["DrDr_alpha"]
    )
    lap_alpha = finite("vendor lap_alpha", vendor_geometry["lap_alpha"])
    fit_n = min(8, r.size)
    z = r[:fit_n] ** 2
    degree = min(3, fit_n - 1)
    coeff = np.polynomial.polynomial.polyfit(z, alpha[:fit_n], deg=degree)
    fit = np.polynomial.polynomial.polyval(z, coeff)
    residual = alpha[:fit_n] - fit
    cells = []
    for i in range(min(8, r.size)):
        cells.append({
            "i": i,
            "r": float(r[i]),
            "alpha": float(alpha[i]),
            "D_r_alpha": float(alpha_r[i]),
            "vendor_DrDr_alpha": float(alpha_hessian_rr[i]),
            "vendor_lap_alpha": float(lap_alpha[i]),
            "fit_alpha_of_r2": float(fit[i]),
            "alpha_minus_even_fit": float(residual[i]) if i < fit_n else None,
        })
    return {
        "positive": bool(np.min(alpha) > 0.0),
        "alpha_min": float(np.min(alpha)),
        "alpha_max": float(np.max(alpha)),
        "even_series_fit": {
            "samples": fit_n,
            "degree_in_r_squared": degree,
            "coefficients_ascending": [float(x) for x in coeff],
            "alpha0_extrapolated": float(coeff[0]),
            "max_abs_fit_residual": float(np.max(np.abs(residual))),
            "rms_fit_residual": float(np.sqrt(np.mean(residual * residual))),
        },
        "centre_cells": cells,
        "fixed_radius_samples": {
            "alpha": tangent.fixed_radius_samples(alpha, r),
            "D_r_alpha": tangent.fixed_radius_samples(alpha_r, r),
            "vendor_DrDr_alpha": tangent.fixed_radius_samples(alpha_hessian_rr, r),
            "vendor_lap_alpha": tangent.fixed_radius_samples(lap_alpha, r),
        },
        "vendor_operator_note": (
            "DrDr_alpha and lap_alpha are returned by the pinned vendor geometry_terms; "
            "D_r_alpha uses the production cell_derivative_fourth parity=1 operator."
        ),
    }


def flattened_centre_alpha(
    alpha: np.ndarray, r: np.ndarray, flat_radius: float, taper_radius: float
) -> np.ndarray:
    """Counterfactual: constant alpha[0] near centre, smoothstep back to CMC alpha."""
    if not (
        math.isfinite(flat_radius)
        and math.isfinite(taper_radius)
        and 0.0 < flat_radius < taper_radius
    ):
        raise ValueError("require 0 < flat_radius < taper_radius")
    x = np.clip((r - flat_radius) / (taper_radius - flat_radius), 0.0, 1.0)
    blend = x * x * (3.0 - 2.0 * x)
    out = float(alpha[0]) + blend * (alpha - float(alpha[0]))
    out = finite("flattened lapse counterfactual", out)
    if np.any(out <= 0.0):
        raise FloatingPointError("flattened lapse counterfactual is not positive")
    return out


def with_lapse_profile(
    state: ProductionState, alpha: np.ndarray
) -> ProductionState:
    geom = state.geometry.copy()
    geom.alpha = finite("assigned lapse profile", alpha).copy()
    # Preserve the CMC branch's zero shift/damping variables.
    geom.beta.fill(0.0)
    geom.B.fill(0.0)
    return ProductionState(
        grid=state.grid,
        geometry=geom,
        scalars=state.scalars,
        matter=state.matter,
        t=state.t,
        tau=state.tau,
        e_folds=state.e_folds,
        electromagnetic=state.electromagnetic,
    )


def directional_state_with_fixed_lapse(
    state: ProductionState, rates: dict[str, Any], delta: float
) -> ProductionState:
    """State +/- delta*RHS with the arm's lapse profile held unchanged."""
    geom = state.geometry.copy()
    for name in tangent.GEOMETRY_RATE_NAMES:
        values = (
            np.asarray(getattr(state.geometry, name), dtype=np.float64)
            + delta * np.asarray(rates["geometry"][name], dtype=np.float64)
        )
        setattr(geom, name, finite(f"directional geometry {name}", values))
    geom.alpha = np.asarray(state.geometry.alpha, dtype=np.float64).copy()
    geom.beta.fill(0.0)
    geom.B.fill(0.0)

    scalars = ScalarFields(*(
        finite(
            f"directional scalar {name}",
            np.asarray(getattr(state.scalars, name), dtype=np.float64)
            + delta * np.asarray(rates["scalars"][name], dtype=np.float64),
        )
        for name in tangent.SCALAR_NAMES
    ))
    species = {}
    for species_name in ("dark_matter", "baryons", "radiation"):
        base = getattr(state.matter, species_name)
        rate = rates["matter"][species_name]
        species[species_name] = ConservedSpecies(*(
            finite(
                f"directional matter {species_name}.{name}",
                np.asarray(getattr(base, name), dtype=np.float64)
                + delta * np.asarray(rate[name], dtype=np.float64),
            )
            for name in ("rest", "energy_t", "momentum_r")
        ))
    matter = V55MatterState(
        species["dark_matter"], species["baryons"], species["radiation"]
    )
    return ProductionState(
        grid=state.grid,
        geometry=geom,
        scalars=scalars,
        matter=matter,
        t=state.t,
        tau=state.tau,
        e_folds=state.e_folds,
        electromagnetic=None,
    )


def tangent_for_fixed_profile(
    state: ProductionState, rates: dict[str, Any], epsilon: float
) -> dict[str, Any]:
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    plus = tangent.observe_constraint(
        directional_state_with_fixed_lapse(state, rates, epsilon)
    )
    minus = tangent.observe_constraint(
        directional_state_with_fixed_lapse(state, rates, -epsilon)
    )
    denom = 2.0 * epsilon
    mdot = (plus["M_total"] - minus["M_total"]) / denom
    mdot_geom = (plus["M_geom"] - minus["M_geom"]) / denom
    geom_dots = {
        key: (
            plus["M_geom_pieces"][key] - minus["M_geom_pieces"][key]
        ) / denom
        for key in tangent.GEOMETRIC_M_NAMES
    }
    j_dots = {
        key: (plus["J_pieces"][key] - minus["J_pieces"][key]) / denom
        for key in tangent.SECTOR_NAMES
    }
    reconstructed = mdot_geom.copy()
    for value in j_dots.values():
        reconstructed = reconstructed - 8.0 * math.pi * value
    closure = mdot - reconstructed
    operand_scale = (
        np.abs(plus["M_total"]) + np.abs(minus["M_total"])
        + np.abs(plus["M_geom"]) + np.abs(minus["M_geom"])
    )
    for key in tangent.SECTOR_NAMES:
        operand_scale = operand_scale + 8.0 * math.pi * (
            np.abs(plus["J_pieces"][key]) + np.abs(minus["J_pieces"][key])
        )
    bound = (128.0 * MACHINE_EPS / denom) * operand_scale
    closure_ok = bool(np.all(np.abs(closure) <= bound))
    r = np.asarray(state.grid.centers, dtype=np.float64)
    return {
        "epsilon": float(epsilon),
        "Mdot": mdot,
        "Mdot_geom": mdot_geom,
        "Mdot_geom_terms": geom_dots,
        "Mdot_current_terms": {
            key: -8.0 * math.pi * value for key, value in j_dots.items()
        },
        "closure": {
            "ok": closure_ok,
            "max_abs": float(np.max(np.abs(closure))),
            "max_bound": float(np.max(bound)),
        },
        "witnesses": {
            "Mdot": compact_profile(mdot, r),
            "Mdot_geom": compact_profile(mdot_geom, r),
            "Mdot_geom_terms": {
                key: compact_profile(value, r) for key, value in geom_dots.items()
            },
            "Mdot_current_terms": {
                key: compact_profile(-8.0 * math.pi * value, r)
                for key, value in j_dots.items()
            },
        },
    }


def kdot_channel_summary(
    state: ProductionState, rates: dict[str, Any]
) -> dict[str, Any]:
    grid, geom = state.grid, state.geometry
    terms = rates["geometry_stage_terms"]
    _, vacuum, _ = adapter.vendor_modules()
    l3_vacuum = vacuum.primary_l3_rhs(grid, geom)
    k_l2 = finite("K L2 RHS", terms["primary_l2"]["K"])
    k_l3_vac = finite("K vacuum L3 RHS", l3_vacuum["K"])
    k_l3_with_matter = finite("K L3 RHS with matter", terms["primary_l3"]["K"])
    k_matter = finite("K matter increment", k_l3_with_matter - k_l3_vac)
    channels = {
        "L2_minus_laplacian_alpha": k_l2,
        "vacuum_L3_quadratic_shift_and_curvature": k_l3_vac,
        "matter_increment": k_matter,
    }
    channels["total"] = k_l2 + k_l3_with_matter
    r = np.asarray(grid.centers, dtype=np.float64)
    out = {}
    for name, values in channels.items():
        derivative = finite(
            f"D_r Kdot channel {name}",
            grid.cell_derivative_fourth(values, parity=1),
        )
        mapped = -(2.0 / 3.0) * derivative
        out[name] = {
            "Kdot": compact_profile(values, r),
            "D_r_Kdot": compact_profile(derivative, r),
            "minus_two_thirds_D_r_Kdot": compact_profile(mapped, r),
        }
    return out


def rhs_arm(
    kernel: V55ProductionKernel,
    state: ProductionState,
    epsilons: tuple[float, float],
) -> dict[str, Any]:
    original_digest = tangent.state_digest(state)
    rates = tangent.rhs_direction(kernel, state)
    tangent_rows = [
        tangent_for_fixed_profile(state, rates, eps) for eps in epsilons
    ]
    r = np.asarray(state.grid.centers, dtype=np.float64)
    diff = tangent_rows[0]["Mdot"] - tangent_rows[1]["Mdot"]
    l2split = rates["Aa_l2_split"]
    output = {
        "state_sha256": original_digest,
        "alpha_profile_sha256": hashlib.sha256(
            np.ascontiguousarray(state.geometry.alpha, dtype=np.float64).tobytes()
        ).hexdigest(),
        "initial_constraints": {
            key: compact_profile(value, r)
            for key, value in tangent.observe_constraint(state).items()
            if isinstance(value, np.ndarray) and value.shape == r.shape
            and key in ("M_geom", "M_total", "J_total")
        },
        "Mdot_by_epsilon": [
            {
                "epsilon": row["epsilon"],
                "Mdot": row["witnesses"]["Mdot"],
                "Mdot_geom": row["witnesses"]["Mdot_geom"],
                "Mdot_geom_terms": row["witnesses"]["Mdot_geom_terms"],
                "Mdot_current_terms": row["witnesses"]["Mdot_current_terms"],
                "closure": row["closure"],
            }
            for row in tangent_rows
        ],
        "epsilon_sensitivity": compact_profile(diff, r),
        "kdot_channels": kdot_channel_summary(state, rates),
        "Aa_channels": {
            "Aa_dot": compact_profile(rates["geometry"]["Aa"], r),
            "lapse_hessian_TF": compact_profile(l2split["lapse_hessian"], r),
            "TF_Ricci": compact_profile(l2split["tracefree_ricci"], r),
            "L2_split_closure": l2split["closure"],
        },
        "RHS_direction_note": (
            "RHS rates are recomputed using this arm's lapse profile. The +/- "
            "directional states keep that same profile fixed; this estimates the "
            "constraint tangent along the selected RHS vector, without a lapse "
            "re-solve at the perturbation endpoints."
        ),
    }
    if tangent.state_digest(state) != original_digest:
        raise RuntimeError("RHS observer mutated its frozen input state")
    return output



def operator_compatible_cmc_arm(
    kernel: V55ProductionKernel,
    state: ProductionState,
    epsilons: tuple[float, float],
) -> dict[str, Any]:
    """Diagnostic-only CMC solve using the exact production D1/D2 operators.

    This copied solver leaves engine/cmc_gauge.py untouched. It solves the
    CMC equation on cells 0..n-2 with the same parity-aware fourth-order
    operators used by the production K RHS, includes the actual K-L3
    dissipation term, and retains alpha[-1]=1 as outer normalization.
    """
    try:
        grid, geom = state.grid, state.geometry
        r = finite("operator-compatible lapse radii", grid.centers)
        n = int(r.size)
        if n < 5:
            raise ValueError("operator-compatible lapse solve requires at least five cells")
        if float(np.max(np.abs(np.asarray(geom.beta, dtype=np.float64)))) != 0.0:
            raise ValueError("operator-compatible lapse probe expects the CMC zero-shift branch")

        # Use the existing production CMC projection's target Kdot on this
        # frozen state, without mutating the input or changing production.
        _, target_kdot = solve_cmc_lapse(
            grid, geom.copy(), state.scalars, state.matter
        )

        a = finite("operator-compatible metric a", geom.a)
        b = finite("operator-compatible metric b", geom.b)
        X = finite("operator-compatible metric X", geom.X)
        Aa = finite("operator-compatible metric Aa", geom.Aa)
        K = finite("operator-compatible metric K", geom.K)
        inv = finite("operator-compatible inverse physical radial metric", X * X / a)
        if np.any(inv <= 0.0):
            raise FloatingPointError("operator-compatible inverse radial metric must be positive")

        d1 = lambda values: finite(
            "operator-compatible D1", grid.cell_derivative_fourth(values, parity=1)
        )
        d2 = lambda values: finite(
            "operator-compatible D2", grid.cell_second_derivative_fourth(values, parity=1)
        )
        c = finite(
            "operator-compatible radial Laplacian coefficient",
            -0.5 * d1(a) / a + d1(b) / b - d1(X) / X + 2.0 / r,
        )

        total = total_matter_projection(grid, geom, state.scalars, state.matter)
        q = finite(
            "operator-compatible lapse potential",
            1.5 * Aa * Aa + K * K / 3.0
            + 4.0 * math.pi * (
                np.asarray(total["rho"], dtype=np.float64)
                + np.asarray(total["pr"], dtype=np.float64)
                + 2.0 * np.asarray(total["pt"], dtype=np.float64)
            ),
        )
        m = finite("operator-compatible scaled lapse potential", q / inv)

        # At zero shift the pinned vendor has K_L3=alpha*(Aa^2+2 Ab^2+
        # K^2/3)+KO(K). Extract KO from the actual vendor RHS, not a guessed
        # dissipation stencil, so the copied equation targets the full K RHS.
        _, vacuum, _ = adapter.vendor_modules()
        l3_vacuum = vacuum.primary_l3_rhs(grid, geom)
        quadratic_k = finite(
            "operator-compatible vacuum quadratic K source",
            Aa * Aa + 2.0 * (-0.5 * Aa) ** 2 + K * K / 3.0,
        )
        ko_k = finite(
            "operator-compatible K dissipation source",
            np.asarray(l3_vacuum["K"], dtype=np.float64)
            - np.asarray(geom.alpha, dtype=np.float64) * quadratic_k,
        )
        # Enforce Kdot=target_kdot: D2_h(alpha)+c_h*D1_h(alpha)-m*alpha
        # = (KO_K-target_kdot)/inv. The last row is outer normalization.
        rhs = finite(
            "operator-compatible lapse equation RHS",
            (ko_k - float(target_kdot)) / inv,
        )

        # Assemble D1 and D2 matrices independently through the actual
        # parity-aware stencil action on each basis vector. Keeping the
        # operators separate lets the roundoff envelope retain their raw
        # operands even when their outputs nearly cancel in the CMC equation.
        d1_matrix = np.zeros((n, n), dtype=np.float64)
        d2_matrix = np.zeros((n, n), dtype=np.float64)
        for j in range(n):
            basis = np.zeros(n, dtype=np.float64)
            basis[j] = 1.0
            d1_matrix[: n - 1, j] = d1(basis)[: n - 1]
            d2_matrix[: n - 1, j] = d2(basis)[: n - 1]
        matrix = d2_matrix + c[:, None] * d1_matrix - np.diag(m)
        matrix[n - 1, :] = 0.0
        matrix[n - 1, n - 1] = 1.0
        rhs_vector = rhs.copy()
        rhs_vector[n - 1] = 1.0
        alpha = finite(
            "operator-compatible lapse solution",
            np.linalg.solve(matrix, rhs_vector),
        )
        if np.any(alpha <= 0.0):
            raise ValueError(
                "operator-compatible lapse solution is not positive: "
                f"min_alpha={float(np.min(alpha)):.17e}"
            )

        d1_alpha = d1(alpha)
        d2_alpha = d2(alpha)
        lhs = finite(
            "operator-compatible lapse equation LHS",
            d2_alpha + c * d1_alpha - m * alpha,
        )
        operator_residual = lhs - rhs
        operator_residual[n - 1] = alpha[n - 1] - 1.0
        linear_residual = matrix @ alpha - rhs_vector
        raw_row_scale = (
            np.abs(d2_matrix) @ np.abs(alpha)
            + np.abs(c) * (np.abs(d1_matrix) @ np.abs(alpha))
            + np.abs(m * alpha) + np.abs(rhs)
        )
        raw_row_scale[n - 1] = abs(alpha[n - 1]) + 1.0
        row_bound = 512.0 * MACHINE_EPS * raw_row_scale
        linear_row_scale = np.abs(matrix) @ np.abs(alpha) + np.abs(rhs_vector)
        linear_bound = 512.0 * MACHINE_EPS * linear_row_scale
        operator_ratios = np.divide(
            np.abs(operator_residual), row_bound,
            out=np.zeros_like(operator_residual), where=row_bound > 0.0,
        )
        linear_ratios = np.divide(
            np.abs(linear_residual), linear_bound,
            out=np.zeros_like(linear_residual), where=linear_bound > 0.0,
        )
        i_ratio = int(np.argmax(operator_ratios))
        i_abs = int(np.argmax(np.abs(operator_residual)))
        i_linear_ratio = int(np.argmax(linear_ratios))
        if (
            not np.all(np.abs(operator_residual) <= row_bound)
            or not np.all(np.abs(linear_residual) <= linear_bound)
        ):
            raise FloatingPointError(
                "operator-compatible lapse residual exceeds the raw-operand "
                f"roundoff envelope; operator cell={i_ratio}, "
                f"linear cell={i_linear_ratio}"
            )

        diagnostic_state = with_lapse_profile(state, alpha)
        rates = tangent.rhs_direction(kernel, diagnostic_state)
        actual_kdot = finite(
            "operator-compatible production Kdot", rates["geometry"]["K"]
        )
        kdot_error = finite(
            "operator-compatible Kdot minus selected target",
            actual_kdot - float(target_kdot),
        )
        # Independent gate on the *actual production K RHS*, with a bound
        # formed from the raw D1/D2 operator operands rather than their
        # cancellation-prone outputs.
        lapse_operator_operand_scale = (
            np.abs(d2_matrix) @ np.abs(alpha)
            + np.abs(c) * (np.abs(d1_matrix) @ np.abs(alpha))
        )
        kdot_operand_scale = finite(
            "operator-compatible Kdot operand scale",
            np.abs(inv) * lapse_operator_operand_scale
            + np.abs(alpha * q) + np.abs(ko_k) + abs(float(target_kdot)),
        )
        kdot_bound = 1024.0 * MACHINE_EPS * kdot_operand_scale
        kdot_ratios = np.divide(
            np.abs(kdot_error), kdot_bound,
            out=np.zeros_like(kdot_error), where=kdot_bound > 0.0,
        )
        interior_error = kdot_error[: n - 1]
        i_kdot = int(np.argmax(np.abs(interior_error)))
        i_kdot_ratio = int(np.argmax(kdot_ratios[: n - 1]))
        kdot_closure_ok = bool(
            np.all(np.abs(interior_error) <= kdot_bound[: n - 1])
        )
        if not kdot_closure_ok:
            raise FloatingPointError(
                "actual production Kdot does not match the selected CMC target "
                f"within the raw-operand bound: cell={i_kdot_ratio}, "
                f"ratio={float(kdot_ratios[i_kdot_ratio]):.17e}"
            )
        original_digest = tangent.state_digest(state)
        arm = rhs_arm(kernel, diagnostic_state, epsilons)
        if tangent.state_digest(state) != original_digest:
            raise RuntimeError("operator-compatible lapse probe mutated its frozen input")

        return {
            "status": "completed_diagnostic_only",
            "production_source_lapse_solver_changed": False,
            "production_state_mutated": False,
            "formula": (
                "D2_h(alpha)+c_h*D1_h(alpha)-(Q/inv)*alpha "
                "=(KO_K-target_Kdot)/inv on cells 0..n-2; alpha[n-1]=1. "
                "D1_h and D2_h are the pinned fourth-order parity-aware "
                "operators; KO_K is extracted from the pinned vendor primary-L3 K RHS."
            ),
            "target_Kdot_from_existing_CMC_projection": float(target_kdot),
            "matrix_residual": {
                "operator_max_abs": float(np.max(np.abs(operator_residual))),
                "operator_max_abs_cell": i_abs,
                "operator_max_abs_cell_radius": float(r[i_abs]),
                "operator_max_local_scaled_ratio": float(np.max(operator_ratios)),
                "operator_worst_ratio_cell": i_ratio,
                "operator_worst_ratio_radius": float(r[i_ratio]),
                "linear_solve_max_abs": float(np.max(np.abs(linear_residual))),
                "linear_solve_max_local_scaled_ratio": float(np.max(linear_ratios)),
                "linear_solve_worst_ratio_cell": i_linear_ratio,
                "bound_all_rows_pass": bool(
                    np.all(np.abs(operator_residual) <= row_bound)
                    and np.all(np.abs(linear_residual) <= linear_bound)
                ),
                "pde_rows_count": int(n - 1),
                "outer_normalization_residual": float(operator_residual[n - 1]),
                "roundoff_scale": "512*eps times the sum of absolute D2/D1 matrix operands, potential term, and RHS; linear solve checked separately",
            },
            "lapse": {
                "minimum": float(np.min(alpha)),
                "maximum": float(np.max(alpha)),
                "max_abs_change_vs_production": float(np.max(np.abs(alpha - geom.alpha))),
                "profile": compact_profile(alpha, r),
                "change_vs_production": compact_profile(alpha - geom.alpha, r),
            },
            "actual_production_Kdot": compact_profile(actual_kdot, r),
            "Kdot_minus_target": {
                "closure_ok_in_PDE_rows": kdot_closure_ok,
                "max_abs_in_PDE_rows": float(np.max(np.abs(interior_error))),
                "worst_cell_in_PDE_rows": i_kdot,
                "worst_radius_in_PDE_rows": float(r[i_kdot]),
                "max_local_raw_operand_scaled_ratio": float(np.max(kdot_ratios[: n - 1])),
                "worst_scaled_ratio_cell": i_kdot_ratio,
                "worst_scaled_ratio_radius": float(r[i_kdot_ratio]),
                "max_roundoff_bound_in_PDE_rows": float(np.max(kdot_bound[: n - 1])),
                "profile": compact_profile(kdot_error, r),
            },
            "M_tangent": arm,
            "scientific_admission": False,
        }
    except Exception as exc:
        return {
            "status": "failed_closed",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "scientific_admission": False,
        }


def audit_resolution(
    kernel: V55ProductionKernel,
    resolution: int,
    r_max: float,
    epsilons: tuple[float, float],
    flat_radius: float,
    taper_radius: float,
) -> dict[str, Any]:
    state = kernel.initialize(
        resolution=resolution, r_max=r_max, include_radiation=True
    )
    initial_digest = tangent.state_digest(state)
    base_alpha = np.asarray(state.geometry.alpha, dtype=np.float64).copy()
    base_lapse_series = lapse_series(state)
    cmc = c_m_c_residual(state)

    flat_alpha = flattened_centre_alpha(
        base_alpha, np.asarray(state.grid.centers, dtype=np.float64),
        flat_radius, taper_radius,
    )
    flat_state = with_lapse_profile(state, flat_alpha)
    base_arm = rhs_arm(kernel, state, epsilons)
    flat_arm = rhs_arm(kernel, flat_state, epsilons)
    operator_compatible_arm = operator_compatible_cmc_arm(kernel, state, epsilons)
    r = np.asarray(state.grid.centers, dtype=np.float64)
    alpha_delta = flat_alpha - base_alpha
    output = {
        "resolution": int(resolution),
        "r_max": float(r_max),
        "dr": float(state.grid.dr),
        "r_cell_0": float(r[0]),
        "r_cell_1": float(r[1]),
        "initial_state_sha256": initial_digest,
        "initial_CMC_lapse": {
            "series": base_lapse_series,
            "CMC_solver_residual": cmc,
        },
        "counterfactual_lapse_profile": {
            "definition": (
                f"alpha_flat=alpha_CMC[0] for r <= {flat_radius}; cubic smoothstep "
                f"blend back to alpha_CMC over {flat_radius} < r < {taper_radius}; "
                f"alpha_flat=alpha_CMC for r >= {taper_radius}."
            ),
            "same_metric_scalar_matter_fields": True,
            "is_a_CMC_solution": False,
            "alpha_delta_max_abs": float(np.max(np.abs(alpha_delta))),
            "alpha_delta_rms": float(np.sqrt(np.mean(alpha_delta * alpha_delta))),
            "alpha_delta_profile": compact_profile(alpha_delta, r),
        },
        "arms": {
            "production_initializer_CMC_alpha": base_arm,
            "counterfactual_flattened_centre_alpha": flat_arm,
            "operator_compatible_fourth_order_cmc_lapse": operator_compatible_arm,
        },
        "input_mutation_check": {
            "production_state_unchanged": tangent.state_digest(state) == initial_digest
        },
        "scientific_admission": False,
    }
    if not output["input_mutation_check"]["production_state_unchanged"]:
        raise RuntimeError("production initial state was mutated by counterfactual")
    return output


def run(args: argparse.Namespace) -> dict[str, Any]:
    provenance = tangent.provenance()
    if provenance["vendor_commit"] != EXPECTED_VENDOR:
        raise RuntimeError("vendor revision differs from the declared expected pin")
    if len(args.epsilons) != 2 or any(not math.isfinite(e) or e <= 0 for e in args.epsilons):
        raise ValueError("provide exactly two positive finite directional epsilons")
    if not args.resolutions or any(n < 8 for n in args.resolutions):
        raise ValueError("resolutions must each be at least eight")
    if not math.isfinite(args.r_max) or args.r_max <= 0:
        raise ValueError("r_max must be finite and positive")
    kernel = V55ProductionKernel()
    rows = [
        audit_resolution(
            kernel, int(n), float(args.r_max),
            (float(args.epsilons[0]), float(args.epsilons[1])),
            float(args.flat_radius), float(args.taper_radius),
        )
        for n in args.resolutions
    ]
    baseline_closure = all(
        all(item["closure"]["ok"] for item in row["arms"]["production_initializer_CMC_alpha"]["Mdot_by_epsilon"])
        for row in rows
    )
    flattened_closure = all(
        all(item["closure"]["ok"] for item in row["arms"]["counterfactual_flattened_centre_alpha"]["Mdot_by_epsilon"])
        for row in rows
    )
    operator_compatible_lapse_constructed = all(
        row["arms"]["operator_compatible_fourth_order_cmc_lapse"].get("status") == "completed_diagnostic_only"
        for row in rows
    )
    return {
        "schema": SCHEMA,
        "status": "completed_diagnostic_only",
        "provenance": provenance,
        "configuration": {
            "resolutions": [int(n) for n in args.resolutions],
            "r_max": float(args.r_max),
            "directional_epsilons": [float(x) for x in args.epsilons],
            "flat_radius": float(args.flat_radius),
            "taper_radius": float(args.taper_radius),
            "initialization": "pinned repository production initializer; radiation enabled",
        },
        "experiment_scope": {
            "hybrid_state_reproduced": False,
            "reason": (
                "The user-supplied hybrid B/Newton result has no pinned state array, "
                "source patch, run configuration, or artifact hash in the repository. "
                "This audit therefore runs the exact pinned production-initializer "
                "baseline and a lapse-profile counterfactual only."
            ),
            "counterfactual_is_production_change": False,
            "evolution_steps_advanced": 0,
        },
        "closures": {
            "production_CMC_alpha_RHS_tangent": baseline_closure,
            "flattened_alpha_RHS_tangent": flattened_closure,
            "operator_compatible_fourth_order_cmc_lapse_constructed": operator_compatible_lapse_constructed,
        },
        "resolutions": rows,
        "scientific_admission": False,
        "interpretation_limits": [
            "This is a lapse-sensitivity counterfactual on the pinned production-initializer state, not a reproduction of the supplied hybrid state.",
            "Flattening alpha changes lapse values and derivatives within the declared region and breaks the CMC solve there; it tests RHS sensitivity, not a physically admissible gauge replacement.",
            "A changed Mdot under the counterfactual demonstrates lapse-profile sensitivity, not a faulty CMC solve or a propagation-law defect.",
            "Cell 0 is at r=dr/2 and changes physical radius under refinement; fixed-radius samples are reported separately.",
            "The independent, normalization-reviewed full constraint-propagation prediction remains the scientific acceptance gate.",
            "No production source, gauge, stencil, coefficient, boundary rule, or evolution equation is changed.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", nargs="+", type=int, default=list(DEFAULT_RESOLUTIONS))
    parser.add_argument("--r-max", type=float, default=80.0)
    parser.add_argument("--epsilons", nargs=2, type=float, default=list(DEFAULT_EPSILONS))
    parser.add_argument("--flat-radius", type=float, default=2.0)
    parser.add_argument("--taper-radius", type=float, default=5.0)
    parser.add_argument(
        "--output", type=Path,
        default=Path("results/cmc-lapse-centre-sensitivity/report.json"),
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        report = run(args)
        code = 0 if (
            report["closures"]["production_CMC_alpha_RHS_tangent"]
            and report["closures"]["flattened_alpha_RHS_tangent"]
            and report["closures"]["operator_compatible_fourth_order_cmc_lapse_constructed"]
        ) else 2
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
    print(json.dumps({
        "schema": report.get("schema"),
        "status": report.get("status"),
        "provenance": report.get("provenance"),
        "closures": report.get("closures"),
        "experiment_scope": report.get("experiment_scope"),
        "error_type": report.get("error_type"),
        "error": report.get("error"),
        "scientific_admission": False,
    }, indent=2, sort_keys=True), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
