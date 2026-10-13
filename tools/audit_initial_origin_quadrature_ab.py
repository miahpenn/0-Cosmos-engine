#!/usr/bin/env python3
"""Current-source t=0 A/B: production origin integration vs even-series quadrature.

The production trapezoid arm is rebuilt independently and must match the real
kernel initializer to a tight numerical bound before the alternative arm is
interpreted. Neither arm advances time. Only observer/workflow/test files may
differ from the exact pinned main baseline.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields, cosmos_potential
from engine import v55_pirk_adapter as adapter
from scipy.linalg import solve_banded
from engine.cmc_gauge import solve_cmc_lapse
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import (
    RHO_DM_PRESENT,
    RHO_B_PRESENT,
    RHO_R_PRESENT,
    initialize_from_archive,
    metric_slice_from_q,
    total_matter_projection,
)
from engine.matter_system import Species, validate_species_admissibility
from tools.audit_initial_centre_series import (
    RESOLUTIONS,
    R_MAX,
    SOURCE_BASELINE,
    SOURCE_BLOBS,
    VENDOR_COMMIT,
    digest_state,
    require_finite_profiles,
    sha256_file,
    stats,
)

ORIGIN_WEIGHTS = np.asarray([401.0 / 720.0, -31.0 / 480.0, 11.0 / 1440.0])
MACHINE_EPS = np.finfo(np.float64).eps


def verify_origin_moments() -> dict:
    """Assert the origin half-cell rule integrates even 1, r², r⁴ exactly."""
    grid_unit_samples = np.asarray([0.5, 1.5, 2.5], dtype=float)
    records = []
    for power in (0, 2, 4):
        actual = float(ORIGIN_WEIGHTS @ (grid_unit_samples ** power))
        expected = float((0.5 ** (power + 1)) / (power + 1))
        error = abs(actual - expected)
        bound = 64.0 * MACHINE_EPS * max(abs(actual), abs(expected), 1.0)
        if error > bound:
            raise RuntimeError(f"origin half-cell moment r^{power} failed: {error}")
        records.append({
            "power": power,
            "unit_dr_moment": actual,
            "exact": expected,
            "absolute_error": error,
            "pass": bool(error <= float(bound)),
        })
    return {"moments": records, "passed": all(x["pass"] for x in records)}


def integrate_origin_even_series(src: np.ndarray, dr: float) -> np.ndarray:
    """PR #18's origin half-cell plus cubic interval-integration counterfactual."""
    values = np.asarray(src, dtype=float)
    h = float(dr)
    if values.ndim != 1 or values.size < 4:
        raise ValueError("origin quadrature requires at least four source samples")
    if not np.isfinite(h) or h <= 0.0 or not np.all(np.isfinite(values)):
        raise ValueError("origin quadrature requires finite samples and positive dr")

    half_cell = h * float(ORIGIN_WEIGHTS @ values[:3])
    increments = np.empty(values.size - 1, dtype=float)
    # Interval [r_0,r_1] with an even source ghost at the origin.
    increments[0] = h * (
        12.0 * values[0] + 13.0 * values[1] - values[2]
    ) / 24.0
    # Interior intervals integrated from the cubic Lagrange interpolant.
    for i in range(2, values.size - 1):
        increments[i - 1] = h * (
            -values[i - 2] + 13.0 * values[i - 1]
            + 13.0 * values[i] - values[i + 1]
        ) / 24.0
    # Final interval uses the one-sided cubic interpolant.
    increments[-1] = h * (
        values[-4] - 5.0 * values[-3]
        + 19.0 * values[-2] + 9.0 * values[-1]
    ) / 24.0
    integral = np.cumsum(np.concatenate(([half_cell], increments)))
    if not np.all(np.isfinite(integral)):
        raise FloatingPointError("origin quadrature generated non-finite integral")
    return integral


def build_initial_state(kernel, grid, *, quadrature: str) -> ProductionState:
    """Same current-main initial equations, varying only radial integration."""
    if quadrature not in ("production_trapezoid", "even_series"):
        raise ValueError(f"unknown quadrature arm: {quadrature}")
    r = np.asarray(grid.centers)
    n = len(r)
    dr = float(grid.dr)
    include_radiation = True

    phi0 = 0.179055
    pi0 = 0.00341913
    dm0 = RHO_DM_PRESENT
    b0 = RHO_B_PRESENT
    rho_r = RHO_R_PRESENT if include_radiation else 0.0

    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.asarray([phi0]))[0])
    rho0 += dm0 + b0 + rho_r
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))

    amplitude, width, D_amplitude = 0.01, 7.0, 1.0e-10
    S = amplitude * np.exp(-(r / width) ** 2)
    D = D_amplitude * np.exp(-(r / width) ** 2)
    PS = np.zeros_like(r)
    PD = D.copy()
    phi = np.full(n, phi0)
    Pi = np.full(n, pi0)
    Kt = np.full(n, -H0)

    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    Sp = d1(S, 1)
    Dp = d1(D, 1)
    Ktp = d1(Kt, 1)
    j = -(PD * Dp)
    Bmass = np.ones_like(r)

    # Deliberately keep the original six source iterations/order of operations.
    for _ in range(6):
        rhoS = 0.5 * (PS**2 + Bmass * Sp**2) + 0.5 * S**2
        rhoD = 0.5 * (PD**2 + Bmass * Dp**2) - 0.5 * D**2
        rhoPhi = (0.5 * Pi**2 + cosmos_potential(phi)) / (8.0 * math.pi)
        rhoMatter = (dm0 + b0 + rho_r) / (8.0 * math.pi)
        rho = rhoS + rhoD + rhoPhi + rhoMatter
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
        src = (
            1.0 - 8.0 * math.pi * r * r * rho
            + 2.0 * r * r * Kr * Kt
            + r * r * Kt * Kt
        )
        if quadrature == "production_trapezoid":
            integ = np.cumsum(np.r_[
                0.5 * dr * src[0],
                0.5 * dr * (src[1:] + src[:-1]),
            ])
        else:
            integ = integrate_origin_even_series(src, dr)
        Bmass = integ / r

    if np.any(~np.isfinite(Bmass)) or np.any(Bmass <= 0.0):
        raise FloatingPointError(f"invalid B(r) in {quadrature} arm")

    A = 1.0 / np.sqrt(Bmass)
    Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
    K = Kr + 2.0 * Kt
    Aa = Kr - K / 3.0
    a = A ** (4.0 / 3.0)
    b = A ** (-2.0 / 3.0)
    X = A ** (-1.0 / 3.0)
    Lambda = (
        d1(a, 1) / (2.0 * a * a)
        - d1(b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )

    _, vacuum, _ = adapter.vendor_modules()
    geometry = vacuum.flat_state(grid)
    geometry.a[:] = a
    geometry.b[:] = b
    geometry.X[:] = X
    geometry.alpha[:] = 1.0
    geometry.beta[:] = 0.0
    geometry.Aa[:] = Aa
    geometry.K[:] = K
    geometry.Lambda[:] = Lambda
    geometry.B[:] = 0.75 * Lambda
    scalars = ScalarFields(S, PS, D, PD, phi, Pi)
    matter = initialize_from_archive(grid, geometry, include_radiation=include_radiation)
    validate_species_admissibility(
        metric_slice_from_q(grid, geometry), matter.radiation,
        Species.RADIATION, stage=f"observer-{quadrature}-initialization",
    )

    # The normal production initialization then solves CMC and zeros beta/B.
    geometry.alpha = kernel._solve_lapse(grid, geometry, scalars, matter)[0]
    geometry.beta.fill(0.0)
    geometry.B.fill(0.0)
    return ProductionState(grid=grid, geometry=geometry, scalars=scalars, matter=matter)


def solve_even_reflection_centre_pde_lapse(kernel, state) -> tuple[np.ndarray, dict]:
    """Solve the same frozen CMC target/operator with a second-order even-ghost centre row."""
    grid = state.grid
    geometry_input = state.geometry.copy()
    # Current production initialize() supplies alpha=1 to the CMC solve.
    geometry_input.alpha[:] = 1.0
    alpha_current, kdot = solve_cmc_lapse(
        grid, geometry_input, state.scalars, state.matter
    )
    current_match = float(np.max(np.abs(
        np.asarray(alpha_current) - np.asarray(state.geometry.alpha)
    )))
    if current_match > 2.0e-13:
        raise RuntimeError(
            f"production CMC reconstruction differs from initialized state: {current_match:.3e}"
        )

    total = total_matter_projection(
        metric_slice_from_q(grid, geometry_input), state.matter
    )
    a = np.asarray(geometry_input.a, dtype=float)
    b = np.asarray(geometry_input.b, dtype=float)
    X = np.asarray(geometry_input.X, dtype=float)
    r = np.asarray(grid.centers, dtype=float)
    inv = X * X / a
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    ap, bp, Xp = d1(a, 1), d1(b, 1), d1(X, 1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r
    Q = (
        1.5 * geometry_input.Aa * geometry_input.Aa
        + geometry_input.K * geometry_input.K / 3.0
        + 4.0 * math.pi * (
            total["rho"] + total["pr"] + 2.0 * total["pt"]
        )
    )
    m = Q / inv
    rhs = -float(kdot) / inv
    n, h = grid.n, float(grid.dr)
    band = np.zeros((3, n), dtype=float)
    bvec = np.zeros(n, dtype=float)

    # Same second-order CMC interior and outer normalization as production;
    # only the first row changes from alpha[0]=alpha[1] to the even-ghost PDE.
    ci0 = c[0] / (2.0 * h)
    band[1, 0] = -1.0 / h**2 - ci0 - m[0]
    band[0, 1] = 1.0 / h**2 + ci0
    bvec[0] = rhs[0]
    for i in range(1, n - 1):
        ci = c[i] / (2.0 * h)
        band[2, i - 1] = 1.0 / h**2 - ci
        band[1, i] = -2.0 / h**2 - m[i]
        band[0, i + 1] = 1.0 / h**2 + ci
        bvec[i] = rhs[i]
    band[1, n - 1] = 1.0
    bvec[n - 1] = 1.0

    alpha_even = solve_banded((1, 1), band, bvec, check_finite=False)
    if not np.all(np.isfinite(alpha_even)) or float(np.min(alpha_even)) <= 0.0:
        raise RuntimeError("even-reflection CMC row produced invalid lapse")
    row_residual = float(np.max(np.abs(
        band[1] * alpha_even
        + np.r_[band[0, 1:] * alpha_even[1:], 0.0]
        + np.r_[0.0, band[2, :-1] * alpha_even[:-1]]
        - bvec
    )))
    return np.asarray(alpha_even, dtype=float), {
        "incoming_target_kdot": float(kdot),
        "production_lapse_reconstruction_max_abs": current_match,
        "current_center_row_alpha0_minus_alpha1": float(alpha_current[0] - alpha_current[1]),
        "even_pde_center_row_alpha0_minus_alpha1": float(alpha_even[0] - alpha_even[1]),
        "production_alpha_min": float(np.min(alpha_current)),
        "even_pde_alpha_min": float(np.min(alpha_even)),
        "production_alpha_max": float(np.max(alpha_current)),
        "even_pde_alpha_max": float(np.max(alpha_even)),
        "max_abs_lapse_difference": float(np.max(np.abs(alpha_even - alpha_current))),
        "even_pde_row_matrix_residual_max_abs": row_residual,
        "outer_alpha_production": float(alpha_current[-1]),
        "outer_alpha_even_pde": float(alpha_even[-1]),
        "fixed_non_lapse_state": True,
        "target_held_fixed_from_alpha_one_input": True,
    }


def measure_state(state) -> tuple[dict, dict[str, np.ndarray]]:
    grid, g = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=float)
    _, vacuum, _ = adapter.vendor_modules()
    geom = vacuum.geometry_terms(grid, g)
    l2 = vacuum.primary_l2_rhs(grid, g)
    l3_vac = vacuum.primary_l3_rhs(grid, g)
    l3_total = adapter.primary_l3_with_matter(grid, g, state.scalars, state.matter)

    hessian = -(np.asarray(geom["DrDr_alpha"]) - np.asarray(geom["lap_alpha"]) / 3.0)
    ricci = np.asarray(g.alpha) * (np.asarray(geom["Rr_r"]) - np.asarray(geom["R"]) / 3.0)
    l2_aa = np.asarray(l2["Aa"], dtype=float)
    l3_aa = np.asarray(l3_total["Aa"], dtype=float)
    full_aa = l2_aa + l3_aa
    stress = assemble_total_stress_energy(grid, g, state.scalars, state.matter)
    raw = vacuum.constraints(grid, g)
    H = np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(stress.rho)
    M = np.asarray(raw["momentum"]) - 8.0 * math.pi * np.asarray(stress.j)
    C = np.asarray(raw["connection"])
    det = np.asarray(raw["determinant"])
    r2 = r * r
    profiles = {
        "r": r,
        "geometry__a": np.asarray(g.a),
        "geometry__b": np.asarray(g.b),
        "geometry__X": np.asarray(g.X),
        "geometry__alpha": np.asarray(g.alpha),
        "geometry__Aa": np.asarray(g.Aa),
        "geometry__K": np.asarray(g.K),
        "geometry__Lambda": np.asarray(g.Lambda),
        "R": np.asarray(geom["R"]),
        "Rr_r": np.asarray(geom["Rr_r"]),
        "constraint_H": H,
        "constraint_M": M,
        "constraint_C_Lambda": C,
        "constraint_determinant": det,
        "lapse_hessian_Aa": hessian,
        "ricci_tracefree_Aa": ricci,
        "primary_L2_Aa": l2_aa,
        "primary_L3_vacuum_Aa": np.asarray(l3_vac["Aa"]),
        "primary_L3_with_matter_Aa": l3_aa,
        "full_Aa_RHS": full_aa,
        "metric_ratio_over_r2": (g.a / g.b - 1.0) / r2,
        "Lambda_over_r": g.Lambda / r,
        "tracefree_Ricci_over_r2": ricci / r2,
        "full_Aa_RHS_over_r2": full_aa / r2,
    }
    require_finite_profiles(profiles)
    closure = l2_aa - (hessian + ricci)
    source_closure = full_aa - (hessian + ricci + l3_aa)
    scale = np.abs(l2_aa) + np.abs(hessian) + np.abs(ricci)
    bound = 64.0 * MACHINE_EPS * np.maximum(scale, np.finfo(float).tiny)
    if np.any(np.abs(closure) > bound):
        raise RuntimeError("L2 source decomposition failed operand-scaled roundoff")
    if float(np.max(np.abs(source_closure))) > 1.0e-14:
        raise RuntimeError("full Aa source decomposition failed roundoff check")

    stats_out = {
        "state_digest": digest_state(state),
        "time": float(state.t),
        "history_rows": len(state.history),
        "initial_constraints": {
            "H_cell_0": float(H[0]),
            "H_full_absmax": float(np.max(np.abs(H))),
            "H_r_ge_5_absmax": float(np.max(np.abs(H[r >= 5.0]))),
            "M_cell_0": float(M[0]),
            "M_full_absmax": float(np.max(np.abs(M))),
            "M_r_ge_5_absmax": float(np.max(np.abs(M[r >= 5.0]))),
            "C_cell_0": float(C[0]),
            "C_full_absmax": float(np.max(np.abs(C))),
            "C_r_ge_5_absmax": float(np.max(np.abs(C[r >= 5.0]))),
            "determinant_full_absmax": float(np.max(np.abs(det))),
        },
        "first_six_cells": [{
            "i": int(i),
            "r": float(r[i]),
            "metric_ratio_over_r2": float(profiles["metric_ratio_over_r2"][i]),
            "Lambda_over_r": float(profiles["Lambda_over_r"][i]),
            "tracefree_Ricci_Aa": float(ricci[i]),
            "tracefree_Ricci_over_r2": float((ricci/r2)[i]),
            "lapse_hessian_Aa": float(hessian[i]),
            "primary_L2_Aa": float(l2_aa[i]),
            "primary_L3_with_matter_Aa": float(l3_aa[i]),
            "full_Aa_RHS": float(full_aa[i]),
            "full_Aa_RHS_over_r2": float((full_aa/r2)[i]),
        } for i in range(min(6, grid.n))],
        "source_closure": {
            "l2_max_abs": float(np.max(np.abs(closure))),
            "l2_max_operand_scaled": float(np.max(np.abs(closure) / np.maximum(bound, np.finfo(float).tiny))),
            "full_max_abs": float(np.max(np.abs(source_closure))),
        },
    }
    return stats_out, profiles


def compare_fields(expected, actual) -> dict:
    differences = {}
    names = ("a", "b", "X", "alpha", "Aa", "K", "Lambda", "beta", "B")
    for name in names:
        lhs = np.asarray(getattr(expected.geometry, name), dtype=float)
        rhs = np.asarray(getattr(actual.geometry, name), dtype=float)
        differences[f"geometry.{name}"] = float(np.max(np.abs(lhs-rhs)))
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        lhs = np.asarray(getattr(expected.scalars, name), dtype=float)
        rhs = np.asarray(getattr(actual.scalars, name), dtype=float)
        differences[f"scalar.{name}"] = float(np.max(np.abs(lhs-rhs)))
    for sp in ("dark_matter", "baryons", "radiation"):
        for field in ("rest", "energy_t", "momentum_r"):
            lhs = np.asarray(getattr(getattr(expected.matter, sp), field), dtype=float)
            rhs = np.asarray(getattr(getattr(actual.matter, sp), field), dtype=float)
            differences[f"matter.{sp}.{field}"] = float(np.max(np.abs(lhs-rhs)))
    return differences


def main() -> None:
    output = Path("artifacts/initial-origin-quadrature-ab")
    output.mkdir(parents=True, exist_ok=True)
    moment_gate = verify_origin_moments()
    kernel = V55ProductionKernel()
    records = []
    profile_archive = {}
    for n in RESOLUTIONS:
        production = kernel.initialize(
            resolution=n, r_max=R_MAX, amplitude=0.01, width=7.0,
            D_amplitude=1.0e-10, include_radiation=True,
        )
        if production.t != 0.0 or production.history:
            raise RuntimeError("production initializer unexpectedly advanced time or wrote history")
        before = digest_state(production)
        grid = production.grid
        reconstructed_trap = build_initial_state(kernel, grid, quadrature="production_trapezoid")
        trap_diffs = compare_fields(production, reconstructed_trap)
        trap_max = max(trap_diffs.values())
        trap_scale = max(
            1.0,
            max(float(np.max(np.abs(getattr(production.geometry, name)))) for name in
                ("a", "b", "X", "alpha", "Aa", "K", "Lambda", "beta", "B"))
        )
        if trap_max > 2.0e-13 * trap_scale:
            raise RuntimeError(f"independent trapezoid reconstruction differs from production initializer at N={n}: {trap_max}")

        candidate = build_initial_state(kernel, grid, quadrature="even_series")
        if candidate.t != 0.0 or candidate.history:
            raise RuntimeError("counterfactual initialization advanced time or wrote history")
        baseline_report, baseline_profiles = measure_state(production)
        candidate_report, candidate_profiles = measure_state(candidate)
        if digest_state(production) != before or baseline_report["state_digest"] != before:
            raise RuntimeError("baseline observer mutated the production initial state")

        # Apply PR #18's centre-row-only alternative to this same quadrature arm.
        # All non-lapse fields, conservative matter, and the incoming target stay fixed.
        alpha_even, row_meta = solve_even_reflection_centre_pde_lapse(kernel, candidate)
        centre_row_state = ProductionState(
            grid=grid,
            geometry=candidate.geometry.copy(),
            scalars=candidate.scalars,
            matter=candidate.matter,
        )
        centre_row_state.geometry.alpha = alpha_even.copy()
        centre_row_state.geometry.beta.fill(0.0)
        centre_row_state.geometry.B.fill(0.0)
        centre_row_report, centre_row_profiles = measure_state(centre_row_state)
        if centre_row_state.t != 0.0 or centre_row_state.history:
            raise RuntimeError("centre-row counterfactual unexpectedly advanced time")
        for field in ("a", "b", "X", "Aa", "K", "Lambda", "beta", "B"):
            if not np.array_equal(
                np.asarray(getattr(candidate.geometry, field)),
                np.asarray(getattr(centre_row_state.geometry, field)),
            ):
                raise RuntimeError(f"centre-row counterfactual changed non-lapse geometry: {field}")

        # Save every radial profile from all arms for independent signed review.
        for name, values in baseline_profiles.items():
            profile_archive[f"N{n}__production__{name}"] = values
        for name, values in candidate_profiles.items():
            profile_archive[f"N{n}__even_series__{name}"] = values
        for name, values in centre_row_profiles.items():
            profile_archive[f"N{n}__even_series_centre_pde__{name}"] = values

        records.append({
            "resolution": n,
            "r_max": R_MAX,
            "production_state_digest": before,
            "baseline_reconstruction": {
                "max_abs_field_difference": trap_max,
                "field_differences": trap_diffs,
                "pass": bool(trap_max <= 2.0e-13 * trap_scale),
            },
            "production_trapezoid": baseline_report,
            "origin_even_series_counterfactual": candidate_report,
            "origin_even_series_centre_pde_lapse_counterfactual": centre_row_report,
            "centre_row_counterfactual": row_meta,
            "delta_counterfactual_minus_production": {
                "H_cell_0": candidate_report["initial_constraints"]["H_cell_0"] - baseline_report["initial_constraints"]["H_cell_0"],
                "M_cell_0": candidate_report["initial_constraints"]["M_cell_0"] - baseline_report["initial_constraints"]["M_cell_0"],
                "C_cell_0": candidate_report["initial_constraints"]["C_cell_0"] - baseline_report["initial_constraints"]["C_cell_0"],
                "metric_ratio_over_r2_cell_0": candidate_report["first_six_cells"][0]["metric_ratio_over_r2"] - baseline_report["first_six_cells"][0]["metric_ratio_over_r2"],
                "Lambda_over_r_cell_0": candidate_report["first_six_cells"][0]["Lambda_over_r"] - baseline_report["first_six_cells"][0]["Lambda_over_r"],
                "full_Aa_RHS_cell_0": candidate_report["first_six_cells"][0]["full_Aa_RHS"] - baseline_report["first_six_cells"][0]["full_Aa_RHS"],
                "full_Aa_RHS_over_r2_cell_0": candidate_report["first_six_cells"][0]["full_Aa_RHS_over_r2"] - baseline_report["first_six_cells"][0]["full_Aa_RHS_over_r2"],
                "H_full_absmax": candidate_report["initial_constraints"]["H_full_absmax"] - baseline_report["initial_constraints"]["H_full_absmax"],
            },
            "delta_centre_pde_lapse_minus_even_series_lapse": {
                "alpha0": float(alpha_even[0] - candidate.geometry.alpha[0]),
                "Aa_rhs_cell0": centre_row_report["first_six_cells"][0]["full_Aa_RHS"] - candidate_report["first_six_cells"][0]["full_Aa_RHS"],
                "Aa_rhs_over_r2_cell0": centre_row_report["first_six_cells"][0]["full_Aa_RHS_over_r2"] - candidate_report["first_six_cells"][0]["full_Aa_RHS_over_r2"],
                "lapse_hessian_cell0": centre_row_report["first_six_cells"][0]["lapse_hessian_Aa"] - candidate_report["first_six_cells"][0]["lapse_hessian_Aa"],
                "tracefree_Ricci_cell0": centre_row_report["first_six_cells"][0]["tracefree_Ricci_Aa"] - candidate_report["first_six_cells"][0]["tracefree_Ricci_Aa"],
                "H_cell0": centre_row_report["initial_constraints"]["H_cell_0"] - candidate_report["initial_constraints"]["H_cell_0"],
                "M_cell0": centre_row_report["initial_constraints"]["M_cell_0"] - candidate_report["initial_constraints"]["M_cell_0"],
                "C_cell0": centre_row_report["initial_constraints"]["C_cell_0"] - candidate_report["initial_constraints"]["C_cell_0"],
            },
        })

    npz = output / "initial_origin_quadrature_ab_profiles.npz"
    np.savez_compressed(npz, **profile_archive)
    report = {
        "schema": "v55_initial_origin_quadrature_ab_v1",
        "status": "completed_initialization_only",
        "diagnostic_only": True,
        "scientific_admission": False,
        "evolution_steps": 0,
        "source_baseline": SOURCE_BASELINE,
        "source_blobs": SOURCE_BLOBS,
        "vendor_commit": VENDOR_COMMIT,
        "moment_gate": moment_gate,
        "resolutions": records,
        "interpretation_limits": [
            "No time step or evolution was performed. This is only a current-source initialized-state comparison.",
            "The production trapezoid arm must reproduce the current production initializer to roundoff before the candidate is interpreted.",
            "The even-series arm changes only radial integration in the initializer and retains the current source equations, parameters, matter initialization and CMC solve.",
            "No production physics, initializer, gauge, source, vendor, boundary or campaign settings are changed.",
            "A better first-cell centre-series/H constraint is not enough for promotion: fixed-coordinate/fixed-areal-radius Misner-Sharp current/flux/work balance remains a separate gate.",
        ],
    }
    report_path = output / "report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    manifest = output / "SHA256SUMS.txt"
    manifest.write_text(
        f"{sha256_file(report_path)}  {report_path.name}\n"
        f"{sha256_file(npz)}  {npz.name}\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "status": report["status"],
        "scientific_admission": report["scientific_admission"],
        "evolution_steps": report["evolution_steps"],
        "moment_gate": moment_gate,
        "resolutions": [{
            "N": x["resolution"],
            "trap_reconstruction_max": x["baseline_reconstruction"]["max_abs_field_difference"],
            "baseline_H0": x["production_trapezoid"]["initial_constraints"]["H_cell_0"],
            "candidate_H0": x["origin_even_series_counterfactual"]["initial_constraints"]["H_cell_0"],
            "baseline_metric_ratio0": x["production_trapezoid"]["first_six_cells"][0]["metric_ratio_over_r2"],
            "candidate_metric_ratio0": x["origin_even_series_counterfactual"]["first_six_cells"][0]["metric_ratio_over_r2"],
            "baseline_AaDotOverR2_0": x["production_trapezoid"]["first_six_cells"][0]["full_Aa_RHS_over_r2"],
            "candidate_AaDotOverR2_0": x["origin_even_series_counterfactual"]["first_six_cells"][0]["full_Aa_RHS_over_r2"],
            "baseline_C0": x["production_trapezoid"]["initial_constraints"]["C_cell_0"],
            "candidate_C0": x["origin_even_series_counterfactual"]["initial_constraints"]["C_cell_0"],
            "candidate_AaDotOverR20_even_pde_centre_row": x["origin_even_series_centre_pde_lapse_counterfactual"]["first_six_cells"][0]["full_Aa_RHS_over_r2"],
            "candidate_Hessian_Aa0_even_pde_centre_row": x["origin_even_series_centre_pde_lapse_counterfactual"]["first_six_cells"][0]["lapse_hessian_Aa"],
            "candidate_alpha0_minus_alpha1_even_pde_centre_row": x["centre_row_counterfactual"]["even_pde_center_row_alpha0_minus_alpha1"],
            "centre_row_max_alpha_difference": x["centre_row_counterfactual"]["max_abs_lapse_difference"],
        } for x in records],
        "hashes": {
            report_path.name: sha256_file(report_path),
            npz.name: sha256_file(npz),
        },
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
