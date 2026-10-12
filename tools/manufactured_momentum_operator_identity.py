#!/usr/bin/env python3
"""Diagnostic-only analytic consistency audit for the production momentum constraint.

A smooth regular profile is sampled on the exact pinned cell grid. Each signed
momentum-constraint term is evaluated both with the production discrete
operator and with independent analytic derivatives. The manufactured state is
NOT claimed to solve the Einstein-matter constraints or the evolution system;
the purpose is to measure operator consistency and centre-versus-annulus error.
No evolution step and no production-source edit is performed.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path("source").resolve()))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_pirk_adapter import vendor_modules

SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 40.0
PROFILE_SCALE = 6.0
ANNULUS = (3.0, 6.0)
CENTER_CELLS = (0, 1, 2, 3)
SAMPLE_RADII = (3.0, 4.0, 5.0)

# Explicitly declared manufactured amplitudes. These are test profiles, not
# fitted evolution coefficients or proposed production parameters.
PROFILE = {
    "b_amplitude": 0.020,
    "X_amplitude": 0.015,
    "alpha_amplitude": 0.020,
    "Aa_amplitude": 0.020,
    "K_offset": 0.010,
    "K_amplitude": 0.015,
    "S_offset": 0.200,
    "S_amplitude": 0.050,
    "PS_offset": 0.015,
    "PS_amplitude": 0.008,
    "D_offset": 0.040,
    "D_amplitude": 0.015,
    "PD_offset": 0.003,
    "PD_amplitude": 0.002,
    "phi_offset": 0.100,
    "phi_amplitude": 0.030,
    "Pi_offset": 0.004,
    "Pi_amplitude": 0.001,
}


def shape_and_derivative(r):
    """g(r)=s/(1+s), s=(r/L)^2, and its exact radial derivative."""
    s = (np.asarray(r, dtype=float) / PROFILE_SCALE) ** 2
    g = s / (1.0 + s)
    gp = (2.0 * np.asarray(r, dtype=float) / PROFILE_SCALE**2) / (1.0 + s) ** 2
    return g, gp


def manufactured_profiles(r):
    r = np.asarray(r, dtype=float)
    g, gp = shape_and_derivative(r)
    p = PROFILE

    b = 1.0 + p["b_amplitude"] * g
    bp = p["b_amplitude"] * gp
    a = b ** -2
    ap = -2.0 * b ** -3 * bp
    X = 1.0 + p["X_amplitude"] * g
    Xp = p["X_amplitude"] * gp
    alpha = 1.0 + p["alpha_amplitude"] * g

    Aa = p["Aa_amplitude"] * g
    Aap = p["Aa_amplitude"] * gp
    K = p["K_offset"] + p["K_amplitude"] * g
    Kp = p["K_amplitude"] * gp

    S = p["S_offset"] + p["S_amplitude"] * g
    Sp = p["S_amplitude"] * gp
    PS = p["PS_offset"] + p["PS_amplitude"] * g
    D = p["D_offset"] + p["D_amplitude"] * g
    Dp = p["D_amplitude"] * gp
    PD = p["PD_offset"] + p["PD_amplitude"] * g
    phi = p["phi_offset"] + p["phi_amplitude"] * g
    phip = p["phi_amplitude"] * gp
    Pi = p["Pi_offset"] + p["Pi_amplitude"] * g

    # The continuum conformal connection for this analytic metric.
    # Since the grid is cell-centred, r is strictly positive here.
    Lambda = ap / (2.0 * a**2) - bp / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a)

    values = {
        "a": a, "ap": ap, "b": b, "bp": bp, "X": X, "Xp": Xp,
        "alpha": alpha, "Aa": Aa, "Aap": Aap, "K": K, "Kp": Kp,
        "S": S, "Sp": Sp, "PS": PS,
        "D": D, "Dp": Dp, "PD": PD,
        "phi": phi, "phip": phip, "Pi": Pi, "Lambda": Lambda,
    }
    for name, value in values.items():
        if not np.all(np.isfinite(value)):
            raise FloatingPointError(f"non-finite manufactured profile {name}")
    return values


def rms(values):
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def maxabs(values):
    return float(np.max(np.abs(np.asarray(values, dtype=float))))


def observed_order(errors, spacings):
    errors = np.asarray(errors, dtype=float)
    spacings = np.asarray(spacings, dtype=float)
    if np.any(errors <= 10.0 * np.finfo(float).eps):
        return None
    slope = float(np.polyfit(np.log(spacings), np.log(errors), 1)[0])
    return slope


def audit_resolution(kernel, n):
    state = kernel.initialize(
        resolution=n, r_max=R_MAX, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    grid = state.grid
    r = np.asarray(grid.centers, dtype=float)
    exact = manufactured_profiles(r)
    geom = state.geometry.copy()

    for name in ("a", "b", "X", "alpha", "Aa", "K", "Lambda"):
        setattr(geom, name, np.asarray(exact[name], dtype=float).copy())
    geom.beta = np.zeros_like(r)
    geom.B = np.zeros_like(r)
    scalars = ScalarFields(
        S=exact["S"].copy(), PS=exact["PS"].copy(),
        D=exact["D"].copy(), PD=exact["PD"].copy(),
        phi=exact["phi"].copy(), Pi=exact["Pi"].copy(),
    )

    _, vacuum, _ = vendor_modules()
    raw_constraints = vacuum.constraints(grid, geom)
    total = assemble_total_stress_energy(grid, geom, scalars, state.matter)

    # Discrete signed terms, independently listed rather than using the
    # bundled momentum constraint as the only reported value.
    d1 = lambda x, parity=1: np.asarray(
        grid.cell_derivative_fourth(x, parity=parity), dtype=float
    )
    Ab = -0.5 * np.asarray(geom.Aa)
    chi_prime_discrete = -0.5 * d1(geom.X) / np.asarray(geom.X)
    discrete_terms = {
        "Aa_prime": d1(geom.Aa),
        "minus_two_thirds_K_prime": -(2.0 / 3.0) * d1(geom.K),
        "six_Aa_chi_prime": 6.0 * np.asarray(geom.Aa) * chi_prime_discrete,
        "radial_tensor": (np.asarray(geom.Aa) - Ab) * (
            2.0 / r + d1(geom.b) / np.asarray(geom.b)
        ),
        "minus_eight_pi_current": -8.0 * math.pi * np.asarray(total.j),
    }

    # Analytic counterparts at the exact same cell centres.
    chi_prime_exact = -0.5 * exact["Xp"] / exact["X"]
    j_exact = (
        -exact["PS"] * exact["Sp"]
        -exact["PD"] * exact["Dp"]
        -exact["Pi"] * exact["phip"] / (8.0 * math.pi)
    )
    continuum_terms = {
        "Aa_prime": exact["Aap"],
        "minus_two_thirds_K_prime": -(2.0 / 3.0) * exact["Kp"],
        "six_Aa_chi_prime": 6.0 * exact["Aa"] * chi_prime_exact,
        "radial_tensor": 1.5 * exact["Aa"] * (
            2.0 / r + exact["bp"] / exact["b"]
        ),
        "minus_eight_pi_current": -8.0 * math.pi * j_exact,
    }

    M_discrete_rebuilt = sum(discrete_terms.values(), np.zeros_like(r))
    M_production = np.asarray(raw_constraints["momentum"]) - 8.0 * math.pi * np.asarray(total.j)
    M_continuum = sum(continuum_terms.values(), np.zeros_like(r))
    connection_residual = np.asarray(raw_constraints["connection"])
    determinant_residual = np.asarray(raw_constraints["determinant"])
    fluid_current = np.asarray(total.fluid_j)

    if maxabs(M_discrete_rebuilt - M_production) > 2.0e-12:
        raise AssertionError("signed discrete momentum-term reconstruction did not close")
    if maxabs(fluid_current) > 1.0e-12:
        raise AssertionError("static fluid control unexpectedly carries radial current")
    if maxabs(determinant_residual) > 2.0e-12:
        raise AssertionError("manufactured metric does not satisfy a*b^2=1")
    if min(float(np.min(getattr(geom, q))) for q in ("a", "b", "X", "alpha")) <= 0.0:
        raise AssertionError("manufactured metric/lapse is not positive")

    annulus = (r >= ANNULUS[0]) & (r <= ANNULUS[1])
    if not np.any(annulus):
        raise AssertionError(f"empty fixed physical annulus for N={n}")

    term_summary = {}
    for name in discrete_terms:
        err = discrete_terms[name] - continuum_terms[name]
        term_summary[name] = {
            "annulus_rms_discretization_error": rms(err[annulus]),
            "annulus_max_abs_discretization_error": maxabs(err[annulus]),
            "first_center_cells": [
                {
                    "i": i, "r": float(r[i]),
                    "discrete": float(discrete_terms[name][i]),
                    "analytic": float(continuum_terms[name][i]),
                    "signed_error": float(err[i]),
                }
                for i in CENTER_CELLS if i < n
            ],
        }

    m_error = M_production - M_continuum
    c_error = connection_residual  # continuum connection constraint is zero
    sample_rows = {}
    for target in SAMPLE_RADII:
        i = int(np.argmin(np.abs(r - target)))
        sample_rows[str(target)] = {
            "target_radius": target, "i": i, "actual_radius": float(r[i]),
            "M_discrete": float(M_production[i]),
            "M_continuum": float(M_continuum[i]),
            "M_signed_error": float(m_error[i]),
            "C_lambda_discrete": float(connection_residual[i]),
        }

    report = {
        "resolution": n, "r_max": R_MAX, "dr": float(grid.dr),
        "source_commit": SOURCE_COMMIT, "profile_scale": PROFILE_SCALE,
        "profile_coefficients": PROFILE,
        "interpretation_limit": (
            "Manufactured smooth regular profiles measure spatial-operator consistency. "
            "They are not asserted to solve the continuum constraints or full evolution equations."
        ),
        "discrete_signed_momentum_terms": term_summary,
        "centre_cells": [
            {
                "i": i, "r": float(r[i]),
                "M_discrete": float(M_production[i]),
                "M_continuum": float(M_continuum[i]),
                "M_signed_error": float(m_error[i]),
                "C_lambda_discrete": float(connection_residual[i]),
                "determinant_residual": float(determinant_residual[i]),
            }
            for i in CENTER_CELLS if i < n
        ],
        "fixed_radius_samples": sample_rows,
        "fixed_annulus": {
            "r_min": ANNULUS[0], "r_max": ANNULUS[1],
            "cell_count": int(np.count_nonzero(annulus)),
            "M_discrete_rms": rms(M_production[annulus]),
            "M_continuum_rms": rms(M_continuum[annulus]),
            "M_rms_discretization_error": rms(m_error[annulus]),
            "M_max_abs_discretization_error": maxabs(m_error[annulus]),
            "C_lambda_rms_discretization_error": rms(c_error[annulus]),
            "C_lambda_max_abs_discretization_error": maxabs(c_error[annulus]),
            "continuum_M_nontrivial": bool(maxabs(M_continuum[annulus]) > 1.0e-8),
        },
        "closure_checks": {
            "max_abs_discrete_rebuild_minus_production_M": maxabs(M_discrete_rebuilt - M_production),
            "max_abs_fluid_current": maxabs(fluid_current),
            "max_abs_determinant_residual": maxabs(determinant_residual),
            "all_finite": bool(
                np.all(np.isfinite(M_production))
                and np.all(np.isfinite(M_continuum))
                and np.all(np.isfinite(connection_residual))
            ),
        },
    }
    if not report["closure_checks"]["all_finite"]:
        raise FloatingPointError(f"non-finite constraint values at N={n}")
    if not report["fixed_annulus"]["continuum_M_nontrivial"]:
        raise AssertionError(f"manufactured constraint profile is too trivial at N={n}")
    return report


def main():
    kernel = V55ProductionKernel()
    reports = [audit_resolution(kernel, n) for n in RESOLUTIONS]
    spacings = [x["dr"] for x in reports]
    annulus_errors = {
        "M": [x["fixed_annulus"]["M_rms_discretization_error"] for x in reports],
        "C_lambda": [x["fixed_annulus"]["C_lambda_rms_discretization_error"] for x in reports],
    }
    orders = {
        name: observed_order(values, spacings)
        for name, values in annulus_errors.items()
    }
    result = {
        "schema": "v55_manufactured_momentum_operator_identity_v1",
        "status": "completed",
        "diagnostic_only": True,
        "source_commit": SOURCE_COMMIT,
        "resolutions": list(RESOLUTIONS),
        "annulus_observed_order_descriptive_only": orders,
        "resolution_reports": reports,
        "guardrails": [
            "Pinned production source is imported read-only at the recorded commit.",
            "No production RHS, parameter, boundary condition, or initial-data routine is modified.",
            "No evolution step is taken; only the analytic profile and the production operators are evaluated.",
            "A successful workflow means the diagnostic ran and its algebraic closures passed, not scientific admission of the production machine.",
        ],
    }
    out = Path("results/manufactured-momentum-operator-identity")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    summary = {
        "status": result["status"],
        "source_commit": SOURCE_COMMIT,
        "observed_order_descriptive_only": orders,
        "rows": [
            {
                "N": x["resolution"], "dr": x["dr"],
                "annulus_cells": x["fixed_annulus"]["cell_count"],
                "M_annulus_error_rms": x["fixed_annulus"]["M_rms_discretization_error"],
                "C_lambda_annulus_error_rms": x["fixed_annulus"]["C_lambda_rms_discretization_error"],
                "center_M_error_i0": x["centre_cells"][0]["M_signed_error"],
                "center_M_error_i1": x["centre_cells"][1]["M_signed_error"],
            }
            for x in reports
        ],
    }
    print("MANUFACTURED_MOMENTUM_OPERATOR_IDENTITY " + json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
