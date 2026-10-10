"""Curved-space scalar energy-balance diagnostic; measurement-only.

Evaluates the production scalar RHS and the production geometry RHS on one
post-regularity-projection initial slice.  For each canonical scalar, it
compares the chain-rule time derivative of the code's Eulerian energy density
against the curved radial flux divergence plus lapse, extrinsic-curvature,
inverse-metric, and explicit scalar-source terms.

This is a diagnostic identity check, not a production change or physics
admission.  In particular, the initial regularity-projection jump is recorded
separately and is not included in the instantaneous balance residual.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.hamiltonian_evolution_increment_decomposition import hamiltonian_residual
from engine.production_kernel import V55ProductionKernel, BETA_DM
from engine.scalar_system import (
    cosmos_potential,
    cosmos_potential_prime,
    scalar_rhs_arrays,
)
from engine.v55_matter import dm_density, metric_slice_from_q

KAPPA = 8.0 * math.pi
SETTINGS = {
    "resolution": 160,
    "r_max": 40.0,
    "amplitude": 0.01,
    "width": 7.0,
    "D_amplitude": 1.0e-10,
    "include_radiation": True,
    "max_iter": 12,
    "tol": 1.0e-13,
    "max_cond": 1.0e12,
}
FIELD_DEFS = (
    ("S", "S", "PS", 1.0),
    ("D", "D", "PD", 1.0),
    ("phi", "phi", "Pi", KAPPA),
)


def first(values):
    return [float(v) for v in np.asarray(values)[:5]]


def _d1(grid, values, parity):
    return grid.cell_derivative_fourth(np.asarray(values), parity=parity)


def _hamiltonian(state):
    return hamiltonian_residual(state)


def _energy_terms(grid, geometry, fields, srhs, geometry_rhs, source_phi):
    """Return cellwise chain-rule and curved-balance terms for S, D, and phi.

    The shift is zero on the production CMC branch.  With gamma^rr=X^2/a,
    sqrt(gamma)=sqrt(a)*b*r^2/X^3, and j_r=-Pi*partial_r(phi)/scale,
    the radial flux is sqrt(gamma)*alpha*gamma^rr*Pi*partial_r(phi).
    """
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    alpha = np.asarray(geometry.alpha, dtype=float)
    invr = X * X / a
    sqrt_gamma = np.sqrt(a) * b * r * r / (X ** 3)
    alphap = _d1(grid, alpha, 1)
    K = np.asarray(geometry.K, dtype=float)
    Aa = np.asarray(geometry.Aa, dtype=float)
    kr_r = K / 3.0 + Aa

    adot = np.asarray(geometry_rhs["explicit"]["a"], dtype=float)
    Xdot = np.asarray(geometry_rhs["explicit"]["X"], dtype=float)
    invrdot = 2.0 * X * Xdot / a - X * X * adot / (a * a)
    invrdot_adm = 2.0 * alpha * invr * kr_r
    invrdot_error = invrdot - invrdot_adm

    result = {}
    for name, f_name, p_name, scale in FIELD_DEFS:
        f = np.asarray(getattr(fields, f_name), dtype=float)
        p = np.asarray(getattr(fields, p_name), dtype=float)
        ft = np.asarray(getattr(srhs, f_name), dtype=float)
        pt = np.asarray(getattr(srhs, p_name), dtype=float)
        fp = _d1(grid, f, 1)
        ftp = _d1(grid, ft, 1)

        if name == "S":
            V = 0.5 * f * f
            Vp = f
            source = np.zeros_like(f)
        elif name == "D":
            V = -0.5 * f * f
            Vp = -f
            source = np.zeros_like(f)
        else:
            V = cosmos_potential(f)
            Vp = cosmos_potential_prime(f)
            source = np.asarray(source_phi, dtype=float)

        rho = (0.5 * p * p + 0.5 * invr * fp * fp + V) / scale
        # Direct chain rule using the production scalar RHS and geometry RHS.
        chain = (
            p * pt
            + invr * fp * ftp
            + 0.5 * invrdot * fp * fp
            + Vp * ft
        ) / scale

        flux = sqrt_gamma * alpha * invr * p * fp
        flux_div = _d1(grid, flux, -1) / sqrt_gamma
        lapse_source = invr * alphap * p * fp
        trace_K_source = alpha * K * p * p
        inverse_metric_source = 0.5 * invrdot * fp * fp
        scalar_source = p * source
        balance = (
            flux_div
            + lapse_source
            + trace_K_source
            + inverse_metric_source
            + scalar_source
        ) / scale
        residual = chain - balance

        # Stress projections allow an independent reconstruction of the
        # geometric contraction K_ij S^ij used in the covariant identity.
        pr = (0.5 * p * p + 0.5 * invr * fp * fp - V) / scale
        ptan = (0.5 * p * p - 0.5 * invr * fp * fp - V) / scale
        KijSij = kr_r * pr + 2.0 * (K / 3.0 - Aa / 2.0) * ptan

        n = len(r)
        bulk_end = max(5, n - 5)
        result[name] = {
            "energy_density_first_five": first(rho),
            "chain_rule_rate_first_five": first(chain),
            "flux_divergence_first_five": first(flux_div / scale),
            "lapse_gradient_term_first_five": first(lapse_source / scale),
            "trace_K_term_first_five": first(trace_K_source / scale),
            "inverse_metric_term_first_five": first(inverse_metric_source / scale),
            "explicit_scalar_source_first_five": first(scalar_source / scale),
            "balance_rate_first_five": first(balance),
            "residual_first_five": first(residual),
            "max_abs_residual_cells_0_to_4": float(np.max(np.abs(residual[:min(5,n)]))),
            "max_abs_residual_bulk_cells_5_to_n_minus_5": float(
                np.max(np.abs(residual[5:bulk_end])) if bulk_end > 5 else 0.0
            ),
            "max_abs_residual_outer_five_cells": float(np.max(np.abs(residual[max(0,n-5):]))),
            "max_abs_residual_all_cells": float(np.max(np.abs(residual))),
            "max_abs_KijSij": float(np.max(np.abs(KijSij))),
            "all_finite": bool(
                np.all(np.isfinite(rho)) and np.all(np.isfinite(chain))
                and np.all(np.isfinite(balance)) and np.all(np.isfinite(residual))
            ),
        }

    result["_geometry"] = {
        "inverse_metric_rate_actual_first_five": first(invrdot),
        "inverse_metric_rate_adm_first_five": first(invrdot_adm),
        "inverse_metric_rate_closure_error_first_five": first(invrdot_error),
        "max_abs_inverse_metric_rate_closure_error_bulk": float(
            np.max(np.abs(invrdot_error[5:max(5,len(r)-5)]))
        ),
        "max_abs_inverse_metric_rate_closure_error_all_cells": float(
            np.max(np.abs(invrdot_error))
        ),
        "all_finite": bool(np.all(np.isfinite(invrdot_error))),
    }
    return result


def run_case():
    kernel = V55ProductionKernel()
    state, B, residual_history, solve = discrete_consistent_state(
        kernel,
        resolution=SETTINGS["resolution"],
        r_max=SETTINGS["r_max"],
        amplitude=SETTINGS["amplitude"],
        width=SETTINGS["width"],
        D_amplitude=SETTINGS["D_amplitude"],
        include_radiation=SETTINGS["include_radiation"],
        max_iter=SETTINGS["max_iter"],
        tol=SETTINGS["tol"],
        max_cond=SETTINGS["max_cond"],
        return_info=True,
    )

    H_before_projection = _hamiltonian(state)
    pre_geometry = state.geometry.copy()
    state.geometry = kernel._enforce_center_regularity(
        state.grid, state.geometry.copy()
    )
    H_after_projection = _hamiltonian(state)

    # Projection is now complete. Re-solve the gauge, without mixing this
    # gauge update into the separately recorded projection jump.
    state.geometry.alpha = np.asarray(
        kernel._solve_lapse(
            state.grid, state.geometry, state.scalars, state.matter
        )[0]
    ).copy()
    state.geometry.beta.fill(0.0)
    state.geometry.B.fill(0.0)

    metric = metric_slice_from_q(state.grid, state.geometry)
    rho_dm = dm_density(metric, state.matter)
    srhs = scalar_rhs_arrays(
        state.grid, state.geometry, state.scalars,
        beta_dm=BETA_DM, rho_dm=rho_dm,
    )
    geometry_rhs = adapter.geometry_stage_terms(
        state.grid, state.geometry, state.scalars, state.matter,
        lambda_m=2.0,
    )
    report_terms = _energy_terms(
        state.grid, state.geometry, state.scalars, srhs, geometry_rhs,
        source_phi=BETA_DM * rho_dm,
    )

    projection_delta = state.geometry.a - pre_geometry.a
    projection_delta_b = state.geometry.b - pre_geometry.b
    projection_delta_X = state.geometry.X - pre_geometry.X
    finite = all(
        entry.get("all_finite", False)
        for name, entry in report_terms.items()
    )

    report = {
        "schema": "curved_scalar_energy_balance_v1",
        "kind": "diagnostic_only_instantaneous_curved_scalar_balance",
        "source_commit": os.environ.get("GITHUB_SHA", "unavailable"),
        "settings": SETTINGS,
        "admission": "NOT_ADMITTED_DIAGNOSTIC_ONLY",
        "production_equations_changed": False,
        "production_defaults_changed": False,
        "physical_parameters_changed": False,
        "initial_solver": {
            "iterations": int(solve.get("iterations", -1)),
            "final_max_residual": float(solve.get("final_max_residual", math.nan)),
            "max_condition_number": float(solve.get("max_condition_number", math.nan)),
            "residual_history": [float(x) for x in residual_history],
            "B_min": float(np.min(B)),
            "B_max": float(np.max(B)),
        },
        "projection_separated": {
            "max_abs_delta_a": float(np.max(np.abs(projection_delta))),
            "max_abs_delta_b": float(np.max(np.abs(projection_delta_b))),
            "max_abs_delta_X": float(np.max(np.abs(projection_delta_X))),
            "H_first_five_before_projection": first(H_before_projection),
            "H_first_five_after_projection": first(H_after_projection),
            "delta_H_first_five_due_to_projection": first(
                H_after_projection - H_before_projection
            ),
            "max_abs_delta_H_all_cells_due_to_projection": float(
                np.max(np.abs(H_after_projection - H_before_projection))
            ),
        },
        "post_projection_state": {
            "t": float(state.t),
            "alpha_first_five": first(state.geometry.alpha),
            "H_first_five": first(_hamiltonian(state)),
            "max_abs_H_all_cells": float(np.max(np.abs(_hamiltonian(state)))),
        },
        "balance": report_terms,
        "summary": {
            "all_finite": bool(finite),
            "admission": "NOT_ADMITTED_DIAGNOSTIC_ONLY",
            "interpretation_guard": (
                "Instantaneous balance residuals are numerical diagnostics, not "
                "evidence for a physical event. Projection and post-projection "
                "evolution-rate effects are separate."
            ),
        },
    }
    out = Path("runs/curved-scalar-energy-balance")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    print(json.dumps({
        "schema": report["schema"],
        "source_commit": report["source_commit"],
        "settings": report["settings"],
        "projection_separated": report["projection_separated"],
        "post_projection_state": report["post_projection_state"],
        "balance": {
            name: {
                "max_abs_residual_cells_0_to_4": values.get("max_abs_residual_cells_0_to_4"),
                "max_abs_residual_bulk_cells_5_to_n_minus_5": values.get("max_abs_residual_bulk_cells_5_to_n_minus_5"),
                "max_abs_residual_all_cells": values.get("max_abs_residual_all_cells"),
                "all_finite": values.get("all_finite"),
            }
            for name, values in report_terms.items()
        },
        "summary": report["summary"],
    }, indent=2, allow_nan=False), flush=True)
    if not finite:
        raise RuntimeError("Non-finite values in curved scalar energy balance")
    return report


if __name__ == "__main__":
    run_case()
