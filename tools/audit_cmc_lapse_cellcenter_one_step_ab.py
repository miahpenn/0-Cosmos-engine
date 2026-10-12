#!/usr/bin/env python3
"""Matched one-step A/B for the CMC lapse centre-row candidate.

The old tridiagonal CMC solver is retained here as an audit-only reference
function. Both runs start with identical non-lapse geometry, scalar fields and
conservative matter. The only intended difference is the lapse row/operator.
No long trajectory or source/physics change is part of this audit.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys
import numpy as np
from scipy.linalg import solve_banded

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.audit_uniform_radial_quadrature_initial_profile import candidate_uniform_initial_state
from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse, target_kdot
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

RESOLUTIONS = (40, 60, 80, 120, 160)
R_MAX = 80.0
CFL = 0.03
OUT = Path("results/cmc-lapse-cellcenter-one-step-ab")


def old_tridiagonal_cmc_lapse(grid, geometry, scalars, matter, outer_frac=0.20):
    """Exact audit copy of the pre-candidate CMC matrix row construction."""
    del outer_frac
    metric_total = assemble_total_stress_energy(grid, geometry, scalars, matter)
    a = np.asarray(geometry.a)
    b = np.asarray(geometry.b)
    X = np.asarray(geometry.X)
    r = np.asarray(grid.centers)
    inv = X * X / a

    d1 = lambda values, parity: grid.cell_derivative_fourth(values, parity=parity)
    ap = d1(a, 1)
    bp = d1(b, 1)
    Xp = d1(X, 1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r

    Q = (
        1.5 * geometry.Aa * geometry.Aa
        + geometry.K * geometry.K / 3.0
        + 4.0 * math.pi * (
            metric_total.rho + metric_total.pr + 2.0 * metric_total.pt
        )
    )
    kdot = target_kdot(grid, geometry, scalars, matter)
    m = Q / inv
    rhs = -kdot / inv

    n = grid.n
    h = grid.dr
    band = np.zeros((3, n), dtype=float)
    bvec = np.zeros(n, dtype=float)

    # Historical row: on a cell-centred grid this incorrectly equates alpha
    # at r=dr/2 and r=3dr/2, while not enforcing the PDE at cell 0.
    band[1, 0] = 1.0
    band[0, 1] = -1.0
    for i in range(1, n - 1):
        ci = c[i] / (2.0 * h)
        band[2, i - 1] = 1.0 / h**2 - ci
        band[1, i] = -2.0 / h**2 - m[i]
        band[0, i + 1] = 1.0 / h**2 + ci
        bvec[i] = rhs[i]
    band[1, n - 1] = 1.0
    bvec[n - 1] = 1.0

    alpha = solve_banded((1, 1), band, bvec, check_finite=False)
    if not np.all(np.isfinite(alpha)) or float(np.min(alpha)) <= 0.0:
        raise ValueError("old CMC reference produced invalid/nonpositive lapse")
    return alpha, kdot


def clone_state(template: ProductionState, alpha: np.ndarray) -> ProductionState:
    g = template.geometry.copy()
    g.alpha = np.asarray(alpha, dtype=float).copy()
    g.beta.fill(0.0)
    g.B.fill(0.0)
    return ProductionState(
        grid=template.grid,
        geometry=g,
        scalars=template.scalars.copy(),
        matter=template.matter.copy(),
        t=template.t,
        tau=template.tau,
        e_folds=template.e_folds,
    )


def _region(values: np.ndarray, r: np.ndarray, mask: np.ndarray) -> dict:
    ids = np.flatnonzero(mask & np.isfinite(values))
    if ids.size == 0:
        return {"count": 0, "max_abs": None, "rms": None, "peak_r": None, "peak_signed": None}
    peak = int(ids[np.argmax(np.abs(values[ids]))])
    return {
        "count": int(ids.size),
        "max_abs": float(abs(values[peak])),
        "rms": float(np.sqrt(np.mean(values[ids] ** 2))),
        "peak_r": float(r[peak]),
        "peak_signed": float(values[peak]),
    }


def summarize_state(state: ProductionState, kernel: V55ProductionKernel) -> dict:
    grid, g = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=float)
    _, vacuum, _ = adapter.vendor_modules()
    total = assemble_total_stress_energy(grid, g, state.scalars, state.matter)
    raw = vacuum.constraints(grid, g)
    H = np.asarray(raw["hamiltonian"], dtype=float) - 16.0 * math.pi * np.asarray(total.rho, dtype=float)
    M = np.asarray(raw["momentum"], dtype=float) - 8.0 * math.pi * np.asarray(total.j, dtype=float)
    C = np.asarray(raw["connection"], dtype=float)

    stages = adapter.geometry_stage_terms(grid, g, state.scalars, state.matter, lambda_m=2.0)
    Aa_rhs = np.asarray(stages["primary_l2"]["Aa"]) + np.asarray(stages["primary_l3"]["Aa"])
    idx = np.arange(grid.n)
    masks = {
        "full_domain": np.ones(grid.n, dtype=bool),
        "first_two_cells": idx < 2,
        "first_three_cells": idx < 3,
        "admitted_inner_i_ge_2_r_lt_5": (idx >= 2) & (r < 5.0),
        "fixed_annulus_2_le_r_lt_5": (r >= 2.0) & (r < 5.0),
        "outer_r_ge_0p8_rmax": r >= 0.8 * grid.r_max,
    }
    constraint_norms = {
        name: {region: _region(values, r, mask) for region, mask in masks.items()}
        for name, values in (("H", H), ("M", M), ("C_lambda", C))
    }

    # Check what the selected solver would return for this exact frozen state
    # and compare its fourth-order PDE residual with the production K RHS.
    alpha_check, kdot_check = kernel._solve_lapse(grid, g.copy(), state.scalars, state.matter)
    gcheck = g.copy()
    gcheck.alpha = alpha_check
    inv = np.asarray(gcheck.X) ** 2 / np.asarray(gcheck.a)
    a_r = grid.cell_derivative_fourth(gcheck.a, parity=1)
    b_r = grid.cell_derivative_fourth(gcheck.b, parity=1)
    X_r = grid.cell_derivative_fourth(gcheck.X, parity=1)
    alpha_r = grid.cell_derivative_fourth(gcheck.alpha, parity=1)
    alpha_rr = grid.cell_second_derivative_fourth(gcheck.alpha, parity=1)
    c = -0.5 * a_r / gcheck.a + b_r / gcheck.b - X_r / gcheck.X + 2.0 / r
    total_check = assemble_total_stress_energy(grid, gcheck, state.scalars, state.matter)
    Q = (
        1.5 * np.asarray(gcheck.Aa) ** 2
        + np.asarray(gcheck.K) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total_check.rho) + np.asarray(total_check.pr)
            + 2.0 * np.asarray(total_check.pt)
        )
    )
    row_residual = alpha_rr + c * alpha_r - (Q / inv) * gcheck.alpha + kdot_check / inv
    stages_check = adapter.geometry_stage_terms(
        grid, gcheck, state.scalars, state.matter, lambda_m=2.0
    )
    k_rhs_check = np.asarray(stages_check["primary_l2"]["K"]) + np.asarray(stages_check["primary_l3"]["K"])
    pde_mask = idx < grid.n - 1

    # Kernel diagnostics provide the invariant/worldtube quantities already
    # defined by production. No derived current is fitted from one step.
    obs = kernel.diagnostics(state, profiles=False)
    return {
        "t": float(state.t),
        "tau": float(state.tau),
        "dtau_over_dt": float(obs["tau_rate"]),
        "alpha_min": float(np.min(g.alpha)),
        "alpha_max": float(np.max(g.alpha)),
        "alpha0_minus_alpha1": float(g.alpha[0] - g.alpha[1]),
        "constraints": constraint_norms,
        "centre_cells": {
            str(i): {
                "r": float(r[i]),
                "H": float(H[i]),
                "M": float(M[i]),
                "C_lambda": float(C[i]),
                "Aa": float(g.Aa[i]),
                "Aa_over_r2": float(g.Aa[i] / r[i]**2),
                "Aa_rhs": float(Aa_rhs[i]),
                "Aa_rhs_over_r2": float(Aa_rhs[i] / r[i]**2),
                "alpha": float(g.alpha[i]),
            }
            for i in range(min(3, grid.n))
        },
        "cmc_operator": {
            "max_abs_row_residual_i0_to_nminus2": float(np.max(np.abs(row_residual[pde_mask]))),
            "max_abs_production_Kdot_minus_target_i0_to_nminus2": float(
                np.max(np.abs((k_rhs_check - kdot_check)[pde_mask]))
            ),
            "centre_row_residual": float(row_residual[0]),
            "first_interior_row_residual": float(row_residual[1]),
        },
        "worldtube_r10": {
            "M_MS": float(obs["M_MS"]),
            "R_sigma": float(obs["R_sigma"]),
            "H_eff": float(obs["H_eff"]),
            "flux_T": float(obs["flux_T"]),
            "work_pR": float(obs["work_pR"]),
            "coordinate_energy_flux": float(obs["coordinate_energy_flux"]),
            "Rdot_sigma": float(obs["Rdot_sigma"]),
        },
        "all_finite": bool(
            np.all(np.isfinite(H)) and np.all(np.isfinite(M))
            and np.all(np.isfinite(C)) and np.all(np.isfinite(Aa_rhs))
            and np.all(np.isfinite(alpha_check))
        ),
        "evolution_advanced": state.t > 0.0,
    }


def one_case(n: int) -> dict:
    candidate_kernel = V55ProductionKernel()
    grid_ops, _, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    # This helper fixes the uniform-quadrature geometry and all fields. We then
    # restore the common pre-lapse state (alpha=1) before solving both alphas.
    template = candidate_uniform_initial_state(candidate_kernel, grid)
    g_template = template.geometry.copy()
    g_template.alpha.fill(1.0)
    g_template.beta.fill(0.0)
    g_template.B.fill(0.0)

    alpha_old, target_old = old_tridiagonal_cmc_lapse(
        grid, g_template.copy(), template.scalars, template.matter
    )
    alpha_new, target_new = solve_cmc_lapse(
        grid, g_template.copy(), template.scalars, template.matter
    )
    old_state = clone_state(template, alpha_old)
    new_state = clone_state(template, alpha_new)

    old_kernel = V55ProductionKernel()
    # Instance-level replacement makes every old/reference stage use the
    # historical CMC matrix row, while the candidate stages use the new row.
    old_kernel._solve_lapse = old_tridiagonal_cmc_lapse

    dt = CFL * float(grid.dr)
    case = {
        "N": n,
        "r_max": R_MAX,
        "dt": dt,
        "CFL": CFL,
        "alpha_target_from_common_alpha1_slice": {
            "old": float(target_old),
            "candidate": float(target_new),
            "difference": float(target_new - target_old),
        },
        "initial_nonlapse_state_shared": True,
        "old": {},
        "candidate": {},
        "status": "started",
    }
    try:
        case["old"]["initial"] = summarize_state(old_state, old_kernel)
        case["candidate"]["initial"] = summarize_state(new_state, candidate_kernel)
        old_final = old_kernel.step(old_state, dt)
        case["old"]["final"] = summarize_state(old_final, old_kernel)
        new_final = candidate_kernel.step(new_state, dt)
        case["candidate"]["final"] = summarize_state(new_final, candidate_kernel)
        case["candidate_vs_old_final"] = {
            "H_full_absmax_delta": (
                case["candidate"]["final"]["constraints"]["H"]["full_domain"]["max_abs"]
                - case["old"]["final"]["constraints"]["H"]["full_domain"]["max_abs"]
            ),
            "M_full_absmax_delta": (
                case["candidate"]["final"]["constraints"]["M"]["full_domain"]["max_abs"]
                - case["old"]["final"]["constraints"]["M"]["full_domain"]["max_abs"]
            ),
            "C_lambda_full_absmax_delta": (
                case["candidate"]["final"]["constraints"]["C_lambda"]["full_domain"]["max_abs"]
                - case["old"]["final"]["constraints"]["C_lambda"]["full_domain"]["max_abs"]
            ),
            "Aa_rhs_r2_cell0_delta": (
                case["candidate"]["final"]["centre_cells"]["0"]["Aa_rhs_over_r2"]
                - case["old"]["final"]["centre_cells"]["0"]["Aa_rhs_over_r2"]
            ),
        }
        case["status"] = "completed"
    except Exception as exc:
        case["status"] = "failed"
        case["error"] = f"{type(exc).__name__}: {exc}"
    return case


def main() -> None:
    cases = [one_case(n) for n in RESOLUTIONS]
    report = {
        "schema": "cmc_lapse_cellcenter_one_step_ab_v1",
        "status": "completed" if all(c["status"] == "completed" for c in cases) else "partial_failure",
        "diagnostic_only": True,
        "source_under_test": "historical second-order CMC centre row vs parity-aware fourth-order CMC row",
        "pinned_geometry_vendor": "d6052d605673ce9d82cdc99b0df79855fa2ea215",
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "dt_rule": "dt=0.03*dr; exactly one ordinary production step per branch and resolution",
        "cases": cases,
        "guardrails": [
            "Same uniform-quadrature non-lapse geometry, scalar state and conservative matter in A/B.",
            "Only the CMC lapse solver row/operator differs; both variants use the same production physics sources.",
            "No physical source, parameter, threshold, damping, floor, clamp or boundary law is added.",
            "N=40/60/80/120/160 only; high-resolution radiation failure remains a separate gate.",
            "A one-step pass is not long-time or full-machine admission."
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print("CMC_LAPSE_CELLCENTER_ONE_STEP_AB " + json.dumps({
        "status": report["status"],
        "cases": [{
            "N": c["N"], "status": c["status"], "error": c.get("error"),
            "old_cmc_residual": c.get("old", {}).get("final", {}).get("cmc_operator", {}).get("max_abs_row_residual_i0_to_nminus2"),
            "candidate_cmc_residual": c.get("candidate", {}).get("final", {}).get("cmc_operator", {}).get("max_abs_row_residual_i0_to_nminus2"),
            "old_H_full": c.get("old", {}).get("final", {}).get("constraints", {}).get("H", {}).get("full_domain", {}).get("max_abs"),
            "candidate_H_full": c.get("candidate", {}).get("final", {}).get("constraints", {}).get("H", {}).get("full_domain", {}).get("max_abs"),
            "old_M_inner": c.get("old", {}).get("final", {}).get("constraints", {}).get("M", {}).get("admitted_inner_i_ge_2_r_lt_5", {}).get("max_abs"),
            "candidate_M_inner": c.get("candidate", {}).get("final", {}).get("constraints", {}).get("M", {}).get("admitted_inner_i_ge_2_r_lt_5", {}).get("max_abs"),
            "old_Aa_rhs_r2_cell0": c.get("old", {}).get("final", {}).get("centre_cells", {}).get("0", {}).get("Aa_rhs_over_r2"),
            "candidate_Aa_rhs_r2_cell0": c.get("candidate", {}).get("final", {}).get("centre_cells", {}).get("0", {}).get("Aa_rhs_over_r2"),
        } for c in cases],
    }, sort_keys=True))
    if report["status"] != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
