#!/usr/bin/env python3
"""Three-arm one-step test separating CMC centre row from stencil order.

Arms:
  A: historical second-order interior plus alpha[0] == alpha[1];
  B: second-order interior plus a second-order even-parity PDE at cell 0;
  C: parity-aware fourth-order operator on all CMC PDE rows.

All arms share the same non-lapse initial data and physical sources. Only the
CMC lapse solve differs. No long trajectory or production physics edit occurs.
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
from tools.audit_cmc_lapse_cellcenter_one_step_ab import (
    CFL, R_MAX, clone_state, old_tridiagonal_cmc_lapse, summarize_state,
)
from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse, target_kdot
from engine.production_kernel import V55ProductionKernel, ProductionState
from engine.stress_energy import assemble_total_stress_energy

RESOLUTIONS = (80, 120, 160)
OUT = Path("results/cmc-lapse-three-arm-row-order")


def second_order_center_pde_lapse(grid, geometry, scalars, matter, outer_frac=0.20):
    """Historical second-order CMC interior, but PDE row at first cell.

    With even reflection on cell-centred values:
      alpha_r[0]  = (alpha[1]-alpha[0])/(2h)
      alpha_rr[0] = (alpha[1]-alpha[0])/h^2
    Thus the first row applies the same second-order discretization of the
    lapse equation rather than imposing alpha[0] == alpha[1].
    """
    del outer_frac
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    r = np.asarray(grid.centers, dtype=float)
    inv = X * X / a

    d1 = lambda values, parity: np.asarray(
        grid.cell_derivative_fourth(values, parity=parity), dtype=float
    )
    ap = d1(a, 1)
    bp = d1(b, 1)
    Xp = d1(X, 1)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r

    total = assemble_total_stress_energy(grid, geometry, scalars, matter)
    Q = (
        1.5 * np.asarray(geometry.Aa) ** 2
        + np.asarray(geometry.K) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total.rho) + np.asarray(total.pr)
            + 2.0 * np.asarray(total.pt)
        )
    )
    target = target_kdot(grid, geometry, scalars, matter)
    m = Q / inv
    rhs = -target / inv

    n = grid.n
    h = float(grid.dr)
    band = np.zeros((3, n), dtype=float)
    bvec = np.zeros(n, dtype=float)

    # Corrected second-order centre PDE row. The outer bandwidth convention
    # is scipy.solve_banded((1,1), ...): band[u + row-col, col].
    band[1, 0] = -1.0 / h**2 - c[0] / (2.0 * h) - m[0]
    band[0, 1] = 1.0 / h**2 + c[0] / (2.0 * h)
    bvec[0] = rhs[0]

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
        raise ValueError("second-order centre-PDE CMC lapse is invalid/nonpositive")
    return alpha, target


def arm_summary(state: ProductionState, kernel: V55ProductionKernel) -> dict:
    summary = summarize_state(state, kernel)
    c = state.grid.centers
    cell0 = summary["centre_cells"]["0"]
    cell1 = summary["centre_cells"]["1"]
    cell2 = summary["centre_cells"]["2"]
    return {
        "t": summary["t"],
        "tau": summary["tau"],
        "alpha_min": summary["alpha_min"],
        "alpha_max": summary["alpha_max"],
        "alpha0_minus_alpha1": cell0["alpha"] - cell1["alpha"],
        "H_full": summary["constraints"]["H"]["full_domain"],
        "H_annulus": summary["constraints"]["H"]["fixed_annulus_2_le_r_lt_5"],
        "M_full": summary["constraints"]["M"]["full_domain"],
        "M_inner": summary["constraints"]["M"]["admitted_inner_i_ge_2_r_lt_5"],
        "M_annulus": summary["constraints"]["M"]["fixed_annulus_2_le_r_lt_5"],
        "C_full": summary["constraints"]["C_lambda"]["full_domain"],
        "C_annulus": summary["constraints"]["C_lambda"]["fixed_annulus_2_le_r_lt_5"],
        "centre_cells": {
            "0": {k: cell0[k] for k in ("r", "H", "M", "C_lambda", "alpha", "Aa_rhs_over_r2")},
            "1": {k: cell1[k] for k in ("r", "H", "M", "C_lambda", "alpha", "Aa_rhs_over_r2")},
            "2": {k: cell2[k] for k in ("r", "H", "M", "C_lambda", "alpha", "Aa_rhs_over_r2")},
        },
        "worldtube_r10": summary["worldtube_r10"],
        "finite": summary["all_finite"],
    }


def run_case(n: int) -> dict:
    common_kernel = V55ProductionKernel()
    grid_ops, _, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(n, R_MAX)
    template = candidate_uniform_initial_state(common_kernel, grid)
    g0 = template.geometry.copy()
    g0.alpha.fill(1.0)
    g0.beta.fill(0.0)
    g0.B.fill(0.0)

    alpha_a, target_a = old_tridiagonal_cmc_lapse(
        grid, g0.copy(), template.scalars, template.matter
    )
    alpha_b, target_b = second_order_center_pde_lapse(
        grid, g0.copy(), template.scalars, template.matter
    )
    alpha_c, target_c = solve_cmc_lapse(
        grid, g0.copy(), template.scalars, template.matter
    )
    states = {
        "A_old_second_order_equality_row": clone_state(template, alpha_a),
        "B_second_order_centre_PDE_row": clone_state(template, alpha_b),
        "C_fourth_order_production_matched": clone_state(template, alpha_c),
    }
    kernels = {
        "A_old_second_order_equality_row": V55ProductionKernel(),
        "B_second_order_centre_PDE_row": V55ProductionKernel(),
        "C_fourth_order_production_matched": V55ProductionKernel(),
    }
    kernels["A_old_second_order_equality_row"]._solve_lapse = old_tridiagonal_cmc_lapse
    kernels["B_second_order_centre_PDE_row"]._solve_lapse = second_order_center_pde_lapse

    dt = CFL * float(grid.dr)
    row = {
        "N": n,
        "dr": float(grid.dr),
        "dt": dt,
        "targets_from_common_alpha1_slice": {
            "A": float(target_a), "B": float(target_b), "C": float(target_c),
            "max_abs_difference": float(max(abs(target_a-target_b), abs(target_a-target_c), abs(target_b-target_c))),
        },
        "identical_nonlapse_initial_state": True,
        "arms": {},
        "pairwise_final_delta": {},
        "status": "started",
    }
    try:
        for name, state in states.items():
            row["arms"][name] = {"initial": arm_summary(state, kernels[name])}
        final_states = {}
        for name, state in states.items():
            final_states[name] = kernels[name].step(state, dt)
            row["arms"][name]["final"] = arm_summary(final_states[name], kernels[name])

        names = list(states)
        for i, left in enumerate(names):
            for right in names[i+1:]:
                a = row["arms"][left]["final"]
                b = row["arms"][right]["final"]
                row["pairwise_final_delta"][f"{right}_minus_{left}"] = {
                    "full_H_absmax": b["H_full"]["max_abs"] - a["H_full"]["max_abs"],
                    "annulus_H_absmax": b["H_annulus"]["max_abs"] - a["H_annulus"]["max_abs"],
                    "full_M_absmax": b["M_full"]["max_abs"] - a["M_full"]["max_abs"],
                    "inner_M_absmax": (
                        None if a["M_inner"]["max_abs"] is None or b["M_inner"]["max_abs"] is None
                        else b["M_inner"]["max_abs"] - a["M_inner"]["max_abs"]
                    ),
                    "annulus_M_absmax": b["M_annulus"]["max_abs"] - a["M_annulus"]["max_abs"],
                    "annulus_C_absmax": b["C_annulus"]["max_abs"] - a["C_annulus"]["max_abs"],
                    "Aa_rhs_r2_cell0": (
                        b["centre_cells"]["0"]["Aa_rhs_over_r2"]
                        - a["centre_cells"]["0"]["Aa_rhs_over_r2"]
                    ),
                    "cell0_H": b["centre_cells"]["0"]["H"] - a["centre_cells"]["0"]["H"],
                    "cell1_H": b["centre_cells"]["1"]["H"] - a["centre_cells"]["1"]["H"],
                }
        row["status"] = "completed"
    except Exception as exc:
        row["status"] = "failed"
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def main() -> None:
    cases = [run_case(n) for n in RESOLUTIONS]
    report = {
        "schema": "cmc_lapse_three_arm_row_order_v1",
        "status": "completed" if all(c["status"] == "completed" for c in cases) else "partial_failure",
        "diagnostic_only": True,
        "resolutions": list(RESOLUTIONS),
        "r_max": R_MAX,
        "dt_rule": "dt=0.03*dr; one production step per arm",
        "arm_definitions": {
            "A_old_second_order_equality_row": "second-order interior, alpha[0]==alpha[1] row",
            "B_second_order_centre_PDE_row": "same second-order interior, proper even-reflection second-order PDE row at cell 0",
            "C_fourth_order_production_matched": "fourth-order parity-aware PDE rows matching production K RHS, with outer alpha=1 normalization",
        },
        "cases": cases,
        "guardrails": [
            "Same initial geometry except lapse, same scalar fields, same conservative matter and same physics sources.",
            "Target Kdot computed from identical pre-lapse state is reported for all arms.",
            "No physical source, coefficient, radiation rule, damping, floor, clamp or projection is changed.",
            "No long trajectory is advanced; N=80/120/160 are within previously admitted radiation resolutions.",
            "This separates centre-row treatment from interior CMC stencil order but is still a one-step diagnostic only.",
        ],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    printable = []
    for c in cases:
        printable.append({
            "N": c["N"], "status": c["status"], "error": c.get("error"),
            "A_full_H": c.get("arms", {}).get("A_old_second_order_equality_row", {}).get("final", {}).get("H_full", {}).get("max_abs"),
            "B_full_H": c.get("arms", {}).get("B_second_order_centre_PDE_row", {}).get("final", {}).get("H_full", {}).get("max_abs"),
            "C_full_H": c.get("arms", {}).get("C_fourth_order_production_matched", {}).get("final", {}).get("H_full", {}).get("max_abs"),
            "A_annulus_M": c.get("arms", {}).get("A_old_second_order_equality_row", {}).get("final", {}).get("M_annulus", {}).get("max_abs"),
            "B_annulus_M": c.get("arms", {}).get("B_second_order_centre_PDE_row", {}).get("final", {}).get("M_annulus", {}).get("max_abs"),
            "C_annulus_M": c.get("arms", {}).get("C_fourth_order_production_matched", {}).get("final", {}).get("M_annulus", {}).get("max_abs"),
            "pairwise": c.get("pairwise_final_delta", {}),
        })
    print("CMC_LAPSE_THREE_ARM_ROW_ORDER " + json.dumps({"status": report["status"], "cases": printable}, sort_keys=True))
    if report["status"] != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
