#!/usr/bin/env python3
"""No-evolution CMC lapse row A/B on the archived radiation failure state.

The physical geometry, matter, scalar fields, and target Kdot remain frozen.
Only the elliptic lapse linear system's centre and outer row placement changes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.cmc_gauge import solve_cmc_lapse, target_kdot
from engine.v55_matter import total_matter_projection
from tools.audit_cmc_lapse_fixed_point import (
    find_artifact_dir,
    load_terminal_state,
    physical_state_digest,
    proper_volume_weights,
    verify_artifact,
    weighted_mean,
)

EXPECTED_MAIN = "c01f2467347e8aebd57d220374fd7bc4e93045bc"
EXPECTED_VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_INPUT_HEAD = "a854c3d6eafde1a1bae8066b7cc212bac7511ffd"
EXPECTED_INPUT_RUN = 38337821862


def array_digest(*arrays: np.ndarray) -> str:
    digest = hashlib.sha256()
    for index, value in enumerate(arrays):
        arr = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
        digest.update(str(index).encode("ascii") + b"\0")
        digest.update(str(arr.shape).encode("ascii") + b"\0")
        digest.update(arr.tobytes())
    return digest.hexdigest()


def profile_stats(values: np.ndarray, start: int = 0, stop: int | None = None) -> dict:
    values = np.asarray(values, dtype=float)
    if stop is None:
        stop = values.size
    if not (0 <= start < stop <= values.size):
        raise ValueError(f"bad statistics slice {start}:{stop} for length {values.size}")
    v = values[start:stop]
    if not np.all(np.isfinite(v)):
        raise RuntimeError("nonfinite value in reported CMC counterfactual profile")
    i_local = int(np.argmax(np.abs(v)))
    return {
        "cell_start_inclusive": int(start),
        "cell_stop_exclusive": int(stop),
        "count": int(v.size),
        "max_abs": float(abs(v[i_local])),
        "max_signed": float(v[i_local]),
        "max_abs_cell": int(start + i_local),
        "rms": float(np.sqrt(np.mean(v * v))),
    }


def even_series(alpha: np.ndarray, r: np.ndarray) -> dict:
    count = min(8, len(alpha))
    z = np.asarray(r[:count], dtype=float) ** 2
    y = np.asarray(alpha[:count], dtype=float)
    coeff = np.polynomial.polynomial.polyfit(z, y, deg=min(3, count - 1))
    residual = y - np.polynomial.polynomial.polyval(z, coeff)
    return {
        "coefficients_alpha_of_r2_ascending": [float(x) for x in coeff],
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rms_residual": float(np.sqrt(np.mean(residual**2))),
        "cell0_minus_fit": float(residual[0]),
        "alpha0_extrapolated": float(coeff[0]),
    }


def _operator_row(grid, order: int, cell: int = 0) -> np.ndarray:
    n = int(grid.n)
    row = np.zeros(n, dtype=float)
    for j in range(n):
        basis = np.zeros(n, dtype=float)
        basis[j] = 1.0
        if order == 1:
            row[j] = float(grid.cell_derivative_fourth(basis, parity=1)[cell])
        elif order == 2:
            row[j] = float(grid.cell_second_derivative_fourth(basis, parity=1)[cell])
        else:
            raise ValueError("operator order must be 1 or 2")
    return row


def assemble_rows(grid, geometry, scalars, matter, target: float,
                  center_row: str, outer_row: str) -> tuple[np.ndarray, dict]:
    """Build a frozen-target linear system; returns matrix and row semantics."""
    n = int(grid.n)
    h = float(grid.dr)
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    inv = X * X / a
    ap = np.asarray(grid.cell_derivative_fourth(a, parity=1), dtype=float)
    bp = np.asarray(grid.cell_derivative_fourth(b, parity=1), dtype=float)
    Xp = np.asarray(grid.cell_derivative_fourth(X, parity=1), dtype=float)
    c = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r
    total = total_matter_projection(grid, geometry, scalars, matter)
    rho, pr, pt = (np.asarray(total[key], dtype=float) for key in ("rho", "pr", "pt"))
    Q = 1.5 * np.asarray(geometry.Aa, dtype=float)**2
    Q += np.asarray(geometry.K, dtype=float)**2 / 3.0
    Q += 4.0 * math.pi * (rho + pr + 2.0 * pt)
    m = Q / inv
    pde_rhs = -float(target) / inv
    A = np.zeros((n, n), dtype=float)
    v = np.zeros(n, dtype=float)

    # Keep all interior rows exactly as the existing second-order solver uses them.
    for i in range(1, n - 1):
        ci = c[i] / (2.0 * h)
        A[i, i - 1] = 1.0 / h**2 - ci
        A[i, i] = -2.0 / h**2 - m[i]
        A[i, i + 1] = 1.0 / h**2 + ci
        v[i] = pde_rhs[i]

    if center_row == "current_equal_adjacent":
        A[0, 0] = 1.0
        A[0, 1] = -1.0
    elif center_row == "second_order_even_ghost_pde":
        A[0, 0] = -1.0 / h**2 - c[0] / (2.0 * h) - m[0]
        A[0, 1] = 1.0 / h**2 + c[0] / (2.0 * h)
        v[0] = pde_rhs[0]
    elif center_row == "vendor_fourth_order_even_pde":
        A[0, :] = _operator_row(grid, 2, 0) + c[0] * _operator_row(grid, 1, 0)
        A[0, 0] -= m[0]
        v[0] = pde_rhs[0]
    else:
        raise ValueError(f"unknown center-row option {center_row}")

    if outer_row == "last_cell_center_alpha_one":
        A[-1, :] = 0.0
        A[-1, -1] = 1.0
        v[-1] = 1.0
    elif outer_row == "face_alpha_one_with_last_cell_pde":
        # alpha(r_max)=1 on the face; linearly reflect across the face:
        # alpha_N = 2-alpha_{N-1}. Enforce the PDE at center N-1.
        A[-1, :] = 0.0
        A[-1, -2] = 1.0 / h**2 - c[-1] / (2.0 * h)
        A[-1, -1] = -3.0 / h**2 - c[-1] / (2.0 * h) - m[-1]
        v[-1] = pde_rhs[-1] - 2.0 / h**2 - c[-1] / h
    else:
        raise ValueError(f"unknown outer-row option {outer_row}")

    alpha = np.linalg.solve(A, v)
    if not np.all(np.isfinite(alpha)) or float(np.min(alpha)) <= 0.0:
        raise RuntimeError(f"row system produced invalid lapse: {center_row}/{outer_row}")
    return alpha, {
        "center_row": center_row,
        "outer_row": outer_row,
        "matrix": A,
        "rhs_vector": v,
        "c": c,
        "inv": inv,
        "Q": Q,
        "pde_rhs": pde_rhs,
    }


def evaluate_arm(state, alpha: np.ndarray, row_info: dict, target: float,
                 vacuum, baseline_alpha: np.ndarray) -> tuple[dict, dict[str, np.ndarray]]:
    grid = state.grid
    geom = state.geometry.copy()
    geom.alpha = np.asarray(alpha, dtype=float).copy()
    r = np.asarray(grid.centers, dtype=float)
    h = float(grid.dr)
    c = row_info["c"]
    inv = row_info["inv"]
    Q = row_info["Q"]
    A = row_info["matrix"]
    v = row_info["rhs_vector"]

    # The solve's second-order PDE residual, independently evaluated on cell centers.
    d1 = np.zeros_like(alpha)
    d2 = np.zeros_like(alpha)
    d1[1:-1] = (alpha[2:] - alpha[:-2]) / (2.0 * h)
    d2[1:-1] = (alpha[2:] - 2.0 * alpha[1:-1] + alpha[:-2]) / h**2
    interior_pde = inv * d2 - inv * 0.0 + inv * c * d1 - Q * alpha + float(target)

    if row_info["center_row"] == "current_equal_adjacent":
        center_condition = float(alpha[0] - alpha[1])
        center_pde = None
        center_vendor_pde = None
    elif row_info["center_row"] == "second_order_even_ghost_pde":
        center_condition = None
        center_lap = (alpha[1] - alpha[0]) / h**2 + c[0] * (alpha[1] - alpha[0]) / (2.0*h)
        center_pde = float(inv[0] * center_lap - Q[0]*alpha[0] + target)
        center_vendor_pde = None
    else:
        center_condition = None
        center_pde = None
        d1row = _operator_row(grid, 1, 0)
        d2row = _operator_row(grid, 2, 0)
        center_vendor_pde = float(inv[0] * (d2row @ alpha + c[0]*(d1row @ alpha)) - Q[0]*alpha[0] + target)

    if row_info["outer_row"] == "last_cell_center_alpha_one":
        outer_condition = float(alpha[-1] - 1.0)
        outer_face_alpha = None
        outer_pde = None
        outer_vendor_k_pde = None
    else:
        ghost = 2.0 - alpha[-1]
        outer_lap = (alpha[-2] - 2.0*alpha[-1] + ghost)/h**2 + c[-1]*(ghost-alpha[-2])/(2.0*h)
        outer_condition = None
        outer_face_alpha = float(0.5 * (alpha[-1] + ghost))
        outer_pde = float(inv[-1]*outer_lap - Q[-1]*alpha[-1] + target)
        outer_vendor_k_pde = None

    l2 = vacuum.primary_l2_rhs(grid, geom)
    l3_vac = vacuum.primary_l3_rhs(grid, geom)
    l3 = adapter.primary_l3_with_matter(grid, geom, state.scalars, state.matter)
    k_rhs = np.asarray(l2["K"], dtype=float) + np.asarray(l3["K"], dtype=float)
    residual = k_rhs - float(target)
    # Production's fourth-order lapse operator compared against this arm's solve PDE.
    alpha_d1 = grid.cell_derivative_fourth(geom.alpha, parity=1)
    alpha_d2 = grid.cell_second_derivative_fourth(geom.alpha, parity=1)
    prod_lapse_minus_solve = np.asarray(l2["K"], dtype=float) + inv * (
        d2 + c*d1
    )
    refreshed_target = float(target_kdot(grid, geom, state.scalars, state.matter))

    radiation = state.matter.radiation
    E = np.asarray(radiation.energy_t, dtype=float)
    Ur = np.asarray(radiation.momentum_r, dtype=float)
    cone_margin = E - np.sqrt(np.asarray(geom.X)**2 / np.asarray(geom.a))*np.abs(Ur)
    with np.errstate(divide="ignore", invalid="ignore"):
        cone_ratio = np.sqrt(np.asarray(geom.X)**2 / np.asarray(geom.a))*np.abs(Ur)/E
    stats = {
        "name": row_info["center_row"] + " + " + row_info["outer_row"],
        "center_row": row_info["center_row"],
        "outer_row": row_info["outer_row"],
        "positive_finite_alpha": bool(np.all(np.isfinite(alpha)) and np.min(alpha) > 0.0),
        "alpha_min": float(np.min(alpha)),
        "alpha_max": float(np.max(alpha)),
        "alpha_max_abs_from_baseline": float(np.max(np.abs(alpha - baseline_alpha))),
        "alpha_center_minus_adjacent": float(alpha[0] - alpha[1]),
        "alpha_even_series": even_series(alpha, r),
        "fixed_pre_solve_target": float(target),
        "target_recomputed_on_candidate_alpha": refreshed_target,
        "refreshed_target_minus_fixed": float(refreshed_target - target),
        "solver_row_closures": {
            "max_abs_matrix_residual": float(np.max(np.abs(A @ alpha - v))),
            "center_condition_residual": center_condition,
            "center_second_order_pde_residual": center_pde,
            "center_vendor_fourth_order_pde_residual": center_vendor_pde,
            "outer_condition_residual": outer_condition,
            "outer_face_alpha_reconstructed": outer_face_alpha,
            "outer_second_order_pde_residual": outer_pde,
            "interior_second_order_pde_max_abs": float(np.max(np.abs(interior_pde[1:-1]))),
        },
        "production_K_rhs_minus_target": {
            "full_domain": profile_stats(residual),
            "strict_interior_i3_to_i78": profile_stats(residual, 3, grid.n-1),
            "center_cell": float(residual[0]),
            "outer_cell": float(residual[-1]),
            "proper_volume_weighted_mean_rhs": weighted_mean(k_rhs, proper_volume_weights(grid, geom)),
            "proper_volume_weighted_mean_residual": weighted_mean(residual, proper_volume_weights(grid, geom)),
        },
        "production_lapse_operator_minus_solve_operator": {
            "strict_interior_i3_to_i78": profile_stats(prod_lapse_minus_solve, 3, grid.n-1),
        },
        "radiation_cone": {
            "min_margin": float(np.min(cone_margin)),
            "min_margin_cell": int(np.argmin(cone_margin)),
            "max_ratio": float(np.nanmax(cone_ratio)),
            "max_ratio_cell": int(np.nanargmax(cone_ratio)),
            "margin_digest": array_digest(cone_margin),
        },
        "outer_cell_terms": {
            "r": float(r[-1]),
            "alpha": float(alpha[-1]),
            "vendor_l2_K": float(l2["K"][-1]),
            "vacuum_l3_K": float(l3_vac["K"][-1]),
            "matter_l3_K_increment": float(l3["K"][-1] - l3_vac["K"][-1]),
            "K_ko": float(vacuum.ko_dissipation(grid, geom.K, 1)[-1]),
            "Q_alpha": float(Q[-1] * alpha[-1]),
            "production_K_rhs_minus_target": float(residual[-1]),
        },
    }
    profiles = {
        "alpha": np.asarray(alpha, dtype=float),
        "production_K_rhs": k_rhs,
        "production_K_rhs_minus_target": residual,
        "interior_pde_residual": interior_pde,
        "production_lapse_operator_minus_solve_operator": prod_lapse_minus_solve,
        "radiation_cone_margin": cone_margin,
        "radiation_cone_ratio": cone_ratio,
    }
    return stats, profiles


def run_audit(input_root: Path, output_dir: Path) -> dict:
    artifact_dir = find_artifact_dir(input_root)
    verified = verify_artifact(artifact_dir)
    kernel, state = load_terminal_state(artifact_dir, verified["capture"])
    if kernel is None:
        raise RuntimeError("kernel object unavailable")
    if int(state.grid.n) != 80 or abs(float(state.t) - 47.15300258518951) > 2e-14:
        raise RuntimeError("unexpected resolution/time for frozen-state A/B")
    before_digest = physical_state_digest(state)
    grid = state.grid
    weights = proper_volume_weights(grid, state.geometry)
    _, vacuum, _ = adapter.vendor_modules()
    original_alpha = np.asarray(state.geometry.alpha, dtype=float).copy()
    target = float(target_kdot(grid, state.geometry, state.scalars, state.matter))
    independent_rhs, _, _, _ = __import__(
        "tools.audit_cmc_lapse_fixed_point", fromlist=["full_k_rhs"]
    ).full_k_rhs(grid, state.geometry, state.scalars, state.matter, vacuum)
    initial_mean = weighted_mean(independent_rhs, weights)
    target_roundoff = abs(target - initial_mean)
    if target_roundoff > 128 * np.finfo(float).eps * max(abs(target), abs(initial_mean), np.finfo(float).tiny):
        raise RuntimeError("fixed target differs from the independent incoming K-RHS mean")

    baseline_alpha, baseline_target = solve_cmc_lapse(
        grid, state.geometry.copy(), state.scalars, state.matter
    )
    if abs(float(baseline_target) - target) > 128*np.finfo(float).eps*max(abs(target), np.finfo(float).tiny):
        raise RuntimeError("production solver target differs from the frozen target")
    baseline_alpha = np.asarray(baseline_alpha, dtype=float).copy()

    definitions = [
        ("current_equal_adjacent", "last_cell_center_alpha_one"),
        ("second_order_even_ghost_pde", "last_cell_center_alpha_one"),
        ("vendor_fourth_order_even_pde", "last_cell_center_alpha_one"),
        ("current_equal_adjacent", "face_alpha_one_with_last_cell_pde"),
        ("second_order_even_ghost_pde", "face_alpha_one_with_last_cell_pde"),
        ("vendor_fourth_order_even_pde", "face_alpha_one_with_last_cell_pde"),
    ]
    arms = []
    profiles = {}
    cone_reference = None
    for center_row, outer_row in definitions:
        name = center_row + "__" + outer_row
        if center_row == "current_equal_adjacent" and outer_row == "last_cell_center_alpha_one":
            alpha = baseline_alpha.copy()
            # Make the diagnostic row matrices independently available for the residual closure.
            candidate, row_info = assemble_rows(
                grid, state.geometry, state.scalars, state.matter, target,
                center_row, outer_row
            )
            matrix_delta = float(np.max(np.abs(candidate - alpha)))
            if matrix_delta > 1e-8:
                raise RuntimeError(f"reassembled current CMC system disagrees with production solve: {matrix_delta}")
            alpha = candidate
        else:
            alpha, row_info = assemble_rows(
                grid, state.geometry, state.scalars, state.matter, target,
                center_row, outer_row
            )
        stats, arrays = evaluate_arm(state, alpha, row_info, target, vacuum, baseline_alpha)
        stats["arm_id"] = name
        if not stats["positive_finite_alpha"]:
            raise RuntimeError(f"arm produced invalid lapse: {name}")
        if cone_reference is None:
            cone_reference = stats["radiation_cone"]["margin_digest"]
        if stats["radiation_cone"]["margin_digest"] != cone_reference:
            raise RuntimeError("radiation cone margin changed despite frozen non-lapse geometry")
        for profile_name, values in arrays.items():
            profiles[f"{name}__{profile_name}"] = values
        arms.append(stats)

    after_digest = physical_state_digest(state)
    if after_digest != before_digest:
        raise RuntimeError("counterfactual altered the frozen physical state")
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles_path = output_dir / "cmc_row_counterfactual_profiles.npz"
    np.savez_compressed(profiles_path, **profiles)
    report = {
        "schema": "cmc_lapse_terminal_row_counterfactual_v1",
        "status": "completed_diagnostic_only",
        "diagnostic_only": True,
        "input_artifact": {
            "run_id": EXPECTED_INPUT_RUN,
            "head": EXPECTED_INPUT_HEAD,
            "artifact_checksums": verified["checksums"],
        },
        "source_pins": {
            "main": EXPECTED_MAIN,
            "vendor": EXPECTED_VENDOR,
        },
        "state": {
            "resolution": int(grid.n),
            "r_max": 80.0,
            "t": float(state.t),
            "tau": float(state.tau),
            "physical_digest_before": before_digest,
            "physical_digest_after": after_digest,
            "physical_state_unchanged": before_digest == after_digest,
        },
        "fixed_target": {
            "kdot": target,
            "independent_incoming_rhs_weighted_mean": initial_mean,
            "difference": target_roundoff,
            "held_constant_across_all_arms": True,
        },
        "radiation_cone_invariance": {
            "margin_digest_same_across_arms": True,
            "margin_digest": cone_reference,
            "reason": "only alpha solve rows change; a, X, and radiation conserved variables are frozen",
        },
        "arms": arms,
        "interpretation_limits": [
            "This is a no-evolution counterfactual on the archived N=80 terminal radiation-failure state.",
            "The target Kdot is held fixed from the incoming state so row choices can be compared without conflating target refresh.",
            "The face-Dirichlet arm imposes alpha=1 at r_max by ghost reflection and uses the last cell-centred PDE; it is a diagnostic boundary alternative, not a promoted physical boundary condition.",
            "Even-ghost and vendor-fourth-order centre rows are operator counterfactuals, not proof of continuum regularity.",
            "Since a, X and radiation conserved fields are fixed, this audit cannot change the radiation cone margin or itself admit the trajectory.",
        ],
    }
    (output_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    hashes = []
    for path in (profiles_path, output_dir / "report.json"):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        hashes.append(f"{digest}  {path.name}\n")
    (output_dir / "SHA256SUMS.txt").write_text("".join(hashes), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = run_audit(args.input_dir, args.output_dir)
    print(json.dumps({
        "status": report["status"],
        "t": report["state"]["t"],
        "physical_state_unchanged": report["state"]["physical_state_unchanged"],
        "fixed_target": report["fixed_target"],
        "arms": [{
            "arm_id": arm["arm_id"],
            "alpha_max_abs_from_baseline": arm["alpha_max_abs_from_baseline"],
            "interior_K_rhs_residual_max": arm["production_K_rhs_minus_target"]["strict_interior_i3_to_i78"]["max_abs"],
            "center_K_rhs_residual": arm["production_K_rhs_minus_target"]["center_cell"],
            "outer_K_rhs_residual": arm["production_K_rhs_minus_target"]["outer_cell"],
            "center_pde_residual": arm["solver_row_closures"]["center_second_order_pde_residual"],
            "center_vendor_pde_residual": arm["solver_row_closures"]["center_vendor_fourth_order_pde_residual"],
            "outer_pde_residual": arm["solver_row_closures"]["outer_second_order_pde_residual"],
            "alpha_min": arm["alpha_min"],
            "alpha_max": arm["alpha_max"],
        } for arm in report["arms"]],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
