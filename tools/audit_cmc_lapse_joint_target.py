#!/usr/bin/env python3
"""Frozen-state joint solve of the production discrete K RHS and CMC target.

This is an operator counterfactual only. It never advances matter or geometry.
The production fourth-order lapse operator is reused, including its exact
cell-centred/one-sided derivative rows. The centre regularity and outer lapse
normalization rows are explicit. A final weighted-mean equation defines the
single target consistently from the full production K RHS.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.audit_cmc_lapse_fixed_point import (
    EXPECTED_CAPTURE_HEAD,
    EXPECTED_CAPTURE_RUN_ID,
    EXPECTED_VENDOR_COMMIT,
    find_artifact_dir,
    full_k_rhs,
    load_terminal_state,
    physical_state_digest,
    profile_stats,
    proper_volume_weights,
    sha256,
    verify_artifact,
    weighted_mean,
)
from engine.cmc_gauge import solve_cmc_lapse
from engine import v55_pirk_adapter as adapter
from engine.v55_matter import total_matter_projection


BASE_MAIN_COMMIT = "c01f2467347e8aebd57d220374fd7bc4e93045bc"
PARENT_AUDIT_HEAD = "2400937f913a98723d93afe95a08c3293241d918"
TERMINAL_TIME = 47.15300258518951
GEOMETRY_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")
RADIATION_FIELDS = ("energy_t", "momentum_r")


def array_sha256(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    return hashlib.sha256(array.tobytes()).hexdigest()


def operator_matrices(grid) -> tuple[np.ndarray, np.ndarray]:
    """Build exact discrete derivative matrices from the pinned grid methods."""
    n = int(grid.n)
    identity = np.eye(n, dtype=float)
    first = np.column_stack([
        grid.cell_derivative_fourth(identity[:, j], parity=1)
        for j in range(n)
    ])
    second = np.column_stack([
        grid.cell_second_derivative_fourth(identity[:, j], parity=1)
        for j in range(n)
    ])
    return first, second


def production_affine_operator(state, vacuum) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """Return A, c with production K_RHS(alpha) = A @ alpha + c.

    The alpha-dependent terms are the pinned vendor fourth-order lapse
    Laplacian and the production L3 coefficient multiplying alpha. The
    intercept is recovered from a direct production RHS evaluation on the
    reference lapse in this frozen diagnostic state, retaining alpha-independent
    shift/KO terms.
    """
    grid, geometry = state.grid, state.geometry
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    inv = X * X / a
    if np.any(inv <= 0.0) or not np.all(np.isfinite(inv)):
        raise RuntimeError("frozen geometry has invalid inverse radial metric")

    ap = grid.cell_derivative_fourth(a, parity=1)
    bp = grid.cell_derivative_fourth(b, parity=1)
    Xp = grid.cell_derivative_fourth(X, parity=1)
    c_radial = -0.5 * ap / a + bp / b - Xp / X + 2.0 / r

    d1, d2 = operator_matrices(grid)
    lapse_laplacian = inv[:, None] * (d2 + c_radial[:, None] * d1)

    total = total_matter_projection(grid, geometry, state.scalars, state.matter)
    alpha_coefficient = (
        1.5 * np.asarray(geometry.Aa, dtype=float) ** 2
        + np.asarray(geometry.K, dtype=float) ** 2 / 3.0
        + 4.0 * math.pi * (
            np.asarray(total["rho"], dtype=float)
            + np.asarray(total["pr"], dtype=float)
            + 2.0 * np.asarray(total["pt"], dtype=float)
        )
    )
    matrix = -lapse_laplacian + np.diag(alpha_coefficient)

    reference_alpha = np.asarray(geometry.alpha, dtype=float).copy()
    reference_rhs = full_k_rhs(
        grid, geometry, state.scalars, state.matter, vacuum
    )[0]
    intercept = reference_rhs - matrix @ reference_alpha
    meta = {
        "max_abs_alpha_coefficient": float(np.max(np.abs(alpha_coefficient))),
        "alpha_coefficient_rms": float(np.sqrt(np.mean(alpha_coefficient**2))),
        "reference_alpha_rhs_max_abs": float(np.max(np.abs(reference_rhs))),
        "reference_lapse_operator_reconstruction_max_abs": float(
            np.max(np.abs(matrix @ reference_alpha + intercept - reference_rhs))
        ),
    }
    return matrix, intercept, reference_rhs, meta


def solve_joint_target_lapse(
    matrix: np.ndarray,
    intercept: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, float, dict]:
    """Solve centre row + production K-RHS interior + outer row + mean target."""
    matrix = np.asarray(matrix, dtype=float)
    intercept = np.asarray(intercept, dtype=float)
    weights = np.asarray(weights, dtype=float)
    n = int(weights.size)
    if matrix.shape != (n, n) or intercept.shape != (n,):
        raise ValueError("operator/offset dimensions do not match weights")
    if n < 5 or not np.all(np.isfinite(matrix)):
        raise ValueError("joint solve requires a finite operator with at least five cells")
    if not np.all(np.isfinite(intercept)) or not np.all(np.isfinite(weights)):
        raise ValueError("joint solve inputs must be finite")
    weight_sum = float(np.sum(weights))
    if not np.isfinite(weight_sum) or weight_sum <= 0.0:
        raise ValueError("joint solve requires positive finite total proper-volume weight")
    normalized = weights / weight_sum

    # Unknown vector is [alpha_0, ..., alpha_(N-1), kdot].
    # Rows: alpha[0] = alpha[1]; production K_RHS = kdot at i=1..N-2;
    # alpha[-1] = 1; and the full-slice proper-volume mean K_RHS = kdot.
    augmented = np.zeros((n + 1, n + 1), dtype=float)
    rhs = np.zeros(n + 1, dtype=float)
    augmented[0, 0] = 1.0
    augmented[0, 1] = -1.0

    for i in range(1, n - 1):
        augmented[i, :n] = matrix[i, :]
        augmented[i, n] = -1.0
        rhs[i] = -intercept[i]

    augmented[n - 1, n - 1] = 1.0
    rhs[n - 1] = 1.0

    augmented[n, :n] = normalized @ matrix
    augmented[n, n] = -1.0
    rhs[n] = -float(normalized @ intercept)

    solution = np.linalg.solve(augmented, rhs)
    alpha = np.asarray(solution[:n], dtype=float)
    target = float(solution[n])
    algebraic_residual = augmented @ solution - rhs
    scale = (
        np.abs(augmented) @ np.abs(solution) + np.abs(rhs)
    )
    bound = 1024.0 * np.finfo(float).eps * np.maximum(
        scale, np.finfo(float).tiny
    )
    ratio = np.abs(algebraic_residual) / bound
    return alpha, target, {
        "augmented_matrix_shape": [int(n + 1), int(n + 1)],
        "condition_number_2": float(np.linalg.cond(augmented)),
        "max_abs_augmented_equation_residual": float(
            np.max(np.abs(algebraic_residual))
        ),
        "max_operand_scaled_roundoff_ratio": float(np.max(ratio)),
        "augmented_roundoff_closure_passed": bool(np.all(np.abs(algebraic_residual) <= bound)),
    }



def solve_joint_target_lapse_equilibrated(
    matrix: np.ndarray,
    intercept: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, float, dict]:
    """Solve the same augmented equations after exact row/column equilibration.

    Equilibration changes only the floating-point representation of the linear
    system, not its equations: R @ A @ C @ y = R @ b, with x = C @ y.
    Acceptance checks are still evaluated on the original unscaled equations.
    """
    matrix = np.asarray(matrix, dtype=float)
    intercept = np.asarray(intercept, dtype=float)
    weights = np.asarray(weights, dtype=float)
    n = int(weights.size)
    if matrix.shape != (n, n) or intercept.shape != (n,):
        raise ValueError("operator/offset dimensions do not match weights")
    if n < 5 or not np.all(np.isfinite(matrix)):
        raise ValueError("joint solve requires a finite operator with at least five cells")
    if not np.all(np.isfinite(intercept)) or not np.all(np.isfinite(weights)):
        raise ValueError("joint solve inputs must be finite")
    weight_sum = float(np.sum(weights))
    if not np.isfinite(weight_sum) or weight_sum <= 0.0:
        raise ValueError("joint solve requires positive finite total proper-volume weight")
    normalized = weights / weight_sum

    augmented = np.zeros((n + 1, n + 1), dtype=float)
    rhs = np.zeros(n + 1, dtype=float)
    augmented[0, 0] = 1.0
    augmented[0, 1] = -1.0
    for i in range(1, n - 1):
        augmented[i, :n] = matrix[i, :]
        augmented[i, n] = -1.0
        rhs[i] = -intercept[i]
    augmented[n - 1, n - 1] = 1.0
    rhs[n - 1] = 1.0
    augmented[n, :n] = normalized @ matrix
    augmented[n, n] = -1.0
    rhs[n] = -float(normalized @ intercept)

    row_max = np.max(np.abs(augmented), axis=1)
    if np.any(~np.isfinite(row_max)) or np.any(row_max <= 0.0):
        raise ValueError("joint solve has a zero or invalid augmented row")
    row_scale = 1.0 / row_max
    row_balanced = row_scale[:, None] * augmented
    rhs_balanced = row_scale * rhs

    col_max = np.max(np.abs(row_balanced), axis=0)
    if np.any(~np.isfinite(col_max)) or np.any(col_max <= 0.0):
        raise ValueError("joint solve has a zero or invalid equilibrated column")
    col_scale = 1.0 / col_max
    balanced = row_balanced * col_scale[None, :]

    y = np.linalg.solve(balanced, rhs_balanced)
    solution = col_scale * y
    alpha = np.asarray(solution[:n], dtype=float)
    target = float(solution[n])

    # Assess both systems, but base equation admission on the original rows.
    raw_residual = augmented @ solution - rhs
    raw_scale = np.abs(augmented) @ np.abs(solution) + np.abs(rhs)
    raw_bound = 1024.0 * np.finfo(float).eps * np.maximum(
        raw_scale, np.finfo(float).tiny
    )
    raw_ratio = np.abs(raw_residual) / raw_bound

    balanced_residual = balanced @ y - rhs_balanced
    balanced_scale = np.abs(balanced) @ np.abs(y) + np.abs(rhs_balanced)
    balanced_bound = 1024.0 * np.finfo(float).eps * np.maximum(
        balanced_scale, np.finfo(float).tiny
    )
    balanced_ratio = np.abs(balanced_residual) / balanced_bound
    return alpha, target, {
        "augmented_matrix_shape": [int(n + 1), int(n + 1)],
        "raw_condition_number_2": float(np.linalg.cond(augmented)),
        "equilibrated_condition_number_2": float(np.linalg.cond(balanced)),
        "row_scale_min": float(np.min(row_scale)),
        "row_scale_max": float(np.max(row_scale)),
        "column_scale_min": float(np.min(col_scale)),
        "column_scale_max": float(np.max(col_scale)),
        "max_abs_original_equation_residual": float(np.max(np.abs(raw_residual))),
        "max_original_operand_scaled_roundoff_ratio": float(np.max(raw_ratio)),
        "original_equations_close": bool(np.all(np.abs(raw_residual) <= raw_bound)),
        "max_abs_equilibrated_equation_residual": float(
            np.max(np.abs(balanced_residual))
        ),
        "max_equilibrated_operand_scaled_roundoff_ratio": float(
            np.max(balanced_ratio)
        ),
        "equilibrated_system_roundoff_closed": bool(
            np.all(np.abs(balanced_residual) <= balanced_bound)
        ),
    }


def solve_joint_target_lapse_sparse(
    matrix: np.ndarray,
    intercept: np.ndarray,
    weights: np.ndarray,
    *,
    equilibrate: bool = False,
) -> tuple[np.ndarray, float, dict]:
    """Solve identical augmented equations using sparse row-local matrix action."""
    matrix = np.asarray(matrix, dtype=float)
    intercept = np.asarray(intercept, dtype=float)
    weights = np.asarray(weights, dtype=float)
    n = int(weights.size)
    if matrix.shape != (n, n) or intercept.shape != (n,):
        raise ValueError("sparse joint solve dimensions do not match")
    if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(intercept)):
        raise ValueError("sparse joint solve inputs must be finite")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise ValueError("sparse joint solve weights must be finite and nonnegative")
    total_weight = float(np.sum(weights))
    if not np.isfinite(total_weight) or total_weight <= 0.0:
        raise ValueError("sparse joint solve requires positive finite total weight")
    normalized = weights / total_weight

    # Same augmented equations; CSR storage/matvec avoids the dense BLAS dot
    # product for each row-local K stencil. The equations and boundary rows
    # are unchanged.
    augmented = np.zeros((n + 1, n + 1), dtype=float)
    rhs = np.zeros(n + 1, dtype=float)
    augmented[0, 0] = 1.0
    augmented[0, 1] = -1.0
    for i in range(1, n - 1):
        augmented[i, :n] = matrix[i, :]
        augmented[i, n] = -1.0
        rhs[i] = -intercept[i]
    augmented[n - 1, n - 1] = 1.0
    rhs[n - 1] = 1.0
    augmented[n, :n] = normalized @ matrix
    augmented[n, n] = -1.0
    rhs[n] = -float(normalized @ intercept)

    sparse = csr_matrix(augmented)
    row_scale = np.ones(n + 1, dtype=float)
    col_scale = np.ones(n + 1, dtype=float)
    solve_matrix = sparse
    solve_rhs = rhs.copy()
    if equilibrate:
        row_max = np.asarray(np.abs(sparse).max(axis=1).toarray()).reshape(-1)
        if np.any(~np.isfinite(row_max)) or np.any(row_max <= 0.0):
            raise ValueError("sparse augmented system has zero or invalid row")
        row_scale = 1.0 / row_max
        row_balanced = csr_matrix(row_scale[:, None] * sparse.toarray())
        solve_rhs = row_scale * rhs
        col_max = np.asarray(np.abs(row_balanced).max(axis=0).toarray()).reshape(-1)
        if np.any(~np.isfinite(col_max)) or np.any(col_max <= 0.0):
            raise ValueError("sparse augmented system has zero or invalid column")
        col_scale = 1.0 / col_max
        solve_matrix = csr_matrix(row_balanced.toarray() * col_scale[None, :])

    scaled_solution = np.asarray(spsolve(solve_matrix, solve_rhs), dtype=float)
    solution = col_scale * scaled_solution if equilibrate else scaled_solution
    alpha = np.asarray(solution[:n], dtype=float)
    target = float(solution[n])

    raw_residual = augmented @ solution - rhs
    raw_scale = np.abs(augmented) @ np.abs(solution) + np.abs(rhs)
    raw_bound = 1024.0 * np.finfo(float).eps * np.maximum(
        raw_scale, np.finfo(float).tiny
    )
    raw_ratio = np.abs(raw_residual) / raw_bound

    csr_residual = np.asarray(sparse @ solution, dtype=float) - rhs
    csr_scale = np.asarray(np.abs(sparse) @ np.abs(solution)).reshape(-1) + np.abs(rhs)
    csr_bound = 1024.0 * np.finfo(float).eps * np.maximum(
        csr_scale, np.finfo(float).tiny
    )
    csr_ratio = np.abs(csr_residual) / csr_bound
    return alpha, target, {
        "method": "scipy sparse CSR solve of the same augmented coefficient matrix",
        "equilibrated": bool(equilibrate),
        "condition_number_2_dense_diagnostic": float(np.linalg.cond(augmented)),
        "nnz": int(sparse.nnz),
        "density": float(sparse.nnz / sparse.shape[0] / sparse.shape[1]),
        "row_scale_min": float(np.min(row_scale)),
        "row_scale_max": float(np.max(row_scale)),
        "column_scale_min": float(np.min(col_scale)),
        "column_scale_max": float(np.max(col_scale)),
        "max_abs_original_equation_residual": float(np.max(np.abs(raw_residual))),
        "max_original_operand_scaled_roundoff_ratio": float(np.max(raw_ratio)),
        "original_equations_close": bool(np.all(np.abs(raw_residual) <= raw_bound)),
        "max_abs_csr_matvec_residual": float(np.max(np.abs(csr_residual))),
        "max_csr_operand_scaled_roundoff_ratio": float(np.max(csr_ratio)),
        "csr_system_roundoff_closed": bool(np.all(np.abs(csr_residual) <= csr_bound)),
    }


def run_audit(input_root: Path, output_dir: Path) -> dict:
    artifact_dir = find_artifact_dir(input_root)
    verified = verify_artifact(artifact_dir)
    kernel, state = load_terminal_state(artifact_dir, verified["capture"])
    if kernel is None:
        raise RuntimeError("kernel object unavailable")
    if int(state.grid.n) != 80:
        raise RuntimeError("this frozen-state discriminator requires N=80")
    if abs(float(state.t) - TERMINAL_TIME) > 2e-14:
        raise RuntimeError("input is not the admitted archived terminal radiation-failure state")

    original_digest = physical_state_digest(state)
    original_alpha = np.asarray(state.geometry.alpha, dtype=float).copy()
    original_matter = {
        name: {
            field: np.asarray(getattr(getattr(state.matter, name), field), dtype=float).copy()
            for field in RADIATION_FIELDS
        }
        for name in ("radiation",)
    }
    original_cone = (
        np.asarray(state.matter.radiation.energy_t, dtype=float)
        - np.sqrt(
            np.asarray(state.geometry.X, dtype=float) ** 2
            / np.asarray(state.geometry.a, dtype=float)
        ) * np.abs(np.asarray(state.matter.radiation.momentum_r, dtype=float))
    )
    original_cone_digest = array_sha256(original_cone)

    _, vacuum, _ = adapter.vendor_modules()
    weights = proper_volume_weights(state.grid, state.geometry)
    accepted_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    accepted_target = weighted_mean(accepted_rhs, weights)
    current_solve_alpha, current_solve_target = solve_cmc_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )
    state.geometry.alpha = np.asarray(current_solve_alpha, dtype=float).copy()
    current_solve_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    current_solve_digest = array_sha256(current_solve_rhs)

    matrix, intercept, reference_recomputed_rhs, operator_meta = production_affine_operator(
        state, vacuum
    )
    # Verify the affine map at a second, positive lapse profile independently
    # of the accepted lapse used to recover its intercept.
    probe_alpha = np.asarray(current_solve_alpha, dtype=float) + (
        0.013 * np.cos(np.asarray(state.grid.centers, dtype=float) / 80.0)
    )
    if float(np.min(probe_alpha)) <= 0.0:
        raise RuntimeError("affine-operator probe unexpectedly has non-positive lapse")
    state.geometry.alpha = probe_alpha.copy()
    direct_probe_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    predicted_probe_rhs = matrix @ probe_alpha + intercept
    probe_scale = (
        np.abs(direct_probe_rhs)
        + np.abs(matrix) @ np.abs(probe_alpha)
        + np.abs(intercept)
    )
    probe_bound = 2048.0 * np.finfo(float).eps * np.maximum(
        probe_scale, np.finfo(float).tiny
    )
    probe_error = direct_probe_rhs - predicted_probe_rhs
    affine_probe_passed = bool(np.all(np.abs(probe_error) <= probe_bound))

    joint_alpha, joint_target, solve_meta = solve_joint_target_lapse(
        matrix, intercept, weights
    )
    state.geometry.alpha = joint_alpha.copy()
    direct_joint_rhs = full_k_rhs(
        state.grid, state.geometry, state.scalars, state.matter, vacuum
    )[0]
    predicted_joint_rhs = matrix @ joint_alpha + intercept
    production_reconstruction_error = direct_joint_rhs - predicted_joint_rhs
    reconstruction_scale = (
        np.abs(direct_joint_rhs)
        + np.abs(matrix) @ np.abs(joint_alpha)
        + np.abs(intercept)
    )
    reconstruction_bound = 2048.0 * np.finfo(float).eps * np.maximum(
        reconstruction_scale, np.finfo(float).tiny
    )
    direct_reconstruction_passed = bool(
        np.all(np.abs(production_reconstruction_error) <= reconstruction_bound)
    )

    joint_residual = direct_joint_rhs - joint_target
    mean_target = weighted_mean(direct_joint_rhs, weights)
    interior = joint_residual[1:-1]
    mean_closure_scale = (
        float(np.sum(weights * np.abs(direct_joint_rhs)))
        + abs(joint_target) * float(np.sum(weights))
    )
    mean_closure_bound = 4096.0 * np.finfo(float).eps * max(
        mean_closure_scale / float(np.sum(weights)), np.finfo(float).tiny
    )
    mean_closure = mean_target - joint_target
    positive_finite = bool(
        np.all(np.isfinite(joint_alpha)) and float(np.min(joint_alpha)) > 0.0
    )
    interior_scale = (
        np.abs(direct_joint_rhs[1:-1])
        + abs(joint_target)
        + np.abs(predicted_joint_rhs[1:-1])
    )
    interior_bound = 4096.0 * np.finfo(float).eps * np.maximum(
        interior_scale, np.finfo(float).tiny
    )
    interior_closure_passed = bool(
        np.all(np.abs(interior) <= interior_bound)
    )

    state.geometry.alpha = original_alpha.copy()
    after_digest = physical_state_digest(state)
    after_cone = (
        np.asarray(state.matter.radiation.energy_t, dtype=float)
        - np.sqrt(
            np.asarray(state.geometry.X, dtype=float) ** 2
            / np.asarray(state.geometry.a, dtype=float)
        ) * np.abs(np.asarray(state.matter.radiation.momentum_r, dtype=float))
    )
    after_cone_digest = array_sha256(after_cone)
    radiation_profiles_unchanged = all(
        np.array_equal(
            original_matter["radiation"][field],
            np.asarray(getattr(state.matter.radiation, field), dtype=float),
        )
        for field in RADIATION_FIELDS
    )
    state_unchanged = (
        original_digest == after_digest
        and original_cone_digest == after_cone_digest
        and radiation_profiles_unchanged
    )

    report = {
        "schema": "cmc_joint_discrete_target_lapse_counterfactual_v1",
        "status": "completed_diagnostic_only",
        "diagnostic_only": True,
        "no_evolution": True,
        "source_pins": {
            "main_base": BASE_MAIN_COMMIT,
            "parent_cmc_audit_head": PARENT_AUDIT_HEAD,
            "production_cmc_blob": "b1e6b3e95ff37702509e1b38ec6401d10e2482c3",
            "vendor": EXPECTED_VENDOR_COMMIT,
            "input_replay_head": EXPECTED_CAPTURE_HEAD,
            "input_replay_run_id": EXPECTED_CAPTURE_RUN_ID,
            "input_manifest_sha256": verified["checksums"],
        },
        "state": {
            "resolution": int(state.grid.n),
            "r_max": float(state.grid.centers[-1] + 0.5 * state.grid.dr),
            "t": float(state.t),
            "tau": float(state.tau),
            "physical_digest_before": original_digest,
            "physical_digest_after": after_digest,
            "non_lapse_physical_state_unchanged": original_digest == after_digest,
            "radiation_cone_digest_before": original_cone_digest,
            "radiation_cone_digest_after": after_cone_digest,
            "radiation_cone_digest_unchanged": original_cone_digest == after_cone_digest,
            "radiation_profiles_unchanged": radiation_profiles_unchanged,
            "all_frozen_state_checks_passed": state_unchanged,
        },
        "operator": {
            "definition": "K_RHS(alpha) = A @ alpha + intercept, where A is reconstructed from the exact pinned fourth-order cell derivative operators plus the actual L3 alpha coefficient; the intercept is recovered from a direct production RHS on the saved accepted lapse.",
            "matrix_reconstruction": operator_meta,
            "affine_probe": {
                "independent_positive_lapse_probe": True,
                "max_abs_direct_minus_affine_rhs": float(np.max(np.abs(probe_error))),
                "max_operand_scaled_roundoff_ratio": float(
                    np.max(np.abs(probe_error) / probe_bound)
                ),
                "passed": affine_probe_passed,
            },
            "joint_solution_direct_rhs_reconstruction": {
                "max_abs_direct_minus_affine_rhs": float(
                    np.max(np.abs(production_reconstruction_error))
                ),
                "max_operand_scaled_roundoff_ratio": float(
                    np.max(np.abs(production_reconstruction_error) / reconstruction_bound)
                ),
                "passed": direct_reconstruction_passed,
            },
        },
        "baseline": {
            "accepted_state_weighted_target": float(accepted_target),
            "accepted_state_rhs_max_abs": float(np.max(np.abs(accepted_rhs))),
            "current_solver_target": float(current_solve_target),
            "current_solver_lapse_min": float(np.min(current_solve_alpha)),
            "current_solver_lapse_max": float(np.max(current_solve_alpha)),
            "current_solver_rhs_residual": {
                "center_cell": float(current_solve_rhs[0] - current_solve_target),
                "strict_interior_i1_to_i78_max_abs": float(
                    np.max(np.abs(current_solve_rhs[1:-1] - current_solve_target))
                ),
                "outer_cell": float(current_solve_rhs[-1] - current_solve_target),
                "full_slice_weighted_mean_residual": float(
                    weighted_mean(current_solve_rhs - current_solve_target, weights)
                ),
            },
            "current_solver_rhs_digest": current_solve_digest,
        },
        "joint_counterfactual": {
            "target_kdot": joint_target,
            "direct_production_weighted_mean_kdot": mean_target,
            "weighted_mean_minus_target": float(mean_closure),
            "weighted_mean_roundoff_bound": float(mean_closure_bound),
            "weighted_mean_closure_passed": bool(abs(mean_closure) <= mean_closure_bound),
            "positive_finite_lapse": positive_finite,
            "lapse_min": float(np.min(joint_alpha)),
            "lapse_max": float(np.max(joint_alpha)),
            "max_abs_lapse_change_from_current_solver": float(
                np.max(np.abs(joint_alpha - current_solve_alpha))
            ),
            "center_alpha_minus_adjacent": float(joint_alpha[0] - joint_alpha[1]),
            "center_production_rhs_minus_target": float(joint_residual[0]),
            "strict_interior_i1_to_i78_max_abs_rhs_residual": float(
                np.max(np.abs(interior))
            ),
            "strict_interior_i1_to_i78_rms_rhs_residual": float(
                np.sqrt(np.mean(interior**2))
            ),
            "outer_production_rhs_minus_target": float(joint_residual[-1]),
            "full_slice_production_rhs_minus_target": profile_stats(
                joint_residual, weights, 0, state.grid.n
            ),
            "solver": solve_meta,
            "interior_equations_close": interior_closure_passed,
            "interior_roundoff_max_ratio": float(
                np.max(np.abs(interior) / interior_bound)
            ),
            "center_row": "alpha[0] = alpha[1], matching current production centre condition",
            "outer_row": "alpha[-1] = 1, matching current production normalization row",
            "target_row": "proper-volume weighted mean of the full production K RHS equals target_kdot",
            "production_rows": "pinned fourth-order production K RHS equals target_kdot at cells i=1..N-2, including the actual one-sided outer-interior stencil row",
        },
        "interpretation_limits": [
            "The solve freezes every non-lapse geometry, scalar, matter, and time field; only the lapse is counterfactually replaced.",
            "The counterfactual is designed to close the production-discrete K equation on interior cells and the existing lapse boundary rows, with target fixed by the full-slice weighted mean.",
            "Centre and outer cells retain boundary-condition residuals; the target mean equation balances their weighted contributions and does not assert the boundary rows are evolution PDE rows.",
            "This diagnostic does not establish physical validity, continuum convergence, radiation admission, or authority to change the production gauge.",
        ],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    profiles_path = output_dir / "cmc_joint_target_profiles.npz"
    np.savez_compressed(
        profiles_path,
        alpha_accepted=original_alpha,
        alpha_current_solver=current_solve_alpha,
        alpha_joint=joint_alpha,
        operator_matrix=matrix,
        operator_intercept=intercept,
        accepted_rhs=accepted_rhs,
        current_solver_rhs=current_solve_rhs,
        joint_rhs_predicted=predicted_joint_rhs,
        joint_rhs_direct=direct_joint_rhs,
        joint_rhs_residual=joint_residual,
        affine_probe_alpha=probe_alpha,
        affine_probe_rhs_direct=direct_probe_rhs,
        affine_probe_rhs_predicted=predicted_probe_rhs,
        affine_probe_residual=probe_error,
        proper_volume_weights=weights,
    )
    report_path = output_dir / "report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (output_dir / "SHA256SUMS.txt").write_text(
        "".join(
            f"{sha256(path)}  {path.name}\n"
            for path in (profiles_path, report_path)
        ),
        encoding="utf-8",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = run_audit(args.input_dir, args.output_dir)
    summary = {
        "status": report["status"],
        "t": report["state"]["t"],
        "all_frozen_state_checks_passed": report["state"]["all_frozen_state_checks_passed"],
        "affine_probe_passed": report["operator"]["affine_probe"]["passed"],
        "direct_rhs_reconstruction_passed": report["operator"]["joint_solution_direct_rhs_reconstruction"]["passed"],
        "positive_finite_lapse": report["joint_counterfactual"]["positive_finite_lapse"],
        "target_kdot": report["joint_counterfactual"]["target_kdot"],
        "weighted_mean_minus_target": report["joint_counterfactual"]["weighted_mean_minus_target"],
        "interior_max_abs_rhs_residual": report["joint_counterfactual"]["strict_interior_i1_to_i78_max_abs_rhs_residual"],
        "center_rhs_residual": report["joint_counterfactual"]["center_production_rhs_minus_target"],
        "outer_rhs_residual": report["joint_counterfactual"]["outer_production_rhs_minus_target"],
        "artifact_output": str(args.output_dir),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    passed = all((
        report["state"]["all_frozen_state_checks_passed"],
        report["operator"]["affine_probe"]["passed"],
        report["operator"]["joint_solution_direct_rhs_reconstruction"]["passed"],
        report["joint_counterfactual"]["solver"]["augmented_roundoff_closure_passed"],
        report["joint_counterfactual"]["interior_equations_close"],
        report["joint_counterfactual"]["weighted_mean_closure_passed"],
        report["joint_counterfactual"]["positive_finite_lapse"],
    ))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
