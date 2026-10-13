import numpy as np
import pytest

from tools.audit_cmc_lapse_joint_target import (
    solve_joint_target_lapse,
    solve_joint_target_lapse_equilibrated,
)


def test_joint_target_lapse_closes_interior_and_full_slice_weighted_mean():
    n = 8
    alpha_expected = np.array([1.0, 1.0, 0.97, 0.94, 0.92, 0.95, 0.98, 1.0])
    target_expected = 0.125
    matrix = np.diag(np.linspace(0.2, 0.9, n))
    matrix[2, 1] = -0.03
    matrix[3, 2] = 0.04
    matrix[5, 4] = -0.02
    intercept = target_expected - matrix @ alpha_expected
    weights = np.linspace(0.5, 2.0, n)

    alpha, target, meta = solve_joint_target_lapse(matrix, intercept, weights)
    rhs = matrix @ alpha + intercept

    assert meta["augmented_roundoff_closure_passed"]
    assert np.max(np.abs(alpha - alpha_expected)) < 2e-12
    assert abs(target - target_expected) < 2e-12
    assert abs(alpha[0] - alpha[1]) < 2e-12
    assert abs(alpha[-1] - 1.0) < 2e-12
    assert np.max(np.abs(rhs[1:-1] - target)) < 2e-12
    assert abs(np.average(rhs, weights=weights) - target) < 2e-12


def test_joint_target_lapse_rejects_invalid_weight_measure():
    n = 6
    matrix = np.eye(n)
    intercept = np.zeros(n)
    with pytest.raises(ValueError, match="positive finite total proper-volume weight"):
        solve_joint_target_lapse(matrix, intercept, np.zeros(n))



def test_equilibrated_joint_target_lapse_solves_same_equations():
    n = 8
    alpha_expected = np.array([1.0, 1.0, 0.97, 0.94, 0.92, 0.95, 0.98, 1.0])
    target_expected = -0.004
    matrix = np.diag(np.linspace(0.2, 0.9, n))
    matrix[2, 1] = -0.03
    matrix[3, 2] = 0.04
    matrix[5, 4] = -0.02
    # Deliberately mix the equation row scales; equilibration must preserve
    # the original equations, target, centre/outer conditions, and mean target.
    row_scale = np.logspace(-3, 3, n)
    matrix = row_scale[:, None] * matrix
    intercept = target_expected - matrix @ alpha_expected
    weights = np.linspace(0.5, 2.0, n)

    alpha, target, meta = solve_joint_target_lapse_equilibrated(
        matrix, intercept, weights
    )
    rhs = matrix @ alpha + intercept

    assert meta["original_equations_close"]
    assert meta["equilibrated_system_roundoff_closed"]
    assert np.max(np.abs(alpha - alpha_expected)) < 2e-8
    assert abs(target - target_expected) < 2e-8
    assert abs(alpha[0] - alpha[1]) < 2e-8
    assert abs(alpha[-1] - 1.0) < 2e-8
    assert np.max(np.abs(rhs[1:-1] - target)) < 2e-8
    assert abs(np.average(rhs, weights=weights) - target) < 2e-8
