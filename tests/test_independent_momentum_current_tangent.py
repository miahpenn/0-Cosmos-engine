import numpy as np
import pytest

from tools.independent_momentum_current_tangent import (
    KAPPA,
    compare_directional_epsilons,
    independent_momentum_residual,
    projected_current_directional_tangent,
    scalar_current_tangent,
)


def test_scalar_current_tangent_matches_centered_directional_difference():
    matrix = np.array([
        [0.0, 1.0, 0.0],
        [-0.5, 0.0, 0.5],
        [0.0, -1.0, 0.0],
    ])
    f = np.array([1.0, 2.0, 4.0])
    p = np.array([0.5, -1.0, 2.0])
    fdot = np.array([0.2, -0.3, 0.7])
    pdot = np.array([-0.4, 0.1, 0.6])
    d1 = lambda x: matrix @ x

    tangent = scalar_current_tangent(d1, f, p, fdot, pdot, normalization=KAPPA)
    eps = 1.0e-6
    j_plus = -(p + eps * pdot) * d1(f + eps * fdot) / KAPPA
    j_minus = -(p - eps * pdot) * d1(f - eps * fdot) / KAPPA
    measured = (j_plus - j_minus) / (2.0 * eps)

    np.testing.assert_allclose(tangent["total"], measured, rtol=1e-9, atol=1e-10)
    np.testing.assert_allclose(
        tangent["total"],
        tangent["momentum_rate_term"] + tangent["field_gradient_rate_term"],
    )


def test_projection_tangent_only_differences_sector_projection():
    base = {"x": np.array([1.0, 2.0]), "p": np.array([3.0, 4.0])}
    rate = {"x": np.array([0.5, -1.0]), "p": np.array([2.0, 0.25])}

    def shifted(delta):
        return {
            "x": base["x"] + delta * rate["x"],
            "p": base["p"] + delta * rate["p"],
        }

    def project(state):
        return -state["p"] * state["x"]

    coarse = projected_current_directional_tangent(project, shifted, epsilon=1e-3)
    fine = projected_current_directional_tangent(project, shifted, epsilon=5e-4)
    expected = -(rate["p"] * base["x"] + base["p"] * rate["x"])
    np.testing.assert_allclose(coarse["Jdot"], expected, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(fine["Jdot"], expected, rtol=1e-10, atol=1e-10)
    convergence = compare_directional_epsilons(coarse, fine)
    assert convergence["max_abs_difference"] < 1e-10
    assert convergence["is_convergence_estimate_only"] is True


def test_independent_residual_sign_and_admission_gates():
    measured = np.array([10.0, -2.0])
    geom = np.array([3.0, 4.0])
    sector_jdots = {
        "S": np.array([0.1, -0.2]),
        "fluid": np.array([-0.05, 0.01]),
    }
    result = independent_momentum_residual(measured, geom, sector_jdots)
    expected = measured - geom + KAPPA * (sector_jdots["S"] + sector_jdots["fluid"])
    np.testing.assert_allclose(result["residual"], expected)
    assert result["scientific_admission"] is False
    assert result["acceptance_oracle_satisfied"] is False


def test_rejects_invalid_epsilon_and_shape_mismatch():
    with pytest.raises(ValueError):
        projected_current_directional_tangent(lambda s: s, lambda d: np.array([d]), epsilon=0)
    with pytest.raises(ValueError):
        independent_momentum_residual(np.zeros(2), np.zeros(3), {"S": np.zeros(2)})
