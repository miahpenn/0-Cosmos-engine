"""Contract tests for the exact discrete geometric-momentum chain rule."""
from __future__ import annotations

import numpy as np

from tools.momentum_constraint_chain_rule import analytic_geometric_momentum_tangent


def _linear_operator(n: int):
    """A fixed linear stencil; the identity must hold for any fixed linear D1."""
    matrix = np.zeros((n, n), dtype=np.float64)
    matrix[0, :3] = (-1.5, 2.0, -0.5)
    for i in range(1, n - 1):
        matrix[i, i - 1] = -0.5
        matrix[i, i + 1] = 0.5
    matrix[-1, -3:] = (0.5, -2.0, 1.5)
    return lambda values: matrix @ np.asarray(values, dtype=np.float64)


def _discrete_momentum(r, Aa, K, b, X, d1):
    Aa = np.asarray(Aa, dtype=np.float64)
    K = np.asarray(K, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    X = np.asarray(X, dtype=np.float64)
    chi_r = -0.5 * d1(X) / X
    return (
        d1(Aa)
        - (2.0 / 3.0) * d1(K)
        + 6.0 * Aa * chi_r
        + 1.5 * Aa * (2.0 / r + d1(b) / b)
    )


def test_exact_chain_rule_matches_two_directional_sizes_and_detects_sign_mutation():
    n = 9
    r = 0.5 + np.arange(n, dtype=np.float64)
    d1 = _linear_operator(n)
    Aa = 1.0e-3 + 2.0e-4 * r + 1.0e-5 * r**2
    K = -2.0e-2 + 3.0e-4 * r + 7.0e-5 * r**2
    b = 1.0 + 2.0e-4 * r + 1.0e-5 * r**2
    X = 1.1 + 4.0e-4 * r + 1.0e-5 * r**2
    dot_Aa = 2.0e-5 + 1.0e-5 * r
    dot_K = 1.0e-4 - 2.0e-6 * r**2
    dot_b = 2.0e-5 * (1.0 + r)
    dot_X = -3.0e-5 + 1.0e-6 * r**2

    analytic = analytic_geometric_momentum_tangent(
        d1, r, Aa, K, b, X, dot_Aa, dot_K, dot_b, dot_X
    )
    differences = []
    for epsilon in (2.0e-3, 1.0e-3):
        plus = _discrete_momentum(
            r, Aa + epsilon * dot_Aa, K + epsilon * dot_K,
            b + epsilon * dot_b, X + epsilon * dot_X, d1
        )
        minus = _discrete_momentum(
            r, Aa - epsilon * dot_Aa, K - epsilon * dot_K,
            b - epsilon * dot_b, X - epsilon * dot_X, d1
        )
        measured = (plus - minus) / (2.0 * epsilon)
        differences.append(measured - analytic["total"])
        assert float(np.max(np.abs(differences[-1]))) < 5.0e-10

    # Deliberately flip one signed term. The same cross-check must reject it.
    wrong_terms = dict(analytic["terms"])
    wrong_terms["minus_2third_D_r_K"] = -wrong_terms["minus_2third_D_r_K"]
    wrong_total = np.zeros_like(r)
    for values in wrong_terms.values():
        wrong_total += values

    fine_difference = (
        _discrete_momentum(
            r, Aa + 1.0e-3 * dot_Aa, K + 1.0e-3 * dot_K,
            b + 1.0e-3 * dot_b, X + 1.0e-3 * dot_X, d1
        )
        - _discrete_momentum(
            r, Aa - 1.0e-3 * dot_Aa, K - 1.0e-3 * dot_K,
            b - 1.0e-3 * dot_b, X - 1.0e-3 * dot_X, d1
        )
    ) / (2.0e-3)
    wrong_residual = float(np.max(np.abs(fine_difference - wrong_total)))
    assert wrong_residual > 1.0e-6
    assert wrong_residual > 100.0 * float(np.max(np.abs(differences[-1])))
