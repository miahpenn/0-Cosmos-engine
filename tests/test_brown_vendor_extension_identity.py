"""Symbolic source check for Brown vs vendor off-constraint C-linear propagation.

The pinned vendor's radial Ricci bracket contains -0.5*Lambda*a'. Holding
the metric fixed and perturbing Lambda by C therefore changes the radial
Ricci eigenvalue by chi*C' + chi*(a'/(2a))*C, while scalar Ricci changes by
chi*C'. This test expands the resulting A_a and M_r tangent symbolically.

This is complementary to the exact-pinned-vendor runtime discriminator. It
uses no fitted coefficients, production evolution, or scientific-admission
thresholds.
"""
from __future__ import annotations

import sympy as sp


def test_vendor_c_linear_response_equals_brown_plus_derived_extension():
    r = sp.symbols("r", real=True)
    alpha = sp.Function("alpha")(r)
    chi = sp.Function("chi")(r)
    q_a = sp.Function("q_a")(r)  # q_a = a'/(2a)
    C = sp.Function("C")(r)      # vendor connection constraint
    Q = sp.Function("Q")(r)      # Q = 2/r + b'/b

    C_r = sp.diff(C, r)
    C_rr = sp.diff(C, r, 2)
    alpha_chi = alpha * chi

    # Delta R^r_r and Delta R follow directly from the two pinned Ricci
    # brackets, including their different undifferentiated-Lambda content.
    delta_Rrr = chi * C_r + chi * q_a * C
    delta_R = chi * C_r
    delta_Aa = sp.simplify(alpha * (delta_Rrr - delta_R / 3))
    expected_delta_Aa = (
        sp.Rational(2, 3) * alpha_chi * C_r
        + alpha_chi * q_a * C
    )
    assert sp.simplify(delta_Aa - expected_delta_Aa) == 0

    # The Aa coefficient in the implemented M_r expression is Gamma_M.
    # Since chip = -X'/(2X) and chi = X^2, 6*chip = -3/2*chi'/chi.
    gamma_M = (
        -sp.Rational(3, 2) * sp.diff(chi, r) / chi
        + sp.Rational(3, 2) * Q
    )
    actual_delta_M = sp.diff(delta_Aa, r) + gamma_M * delta_Aa

    brown_C_terms = (
        sp.Rational(2, 3) * alpha_chi * C_rr
        + (
            sp.Rational(2, 3) * sp.diff(alpha, r) * chi
            - sp.Rational(1, 3) * alpha * sp.diff(chi, r)
            + alpha_chi * Q
        ) * C_r
    )
    extension = (
        alpha_chi * q_a * C_r
        + (
            q_a * (sp.diff(alpha_chi, r) + gamma_M * alpha_chi)
            + alpha_chi * sp.diff(q_a, r)
        ) * C
    )

    # Exact identity: no fitted coefficient and no evaluation at sample data.
    assert sp.simplify(actual_delta_M - brown_C_terms - extension) == 0

    # Deliberate wrong-sign control must not satisfy the same identity.
    assert sp.simplify(actual_delta_M - brown_C_terms + extension) != 0
