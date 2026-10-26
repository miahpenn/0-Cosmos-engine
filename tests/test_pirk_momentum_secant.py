import unittest

import numpy as np

from tools.audit_pirk_momentum_secant import (
    constraint_secant_components,
    d1,
    d2,
    momentum_constraint_geom,
)


class PirkMomentumSecantTests(unittest.TestCase):
    def test_pinned_derivatives_reproduce_regular_polynomials(self):
        dr = 0.2
        r = (np.arange(11, dtype=float) + 0.5) * dr
        even = 2.0 + 3.0 * r**2 + 0.25 * r**4
        odd = r + 0.5 * r**3
        np.testing.assert_allclose(d1(even, dr, 1), 6.0*r + r**3, rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(d2(even, dr, 1), 6.0 + 3.0*r**2, rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(d1(odd, dr, -1), 1.0 + 1.5*r**2, rtol=2e-12, atol=2e-12)

    def test_exact_constraint_secant_includes_nonlinear_cross_terms(self):
        dr = 0.1
        r = (np.arange(61, dtype=float) + 0.5) * dr
        b0 = 1.0 + 0.03*r**2/(1.0+r**2)
        b1 = b0 + 2.0e-4*r**2/(1.0+r**2)
        X0 = 1.0 + 0.02/(1.0+r**2)
        X1 = X0 - 1.0e-4/(1.0+r**2)
        A0 = 2.0e-3*r**2/(1.0+r**2)
        A1 = A0 + 7.0e-5*r**2/(1.0+r**2)
        K0 = 4.0e-3*np.exp(-r**2/8.0)
        K1 = K0 - 3.0e-5*np.exp(-r**2/8.0)
        j0 = 2.0e-5*np.exp(-r**2/10.0)
        j1 = j0 + 3.0e-6*np.exp(-r**2/9.0)
        delta_j = j1-j0

        M0 = momentum_constraint_geom(r, b0, X0, A0, K0, dr) - 8.0*np.pi*j0
        M1 = momentum_constraint_geom(r, b1, X1, A1, K1, dr) - 8.0*np.pi*j1
        linear, cross = constraint_secant_components(
            r, dr, b0, b1, X0, X1, A0, A1, K0, K1, delta_j
        )
        np.testing.assert_allclose(M1-M0, linear+cross, rtol=2e-10, atol=2e-13)
        self.assertGreater(float(np.max(np.abs(cross))), 0.0)


if __name__ == "__main__":
    unittest.main()
