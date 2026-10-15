import math
import unittest

from tools.derived_physics_diagnostics import (
    bianchi_ix_S_tau_derivative,
    bianchi_ix_trace_term,
    directional_kasner_exponents,
    heron_discriminant,
    omega_sigma_from_directional_shear,
    raychaudhuri_bounce_margin,
    scalar_kasner_exponent_square_sum,
)


class ShearDiagnosticsTests(unittest.TestCase):
    def test_omega_sigma_uses_archive_normalization(self):
        self.assertAlmostEqual(
            omega_sigma_from_directional_shear((2.0, -1.0, -1.0), 3.0),
            1.0 / 9.0,
        )

    def test_omega_sigma_rejects_h_zero_and_nontraceless_shear(self):
        with self.assertRaisesRegex(ValueError, "H must be nonzero"):
            omega_sigma_from_directional_shear((0.0, 0.0, 0.0), 0.0)
        with self.assertRaisesRegex(ValueError, "sum to zero"):
            omega_sigma_from_directional_shear((1.0, 0.0, 0.0), 1.0)

    def test_orientation_that_is_dangerous_flips_with_hubble_sign(self):
        omega = 0.30
        s = math.sqrt(omega)  # |H|=1 and sum sigma_i**2=6 Omega_sigma H**2
        expanding_dangerous = directional_kasner_exponents((-2*s, s, s), 1.0)
        expanding_safe = directional_kasner_exponents((2*s, -s, -s), 1.0)
        contracting_dangerous = directional_kasner_exponents((2*s, -s, -s), -1.0)
        contracting_safe = directional_kasner_exponents((-2*s, s, s), -1.0)
        self.assertLess(min(expanding_dangerous), 0.0)
        self.assertGreater(min(expanding_safe), 0.0)
        self.assertLess(min(contracting_dangerous), 0.0)
        self.assertGreater(min(contracting_safe), 0.0)


class BounceDiagnosticTests(unittest.TestCase):
    def test_raychaudhuri_bounce_margin(self):
        self.assertGreater(raychaudhuri_bounce_margin(-3.0, 0.5), 0.0)
        self.assertEqual(raychaudhuri_bounce_margin(-1.0, 0.5), 0.0)
        self.assertLess(raychaudhuri_bounce_margin(0.1, 0.5), 0.0)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            raychaudhuri_bounce_margin(-1.0, -0.1)


class BianchiIXDiagnosticsTests(unittest.TestCase):
    def test_heron_discriminant_triangle_and_nontriangle_cases(self):
        self.assertAlmostEqual(heron_discriminant(3.0, 4.0, 5.0), 576.0)
        self.assertAlmostEqual(heron_discriminant(1.0, 1.0, 3.0), -45.0)
        self.assertEqual(heron_discriminant(1.0, 1.0, 2.0), 0.0)

    def test_stiff_trace_correction_vanishes_even_with_anisotropic_pressures(self):
        self.assertEqual(bianchi_ix_trace_term(1.0, 1.0, (0.0, 2.0, 4.0)), -3.0)
        self.assertEqual(bianchi_ix_trace_term(1.0, 1.0, (1.0, 1.0, 1.0)), 0.0)

    def test_S_slope_contains_trace_corrected_term(self):
        # For a=b=c=1, T=3 and V=1. Isotropic stiff pressure makes trace correction zero.
        self.assertAlmostEqual(
            bianchi_ix_S_tau_derivative((1.0, 1.0, 1.0), 1.0, (1.0, 1.0, 1.0)),
            -1.5,
        )
        # Pressures (0,2,4) sum to 6; for rho=1 the trace correction is -3 and S'=-3.
        self.assertAlmostEqual(
            bianchi_ix_S_tau_derivative((1.0, 1.0, 1.0), 1.0, (0.0, 2.0, 4.0)),
            -3.0,
        )

    def test_scalar_kasner_square_sum_formula_and_vacuum_control(self):
        self.assertEqual(scalar_kasner_exponent_square_sum(0.0, 5.0), 1.0)
        self.assertAlmostEqual(scalar_kasner_exponent_square_sum(0.5, 1.0), 0.75)
        with self.assertRaisesRegex(ValueError, "strictly positive"):
            scalar_kasner_exponent_square_sum(0.5, 0.0)


if __name__ == "__main__":
    unittest.main()
