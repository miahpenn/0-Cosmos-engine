"""Tests for the accepted-state Bianchi reference regression gate."""
from __future__ import annotations

import unittest
import numpy as np

from tools.issue4_bianchi_reference_regression import assert_reference_regression


class BianchiReferenceRegressionTests(unittest.TestCase):
    def _fixture(self):
        n = 60
        radii = (np.arange(n, dtype=np.float64) + 0.5) * (80.0 / n)
        point = np.zeros(n)
        point[:3] = (1.8406e-5, 1.4017e-5, 7.7777e-6)
        point[4] = 2.0498e-6
        fv = np.zeros(n)
        fv[:3] = (-1.2950e-6, 9.6404e-6, 6.9830e-6)
        fv[4] = 1.8089e-6
        tangent_a = np.zeros(n)
        tangent_b = tangent_a.copy()
        tangent_b[0] = 1.38e-13
        predictor_a = np.zeros(n)
        predictor_b = predictor_a.copy()
        predictor_b[0] = 7.11e-15
        return {
            "resolution": n,
            "grid_centers": radii,
            "channels": {
                "eps2e-05": {
                    "measured_tangent": tangent_a,
                    "pointwise": {"predictor": predictor_a, "residual_tangent_minus_predictor": point},
                    "finite_volume": {"predictor": predictor_a, "residual_tangent_minus_predictor": fv},
                },
                "eps1e-05": {
                    "measured_tangent": tangent_b,
                    "pointwise": {"predictor": predictor_b, "residual_tangent_minus_predictor": point},
                    "finite_volume": {"predictor": predictor_b, "residual_tangent_minus_predictor": fv},
                },
            },
        }

    def test_signed_reference_witnesses_pass(self):
        result = assert_reference_regression(self._fixture())
        self.assertTrue(result["passed"])
        self.assertEqual(result["archived_signed_witness_orientation"], "same")

    def test_wrong_residual_sign_is_rejected(self):
        fixture = self._fixture()
        fixture["channels"]["eps2e-05"]["pointwise"]["residual_tangent_minus_predictor"] *= -1.0
        with self.assertRaisesRegex(RuntimeError, "regression mismatch"):
            assert_reference_regression(fixture)


if __name__ == "__main__":
    unittest.main()
