"""Structural tests for the initializer-free Issue #4 manufactured family."""
from __future__ import annotations

import unittest
import numpy as np

from tools.issue4_manufactured_regular_family import (
    PROFILE,
    STATE_ARRAY_NAMES,
    construct_regular_family,
    construct_triplet,
)


class ManufacturedRegularFamilyTests(unittest.TestCase):
    def setUp(self):
        self.r = (np.arange(60, dtype=np.float64) + 0.5) * (40.0 / 60.0)

    def test_complete_canonical_six_sector_state_is_finite_float64(self):
        state = construct_regular_family(self.r)
        self.assertEqual(set(state), set(STATE_ARRAY_NAMES))
        for name, value in state.items():
            self.assertEqual(value.dtype, np.float64, name)
            self.assertEqual(value.shape, self.r.shape, name)
            self.assertTrue(np.all(np.isfinite(value)), name)

    def test_exact_flat_geometry_zero_shift_and_no_em_carrier(self):
        state = construct_regular_family(self.r)
        for name in ("a", "b", "X", "alpha"):
            self.assertTrue(np.all(state[f"geometry.{name}"] == 1.0))
        for name in ("Aa", "K", "Lambda", "beta", "B"):
            self.assertTrue(np.all(state[f"geometry.{name}"] == 0.0))

    def test_family_changes_only_declared_scalar_momenta_and_fluid_velocities(self):
        base = construct_regular_family(self.r, 0.0)
        plus = construct_regular_family(self.r, 1e-5)
        minus = construct_regular_family(self.r, -1e-5)
        changing = {
            "scalars.PS", "scalars.PD", "scalars.Pi",
            *(f"matter.{s}.{f}" for s in
              ("dark_matter", "baryons", "radiation")
              for f in ("energy_t", "momentum_r")),
        }
        for name in STATE_ARRAY_NAMES:
            if name not in changing:
                np.testing.assert_array_equal(plus[name], base[name], err_msg=name)
                np.testing.assert_array_equal(minus[name], base[name], err_msg=name)
        for name in ("scalars.PS", "scalars.PD", "scalars.Pi"):
            np.testing.assert_allclose(plus[name], -minus[name], rtol=0, atol=0)
        for species in ("dark_matter", "baryons", "radiation"):
            np.testing.assert_allclose(
                plus[f"matter.{species}.momentum_r"],
                -minus[f"matter.{species}.momentum_r"], rtol=1e-14, atol=1e-18
            )

    def test_triplet_members_are_independent_arrays_and_epsilon_is_validated(self):
        triplet = construct_triplet(self.r, 2e-5)
        self.assertEqual(set(triplet), {"base", "plus", "minus"})
        for name in STATE_ARRAY_NAMES:
            self.assertIsNot(triplet["base"][name], triplet["plus"][name])
            self.assertIsNot(triplet["base"][name], triplet["minus"][name])
        with self.assertRaises(ValueError):
            construct_triplet(self.r, 0.0)

    def test_regular_profile_parity_and_manifest_coefficients_are_explicit(self):
        self.assertEqual(PROFILE["geometry"], "exact_flat_areal_slice")
        self.assertEqual(PROFILE["family_parameter"], "s")
        # The analytic radial velocity profile is r times an even Gaussian.
        velocity = PROFILE["family_velocity_slopes"]["dark_matter"] * self.r * np.exp(
            -(self.r / PROFILE["fluid_width"])**2
        )
        self.assertTrue(np.all(velocity > 0.0))
        self.assertEqual(PROFILE["radiation_equation_of_state"], "p=rho/3")


if __name__ == "__main__":
    unittest.main()
