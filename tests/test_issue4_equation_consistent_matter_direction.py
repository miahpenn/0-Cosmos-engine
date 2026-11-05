import numpy as np

from tools.issue4_manufactured_regular_family import PROFILE, STATE_ARRAY_NAMES
from tools.issue4_equation_consistent_matter_direction import (
    analytic_matter_rates,
    construct_equation_consistent_family,
    construct_equation_consistent_triplet,
    gaussian_laplacian,
    gaussian_radial_derivative,
)


def grid():
    n, rmax = 60, 40.0
    return (np.arange(n, dtype=np.float64) + 0.5) * rmax / n


def test_exact_gaussian_radial_operators():
    r = grid()
    f = 0.125 * np.exp(-(r / 8.0)**2)
    np.testing.assert_allclose(
        gaussian_laplacian(r, f, 8.0),
        f * (4*r*r/8.0**4 - 6/8.0**2),
        rtol=0, atol=0,
    )
    np.testing.assert_allclose(
        gaussian_radial_derivative(r, f, 8.0), -2*r*f/8.0**2,
        rtol=0, atol=0,
    )


def test_scalar_rates_are_pinned_flat_slice_equation_forms():
    r = grid()
    rates = analytic_matter_rates(r)
    shape = np.exp(-(r / PROFILE["scalar_width"])**2)
    S = PROFILE["scalar_amplitudes"]["S"] * shape
    D = PROFILE["scalar_amplitudes"]["D"] * shape
    phi = PROFILE["scalar_amplitudes"]["phi"] * shape
    lap = lambda f: gaussian_laplacian(r, f, PROFILE["scalar_width"])
    rho_dm = PROFILE["fluid_rest_density_amplitudes"]["dark_matter"] * np.exp(
        -(r / PROFILE["fluid_width"])**2
    )
    Vprime = -8.242522415500654e-5 * np.exp(-phi)
    np.testing.assert_allclose(rates["S_dot"], 0.0, atol=0.0)
    np.testing.assert_allclose(rates["PS_dot"], lap(S)-S, rtol=0, atol=0)
    np.testing.assert_allclose(rates["D_dot"], 0.0, atol=0.0)
    np.testing.assert_allclose(rates["PD_dot"], lap(D)+D, rtol=0, atol=0)
    np.testing.assert_allclose(rates["phi_dot"], 0.0, atol=0.0)
    np.testing.assert_allclose(rates["Pi_dot"], lap(phi)-Vprime-0.04*rho_dm, rtol=0, atol=0)


def test_fluid_continuity_and_euler_rates_at_rest():
    r = grid()
    rates = analytic_matter_rates(r)
    for species in ("dark_matter", "baryons", "radiation"):
        np.testing.assert_allclose(rates[f"{species}_rho_dot"], 0.0, atol=0.0)
    phi = PROFILE["scalar_amplitudes"]["phi"] * np.exp(-(r/8.0)**2)
    phi_r = gaussian_radial_derivative(r, phi, 8.0)
    np.testing.assert_allclose(rates["dark_matter_v_dot"], -0.04*phi_r, rtol=0, atol=0)
    np.testing.assert_allclose(rates["baryons_v_dot"], 0.0, atol=0.0)
    rho = PROFILE["fluid_rest_density_amplitudes"]["radiation"] * np.exp(-(r/10.0)**2)
    dp = gaussian_radial_derivative(r, rho, 10.0)/3.0
    np.testing.assert_allclose(rates["radiation_v_dot"], -dp/(4*rho/3.0), rtol=0, atol=0)


def test_family_is_canonical_finite_and_zero_step():
    r = grid()
    state = construct_equation_consistent_family(r, 0.0)
    assert set(state) == set(STATE_ARRAY_NAMES)
    assert all(value.dtype == np.float64 and value.shape == r.shape for value in state.values())
    assert all(np.all(np.isfinite(value)) for value in state.values())
    assert np.all(state["scalars.PS"] == 0.0)
    assert np.all(state["scalars.PD"] == 0.0)
    assert np.all(state["scalars.Pi"] == 0.0)
    for species in ("dark_matter", "baryons", "radiation"):
        assert np.all(state[f"matter.{species}.momentum_r"] == 0.0)
    assert np.all(state["matter.radiation.rest"] == 0.0)


def test_plus_minus_tangent_is_signed_and_density_continuity_is_zero():
    r = grid()
    eps = 1.0e-5
    triplet = construct_equation_consistent_triplet(r, eps)
    rates = analytic_matter_rates(r)
    for field, rate in (("PS", "PS_dot"), ("PD", "PD_dot"), ("Pi", "Pi_dot")):
        measured = (triplet["plus"][f"scalars.{field}"] -
                    triplet["minus"][f"scalars.{field}"]) / (2*eps)
        np.testing.assert_allclose(measured, rates[rate], rtol=1e-13, atol=1e-15)
    for species in ("dark_matter", "baryons", "radiation"):
        measured = (triplet["plus"][f"matter.{species}.momentum_r"] -
                    triplet["minus"][f"matter.{species}.momentum_r"]) / (2*eps)
        if species == "dark_matter":
            rho = PROFILE["fluid_rest_density_amplitudes"][species] * np.exp(-(r/10.0)**2)
            expected = r*r*rho*rates["dark_matter_v_dot"]
        elif species == "baryons":
            expected = np.zeros_like(r)
        else:
            rho = PROFILE["fluid_rest_density_amplitudes"][species] * np.exp(-(r/10.0)**2)
            expected = r*r*(4*rho/3.0)*rates["radiation_v_dot"]
        np.testing.assert_allclose(measured, expected, rtol=1e-12, atol=1e-15)
        # First-order continuity at the zero-velocity base: density tangent = 0.
        for field in ("rest", "energy_t"):
            diff = (triplet["plus"][f"matter.{species}.{field}"] -
                    triplet["minus"][f"matter.{species}.{field}"]) / (2*eps)
            np.testing.assert_allclose(diff, 0.0, atol=1e-15)
