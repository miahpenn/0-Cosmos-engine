"""Tests for research-only radial EM source-convention identities."""
import numpy as np
import pytest

from tools.research_em_active_source_conversion import radial_lorentz_source_map


def test_lorentz_map_matches_covector_eulerian_sources_and_mixed_reference():
    alpha = np.asarray([0.65, 0.9, 1.2, 1.4])
    beta = np.asarray([-0.3, 0.0, 0.2, -0.15])
    sqrt_gamma = np.asarray([0.0, 0.12, 0.7, 2.1])
    electric_covector_r = np.asarray([0.0, 0.4, -0.85, 1.1])
    rho_e = np.asarray([0.8, -0.1, 0.45, 0.2])
    current_r = np.asarray([0.0, 0.3, -0.7, 0.9])

    source = radial_lorentz_source_map(
        alpha, beta, sqrt_gamma, electric_covector_r, rho_e, current_r
    )

    expected_q_t = (
        -alpha * electric_covector_r * current_r
        + beta * rho_e * electric_covector_r
    )
    expected_q_r = rho_e * electric_covector_r
    expected_eulerian_energy = (
        alpha * sqrt_gamma * electric_covector_r * current_r
    )
    expected_eulerian_momentum = alpha * sqrt_gamma * expected_q_r

    np.testing.assert_allclose(
        source.force_covector_t, expected_q_t, rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        source.force_covector_r, expected_q_r, rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        source.eulerian_energy, expected_eulerian_energy, rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        source.eulerian_momentum_r,
        expected_eulerian_momentum,
        rtol=0.0,
        atol=0.0,
    )

    # The production matter arrays use Eulerian Valencia energy/momentum
    # stores for dust/DM and baryons too, as confirmed by initialize_dust,
    # _valencia_conserved, _valencia_flux and recover_dust. The following is a
    # separate covariant mixed-component reference only, not active storage.
    expected_mixed_energy = (
        -alpha * source.eulerian_energy
        + beta * source.eulerian_momentum_r
    )
    np.testing.assert_allclose(
        source.mixed_reference_energy_t,
        expected_mixed_energy,
        rtol=0.0,
        atol=8.0e-16,
    )
    np.testing.assert_allclose(
        source.mixed_reference_momentum_r,
        source.eulerian_momentum_r,
        rtol=0.0,
        atol=0.0,
    )

    # The independent mixed-reference result follows the covariant identity
    # sqrt(-g) Q_nu, with sqrt(-g)=alpha*sqrt(gamma). This catches sign and
    # lapse/shift errors without claiming mixed storage is used by production.
    np.testing.assert_allclose(
        source.mixed_energy_t,
        alpha * sqrt_gamma * source.force_covector_t,
        rtol=0.0,
        atol=8.0e-16,
    )
    np.testing.assert_allclose(
        source.mixed_momentum_r,
        alpha * sqrt_gamma * source.force_covector_r,
        rtol=0.0,
        atol=0.0,
    )
    assert source.mixed_energy_t[1] == pytest.approx(
        -alpha[1] * source.eulerian_energy[1], abs=1.0e-16
    )


def test_zero_electric_field_has_zero_force_in_every_representation():
    source = radial_lorentz_source_map(
        alpha=[0.8, 1.1],
        beta_r=[0.3, -0.4],
        sqrt_gamma=[0.2, 1.4],
        electric_covector_r=[0.0, 0.0],
        charge_density=[2.0, -3.0],
        spatial_current_r=[-4.0, 5.0],
    )
    for field in (
        source.force_covector_t,
        source.force_covector_r,
        source.eulerian_energy,
        source.eulerian_momentum_r,
        source.mixed_energy_t,
        source.mixed_momentum_r,
    ):
        np.testing.assert_array_equal(field, np.zeros(2))


@pytest.mark.parametrize(
    "overrides",
    [
        {"alpha": [1.0, 0.0]},
        {"alpha": [1.0, np.nan]},
        {"sqrt_gamma": [0.2, -0.1]},
        {"sqrt_gamma": [0.2, np.inf]},
        {"electric_covector_r": [0.2, np.nan]},
        {"charge_density": [0.2, np.inf]},
        {"spatial_current_r": [0.2, np.nan]},
        {"beta_r": [False, True]},
        {"alpha": []},
    ],
)
def test_lorentz_map_fails_closed_on_invalid_arrays(overrides):
    args = {
        "alpha": [1.0, 0.9],
        "beta_r": [0.1, -0.2],
        "sqrt_gamma": [0.2, 0.8],
        "electric_covector_r": [0.3, -0.4],
        "charge_density": [0.5, -0.6],
        "spatial_current_r": [-0.7, 0.8],
    }
    args.update(overrides)
    with pytest.raises(ValueError):
        radial_lorentz_source_map(**args)


def test_lorentz_map_rejects_shape_mismatch_and_nonfinite_products():
    with pytest.raises(ValueError, match="same shape"):
        radial_lorentz_source_map(
            alpha=[1.0, 0.9],
            beta_r=[0.0],
            sqrt_gamma=[0.2, 0.8],
            electric_covector_r=[0.1, 0.2],
            charge_density=[0.3, 0.4],
            spatial_current_r=[0.5, 0.6],
        )

    with pytest.raises(ValueError, match="non-finite"):
        radial_lorentz_source_map(
            alpha=[1.0, 1.0],
            beta_r=[0.0, 0.0],
            sqrt_gamma=[1.0e308, 1.0],
            electric_covector_r=[1.0e308, 0.2],
            charge_density=[0.1, 0.4],
            spatial_current_r=[1.0e308, 0.6],
        )
