"""Tests for the research-only spherical charge-continuity kernel."""
import numpy as np
import pytest

from tools.research_em_charge_continuity import (
    RadialChargeState,
    integrated_outward_current_flux_faces,
    radial_electric_field_faces,
)


def test_zero_face_current_preserves_charge_state_and_shell_charges():
    initial = np.asarray([0.0, 0.15, 0.9, 1.4, 2.25])
    state = RadialChargeState(initial)

    advanced = state.advance(0.125, np.zeros_like(initial))

    np.testing.assert_array_equal(advanced.enclosed_charge, initial)
    np.testing.assert_array_equal(advanced.cell_charge, np.diff(initial))
    assert advanced.enclosed_charge[0] == 0.0


def test_face_flux_update_closes_every_shell_and_global_charge_balance():
    initial = RadialChargeState([0.0, 0.2, 0.75, 1.5, 2.0])
    # This is the sphere-integrated outward current flux, not a pointwise j^r.
    flux = np.asarray([0.0, 0.03, -0.05, 0.2, 0.4])
    dt = 0.125

    updated = initial.advance(dt, flux)

    # Exact finite-volume identity for each shell:
    # Delta(charge_i) = -dt * (I_outer - I_inner).
    np.testing.assert_allclose(
        updated.cell_charge - initial.cell_charge,
        -dt * np.diff(flux),
        rtol=0.0,
        atol=2.0e-16,
    )
    # The total domain charge changes only through the outer boundary.
    assert updated.enclosed_charge[-1] - initial.enclosed_charge[-1] == pytest.approx(
        -dt * flux[-1], rel=0.0, abs=2.0e-16
    )
    assert updated.enclosed_charge[0] == 0.0


def test_eulerian_current_flux_includes_lapse_and_shift_with_declared_sign():
    alpha = np.asarray([1.0, 0.9, 0.8])
    beta = np.asarray([0.0, 0.1, -0.2])
    sqrt_gamma = np.asarray([0.0, 0.25, 1.5])
    rho_e = np.asarray([0.7, 0.6, -0.2])
    current_r = np.asarray([0.4, -0.3, 0.5])

    flux = integrated_outward_current_flux_faces(
        alpha, beta, sqrt_gamma, rho_e, current_r
    )
    expected = 4.0 * np.pi * sqrt_gamma * (alpha * current_r - beta * rho_e)

    np.testing.assert_allclose(flux, expected, rtol=0.0, atol=2.0e-16)
    assert flux[0] == 0.0
    # Positive outward current lowers enclosed charge when used by advance().
    state = RadialChargeState([0.0, 0.5, 1.0])
    advanced = state.advance(0.1, flux)
    assert advanced.enclosed_charge[-1] < state.enclosed_charge[-1]


def test_gauss_reconstruction_uses_areal_radius_and_regular_centre():
    # Radius samples need not be monotone in a generic strong-field slice.
    state = RadialChargeState([0.0, 0.08, -0.12, 0.5])
    areal_radius = np.asarray([0.0, 0.4, 0.25, 1.0])

    field = radial_electric_field_faces(state, areal_radius)

    np.testing.assert_allclose(
        field,
        np.asarray([0.0, 0.08 / 0.16, -0.12 / 0.0625, 0.5]),
        rtol=0.0,
        atol=2.0e-16,
    )


@pytest.mark.parametrize(
    "charge",
    [
        [0.1, 0.2],
        [0.0, np.nan],
        [0.0, np.inf],
        [False, True],
        [[0.0, 0.2]],
    ],
)
def test_charge_state_rejects_malformed_or_nonregular_input(charge):
    with pytest.raises(ValueError):
        RadialChargeState(charge)


@pytest.mark.parametrize("dt", [0.0, -0.1, np.inf, np.nan, True, "bad"])
def test_advance_rejects_invalid_time_steps(dt):
    state = RadialChargeState([0.0, 0.2, 0.7])
    with pytest.raises(ValueError):
        state.advance(dt, [0.0, 0.0, 0.0])


def test_advance_fails_closed_on_nonzero_centre_flux_or_bad_shape():
    state = RadialChargeState([0.0, 0.2, 0.7])
    with pytest.raises(ValueError, match="centre"):
        state.advance(0.1, [1.0e-12, 0.0, 0.0])
    with pytest.raises(ValueError, match="exactly 3"):
        state.advance(0.1, [0.0, 0.0])
    with pytest.raises(ValueError, match="finite"):
        state.advance(0.1, [0.0, np.nan, 0.0])


def test_current_flux_rejects_invalid_metric_factors_and_nonfinite_results():
    base = np.asarray([1.0, 1.0, 1.0])
    zero_center_volume = np.asarray([0.0, 0.2, 0.8])
    rho = np.asarray([0.0, 0.1, 0.1])
    current = np.asarray([0.0, 0.1, 0.1])

    with pytest.raises(ValueError, match="lapse"):
        integrated_outward_current_flux_faces(
            [1.0, 0.0, 1.0], base, zero_center_volume, rho, current
        )
    with pytest.raises(ValueError, match="centre"):
        integrated_outward_current_flux_faces(
            base, base, [0.1, 0.2, 0.8], rho, current
        )
    with pytest.raises(ValueError, match="finite"):
        integrated_outward_current_flux_faces(
            base, base, zero_center_volume, [0.0, np.nan, 0.1], current
        )


def test_gauss_reconstruction_rejects_invalid_areal_radii():
    state = RadialChargeState([0.0, 0.2, 0.7])
    with pytest.raises(ValueError, match="centre"):
        radial_electric_field_faces(state, [0.1, 0.5, 1.0])
    with pytest.raises(ValueError, match="positive"):
        radial_electric_field_faces(state, [0.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="exactly 3"):
        radial_electric_field_faces(state, [0.0, 1.0])
