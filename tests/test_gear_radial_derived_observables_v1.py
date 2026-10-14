"""Unit tests for research-only GEAR derived observables."""
import math

import numpy as np
import pytest

from research.gear_radial_derived_observables_v1 import (
    CORRECTED_ARCHIVE_SOURCE_SHA256,
    LEGACY_HELPER_SOURCE_SHA256,
    balance_residual,
    face_flux_from_centres,
    misner_sharp_face_mass,
    profile_peak_observables,
    rk4_proper_time_increment,
    signed_mode_flux,
)


def test_source_provenance_anchors_are_explicit_and_distinct():
    assert len(CORRECTED_ARCHIVE_SOURCE_SHA256) == 64
    assert len(LEGACY_HELPER_SOURCE_SHA256) == 64
    assert CORRECTED_ARCHIVE_SOURCE_SHA256 != LEGACY_HELPER_SOURCE_SHA256


def test_flat_metric_has_zero_misner_sharp_mass():
    assert misner_sharp_face_mass(10.0, 1.0, 1.0) == pytest.approx(0.0)


def test_flux_preserves_archived_outward_sign_convention():
    r = np.array([1.0, 2.0])
    flux = signed_mode_flux(r, np.ones(2), np.ones(2), -np.ones(2), np.ones(2))
    assert flux[0] == pytest.approx(4.0 * math.pi)
    assert flux[1] == pytest.approx(16.0 * math.pi)


def test_profile_observables_keep_areal_and_proper_radius_separate():
    r = np.array([0.5, 1.5, 2.5])
    a = np.full(3, 2.0)
    result = profile_peak_observables(
        r, a, np.ones(3), np.zeros(3), np.zeros(3),
        np.array([-20.0, -1.0, -0.1]), np.ones(3),
    )
    assert result["proper_radius"] == pytest.approx([1.0, 3.0, 5.0])
    assert result["d_peak_index"] == 0
    assert result["d_peak_areal_radius"] == pytest.approx(0.5)
    assert result["d_peak_proper_radius"] == pytest.approx(1.0)


def test_face_flux_uses_direct_face_gradient_and_same_stage_averages():
    flux = face_flux_from_centres(
        1.0, 0.5, 1.0, 1.0, 1.0, 1.0,
        0.0, -0.5, 0.0, 1.0,
        0.0, -1.0, 0.0, 2.0,
    )
    # Products are (-1)*0.5 and (-2)*1, giving summed product -2.5.
    assert flux == pytest.approx(10.0 * math.pi)


def test_trapezoid_mass_balance_closes_for_constructed_linear_exchange():
    t = np.linspace(0.0, 2.0, 5)
    flux = np.full_like(t, 0.25)
    mass = 5.0 - 0.25 * t
    assert balance_residual(t, mass, flux) == pytest.approx(0.0, abs=1e-14)


def test_rk4_clock_quadrature_uses_the_actual_four_lapse_stages():
    assert rk4_proper_time_increment(0.2, [1.0, 0.9, 0.8, 0.7]) == pytest.approx(
        0.2 * (1 + 2*0.9 + 2*0.8 + 0.7) / 6
    )


@pytest.mark.parametrize(
    "time,mass,flux",
    [
        ([0.0], [1.0], [0.0]),
        ([0.0, 0.0], [1.0, 1.0], [0.0, 0.0]),
        ([0.0, 1.0], [1.0, 1.0], [0.0]),
        ([0.0, float("nan")], [1.0, 1.0], [0.0, 0.0]),
    ],
)
def test_balance_ledger_rejects_invalid_time_series(time, mass, flux):
    with pytest.raises(ValueError):
        balance_residual(time, mass, flux)
