"""Contracts for the representation-aware proper-radial observer."""
import numpy as np
import pytest

from research.gear_proper_radial_v1 import (
    archive_proper_radius,
    bssn_proper_radius,
    profile_on_proper_grid,
)


def test_flat_archive_and_flat_bssn_geometry_preserve_areal_coordinate():
    radius = np.array([0.5, 1.5, 2.5, 3.5])
    archive_ell = archive_proper_radius(radius, np.ones_like(radius))
    bssn_ell = bssn_proper_radius(radius, np.ones_like(radius), np.ones_like(radius))
    assert archive_ell == pytest.approx(radius)
    assert bssn_ell == pytest.approx(radius)


def test_constant_archive_line_factor_scales_proper_distance():
    radius = np.array([0.5, 1.5, 2.5])
    ell = archive_proper_radius(radius, np.full_like(radius, 2.0))
    assert ell == pytest.approx(2.0 * radius)


def test_bssn_line_factor_uses_sqrt_metric_a_over_X_not_metric_a_directly():
    radius = np.array([0.5, 1.5, 2.5])
    metric_a = np.full_like(radius, 4.0)
    X = np.full_like(radius, 2.0)
    bssn_ell = bssn_proper_radius(radius, metric_a, X)
    archive_ell = archive_proper_radius(radius, metric_a)
    assert bssn_ell == pytest.approx(radius)
    assert archive_ell == pytest.approx(4.0 * radius)
    assert not np.allclose(bssn_ell, archive_ell)


def test_nonconstant_line_factor_uses_midpoint_plus_trapezoid_quadrature():
    radius = np.array([0.5, 1.5, 2.5])
    factor = np.array([2.0, 3.0, 4.0])
    ell = archive_proper_radius(radius, factor)
    assert ell == pytest.approx([1.0, 3.5, 7.0])


def test_profile_remap_interpolates_without_changing_source_samples():
    ell = np.array([0.5, 1.5, 2.5])
    profile = np.array([2.0, 4.0, 6.0])
    target = np.array([0.5, 1.0, 2.5])
    remapped = profile_on_proper_grid(ell, profile, target)
    assert remapped == pytest.approx([2.0, 3.0, 6.0])
    assert profile == pytest.approx([2.0, 4.0, 6.0])


@pytest.mark.parametrize(
    "radius,factor",
    [
        ([0.5, 0.4], [1.0, 1.0]),
        ([0.0, 1.0], [1.0, 1.0]),
        ([0.5, 1.5], [1.0, 0.0]),
        ([0.5, 1.5], [1.0, np.nan]),
        ([0.5, 1.5], [1.0]),
    ],
)
def test_archive_map_fails_closed_on_invalid_grid_or_metric(radius, factor):
    with pytest.raises(ValueError):
        archive_proper_radius(radius, factor)


@pytest.mark.parametrize(
    "metric_a,X",
    [
        ([1.0, 1.0], [1.0]),
        ([1.0, -1.0], [1.0, 1.0]),
        ([1.0, np.inf], [1.0, 1.0]),
        ([1.0, 1.0], [1.0, 0.0]),
        ([1.0, 1.0], [1.0, np.nan]),
    ],
)
def test_bssn_map_fails_closed_on_invalid_metric_data(metric_a, X):
    radius = [0.5, 1.5]
    with pytest.raises(ValueError):
        bssn_proper_radius(radius, metric_a, X)


@pytest.mark.parametrize(
    "ell,profile,target",
    [
        ([0.5, 0.4], [1.0, 2.0], [0.45]),
        ([0.5, 1.5], [1.0], [0.5]),
        ([0.5, 1.5], [1.0, 2.0], [0.1]),
        ([0.5, 1.5], [1.0, 2.0], [2.0]),
        ([0.5, np.nan], [1.0, 2.0], [0.5]),
    ],
)
def test_profile_remap_rejects_invalid_or_extrapolated_samples(ell, profile, target):
    with pytest.raises(ValueError):
        profile_on_proper_grid(ell, profile, target)
