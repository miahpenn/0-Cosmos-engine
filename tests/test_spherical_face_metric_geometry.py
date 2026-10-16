"""Regression tests for geometrically located interior finite-volume face metrics."""

from __future__ import annotations

import numpy as np
import pytest

from engine.matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    MetricDerivativeSet,
    Species,
    _face_metric,
    _metric_arrays,
    _species_rhs_components,
)
from engine.valencia import spherical_metric_from_bssn


def _flat_slice(n: int, dr: float) -> BSSNMetricSlice:
    r = (np.arange(n, dtype=float) + 0.5) * dr
    one = np.ones(n, dtype=float)
    zero = np.zeros(n, dtype=float)
    return BSSNMetricSlice(
        r=r,
        a=one.copy(),
        b=one.copy(),
        X=one.copy(),
        alpha=one.copy(),
        beta=zero.copy(),
        Aa=zero.copy(),
        K=zero.copy(),
        Lambda=zero.copy(),
        B=zero.copy(),
    )


@pytest.mark.parametrize("dr", [0.5, 1.0, 4.0 / 3.0])
def test_flat_interior_face_uses_face_radius_for_volume_measure(dr):
    """At r=dr, sqrt(gamma)=dr^2, not the 5/4 cell-centre average."""
    metric = _flat_slice(n=8, dr=dr)
    face = _face_metric(metric, 0)
    cell_metrics = _metric_arrays(metric)
    arithmetic_cell_average = 0.5 * (
        cell_metrics[0].sqrt_gamma + cell_metrics[1].sqrt_gamma
    )

    assert face.sqrt_gamma == pytest.approx(dr * dr, rel=2.0e-15, abs=0.0)
    assert face.gamma_thth == pytest.approx(dr * dr, rel=2.0e-15, abs=0.0)
    assert arithmetic_cell_average == pytest.approx(
        1.25 * dr * dr, rel=2.0e-15, abs=0.0
    )
    assert abs(face.sqrt_gamma - arithmetic_cell_average) > 0.2 * dr * dr


def test_nontrivial_face_metric_is_derived_from_interpolated_primary_fields():
    """Interpolate a,b,X,alpha,beta, then derive every face metric component."""
    metric = BSSNMetricSlice(
        r=np.asarray([0.5, 1.5]),
        a=np.asarray([1.0, 1.4]),
        b=np.asarray([0.8, 1.2]),
        X=np.asarray([0.9, 1.1]),
        alpha=np.asarray([0.95, 0.75]),
        beta=np.asarray([-0.1, 0.2]),
        Aa=np.zeros(2),
        K=np.zeros(2),
        Lambda=np.zeros(2),
        B=np.zeros(2),
    )
    expected = spherical_metric_from_bssn(
        r=1.0,
        a=1.2,
        b=1.0,
        X=1.0,
        alpha=0.85,
        beta=0.05,
    )
    actual = _face_metric(metric, 0)

    for field in (
        "alpha",
        "beta",
        "gamma_rr",
        "gamma_rr_inv",
        "gamma_thth",
        "gamma_thth_inv",
        "sqrt_gamma",
    ):
        assert getattr(actual, field) == pytest.approx(
            getattr(expected, field), rel=2.0e-15, abs=0.0
        )


def test_homogeneous_flat_radiation_momentum_rhs_cancels_in_interior_cells():
    """The exact face volume and spherical geometric source cancel cellwise."""
    dr = 0.5
    metric = _flat_slice(n=8, dr=dr)
    radii = metric.r
    cell_sqrt_gamma = radii * radii
    rho_r = 9.2e-5 * 9.64116e-5
    state = ConservedSpecies(
        rest=np.zeros_like(radii),
        energy_t=cell_sqrt_gamma * rho_r,
        momentum_r=np.zeros_like(radii),
    )
    derivatives = MetricDerivativeSet(
        time={},
        radial={
            "alpha": np.zeros_like(radii),
            "beta": np.zeros_like(radii),
            "rr": np.zeros_like(radii),
            "thth": 2.0 * radii,
        },
    )

    face_flux, _source_e, source_s, _q_e, _q_s = _species_rhs_components(
        metric, derivatives, state, Species.RADIATION
    )
    rhs_momentum = -np.diff(face_flux[:, 2]) / dr + source_s

    # The outermost face is a separate boundary prescription and is excluded.
    # Cells 0..n-2 use geometrically reconstructed interior face metrics.
    assert float(np.max(np.abs(rhs_momentum[:-1]))) <= 5.0e-23
    assert abs(float(rhs_momentum[0])) <= 5.0e-23


def test_interior_face_index_is_checked():
    metric = _flat_slice(n=2, dr=1.0)
    with pytest.raises(IndexError, match="interior face index"):
        _face_metric(metric, -1)
    with pytest.raises(IndexError, match="interior face index"):
        _face_metric(metric, 1)
