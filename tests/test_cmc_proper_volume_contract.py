"""Fail-closed contract tests for the production CMC proper-volume projection."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine import cmc_gauge


def _inputs(*, a=None, b=None, X=None, centers=None, dr=0.5):
    grid = SimpleNamespace(
        centers=np.asarray(
            [0.25, 0.75, 1.25] if centers is None else centers,
            dtype=float,
        ),
        dr=dr,
    )
    geometry = SimpleNamespace(
        a=np.asarray([1.0, 1.2, 1.4] if a is None else a, dtype=float),
        b=np.asarray([1.0, 0.9, 0.8] if b is None else b, dtype=float),
        X=np.asarray([1.0, 1.1, 1.2] if X is None else X, dtype=float),
    )
    return grid, geometry


def test_proper_volume_weights_match_spherical_bssn_measure():
    grid, geometry = _inputs()
    weights, total = cmc_gauge._proper_volume_weights(grid, geometry)
    expected = (
        grid.centers**2
        * np.sqrt(geometry.a)
        * geometry.b
        / geometry.X**3
        * grid.dr
    )
    np.testing.assert_allclose(weights, expected, rtol=0.0, atol=0.0)
    assert total == pytest.approx(float(np.sum(expected)))
    assert np.all(weights > 0.0)


@pytest.mark.parametrize(
    ("component", "index", "bad_value"),
    [
        ("a", 0, 0.0),
        ("a", 1, -1.0),
        ("a", 2, np.nan),
        ("b", 0, -0.1),
        ("b", 1, np.inf),
        ("X", 0, 0.0),
        ("X", 1, -1.0),
        ("X", 2, np.nan),
    ],
)
def test_invalid_metric_component_fails_closed(component, index, bad_value):
    values = {
        "a": [1.0, 1.2, 1.4],
        "b": [1.0, 0.9, 0.8],
        "X": [1.0, 1.1, 1.2],
    }
    values[component][index] = bad_value
    grid, geometry = _inputs(**values)

    with pytest.raises(
        FloatingPointError,
        match=f"metric component {component} must be finite and positive",
    ):
        cmc_gauge._proper_volume_weights(grid, geometry)


def test_finite_inputs_that_overflow_derived_volume_weight_fail_closed():
    grid, geometry = _inputs(
        centers=[0.5, 1.5, 2.5],
        dr=1.0,
        a=[1.0e308, 1.0e308, 1.0e308],
        b=[1.0e308, 1.0e308, 1.0e308],
        X=[1.0e-100, 1.0e-100, 1.0e-100],
    )
    assert np.all(np.isfinite(geometry.a))
    assert np.all(np.isfinite(geometry.b))
    assert np.all(np.isfinite(geometry.X))

    with pytest.raises(
        FloatingPointError,
        match="non-finite or nonpositive cell weight",
    ):
        cmc_gauge._proper_volume_weights(grid, geometry)


def test_invalid_grid_shape_or_spacing_fails_closed():
    grid, geometry = _inputs()
    geometry.b = geometry.b[:-1]
    with pytest.raises(ValueError, match="must match the radial-centre shape"):
        cmc_gauge._proper_volume_weights(grid, geometry)

    grid, geometry = _inputs(dr=0.0)
    with pytest.raises(FloatingPointError, match="spacing must be finite and positive"):
        cmc_gauge._proper_volume_weights(grid, geometry)

    grid, geometry = _inputs(centers=[0.0, 0.75, 1.25])
    with pytest.raises(FloatingPointError, match="centres must be finite and positive"):
        cmc_gauge._proper_volume_weights(grid, geometry)


def test_off_center_grid_is_rejected_at_its_coordinate_scale():
    # The middle centre is displaced by 0.05, well above coordinate-scale
    # roundoff but small enough that it remains inside its nominal cell.
    grid, geometry = _inputs(centers=[0.25, 0.80, 1.25])
    with pytest.raises(
        FloatingPointError,
        match="radial centres must be uniform and cell-centred",
    ):
        cmc_gauge._proper_volume_weights(grid, geometry)


def test_target_kdot_rejects_invalid_metric_before_rhs_consumers(monkeypatch):
    grid, geometry = _inputs(a=[1.0, -0.25, 1.4])

    def unexpected_rhs(*args, **kwargs):
        pytest.fail("CMC RHS consumers must not receive an invalid metric")

    monkeypatch.setattr(cmc_gauge.adapter, "vendor_modules", unexpected_rhs)
    monkeypatch.setattr(cmc_gauge.adapter, "primary_l3_with_matter", unexpected_rhs)

    with pytest.raises(FloatingPointError, match="metric component a must be finite and positive"):
        cmc_gauge.target_kdot(grid, geometry, scalars=None, matter=None)


def test_target_kdot_uses_validated_weights_and_finite_rhs(monkeypatch):
    grid, geometry = _inputs()
    raw_rhs = np.asarray([1.0, -2.0, 4.0])

    class FakeVacuum:
        @staticmethod
        def primary_l2_rhs(grid, geometry):
            return {"K": np.zeros_like(grid.centers)}

    monkeypatch.setattr(
        cmc_gauge.adapter,
        "vendor_modules",
        lambda: (None, FakeVacuum(), None),
    )
    monkeypatch.setattr(
        cmc_gauge.adapter,
        "primary_l3_with_matter",
        lambda *args, **kwargs: {"K": raw_rhs.copy()},
    )

    weights, total = cmc_gauge._proper_volume_weights(grid, geometry)
    expected = float(np.sum(weights * raw_rhs) / total)
    actual = cmc_gauge.target_kdot(grid, geometry, scalars=None, matter=None)
    assert actual == pytest.approx(expected)


def test_target_kdot_rejects_nonfinite_raw_rhs(monkeypatch):
    grid, geometry = _inputs()

    class FakeVacuum:
        @staticmethod
        def primary_l2_rhs(grid, geometry):
            return {"K": np.zeros_like(grid.centers)}

    monkeypatch.setattr(
        cmc_gauge.adapter,
        "vendor_modules",
        lambda: (None, FakeVacuum(), None),
    )
    monkeypatch.setattr(
        cmc_gauge.adapter,
        "primary_l3_with_matter",
        lambda *args, **kwargs: {"K": np.asarray([0.0, np.nan, 1.0])},
    )

    with pytest.raises(FloatingPointError, match="non-finite or mismatched K RHS"):
        cmc_gauge.target_kdot(grid, geometry, scalars=None, matter=None)
