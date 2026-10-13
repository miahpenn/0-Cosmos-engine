"""Read-only finite-volume Gauss contract on frozen production metric slices."""
from functools import lru_cache
from types import SimpleNamespace

import numpy as np
import pytest

from diagnostics.research_em_frozen_metric_gauss import (
    solve_frozen_production_metric_gauss,
)
from engine.production_kernel import V55ProductionKernel


@lru_cache(maxsize=3)
def _initialized_production_state(resolution):
    # Initialization only: no ProductionState.step() and no campaign.
    # R0=2 aligns with a face at r_max=8 for each tested resolution.
    return V55ProductionKernel().initialize(
        resolution=resolution,
        r_max=8.0,
        include_radiation=True,
    )


def test_frozen_production_metric_gauss_closes_on_actual_initialized_geometry():
    state = _initialized_production_state(64)
    original_a = state.geometry.a.copy()
    original_b = state.geometry.b.copy()
    original_X = state.geometry.X.copy()

    result = solve_frozen_production_metric_gauss(
        state, rho0=0.75, source_radius=2.0
    )

    assert state.t == 0.0
    assert state.history == []
    np.testing.assert_array_equal(state.geometry.a, original_a)
    np.testing.assert_array_equal(state.geometry.b, original_b)
    np.testing.assert_array_equal(state.geometry.X, original_X)
    # This is the actual initialized metric, not PR #107's manufactured psi.
    assert not np.allclose(
        original_X,
        1.0 / (1.0 + 0.4 * (state.grid.centers / 2.0) ** 2),
        rtol=0.0,
        atol=1.0e-5,
    )
    assert result.total_charge > 0.0
    assert result.faces.shape == (state.grid.n + 1,)
    assert result.cell_charges.shape == (state.grid.n,)
    assert np.all(np.isfinite(result.cell_charges))
    assert np.all(np.isfinite(result.face_electric_field_r))
    assert np.all(np.isfinite(result.face_electric_field_hat))

    scale = max(1.0, abs(result.total_charge))
    assert result.max_abs_cell_residual <= 128.0 * np.finfo(float).eps * scale
    assert result.face_charges[0] == 0.0
    assert result.face_electric_field_r[0] == 0.0
    assert result.face_electric_field_hat[0] == 0.0
    assert result.face_areal_radius[0] == 0.0

    # Independently re-evaluate the production metric identities at faces.
    centers = state.grid.centers
    for index, radius in enumerate(result.faces):
        def face_value(values):
            if radius < centers[0]:
                return float(values[0])
            if radius > centers[-1]:
                slope = (values[-1] - values[-2]) / (centers[-1] - centers[-2])
                return float(values[-1] + slope * (radius - centers[-1]))
            return float(np.interp(radius, centers, values))

        expected_areal = (
            radius * np.sqrt(face_value(state.geometry.b)) / face_value(state.geometry.X)
        )
        assert result.face_areal_radius[index] == pytest.approx(
            expected_areal, rel=2.0e-13, abs=1.0e-14
        )

    # The integrated constraint and the areal flux form agree at every face.
    face_flux = result.face_sqrt_gamma * result.face_electric_field_r
    np.testing.assert_allclose(
        face_flux, result.face_charges, rtol=2.0e-13, atol=2.0e-13 * scale
    )
    np.testing.assert_allclose(
        result.face_areal_radius[1:] ** 2
        * result.face_electric_field_hat[1:],
        result.face_charges[1:],
        rtol=3.0e-13,
        atol=3.0e-13 * scale,
    )


def test_frozen_production_metric_centre_field_is_regular_and_refines():
    errors = []
    for resolution in (32, 64, 128):
        state = _initialized_production_state(resolution)
        result = solve_frozen_production_metric_gauss(
            state, rho0=0.75, source_radius=2.0
        )
        assert result.face_areal_radius[1] > 0.0
        measured_slope = (
            result.face_electric_field_hat[1] / result.face_areal_radius[1]
        )
        central_metric_slope = (
            4.0
            * np.pi
            * 0.75
            / 3.0
            * np.sqrt(state.geometry.a[0] / state.geometry.b[0])
        )
        assert np.isfinite(measured_slope)
        assert np.isfinite(central_metric_slope)
        assert measured_slope > 0.0
        errors.append(abs(measured_slope - central_metric_slope))

    # The centre-cell field must approach the regular O(R_areal) limit;
    # no fitted coefficient or production correction is introduced.
    assert errors[-1] < errors[0]


def test_frozen_production_gauss_quadrature_and_midpoint_estimate_refine():
    errors = []
    for resolution in (32, 64, 128):
        result = solve_frozen_production_metric_gauss(
            _initialized_production_state(resolution),
            rho0=0.75,
            source_radius=2.0,
        )
        error = abs(result.total_charge - result.midpoint_total_charge)
        assert np.isfinite(error)
        assert error > 0.0
        errors.append(error)

    # The centre-sampled proper-volume estimate converges toward the
    # higher-order integral as the *frozen production grid* is refined.
    assert errors[1] < errors[0]
    assert errors[2] < errors[1]


def test_frozen_production_gauss_operator_detects_a_perturbed_face_field():
    result = solve_frozen_production_metric_gauss(
        _initialized_production_state(64), rho0=0.75, source_radius=2.0
    )
    cell_index = int(round(0.5 / (result.faces[1] - result.faces[0]))) - 1
    altered = result.face_electric_field_r.copy()
    altered[cell_index + 1] += 1.0e-4 * max(
        1.0, abs(altered[cell_index + 1])
    )
    altered_flux = result.face_sqrt_gamma * altered
    residual = (
        altered_flux[cell_index + 1]
        - altered_flux[cell_index]
        - result.cell_charges[cell_index]
    )
    assert abs(residual) > 1.0e-8


def test_frozen_production_gauss_rejects_misaligned_source_edge():
    with pytest.raises(ValueError, match="source_radius must align"):
        solve_frozen_production_metric_gauss(
            _initialized_production_state(32),
            rho0=0.75,
            source_radius=2.1,
        )


def test_frozen_production_gauss_rejects_nonfinite_metric():
    state = _initialized_production_state(32)
    geometry = SimpleNamespace(
        a=state.geometry.a.copy(),
        b=state.geometry.b.copy(),
        X=state.geometry.X.copy(),
    )
    geometry.X[3] = np.nan
    invalid_state = SimpleNamespace(grid=state.grid, geometry=geometry)

    with pytest.raises(ValueError, match="X must be finite"):
        solve_frozen_production_metric_gauss(
            invalid_state, rho0=0.75, source_radius=2.0
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"rho0": np.nan, "source_radius": 2.0},
        {"rho0": 0.75, "source_radius": np.inf},
        {"rho0": 0.75, "source_radius": 0.0},
        {"rho0": True, "source_radius": 2.0},
        {"rho0": 0.75, "source_radius": 2.0, "quadrature_order": 1},
        {"rho0": 0.75, "source_radius": 2.0, "quadrature_order": True},
    ],
)
def test_frozen_production_gauss_rejects_invalid_inputs(kwargs):
    with pytest.raises(ValueError):
        solve_frozen_production_metric_gauss(
            _initialized_production_state(32), **kwargs
        )
