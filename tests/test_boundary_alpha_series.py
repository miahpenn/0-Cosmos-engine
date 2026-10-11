import json

import numpy as np
import pytest

import engine.production_kernel as production_kernel_module
from diagnostics.curvature_consistency.boundary_alpha_series import (
    BoundaryAlphaSeries,
    R80BoundaryFaceSeries,
)
from engine.cmc_gauge import solve_cmc_lapse
from engine.production_kernel import V55ProductionKernel


def _small_state():
    return V55ProductionKernel().initialize(
        resolution=32,
        r_max=16.0,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )


def test_boundary_alpha_series_metadata_interpolation_and_clamps():
    series = BoundaryAlphaSeries()
    assert series.t.size == 49
    assert series.source_run == 38117001920
    assert series.source_index == 159
    assert series.source_radius == 79.75
    assert np.all(np.diff(series.t) > 0)
    assert series( series.t[0] - 1.0) == series.alpha[0]
    assert series(series.t[-1] + 1.0) == series.alpha[-1]
    mid_t = 0.5 * (series.t[10] + series.t[11])
    assert series(mid_t) == pytest.approx(0.5 * (series.alpha[10] + series.alpha[11]))


def test_cmc_default_none_and_override_boundary_contract():
    state = _small_state()
    args = (state.grid, state.geometry, state.scalars, state.matter)
    alpha_default, _ = solve_cmc_lapse(*args)
    alpha_none, _ = solve_cmc_lapse(*args, outer_alpha=None)
    alpha_override, _ = solve_cmc_lapse(*args, outer_alpha=0.987654321)

    np.testing.assert_array_equal(alpha_default, alpha_none)
    assert alpha_default[-1] == 1.0
    assert alpha_override[-1] == 0.987654321
    assert not np.array_equal(alpha_default[:-1], alpha_override[:-1])


def test_cmc_outer_alpha_rejects_nonfinite_value():
    state = _small_state()
    with pytest.raises(FloatingPointError, match="non-finite"):
        solve_cmc_lapse(
            state.grid, state.geometry, state.scalars, state.matter,
            outer_alpha=float("nan"),
        )


def test_production_default_one_step_remains_finite():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32, r_max=16.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    stepped = kernel.step(state, 0.002)
    assert np.all(np.isfinite(stepped.geometry.alpha))
    assert np.all(stepped.geometry.alpha > 0.0)
    assert stepped.geometry.alpha[-1] == 1.0


def test_r80_boundary_face_series_returns_interpolated_metric_inputs(tmp_path):
    path = tmp_path / "profile.jsonl"
    fields = ("a", "b", "X", "alpha", "beta")
    rows = []
    for t, offset in ((0.0, 0.0), (1.0, 2.0)):
        geometry = {name: [1.0, 10.0 + offset + i] for i, name in enumerate(fields)}
        rows.append({
            "schema": "strong_D_evolved_state_r80_capture_v1",
            "t": t,
            "state_snapshot": {"geometry": geometry},
        })
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    series = R80BoundaryFaceSeries(path)
    midpoint = series(0.5)
    assert set(midpoint) == set(fields)
    assert midpoint["a"] == pytest.approx(11.0)
    assert midpoint["alpha"] == pytest.approx(14.0)
    assert series(-1.0)["a"] == pytest.approx(10.0)
    assert series(2.0)["a"] == pytest.approx(12.0)


def test_stage_time_override_and_frozen_matter_boundary_are_separated(monkeypatch):
    baseline = _small_state()
    fixed_boundary = {
        name: float(getattr(baseline.geometry, name)[-1])
        for name in ("a", "b", "X", "alpha", "beta")
    }
    alpha_times = []
    boundary_times = []
    observed = []
    alpha_value = 0.987654321

    def alpha_provider(t):
        alpha_times.append(float(t))
        return alpha_value + 0.001 * float(t)

    def boundary_provider(t):
        boundary_times.append(float(t))
        return dict(fixed_boundary)

    original_species_rhs = production_kernel_module.species_rhs
    original_metric_derivatives = production_kernel_module.geometry_metric_derivatives
    derivative_alpha_boundaries = []

    def recording_species_rhs(metric, *args, **kwargs):
        observed.append({name: float(getattr(metric, name)[-1]) for name in fixed_boundary})
        return original_species_rhs(metric, *args, **kwargs)

    def recording_metric_derivatives(metric, *args, **kwargs):
        derivative_alpha_boundaries.append(float(metric.alpha[-1]))
        return original_metric_derivatives(metric, *args, **kwargs)

    monkeypatch.setattr(production_kernel_module, "species_rhs", recording_species_rhs)
    monkeypatch.setattr(
        production_kernel_module, "geometry_metric_derivatives",
        recording_metric_derivatives,
    )
    kernel = V55ProductionKernel(
        diagnostic_outer_alpha=alpha_provider,
        diagnostic_matter_boundary=boundary_provider,
    )
    state = kernel.initialize(
        resolution=32, r_max=16.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    assert state.geometry.alpha[-1] == alpha_value
    stepped = kernel.step(state, 0.002)
    assert np.all(np.isfinite(stepped.geometry.alpha))

    assert alpha_times == pytest.approx([0.0, 0.0, 0.002, 0.002, 0.002])
    assert boundary_times == pytest.approx([0.0, 0.002])
    assert len(observed) == 6
    for values in observed:
        for name, expected in fixed_boundary.items():
            assert values[name] == expected
    # Matter metric endpoints are restored, but derivatives continue to use
    # the live stage alpha rather than the baseline endpoint value.
    assert derivative_alpha_boundaries == pytest.approx([
        alpha_value, alpha_value + 0.001 * 0.002,
    ])
    assert all(
        abs(value - fixed_boundary["alpha"]) > 1.0e-6
        for value in derivative_alpha_boundaries
    )
    assert stepped.geometry.alpha[-1] == pytest.approx(alpha_value + 0.001 * 0.002)
