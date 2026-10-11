import numpy as np
import pytest

from diagnostics.curvature_consistency.boundary_alpha_series import BoundaryAlphaSeries
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
