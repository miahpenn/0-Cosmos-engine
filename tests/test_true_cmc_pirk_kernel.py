import numpy as np

from engine.true_cmc_pirk_kernel import (
    V55TrueCMCPIRKKernel,
    solve_archive_cmc_lapse,
)
from engine.v55_matter import metric_slice_from_q
from engine.matter_system import _metric_arrays, primitives, Species


def test_true_cmc_kernel_solves_positive_stage_lapse():
    kernel = V55TrueCMCPIRKKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )

    assert np.all(np.isfinite(state.geometry.alpha))
    assert np.all(state.geometry.alpha > 0.0)
    assert np.isclose(state.geometry.alpha[-1], 1.0)
    assert np.allclose(state.geometry.beta, 0.0)
    assert np.allclose(state.geometry.B, 0.0)


def test_true_cmc_kernel_one_step_preserves_stage_domain():
    kernel = V55TrueCMCPIRKKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    stepped = kernel.step(state, 0.002)

    for name in (
        "a", "b", "X", "alpha", "beta",
        "Aa", "K", "Lambda", "B",
    ):
        values = getattr(stepped.geometry, name)
        assert np.all(np.isfinite(values))

    assert np.all(stepped.geometry.a > 0.0)
    assert np.all(stepped.geometry.b > 0.0)
    assert np.all(stepped.geometry.X > 0.0)
    assert np.all(stepped.geometry.alpha > 0.0)
    assert np.allclose(stepped.geometry.beta, 0.0)
    assert np.allclose(stepped.geometry.B, 0.0)

    metric = metric_slice_from_q(stepped.grid, stepped.geometry)
    metrics = _metric_arrays(metric)
    prim = primitives(
        metrics,
        stepped.matter.radiation,
        Species.RADIATION,
    )
    for q in prim:
        assert q.rho >= 0.0
        assert q.pressure >= 0.0
        assert q.gamma_rr * q.v_r * q.v_r < 1.0


def test_true_cmc_lapse_outer_normalization_is_exact():
    kernel = V55TrueCMCPIRKKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    alpha, _ = solve_archive_cmc_lapse(
        state.grid,
        state.geometry,
        state.scalars,
        state.matter,
    )
    assert np.all(np.isfinite(alpha))
    assert np.all(alpha > 0.0)
    assert np.isclose(alpha[-1], 1.0)
