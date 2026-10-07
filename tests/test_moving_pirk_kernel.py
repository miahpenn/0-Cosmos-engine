import numpy as np

from engine.moving_pirk_kernel import V55MovingPIRKKernel


def test_moving_kernel_uses_native_gauge_without_cmc():
    kernel = V55MovingPIRKKernel()

    def forbidden(*args, **kwargs):
        raise AssertionError("moving-gauge kernel must not solve CMC lapse")

    kernel._solve_lapse = forbidden
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )

    assert np.allclose(state.geometry.alpha, 1.0)
    assert np.allclose(state.geometry.beta, 0.0)
    assert np.all(np.isfinite(state.geometry.B))
    assert np.all(np.isfinite(state.scalars.D))

    b0 = state.geometry.B.copy()
    lam0 = state.geometry.Lambda.copy()
    stepped = kernel.step(state, 0.002)

    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        assert np.all(np.isfinite(getattr(stepped.geometry, name)))

    assert np.all(stepped.geometry.a > 0.0)
    assert np.all(stepped.geometry.b > 0.0)
    assert np.all(stepped.geometry.X > 0.0)
    assert np.all(stepped.geometry.alpha > 0.0)

    expected_B = b0 + 0.75 * (stepped.geometry.Lambda - lam0)
    assert np.allclose(stepped.geometry.B, expected_B)


def test_moving_kernel_preserves_radiation_physical_state_on_one_step():
    kernel = V55MovingPIRKKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=16.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    stepped = kernel.step(state, 0.002)

    metric = __import__(
        "engine.v55_matter",
        fromlist=["metric_slice_from_q"],
    ).metric_slice_from_q(stepped.grid, stepped.geometry)

    species = __import__(
        "engine.matter_system",
        fromlist=["Species", "primitives"],
    )
    prim = species.primitives(
        species._metric_arrays(metric),
        stepped.matter.radiation,
        species.Species.RADIATION,
    )
    for q in prim:
        assert q.rho >= 0.0
        assert q.pressure >= 0.0
        assert q.gamma_rr * q.v_r * q.v_r < 1.0
