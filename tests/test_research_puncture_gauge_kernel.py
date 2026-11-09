"""research/puncture_gauge_kernel_v1.py reproduces the vendor damped moving-puncture step in vacuum."""
import dataclasses

import numpy as np
import pytest

from engine import v55_pirk_adapter as adapter


def _vendor_present():
    return adapter.VENDOR.exists() and any(adapter.VENDOR.iterdir())


pytestmark = pytest.mark.skipif(not _vendor_present(), reason="pinned vendor kernel not initialised")


@pytest.mark.parametrize("ko_epsilon", [0.0, 0.1])
def test_vacuum_schwarzschild_matches_the_vendor_step_bit_for_bit(ko_epsilon, monkeypatch):
    import engine.scalar_system as ss
    from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
    from engine.v55_matter import V55MatterState, metric_slice_from_q
    from research.puncture_gauge_kernel_v1 import PunctureGaugeKernel

    monkeypatch.setattr(ss, "V0", 0.0)
    _, vac, moving = adapter.vendor_modules()
    kernel = PunctureGaugeKernel(eta=2.0, ko_epsilon=ko_epsilon)
    state = kernel.initialize(resolution=60, r_max=30.0, include_radiation=False)
    geometry = vac.flat_state(state.grid)
    r = np.asarray(state.grid.centers)
    geometry.X = (1.0 + 1.0 / (2.0 * r)) ** -2
    geometry.alpha.fill(1.0)
    z = np.zeros_like(r)
    tiny = np.full_like(r, 1e-30)
    metrics = _metric_arrays(metric_slice_from_q(state.grid, geometry))
    state = dataclasses.replace(
        state, geometry=geometry,
        scalars=dataclasses.replace(state.scalars, S=z, PS=z, D=z, PD=z, phi=z, Pi=z),
        matter=V55MatterState(initialize_dust(metrics, tiny), initialize_dust(metrics, tiny),
                              initialize_radiation(metrics, tiny), None),
    )
    reference = geometry.copy()
    dt = 0.1 * float(state.grid.dr)
    for _ in range(20):
        state = kernel.step(state, dt)
        reference = moving.damped_moving_puncture_pirk2_step(
            state.grid, reference, dt, eta=2.0, ko_epsilon=ko_epsilon)
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        np.testing.assert_array_equal(getattr(state.geometry, name), getattr(reference, name))
    assert float(np.max(np.abs(state.geometry.beta))) > 1e-3   # the shift is live
    assert state.tau > 0.0


def test_rejects_the_radial_maxwell_sector():
    from research.puncture_gauge_kernel_v1 import PunctureGaugeKernel

    kernel = PunctureGaugeKernel()
    state = kernel.initialize(resolution=24, r_max=20.0)
    state = dataclasses.replace(state, electromagnetic=object())
    with pytest.raises(NotImplementedError):
        kernel.step(state, 0.01)
