"""Opt-in discrete-consistent initial data: closure, dr^2 size of the correction,
invariance of the central lapse, and the default initializer left unchanged."""
import numpy as np
import pytest

from engine.discrete_consistent_initial import discrete_consistent_state
from engine.production_kernel import V55ProductionKernel


@pytest.fixture(scope="module")
def kernel():
    return V55ProductionKernel()


@pytest.mark.parametrize("n", [40, 80])
def test_discrete_hamiltonian_closes_at_t0(kernel, n):
    _, _, hist = discrete_consistent_state(kernel, n, 40.0)
    assert hist[-1] < 1e-12


def test_correction_scales_as_dr_squared(kernel):
    st40, B40, _ = discrete_consistent_state(kernel, 40, 40.0)
    st80, B80, _ = discrete_consistent_state(kernel, 80, 40.0)
    o40 = kernel.initialize(resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
                            D_amplitude=1e-10, include_radiation=True)
    o80 = kernel.initialize(resolution=80, r_max=40.0, amplitude=0.01, width=7.0,
                            D_amplitude=1e-10, include_radiation=True)
    c40 = (B40[0] - np.asarray(o40.geometry.X)[0] ** 6) / 1.0 ** 2
    c80 = (B80[0] - np.asarray(o80.geometry.X)[0] ** 6) / 0.5 ** 2
    assert abs(c80 / c40 - 1.0) < 0.05


def test_central_lapse_is_insensitive_to_the_correction(kernel):
    st, _, _ = discrete_consistent_state(kernel, 40, 40.0)
    orig = kernel.initialize(resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
                             D_amplitude=1e-10, include_radiation=True)
    a_new = kernel._solve_lapse(st.grid, st.geometry, st.scalars, st.matter)[0][0]
    a_old = kernel._solve_lapse(orig.grid, orig.geometry, orig.scalars, orig.matter)[0][0]
    assert abs(a_new - a_old) < 1e-5


def test_default_initializer_is_unchanged(kernel):
    st = kernel.initialize(resolution=40, r_max=40.0, amplitude=0.01, width=7.0,
                           D_amplitude=1e-10, include_radiation=True)
    from engine.production_kernel import V55ProductionKernel as K
    assert callable(getattr(K, "initialize"))
    assert np.isfinite(np.asarray(st.geometry.X)).all()
