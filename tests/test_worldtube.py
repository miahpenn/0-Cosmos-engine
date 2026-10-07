import numpy as np

from engine.worldtube import current_residual


def test_current_residual_excludes_one_sided_time_endpoints():
    t = np.linspace(0.0, 1.0, 6)
    m = t**2
    rhs = 2.0 * t
    residual = current_residual(t, m, rhs)

    assert np.isnan(residual[0])
    assert np.isnan(residual[-1])
    assert np.max(np.abs(residual[1:-1])) < 1.0e-12


def test_current_residual_rejects_shape_mismatch():
    t = np.arange(4, dtype=float)
    m = np.arange(3, dtype=float)
    rhs = np.arange(4, dtype=float)
    try:
        current_residual(t, m, rhs)
    except ValueError:
        return
    raise AssertionError("shape mismatch must be rejected")
