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

def test_current_residual_excludes_stencil_touching_tiny_terminal_step():
    # Simulate accumulated time landing a few picoseconds-in-code-time short
    # of the requested endpoint, followed by the permitted roundoff-sized step.
    regular_t = np.linspace(0.0, 1.0, 6)
    t = np.concatenate([regular_t, [1.0 + 5.0e-12]])
    m = t**2
    rhs = 2.0 * t

    residual = current_residual(t, m, rhs)

    assert np.isnan(residual[-2])
    assert np.isnan(residual[-1])
    assert np.max(np.abs(residual[1:-2])) < 1.0e-12

