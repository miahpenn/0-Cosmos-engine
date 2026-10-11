import numpy as np

from engine.production_kernel import V55ProductionKernel


def _quadratic_extrapolate_first(r, values):
    x = r**2
    target = x[0]
    nodes = x[1:4]
    samples = values[1:4]
    result = 0.0
    for j in range(3):
        weight = 1.0
        for k in range(3):
            if j != k:
                weight *= (target - nodes[k]) / (nodes[j] - nodes[k])
        result += weight * samples[j]
    return float(result)


def _center_regularity_errors(state):
    r = np.asarray(state.grid.centers)
    errors = {}

    metric_q = (1.0 - state.geometry.a / state.geometry.b) / r**2
    aa_q = state.geometry.Aa / r**2
    lambda_q = state.geometry.Lambda / r

    errors["metric"] = abs(
        metric_q[0] - _quadratic_extrapolate_first(r, metric_q)
    )
    errors["Aa"] = abs(
        aa_q[0] - _quadratic_extrapolate_first(r, aa_q)
    )
    errors["Lambda"] = abs(
        lambda_q[0] - _quadratic_extrapolate_first(r, lambda_q)
    )
    errors["determinant"] = float(
        np.max(np.abs(state.geometry.a * state.geometry.b**2 - 1.0))
    )
    return errors


def test_production_initialize_is_center_regular():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32, r_max=16.0, D_amplitude=1.0e-10
    )
    errors = _center_regularity_errors(state)
    assert errors["metric"] < 1.0e-12, errors
    assert errors["Aa"] < 1.0e-12, errors
    assert errors["Lambda"] < 1.0e-12, errors
    assert errors["determinant"] < 1.0e-12, errors


def test_production_step_preserves_center_regularity():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32, r_max=16.0, D_amplitude=1.0e-10
    )
    state = kernel.step(state, 2.0e-3)
    errors = _center_regularity_errors(state)
    assert errors["metric"] < 1.0e-12, errors
    assert errors["Aa"] < 1.0e-12, errors
    assert errors["Lambda"] < 1.0e-12, errors
    assert errors["determinant"] < 1.0e-12, errors
