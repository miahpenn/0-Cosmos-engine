"""Tests for the diagnostic-only accepted-step Hamiltonian decomposition."""
import numpy as np

from engine.hamiltonian_evolution_increment_decomposition import (
    GROUPS,
    decompose_increment,
    hamiltonian_residual,
)
from engine.production_kernel import V55ProductionKernel


def _physical_arrays(state):
    result = {}
    for object_name in ("geometry", "scalars", "matter"):
        obj = getattr(state, object_name)
        for name, value in vars(obj).items():
            if isinstance(value, np.ndarray):
                result[f"{object_name}.{name}"] = value.copy()
            elif hasattr(value, "__dict__"):
                for subname, subvalue in vars(value).items():
                    if isinstance(subvalue, np.ndarray):
                        result[f"{object_name}.{name}.{subname}"] = subvalue.copy()
    return result


def test_shapley_groups_are_explicit_and_include_the_required_evolved_fields():
    labels = [name for name, _, _ in GROUPS]
    assert labels == [
        "metric", "lapse_shift_gauge", "A_a", "K", "Lambda", "matter"
    ]
    fields = {(owner, name) for _, owner, names in GROUPS for name in names}
    for name in ("a", "b", "X", "Aa", "K", "Lambda"):
        assert ("geometry", name) in fields
    assert ("matter", "scalars") in fields
    assert ("matter", "matter") in fields


def test_accepted_step_increment_closes_and_diagnostic_does_not_mutate_states():
    kernel = V55ProductionKernel()
    before = kernel.initialize(
        resolution=16,
        r_max=16.0,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    before_arrays = _physical_arrays(before)
    before_h = hamiltonian_residual(before).copy()

    # One ordinary accepted step only; no long trajectory/campaign here.
    after = kernel.step(before, 0.001)
    before_arrays_after_step = _physical_arrays(before)
    after_arrays_before_witness = _physical_arrays(after)
    witness = decompose_increment(before, after)

    # The evolution input remains an unchanged pre-step state.
    assert before_arrays.keys() == before_arrays_after_step.keys()
    for key, expected in before_arrays.items():
        np.testing.assert_array_equal(before_arrays_after_step[key], expected)

    # Evaluating every hybrid state is strictly observational.
    after_arrays_after_witness = _physical_arrays(after)
    assert after_arrays_before_witness.keys() == after_arrays_after_witness.keys()
    for key, expected in after_arrays_before_witness.items():
        np.testing.assert_array_equal(after_arrays_after_witness[key], expected)

    after_h = hamiltonian_residual(after)
    assert witness["finite"]
    assert witness["closure_max_abs_all_cells"] < 1.0e-9

    cells = witness["cells_0_to_4"]
    for cell in cells:
        index = cell["cell"]
        assert cell["H_before"] == pytest_approx(before_h[index])
        assert cell["H_after"] == pytest_approx(after_h[index])
        assert abs(
            sum(cell["component_delta_H"].values()) - cell["delta_H"]
        ) < 1.0e-9
        assert abs(cell["closure_error"]) < 1.0e-9


def pytest_approx(value):
    # Keep numpy's array comparison separate from this simple scalar tolerance.
    return np.testing.assert_allclose if False else _Approx(float(value))


class _Approx:
    """Tiny local scalar comparator so the test has no implicit tolerance guess."""
    def __init__(self, expected):
        self.expected = expected

    def __eq__(self, actual):
        return bool(np.isclose(actual, self.expected, rtol=0.0, atol=1.0e-12))
