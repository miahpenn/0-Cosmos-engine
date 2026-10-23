"""Unit tests for the research-only stationary evaluator adapter."""
from types import SimpleNamespace

import numpy as np
import pytest

from research.stationary_evaluator_adapter import StationaryEvaluatorAdapter


def test_uses_caller_grid_without_mutating_it():
    r = np.array([0.0, 0.2, 1.1, 3.7, 8.0, 24.0])
    before = r.copy()
    state = StationaryEvaluatorAdapter().evaluate_slice(
        r, tau_local=5.0, t_cosmos=20.0
    )
    np.testing.assert_array_equal(state["r"], before)
    np.testing.assert_array_equal(r, before)
    assert len(state["rho"]) == len(r)
    assert state["model_status"] == "illustrative_profile_not_field_equation_solution"


def test_query_lanes_are_independently_indexed():
    adapter = StationaryEvaluatorAdapter()
    r = np.linspace(0.0, 24.0, 24)
    base = adapter.evaluate_slice(r, tau_local=2.0, t_cosmos=10.0)
    tau_changed = adapter.evaluate_slice(r, tau_local=4.0, t_cosmos=10.0)
    cosmos_changed = adapter.evaluate_slice(r, tau_local=2.0, t_cosmos=12.0)

    # Local profile scaling changes effective radius, not the cosmic shell index.
    assert not np.array_equal(base["effective_r"], tau_changed["effective_r"])
    assert base["shell_center_effective"] == tau_changed["shell_center_effective"]
    # Cosmic query changes shell projection, not local profile scaling.
    np.testing.assert_array_equal(base["effective_r"], cosmos_changed["effective_r"])
    assert base["shell_center_effective"] != cosmos_changed["shell_center_effective"]


def test_metric_slice_adapter_reads_existing_grid_and_lapse():
    r = np.array([0.0, 0.5, 2.0, 7.0])
    alpha = np.array([0.9, 0.92, 0.97, 1.0])
    metric = SimpleNamespace(r=r, alpha=alpha)
    state = StationaryEvaluatorAdapter().evaluate_metric_slice(
        metric, tau_local=1.0, t_cosmos=3.0
    )
    np.testing.assert_array_equal(state["r"], r)
    np.testing.assert_array_equal(state["lapse"], alpha)
    np.testing.assert_array_equal(metric.r, r)
    np.testing.assert_array_equal(metric.alpha, alpha)


def test_adapter_accepts_real_bssn_metric_slice_without_mutating_metric():
    # Exercise the repository's real metric contract rather than only a mock.
    from engine.matter_system import BSSNMetricSlice

    r = np.array([0.0, 0.5, 1.5, 4.0])
    arrays = {
        "a": np.array([1.0, 1.01, 1.02, 1.03]),
        "b": np.array([1.0, 0.99, 0.98, 0.97]),
        "X": np.array([1.0, 0.98, 0.96, 0.94]),
        "alpha": np.array([1.0, 0.97, 0.93, 0.90]),
        "beta": np.array([0.0, 0.01, 0.02, 0.03]),
        "Aa": np.array([0.0, 0.001, 0.002, 0.003]),
        "K": np.array([-0.01, -0.011, -0.012, -0.013]),
        "Lambda": np.array([0.0, 0.001, 0.002, 0.003]),
        "B": np.array([0.0, 0.0, 0.001, 0.001]),
    }
    before = {name: values.copy() for name, values in arrays.items()}
    metric = BSSNMetricSlice(r=r.copy(), **arrays)

    state = StationaryEvaluatorAdapter().evaluate_metric_slice(
        metric, tau_local=2.0, t_cosmos=4.0
    )

    np.testing.assert_array_equal(state["r"], r)
    np.testing.assert_array_equal(state["lapse"], before["alpha"])
    for name, expected in before.items():
        np.testing.assert_array_equal(getattr(metric, name), expected)
    assert state["model_status"] == "illustrative_profile_not_field_equation_solution"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tau_local": float("nan"), "t_cosmos": 1.0},
        {"tau_local": -20.0, "t_cosmos": 1.0},
        {"tau_local": 1.0, "t_cosmos": -6.0},
    ],
)
def test_rejects_invalid_query_domains(kwargs):
    with pytest.raises(ValueError):
        StationaryEvaluatorAdapter().evaluate_slice(np.array([0.0, 1.0]), **kwargs)


def test_rejects_nonmonotone_grid_and_mismatched_lapse():
    adapter = StationaryEvaluatorAdapter()
    with pytest.raises(ValueError, match="strictly increasing"):
        adapter.evaluate_slice([0.0, 2.0, 1.0], tau_local=0.0, t_cosmos=1.0)
    with pytest.raises(ValueError, match="same shape"):
        adapter.evaluate_slice(
            [0.0, 1.0], tau_local=0.0, t_cosmos=1.0, lapse=[1.0]
        )
