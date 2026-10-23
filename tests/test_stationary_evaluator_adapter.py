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
