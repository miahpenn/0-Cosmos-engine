"""Tests for the research-only closed-loop recirculation probe."""
import numpy as np
import pytest

from research.closed_loop_recirculation_probe import StationaryRecirculationProbe


@pytest.mark.parametrize("t_cosmos", [0.0, 5.0, 10.0, 15.0])
def test_weighted_boundary_tail_integral_converges_with_resolution(t_cosmos):
    coarse = StationaryRecirculationProbe(resolution=48).evaluate_slice(
        tau_local=0.0, t_cosmos=t_cosmos
    )
    fine = StationaryRecirculationProbe(resolution=96).evaluate_slice(
        tau_local=0.0, t_cosmos=t_cosmos
    )
    assert coarse["boundary_tail_integral"] >= 0.0
    assert fine["boundary_tail_integral"] >= 0.0
    assert abs(coarse["boundary_tail_integral"] - fine["boundary_tail_integral"]) / max(
        fine["boundary_tail_integral"], 1e-12
    ) < 0.01


def test_late_shell_tail_is_small_with_bounded_absolute_grid_error():
    coarse = StationaryRecirculationProbe(resolution=48).evaluate_slice(
        tau_local=0.0, t_cosmos=20.0
    )
    fine = StationaryRecirculationProbe(resolution=96).evaluate_slice(
        tau_local=0.0, t_cosmos=20.0
    )
    # The relative discrepancy grows when the tail is nearly zero; bound the
    # absolute error instead of misrepresenting this as 1% relative convergence.
    assert fine["boundary_tail_integral"] < 0.01
    assert abs(
        coarse["boundary_tail_integral"] - fine["boundary_tail_integral"]
    ) < 5e-4


def test_recirc_term_scales_with_fixed_feedback_gain():
    base = StationaryRecirculationProbe(resolution=48, feedback_gain=0.15)
    doubled = StationaryRecirculationProbe(resolution=48, feedback_gain=0.30)
    a = base.evaluate_slice(tau_local=0.0, t_cosmos=10.0)
    b = doubled.evaluate_slice(tau_local=0.0, t_cosmos=10.0)
    np.testing.assert_allclose(b["recirculation"], 2.0 * a["recirculation"])
    assert a["boundary_tail_integral"] == pytest.approx(b["boundary_tail_integral"])


def test_grid_and_outputs_are_resolution_sized_and_labelled_as_toy():
    state = StationaryRecirculationProbe(resolution=24).evaluate_slice(
        tau_local=5.0, t_cosmos=10.0
    )
    assert state["r"].shape == (24,)
    assert state["rho"].shape == (24,)
    assert state["lapse_synthetic"].shape == (24,)
    assert state["boundary_measure"] == (
        "spatial_integral_of_shell_density_tail_not_conserved_flux"
    )
    assert state["model_status"] == (
        "illustrative_recirculation_probe_not_field_equation_solution"
    )
    assert "lapse" not in state


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tau_local": float("nan"), "t_cosmos": 1.0},
        {"tau_local": -20.0, "t_cosmos": 1.0},
        {"tau_local": 1.0, "t_cosmos": -6.0},
        {"tau_local": 1.0, "t_cosmos": float("inf")},
    ],
)
def test_invalid_query_domains_fail_closed(kwargs):
    with pytest.raises(ValueError):
        StationaryRecirculationProbe().evaluate_slice(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"resolution": 1},
        {"resolution": 2.5},
        {"resolution": 24.0},
        {"r_min": 24.0, "r_max": 0.0},
        {"core_scale": 0.0},
        {"shell_width": -1.0},
        {"feedback_gain": float("nan")},
    ],
)
def test_invalid_configuration_fails_closed(kwargs):
    with pytest.raises(ValueError):
        StationaryRecirculationProbe(**kwargs)
