"""Regression gates for the diagnostic completed-state radiation budget."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from engine.matter_system import BSSNMetricSlice, ConservedSpecies, _metric_arrays
from engine.radiation_geometry_stage_trace import RadiationStageTrace


def _geometry(n: int, scale: float, alpha: float) -> BSSNMetricSlice:
    r = 0.5 + np.arange(n, dtype=float)
    ones, zeros = np.ones(n), np.zeros(n)
    return BSSNMetricSlice(
        r=r, a=scale * ones.copy(), b=ones.copy(), X=ones.copy(),
        alpha=alpha * ones.copy(), beta=zeros.copy(), Aa=zeros.copy(),
        K=zeros.copy(), Lambda=zeros.copy(), B=zeros.copy(),
    )


def _radiation(energy: np.ndarray, momentum: np.ndarray) -> ConservedSpecies:
    return ConservedSpecies(
        np.zeros_like(energy, dtype=float),
        np.asarray(energy, dtype=float).copy(),
        np.asarray(momentum, dtype=float).copy(),
    )


def _metric_helper(grid, geometry):
    metrics = _metric_arrays(geometry)
    return (
        geometry, metrics,
        np.asarray([metric.sqrt_gamma for metric in metrics], dtype=float),
        np.asarray([metric.gamma_rr_inv for metric in metrics], dtype=float),
    )


def test_completed_state_budget_closes_trapezoid_and_metric_transition(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(RadiationStageTrace, "_patch_production_kernel", lambda self: None)
    trace = RadiationStageTrace(tmp_path)
    trace.geometry_and_metrics = _metric_helper

    n, dt = 8, 0.025
    accepted_geometry = _geometry(n, scale=1.00, alpha=0.92)
    rhs1_geometry = _geometry(n, scale=0.93, alpha=0.81)
    completed_geometry = _geometry(n, scale=0.88, alpha=0.76)
    grid = SimpleNamespace(centers=accepted_geometry.r.copy())
    U0e, U0r = np.linspace(1.0e-4, 1.7e-4, n), np.linspace(0.25e-4, 0.43e-4, n)
    U0 = _radiation(U0e, U0r)
    rhs0e, rhs0r = np.linspace(-1.0e-6, 2.0e-6, n), np.linspace(0.3e-6, -0.7e-6, n)
    rhs1e, rhs1r = np.linspace(1.5e-6, -0.6e-6, n), np.linspace(-0.8e-6, 0.9e-6, n)
    U1 = _radiation(U0e + dt * rhs0e, U0r + dt * rhs0r)
    Uc = _radiation(
        U0e + 0.5 * dt * (rhs0e + rhs1e),
        U0r + 0.5 * dt * (rhs0r + rhs1r),
    )
    outer = range(n - 5, n)
    sources0 = {i: (float(0.17 * rhs0e[i]), float(-0.23 * rhs0r[i])) for i in outer}
    sources1 = {i: (float(-0.11 * rhs1e[i]), float(0.31 * rhs1r[i])) for i in outer}
    rhs0 = {
        "geometry": accepted_geometry, "matter_radiation": U0,
        "rhs_energy": rhs0e, "rhs_momentum": rhs0r, "sources": sources0,
    }
    rhs1 = {
        "geometry": rhs1_geometry, "matter_radiation": U1,
        "rhs_energy": rhs1e, "rhs_momentum": rhs1r, "sources": sources1,
    }
    budget = trace.completed_state_budget(
        grid, completed_geometry, SimpleNamespace(radiation=Uc),
        rhs0, rhs1, dt, t=12.0, step_index=75,
        stage="completed_state_outer_light_boundary",
    )

    assert budget["t_accepted"] == pytest.approx(12.0)
    assert budget["t_completed_candidate"] == pytest.approx(12.0 + dt)
    assert budget["source_capture_missing_cells"] == []
    for key in (
        "max_abs_budget_closure_error",
        "max_abs_predictor_U_E_reconstruction_error",
        "max_abs_predictor_U_r_reconstruction_error",
        "max_abs_completed_U_E_reconstruction_error",
        "max_abs_completed_U_r_reconstruction_error",
        "max_abs_source_flux_split_energy_error",
        "max_abs_source_flux_split_momentum_error",
        "max_abs_predictor_to_completed_closure_error",
    ):
        assert budget[key] == pytest.approx(0.0, abs=1e-18), key

    for row in budget["cells"]:
        assert row["source_capture_missing"] is False
        assert row["U_E_completed_actual"] == pytest.approx(
            row["U_E_completed_reconstructed_from_trapezoid"], abs=1e-18
        )
        assert row["U_r_completed_actual"] == pytest.approx(
            row["U_r_completed_reconstructed_from_trapezoid"], abs=1e-18
        )
        assert row["C_completed_final_metric"] - row["C_accepted"] == pytest.approx(
            row["delta_C_flux_trapezoidal"]
            + row["delta_C_sources_trapezoidal"]
            + row["delta_C_metric_accepted_to_rhs1"]
            + row["delta_C_metric_rhs1_to_completed"], abs=1e-18
        )
        assert row["C_completed_final_metric"] - row["C_predictor_on_rhs1_metric"] == pytest.approx(
            row["delta_C_trapezoid_matter_correction_on_rhs1_metric"]
            + row["delta_C_metric_rhs1_to_completed"], abs=1e-18
        )
    assert budget["cells"][-1]["accepted_metric"]["a"] != budget["cells"][-1]["rhs1_stage_metric"]["a"]
    assert budget["cells"][-1]["rhs1_stage_metric"]["a"] != budget["cells"][-1]["completed_metric"]["a"]


def test_tracer_records_rhs1_rates_and_source_provenance(tmp_path, monkeypatch):
    import engine.production_kernel as pk
    import engine.matter_rhs as mr
    import engine.matter_system as ms

    class Copyable:
        def __init__(self, **kwargs):
            vars(self).update(kwargs)
        def copy(self):
            return Copyable(**{
                name: value.copy() if isinstance(value, np.ndarray) else value
                for name, value in vars(self).items()
            })

    n = 8
    geom = Copyable(**{
        "r": np.arange(n, dtype=float) + 0.5, "a": np.ones(n),
        "b": np.ones(n), "X": np.ones(n), "alpha": np.ones(n),
        "beta": np.zeros(n), "Aa": np.zeros(n), "K": np.zeros(n),
        "Lambda": np.zeros(n), "B": np.zeros(n),
    })
    U = Copyable(rest=np.zeros(n), energy_t=np.ones(n), momentum_r=np.zeros(n))
    rhs_radiation = SimpleNamespace(
        energy_t=np.linspace(0.1, 0.8, n), momentum_r=np.linspace(-0.2, 0.2, n)
    )
    source_map = {i: (0.01 * i, -0.02 * i) for i in range(n - 5, n)}
    def fake_rhs(kernel_self, state, radiation_recovery_metric=None):
        return object(), {"radiation": rhs_radiation}, object()

    monkeypatch.setattr(pk.V55ProductionKernel, "_rhs", fake_rhs)
    trace = RadiationStageTrace(tmp_path)
    trace.record_stage = lambda *args, **kwargs: None
    trace.last_radiation_sources = source_map
    trace.current_step = {
        "t0": 1.0, "tau0": 0.5, "dt": 0.01,
        "solve_count": 0, "rhs_count": 1,
        "rhs0_record": None, "rhs1_record": None, "failing_stage": None,
    }
    state = SimpleNamespace(
        t=1.0, tau=0.5,
        grid=SimpleNamespace(centers=geom.r.copy()), geometry=geom,
        matter=SimpleNamespace(radiation=U),
    )
    try:
        pk.V55ProductionKernel()._rhs(state)
        record = trace.current_step["rhs1_record"]
        assert record is not None
        assert np.array_equal(record["rhs_energy"], rhs_radiation.energy_t)
        assert np.array_equal(record["rhs_momentum"], rhs_radiation.momentum_r)
        assert record["sources"] == source_map
        assert np.array_equal(record["matter_radiation"].energy_t, U.energy_t)
    finally:
        kernel_cls, old_solve, old_rhs, old_step, old_evolve, old_source, old_boundary = trace._restore_hooks
        kernel_cls._solve_lapse = staticmethod(old_solve)
        kernel_cls._rhs = old_rhs
        kernel_cls._apply_outer_light_boundary = staticmethod(old_boundary)
        kernel_cls.step = old_step
        mr.evolve_species = old_evolve
        ms._radiation_source = old_source


def test_trace_instrumentation_preserves_one_step_production_state_bitwise(tmp_path):
    """Observational tracing must not perturb a representative production step."""
    import engine.matter_rhs as mr
    import engine.matter_system as ms
    from engine.production_kernel import V55ProductionKernel

    kernel = V55ProductionKernel()
    state_plain = kernel.initialize(
        resolution=32, r_max=16.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )
    state_traced = kernel.initialize(
        resolution=32, r_max=16.0, amplitude=0.01, width=7.0,
        D_amplitude=1.0e-10, include_radiation=True,
    )

    dt = 0.002
    plain = kernel.step(state_plain, dt)
    trace = RadiationStageTrace(tmp_path)
    try:
        traced = kernel.step(state_traced, dt)
        assert not trace.instrumentation_errors
    finally:
        kernel_cls, old_solve, old_rhs, old_step, old_evolve, old_source, old_boundary = trace._restore_hooks
        kernel_cls._solve_lapse = staticmethod(old_solve)
        kernel_cls._rhs = old_rhs
        kernel_cls._apply_outer_light_boundary = staticmethod(old_boundary)
        kernel_cls.step = old_step
        mr.evolve_species = old_evolve
        ms._radiation_source = old_source

    assert plain.t == traced.t
    assert plain.tau == traced.tau
    assert plain.e_folds == traced.e_folds

    geometry_fields = (
        "a", "b", "X", "alpha", "beta", "Aa", "K",
        "Lambda", "B",
    )
    for name in geometry_fields:
        assert np.array_equal(
            getattr(plain.geometry, name), getattr(traced.geometry, name)
        ), f"traced production geometry differs in {name}"

    for name, value in vars(plain.scalars).items():
        if isinstance(value, np.ndarray):
            assert np.array_equal(value, getattr(traced.scalars, name)), (
                f"traced scalar state differs in {name}"
            )

    for species in ("dark_matter", "baryons", "radiation"):
        a, b = getattr(plain.matter, species), getattr(traced.matter, species)
        for name in ("rest", "energy_t", "momentum_r"):
            assert np.array_equal(getattr(a, name), getattr(b, name)), (
                f"traced {species}.{name} differs"
            )
