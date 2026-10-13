"""E2 gate: one optional normalized EM stress projection reaches shared consumers.

These tests use a manufactured source already in Einstein/CMC consumer units.
They do not select a field normalization, charge carrier/current law, or add
Maxwell state to ProductionState. The default runtime projection remains zero.
"""
import importlib
import math

import numpy as np
import pytest

from engine import cmc_gauge, v55_matter
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_pirk_adapter import (
    lambda_l3_with_matter,
    primary_l3_with_matter,
)


def _fixture_state():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=32.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    return kernel, state


def _pure_radial_e_projection(grid):
    """Manufactured consumer-unit stress (u, -u, +u, 0), not physical data."""
    r = np.asarray(grid.centers, dtype=float)
    u = 1.0e-6 * (1.0 + r / max(float(r[-1]), 1.0))
    return {
        "rho": u.copy(),
        "pr": -u.copy(),
        "pt": u.copy(),
        "j": np.zeros_like(u),
    }


def _assert_projection_equal(got, expected):
    for key in ("rho", "pr", "pt", "j"):
        np.testing.assert_array_equal(got[key], expected[key])


def test_optional_em_projection_is_zero_by_default_and_added_once():
    _kernel, state = _fixture_state()
    baseline = v55_matter.total_matter_projection(
        state.grid, state.geometry, state.scalars, state.matter
    )
    explicit_zero = {
        key: np.zeros_like(baseline[key]) for key in ("rho", "pr", "pt", "j")
    }
    zero_result = v55_matter.total_matter_projection(
        state.grid, state.geometry, state.scalars, state.matter,
        electromagnetic_projection=explicit_zero,
    )
    _assert_projection_equal(zero_result, baseline)
    for key in ("rho", "pr", "pt", "j"):
        np.testing.assert_array_equal(
            zero_result["electromagnetic"][key], explicit_zero[key]
        )

    em = _pure_radial_e_projection(state.grid)
    em_original = {key: value.copy() for key, value in em.items()}
    candidate = v55_matter.total_matter_projection(
        state.grid, state.geometry, state.scalars, state.matter,
        electromagnetic_projection=em,
    )
    for key in ("rho", "pr", "pt", "j"):
        np.testing.assert_array_equal(
            candidate[key], baseline[key] + em[key]
        )
        np.testing.assert_array_equal(candidate["scalar"][key], baseline["scalar"][key])
        np.testing.assert_array_equal(candidate["fluid"][key], baseline["fluid"][key])
        np.testing.assert_array_equal(candidate["electromagnetic"][key], em[key])
        np.testing.assert_array_equal(em[key], em_original[key])

    # The independent diagnostics/boundary object must expose the same sum,
    # with sector components counted exactly once.
    assembled = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter,
        electromagnetic_projection=em,
    )
    assert assembled.finite()
    for key in ("rho", "pr", "pt", "j"):
        np.testing.assert_array_equal(getattr(assembled, key), candidate[key])
        np.testing.assert_array_equal(
            getattr(assembled, f"electromagnetic_{key}"), em[key]
        )


@pytest.mark.parametrize(
    "bad",
    [
        {"rho": np.ones(2), "pr": np.ones(2), "pt": np.ones(2), "j": np.ones(2)},
        {"rho": np.array([np.nan]), "pr": np.array([0.0]),
         "pt": np.array([0.0]), "j": np.array([0.0])},
        {"rho": np.array([0.0]), "pr": np.array([0.0]),
         "pt": np.array([0.0])},
    ],
)
def test_optional_em_projection_fails_closed_on_bad_shape_nonfinite_or_missing_key(bad):
    _kernel, state = _fixture_state()
    with pytest.raises(ValueError):
        v55_matter.total_matter_projection(
            state.grid, state.geometry, state.scalars, state.matter,
            electromagnetic_projection=bad,
        )


def test_em_stress_reaches_actual_geometry_l3_and_cmc_sources(monkeypatch):
    _kernel, state = _fixture_state()
    em = _pure_radial_e_projection(state.grid)
    alpha = np.asarray(state.geometry.alpha, dtype=float)
    a = np.asarray(state.geometry.a, dtype=float)
    b = np.asarray(state.geometry.b, dtype=float)
    X = np.asarray(state.geometry.X, dtype=float)
    inv = X * X / a
    u = em["rho"]

    primary0 = primary_l3_with_matter(
        state.grid, state.geometry, state.scalars, state.matter
    )
    primary1 = primary_l3_with_matter(
        state.grid, state.geometry, state.scalars, state.matter,
        electromagnetic_projection=em,
    )
    np.testing.assert_allclose(
        primary1["Aa"] - primary0["Aa"],
        (32.0 * math.pi / 3.0) * alpha * u,
        rtol=2.0e-10, atol=2.0e-12,
    )
    np.testing.assert_allclose(
        primary1["K"] - primary0["K"],
        8.0 * math.pi * alpha * u,
        rtol=2.0e-10, atol=2.0e-12,
    )

    # Pure radial electric stress has j=0, hence no Lambda-source increment.
    lambda0 = lambda_l3_with_matter(
        state.grid, state.geometry, state.scalars, state.matter
    )
    lambda1 = lambda_l3_with_matter(
        state.grid, state.geometry, state.scalars, state.matter,
        electromagnetic_projection=em,
    )
    np.testing.assert_allclose(lambda1, lambda0, rtol=0.0, atol=2.0e-13)

    weights = (
        np.asarray(state.grid.centers, dtype=float) ** 2
        * np.sqrt(a) * b / (X ** 3)
    )
    weights *= float(state.grid.dr)
    expected_kdot = float(
        np.sum(weights * (8.0 * math.pi * alpha * u)) / np.sum(weights)
    )
    target0 = cmc_gauge.target_kdot(
        state.grid, state.geometry, state.scalars, state.matter
    )
    target1 = cmc_gauge.target_kdot(
        state.grid, state.geometry, state.scalars, state.matter,
        electromagnetic_projection=em,
    )
    np.testing.assert_allclose(
        target1 - target0, expected_kdot, rtol=2.0e-10, atol=2.0e-12
    )

    # Capture the actual CMC matrix and RHS to verify the local potential and
    # global target respond to the same optional canonical source exactly once.
    real_solve_banded = cmc_gauge.solve_banded
    captured = []

    def capture(*args, **kwargs):
        captured.append(
            (np.asarray(args[1], dtype=float).copy(),
             np.asarray(args[2], dtype=float).copy())
        )
        return real_solve_banded(*args, **kwargs)

    monkeypatch.setattr(cmc_gauge, "solve_banded", capture)
    alpha0, kdot0 = cmc_gauge.solve_cmc_lapse(
        state.grid, state.geometry, state.scalars, state.matter
    )
    band0, rhs0 = captured.pop()
    alpha1, kdot1 = cmc_gauge.solve_cmc_lapse(
        state.grid, state.geometry, state.scalars, state.matter,
        electromagnetic_projection=em,
    )
    band1, rhs1 = captured.pop()
    assert np.all(np.isfinite(alpha0)) and np.all(alpha0 > 0.0)
    assert np.all(np.isfinite(alpha1)) and np.all(alpha1 > 0.0)
    np.testing.assert_allclose(kdot1 - kdot0, expected_kdot, rtol=2.0e-10, atol=2.0e-12)
    np.testing.assert_allclose(
        band1[1, 1:-1] - band0[1, 1:-1],
        -8.0 * math.pi * u[1:-1] / inv[1:-1],
        rtol=2.0e-10, atol=2.0e-10,
    )
    np.testing.assert_allclose(
        rhs1[1:-1] - rhs0[1:-1],
        -(kdot1 - kdot0) / inv[1:-1],
        rtol=2.0e-10, atol=2.0e-10,
    )
    np.testing.assert_allclose(rhs1[-1] - rhs0[-1], 0.0, rtol=0.0, atol=0.0)


def test_outer_boundary_consumer_receives_canonical_em_sum(monkeypatch):
    kernel, state = _fixture_state()
    em = _pure_radial_e_projection(state.grid)
    baseline = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    adapter_module = importlib.import_module("engine.v55_pirk_adapter")
    adapter_module.vendor_modules()
    boundary_module = importlib.import_module("bssn_characteristic_boundary")
    captured = {}

    def capture_boundary(grid, geometry, rho, j):
        captured["rho"] = np.asarray(rho).copy()
        captured["j"] = np.asarray(j).copy()

    monkeypatch.setattr(
        boundary_module, "apply_light_constraint_boundary", capture_boundary
    )
    kernel._apply_outer_light_boundary(
        state.grid, state.geometry, state.scalars, state.matter,
        electromagnetic_projection=em,
    )
    np.testing.assert_array_equal(captured["rho"], baseline.rho + em["rho"])
    np.testing.assert_array_equal(captured["j"], baseline.j + em["j"])
