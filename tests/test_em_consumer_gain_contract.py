"""Research preflight for radial-electric stress response in live consumers.

No EM field is added to ProductionState and no physical normalization is chosen.
The manufactured increment below is already in the engine's existing
rho/pr/pt/j consumer units. This pins only the linear response of the active
Einstein/CMC consumer code before a future action/unit bridge is approved.
"""
import math

import numpy as np

from engine import cmc_gauge, v55_matter
from engine.production_kernel import V55ProductionKernel
from engine.v55_pirk_adapter import (
    lambda_l3_with_matter,
    primary_l3_with_matter,
)


_STRESS_INCREMENT_SCALE = 2.5e-5


def _radial_e_increment(grid):
    """Symbolic positive energy profile in consumer units, not Gaussian u_EM."""
    r = np.asarray(grid.centers, dtype=float)
    return _STRESS_INCREMENT_SCALE * (1.0 + r / max(float(r[-1]), 1.0))


def _augment_radial_e_projection(base_projection):
    """Test wrapper: add (u, -u, +u, 0) without mutating baseline arrays."""
    def projected(grid, geometry, scalars, matter):
        baseline = base_projection(grid, geometry, scalars, matter)
        u = _radial_e_increment(grid)
        augmented = dict(baseline)
        augmented["rho"] = np.asarray(baseline["rho"]) + u
        augmented["pr"] = np.asarray(baseline["pr"]) - u
        augmented["pt"] = np.asarray(baseline["pt"]) + u
        # A static pure radial E field has no radial momentum/Poynting flux.
        augmented["j"] = np.asarray(baseline["j"]).copy()
        return augmented
    return projected


def _evaluate_consumers(state, monkeypatch, add_increment):
    """Run actual source consumers and capture the assembled CMC linear system."""
    base_matter_projection = v55_matter.total_matter_projection
    base_cmc_projection = cmc_gauge.total_matter_projection
    original_solve_banded = cmc_gauge.solve_banded
    captured_systems = []

    with monkeypatch.context() as patch:
        if add_increment:
            patch.setattr(
                v55_matter,
                "total_matter_projection",
                _augment_radial_e_projection(base_matter_projection),
            )
            patch.setattr(
                cmc_gauge,
                "total_matter_projection",
                _augment_radial_e_projection(base_cmc_projection),
            )

        def capture_linear_system(*args, **kwargs):
            # scipy.linalg.solve_banded((1,1), band, rhs, ...).
            captured_systems.append(
                (np.asarray(args[1], dtype=float).copy(),
                 np.asarray(args[2], dtype=float).copy())
            )
            return original_solve_banded(*args, **kwargs)

        patch.setattr(cmc_gauge, "solve_banded", capture_linear_system)

        primary = primary_l3_with_matter(
            state.grid, state.geometry, state.scalars, state.matter
        )
        connection = lambda_l3_with_matter(
            state.grid, state.geometry, state.scalars, state.matter
        )
        lapse, target_kdot = cmc_gauge.solve_cmc_lapse(
            state.grid, state.geometry, state.scalars, state.matter
        )
        assert np.all(np.isfinite(lapse)) and np.all(lapse > 0.0)
        assert len(captured_systems) == 1

    return primary, connection, target_kdot, captured_systems[0]


def _assert_live_consumer_gain(state, monkeypatch):
    baseline = _evaluate_consumers(state, monkeypatch, add_increment=False)
    candidate = _evaluate_consumers(state, monkeypatch, add_increment=True)

    primary0, connection0, kdot0, (band0, rhs0) = baseline
    primary1, connection1, kdot1, (band1, rhs1) = candidate

    alpha = np.asarray(state.geometry.alpha, dtype=float)
    a = np.asarray(state.geometry.a, dtype=float)
    b = np.asarray(state.geometry.b, dtype=float)
    X = np.asarray(state.geometry.X, dtype=float)
    assert np.all(np.isfinite(a)) and np.all(a > 0.0)
    assert np.all(np.isfinite(b)) and np.all(b > 0.0)
    assert np.all(np.isfinite(X)) and np.all(X > 0.0)
    inv = X * X / a
    u = _radial_e_increment(state.grid)

    # Derived from the actual active adapter coefficients for the pure radial
    # electric stress: rho=u, pr=-u, pt=+u, j=0. Units are consumer units.
    np.testing.assert_allclose(
        primary1["Aa"] - primary0["Aa"],
        (32.0 * math.pi / 3.0) * alpha * u,
        rtol=2.0e-10,
        atol=2.0e-12,
    )
    np.testing.assert_allclose(
        primary1["K"] - primary0["K"],
        8.0 * math.pi * alpha * u,
        rtol=2.0e-10,
        atol=2.0e-12,
    )
    np.testing.assert_allclose(
        connection1 - connection0, 0.0, rtol=0.0, atol=2.0e-13
    )

    # CMC local source Q receives 4*pi*(rho+pr+2*pt)=8*pi*u.
    # The banded diagonal is -2/dr**2 - Q/inv on interior rows.
    np.testing.assert_allclose(
        band1[1, 1:-1] - band0[1, 1:-1],
        -8.0 * math.pi * u[1:-1] / inv[1:-1],
        rtol=2.0e-10,
        atol=2.0e-10,
    )

    # target_kdot projects the K RHS with the current proper-volume weights;
    # the stress increment changes the source K tendency by 8*pi*alpha*u.
    r = np.asarray(state.grid.centers, dtype=float)
    weights = r**2 * np.sqrt(a) * b / X**3
    assert np.all(np.isfinite(weights)) and np.all(weights > 0.0)
    weights *= float(state.grid.dr)
    expected_delta_kdot = float(
        np.sum(weights * (8.0 * math.pi * alpha * u)) / np.sum(weights)
    )
    np.testing.assert_allclose(
        kdot1 - kdot0, expected_delta_kdot, rtol=2.0e-10, atol=2.0e-12
    )
    # Interior CMC RHS is -target_kdot/inv; the outer row is Dirichlet alpha=1.
    np.testing.assert_allclose(
        rhs1[1:-1] - rhs0[1:-1],
        -(kdot1 - kdot0) / inv[1:-1],
        rtol=2.0e-10,
        atol=2.0e-10,
    )
    np.testing.assert_allclose(rhs1[-1] - rhs0[-1], 0.0, atol=0.0, rtol=0.0)


def test_radial_e_consumer_gains_on_initial_and_accepted_production_slices(
    monkeypatch,
):
    """Pin geometry/CMC consumer response without adding production EM physics."""
    kernel = V55ProductionKernel()
    initial = kernel.initialize(
        resolution=40, r_max=40.0, include_radiation=True
    )
    dt = 0.03 * float(initial.grid.dr)
    accepted = kernel.step(initial, dt)
    assert initial.t == 0.0
    assert accepted.t > initial.t

    _assert_live_consumer_gain(initial, monkeypatch)
    _assert_live_consumer_gain(accepted, monkeypatch)
