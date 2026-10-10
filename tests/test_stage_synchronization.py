"""Deterministic stage-synchronization checks for the V5.5 production step.

No trajectory is interpreted here. The gate checks that predictor RHS calls
consume synchronized geometry/scalar/matter states and that scalar/matter
correctors use the matching two RHS evaluations.
"""
import numpy as np

from engine.matter_system import ConservedSpecies
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState


_GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
_SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
_SPECIES = ("dark_matter", "baryons", "radiation")
_CONSERVED_FIELDS = ("rest", "energy_t", "momentum_r")


def _snapshot(state):
    return {
        "geometry": {
            name: np.array(getattr(state.geometry, name), copy=True)
            for name in _GEOMETRY_FIELDS
        },
        "scalars": {
            name: np.array(getattr(state.scalars, name), copy=True)
            for name in _SCALAR_FIELDS
        },
        "matter": {
            species: {
                name: np.array(
                    getattr(getattr(state.matter, species), name), copy=True
                )
                for name in _CONSERVED_FIELDS
            }
            for species in _SPECIES
        },
    }


def _assert_array_equal(actual, expected, message):
    try:
        np.testing.assert_array_equal(actual, expected)
    except AssertionError as exc:
        raise AssertionError(message) from exc


def _assert_stage_matches_solve_input(stage, solve_call, label):
    # The solver is a one-pass projected CMC solve, not necessarily an
    # idempotent fixed-point map. Compare to the exact input and output from
    # that stage's own solve rather than solving again from the returned lapse.
    # B is a gauge auxiliary that the step intentionally resets to zero after
    # solving the lapse, and is not an input to solve_cmc_lapse.
    for name in _GEOMETRY_FIELDS:
        if name not in ("alpha", "B"):
            _assert_array_equal(
                stage["geometry"][name],
                solve_call["input"]["geometry"][name],
                f"{label} RHS geometry {name} differs from its lapse-solve input",
            )
    for name in _SCALAR_FIELDS:
        _assert_array_equal(
            stage["scalars"][name],
            solve_call["input"]["scalars"][name],
            f"{label} RHS scalar {name} differs from its lapse-solve input",
        )
    for species in _SPECIES:
        for name in _CONSERVED_FIELDS:
            _assert_array_equal(
                stage["matter"][species][name],
                solve_call["input"]["matter"][species][name],
                f"{label} RHS {species}.{name} differs from its lapse-solve input",
            )
    _assert_array_equal(
        stage["geometry"]["alpha"], solve_call["alpha"],
        f"{label} RHS did not use its own CMC solve output",
    )


def test_step_uses_synchronized_predictor_and_matching_heun_rhs(monkeypatch):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32, r_max=16.0, D_amplitude=1.0e-10,
        include_radiation=True,
    )
    dt = 2.0e-3
    initial_scalars = {
        name: np.array(getattr(state.scalars, name), copy=True)
        for name in _SCALAR_FIELDS
    }
    initial_matter = {
        species: {
            name: np.array(
                getattr(getattr(state.matter, species), name), copy=True
            )
            for name in _CONSERVED_FIELDS
        }
        for species in _SPECIES
    }

    original_rhs = kernel._rhs
    original_solve = kernel._solve_lapse
    rhs_calls = []
    solve_calls = []

    def traced_solve(grid, geometry, scalars, matter):
        solve_input = {
            "geometry": {
                name: np.array(getattr(geometry, name), copy=True)
                for name in _GEOMETRY_FIELDS
            },
            "scalars": {
                name: np.array(getattr(scalars, name), copy=True)
                for name in _SCALAR_FIELDS
            },
            "matter": {
                species: {
                    name: np.array(
                        getattr(getattr(matter, species), name), copy=True
                    )
                    for name in _CONSERVED_FIELDS
                }
                for species in _SPECIES
            },
        }
        alpha, kdot = original_solve(grid, geometry, scalars, matter)
        solve_calls.append({
            "input": solve_input,
            "alpha": np.array(alpha, copy=True),
            "kdot": float(kdot),
        })
        return alpha, kdot

    def traced_rhs(stage_state):
        stage_snapshot = _snapshot(stage_state)
        rhs_scalar, rhs_matter, derivatives = original_rhs(stage_state)
        rhs_calls.append((stage_snapshot, rhs_scalar, rhs_matter, derivatives))
        return rhs_scalar, rhs_matter, derivatives

    monkeypatch.setattr(kernel, "_solve_lapse", traced_solve)
    monkeypatch.setattr(kernel, "_rhs", traced_rhs)
    final_state = kernel.step(state, dt)

    assert len(rhs_calls) == 2, (
        "one synchronized initial RHS and one synchronized predictor RHS expected"
    )
    assert len(solve_calls) == 4, (
        "expected initial, predictor-preparation, predictor-RHS, and final lapse solves"
    )
    stage0, rhs0_scalar, rhs0_matter, _ = rhs_calls[0]
    stage1, rhs1_scalar, rhs1_matter, _ = rhs_calls[1]

    # The predictor passed to RHS1 must contain the explicit Euler update
    # generated by RHS0 for every scalar and all conservative fluid variables.
    for name in _SCALAR_FIELDS:
        _assert_array_equal(
            stage1["scalars"][name],
            initial_scalars[name] + dt * getattr(rhs0_scalar, name),
            f"predictor scalar {name} is not synchronized to RHS0",
        )
    for species in _SPECIES:
        for name in _CONSERVED_FIELDS:
            _assert_array_equal(
                stage1["matter"][species][name],
                initial_matter[species][name]
                + dt * getattr(rhs0_matter[species], name),
                f"predictor {species}.{name} is not synchronized to RHS0",
            )

    # RHS0 must consume the current slice and the exact lapse returned by
    # solve call 0. RHS1 consumes the synchronized predictor after call 2;
    # call 1 is intentionally only the pre-primary geometry-stage solve.
    _assert_stage_matches_solve_input(stage0, solve_calls[0], "initial")
    _assert_stage_matches_solve_input(stage1, solve_calls[2], "predictor")
    for label, stage in (("initial", stage0), ("predictor", stage1)):
        _assert_array_equal(
            stage["geometry"]["beta"],
            np.zeros_like(stage["geometry"]["beta"]),
            f"{label} RHS stage did not use the zero-shift CMC branch",
        )

    # The accepted matter/scalar state must use the trapezoidal pairing of
    # exactly those two stage RHS evaluations.
    for name in _SCALAR_FIELDS:
        expected = initial_scalars[name] + 0.5 * dt * (
            getattr(rhs0_scalar, name) + getattr(rhs1_scalar, name)
        )
        _assert_array_equal(
            getattr(final_state.scalars, name), expected,
            f"accepted scalar {name} did not use matching RHS0/RHS1",
        )
    for species in _SPECIES:
        for name in _CONSERVED_FIELDS:
            expected = initial_matter[species][name] + 0.5 * dt * (
                getattr(rhs0_matter[species], name)
                + getattr(rhs1_matter[species], name)
            )
            _assert_array_equal(
                getattr(getattr(final_state.matter, species), name), expected,
                f"accepted {species}.{name} did not use matching RHS0/RHS1",
            )

    # The final lapse must be the fourth solve's output on the accepted
    # geometry/scalar/matter slice. B is reset after that solve by design.
    final_snapshot = _snapshot(final_state)
    _assert_stage_matches_solve_input(final_snapshot, solve_calls[3], "final")
    _assert_array_equal(
        final_snapshot["geometry"]["beta"],
        np.zeros_like(final_snapshot["geometry"]["beta"]),
        "accepted state did not preserve the zero-shift CMC branch",
    )
    _assert_array_equal(
        final_snapshot["geometry"]["B"],
        np.zeros_like(final_snapshot["geometry"]["B"]),
        "accepted state did not reset the auxiliary B gauge variable",
    )

    # The caller's start state is not mutated in place by predictor formation.
    for name in _SCALAR_FIELDS:
        _assert_array_equal(
            getattr(state.scalars, name), initial_scalars[name],
            f"input scalar {name} was mutated during the step",
        )
    for species in _SPECIES:
        for name in _CONSERVED_FIELDS:
            _assert_array_equal(
                getattr(getattr(state.matter, species), name),
                initial_matter[species][name],
                f"input {species}.{name} was mutated during the step",
            )
