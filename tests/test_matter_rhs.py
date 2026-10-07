import numpy as np

import engine.matter_rhs as matter_rhs
from engine.matter_system import (
    ConservedSpecies,
    Species,
)


def test_species_rhs_does_not_validate_unit_bookkeeping_step(monkeypatch):
    """RHS extraction must not reject a hypothetical dt=1 physical update."""
    state = ConservedSpecies(
        np.ones(2),
        np.full(2, 2.0),
        np.zeros(2),
    )
    captured = {}

    def fake_evolve_species(
        metric,
        metric_derivatives,
        state,
        species,
        dt,
        dphi_t=None,
        dphi_r=None,
        beta_dm=-0.04,
        validate_physical_state=True,
    ):
        captured.update(
            dt=dt,
            species=species,
            validate_physical_state=validate_physical_state,
        )
        return ConservedSpecies(
            state.rest + 1.0,
            state.energy_t + 2.0,
            state.momentum_r + 3.0,
        )

    monkeypatch.setattr(matter_rhs, "evolve_species", fake_evolve_species)

    rhs = matter_rhs.species_rhs(
        None,
        None,
        state,
        Species.RADIATION,
    )

    assert captured == {
        "dt": 1.0,
        "species": Species.RADIATION,
        "validate_physical_state": False,
    }
    np.testing.assert_allclose(rhs.rest, 1.0)
    np.testing.assert_allclose(rhs.energy_t, 2.0)
    np.testing.assert_allclose(rhs.momentum_r, 3.0)
