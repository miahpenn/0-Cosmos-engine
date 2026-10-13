"""Regression checks for archive-to-Einstein stress normalization."""

import numpy as np

from engine.matter_system import Species, project_species
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import KAPPA
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import metric_slice_from_q, total_fluid_projection


def test_archive_fluid_stress_is_converted_once_at_einstein_boundary():
    """DM, baryon, and radiation projections share the scalar's 1/KAPPA scale."""
    state = V55ProductionKernel().initialize(
        resolution=40, r_max=16.0, include_radiation=True
    )
    metric = metric_slice_from_q(state.grid, state.geometry)
    species_projections = (
        project_species(metric, state.matter.dark_matter, Species.DARK_MATTER),
        project_species(metric, state.matter.baryons, Species.BARYON),
        project_species(metric, state.matter.radiation, Species.RADIATION),
    )
    component_names = ("rho", "pr", "pt", "j")
    archive_fluid = {
        name: sum(piece[name] for piece in species_projections)
        for name in component_names
    }
    expected_physical = {
        name: archive_fluid[name] / KAPPA for name in component_names
    }

    fluid = total_fluid_projection(metric, state.matter)
    for name in component_names:
        np.testing.assert_allclose(
            fluid[name], expected_physical[name], rtol=2.0e-14, atol=1.0e-24
        )

    # Exercise the exact assembled object consumed by the geometry adapter and
    # constraint diagnostic, not only the helper's return value.
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    for key, field_name in (
        ("rho", "fluid_rho"),
        ("pr", "fluid_pr"),
        ("pt", "fluid_pt"),
        ("j", "fluid_j"),
    ):
        np.testing.assert_allclose(
            getattr(total, field_name),
            expected_physical[key],
            rtol=2.0e-14,
            atol=1.0e-24,
        )
        np.testing.assert_allclose(
            getattr(total, key),
            getattr(total, f"scalar_{key}") + expected_physical[key],
            rtol=2.0e-14,
            atol=1.0e-24,
        )
