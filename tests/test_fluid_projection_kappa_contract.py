"""Independent unit-conversion contract for Eulerian fluid stress projections.

Archive-normalized primitive densities must be converted to Einstein-source
units exactly once at the species projection boundary. The total-fluid
aggregator must sum those already-converted projections without another
KAPPA factor.
"""
import numpy as np
import pytest

from engine.matter_system import (
    BSSNMetricSlice,
    KAPPA,
    Species,
    _metric_arrays,
    initialize_dust,
    initialize_radiation,
    primitives,
    project_species,
)
from engine.v55_matter import V55MatterState, total_fluid_projection


def _metric_slice() -> BSSNMetricSlice:
    return BSSNMetricSlice(
        r=np.asarray([0.25, 0.75, 1.50], dtype=float),
        a=np.asarray([1.03, 1.07, 1.11], dtype=float),
        b=np.asarray([0.98, 1.02, 1.04], dtype=float),
        X=np.asarray([0.96, 0.98, 1.01], dtype=float),
        alpha=np.asarray([0.92, 0.94, 0.96], dtype=float),
        beta=np.asarray([0.006, 0.004, 0.002], dtype=float),
        Aa=np.asarray([0.0010, 0.0007, 0.0004], dtype=float),
        K=np.asarray([0.020, 0.018, 0.015], dtype=float),
        Lambda=np.asarray([0.003, 0.002, 0.001], dtype=float),
        B=np.asarray([0.0, 0.0, 0.0], dtype=float),
    )


def _direct_primitive_projection(metrics, primitive_values):
    """Compute the physical Eulerian projection independently of project_species."""
    raw = {key: np.empty(len(metrics), dtype=float) for key in ("rho", "pr", "pt", "j")}
    for i, (metric, primitive) in enumerate(zip(metrics, primitive_values)):
        w = primitive.lorentz()
        h = primitive.rho + primitive.pressure
        vhat2 = metric.gamma_rr * primitive.v_r * primitive.v_r
        raw["rho"][i] = h * w * w - primitive.pressure
        raw["pr"][i] = h * w * w * vhat2 + primitive.pressure
        raw["pt"][i] = primitive.pressure
        raw["j"][i] = h * w * w * metric.gamma_rr * primitive.v_r
    return {key: values / KAPPA for key, values in raw.items()}, raw


@pytest.fixture
def nontrivial_fluid_bundle():
    metric = _metric_slice()
    metrics = _metric_arrays(metric)
    rho_dm = np.asarray([2.4e-5, 2.1e-5, 1.8e-5], dtype=float)
    rho_baryon = np.asarray([5.5e-6, 4.6e-6, 3.7e-6], dtype=float)
    rho_radiation = np.asarray([8.0e-9, 7.0e-9, 6.0e-9], dtype=float)
    v_dm = np.asarray([0.025, -0.018, 0.011], dtype=float)
    v_baryon = np.asarray([-0.012, 0.008, 0.016], dtype=float)
    v_radiation = np.asarray([0.035, -0.022, 0.013], dtype=float)

    matter = V55MatterState(
        dark_matter=initialize_dust(metrics, rho_dm, v_dm),
        baryons=initialize_dust(metrics, rho_baryon, v_baryon),
        radiation=initialize_radiation(metrics, rho_radiation, v_radiation),
    )
    species_states = {
        Species.DARK_MATTER: matter.dark_matter,
        Species.BARYON: matter.baryons,
        Species.RADIATION: matter.radiation,
    }
    primitives_by_species = {
        species: primitives(metrics, state, species)
        for species, state in species_states.items()
    }
    return metric, metrics, matter, species_states, primitives_by_species


@pytest.mark.parametrize(
    "species",
    [Species.DARK_MATTER, Species.BARYON, Species.RADIATION],
)
def test_species_projection_converts_each_component_exactly_once(
    nontrivial_fluid_bundle, species
):
    metric, metrics, _matter, species_states, primitives_by_species = nontrivial_fluid_bundle
    projected = project_species(metric, species_states[species], species)
    expected, raw = _direct_primitive_projection(metrics, primitives_by_species[species])

    assert set(projected) == {"rho", "pr", "pt", "j"}
    for component in ("rho", "pr", "pt", "j"):
        assert np.all(np.isfinite(projected[component]))
        np.testing.assert_allclose(
            projected[component], expected[component], rtol=3.0e-14, atol=1.0e-25
        )
        # Explicitly check the unit contract independently of the returned
        # projected arrays: KAPPA * physical projection is the raw fluid term.
        np.testing.assert_allclose(
            KAPPA * projected[component], raw[component], rtol=3.0e-14, atol=1.0e-24
        )

    # Ensure nonzero velocity/pressure data makes the current and tangential
    # stress checks meaningful rather than vacuous zero comparisons.
    assert np.any(np.abs(projected["j"]) > 0.0)
    if species is Species.RADIATION:
        assert np.all(projected["pt"] > 0.0)


def test_total_fluid_projection_sums_once_scaled_species_without_a_second_factor(
    nontrivial_fluid_bundle,
):
    metric, metrics, matter, species_states, primitives_by_species = nontrivial_fluid_bundle
    total = total_fluid_projection(metric, matter)
    expected_by_species = {
        species: _direct_primitive_projection(metrics, prims)[0]
        for species, prims in primitives_by_species.items()
    }

    assert set(total) == {"rho", "pr", "pt", "j"}
    for component in ("rho", "pr", "pt", "j"):
        expected_sum = sum(
            expected_by_species[species][component]
            for species in (Species.DARK_MATTER, Species.BARYON, Species.RADIATION)
        )
        np.testing.assert_allclose(total[component], expected_sum, rtol=3.0e-14, atol=1.0e-24)

    # At least one component must materially differ from its raw unconverted
    # total; this catches accidental omission of KAPPA as well as a second
    # aggregate division.
    raw_total_rho = sum(
        _direct_primitive_projection(metrics, primitives_by_species[s])[1]["rho"]
        for s in (Species.DARK_MATTER, Species.BARYON, Species.RADIATION)
    )
    assert np.max(np.abs(total["rho"] - raw_total_rho)) > 0.9 * np.max(np.abs(raw_total_rho))
    np.testing.assert_allclose(total["rho"] * KAPPA, raw_total_rho, rtol=3.0e-14, atol=1.0e-23)
