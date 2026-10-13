"""Pure contracts for the source-pinned joined EM energy ledger."""
import numpy as np
import pytest

from tools.audit_em_total_energy_budget import (
    SOURCE_PINS,
    VENDOR_GITLINK,
    VENDOR_PIN,
    field_energy_functional,
    interpolate_at_areal_radius,
    _reconstruct_energy_rhs_with_roundoff_contract,
)


def test_observer_is_pinned_to_all_ledger_source_boundaries():
    required = {
        "engine/campaign.py",
        "engine/cmc_gauge.py",
        "engine/cosmology_observables.py",
        "engine/invariant_diagnostics.py",
        "engine/matter_rhs.py",
        "engine/matter_system.py",
        "engine/production_contract.py",
        "engine/production_kernel.py",
        "engine/radial_maxwell.py",
        "engine/run_production.py",
        "engine/scalar_system.py",
        "engine/stress_energy.py",
        "engine/v55_initial.py",
        "engine/v55_matter.py",
        "engine/valencia.py",
        "engine/v55_pirk_adapter.py",
        "engine/worldtube.py",
    }
    assert set(SOURCE_PINS) == required
    assert all(len(value) == 40 for value in SOURCE_PINS.values())
    assert VENDOR_GITLINK == "vendor/bb-palatini-unified-r0"
    assert len(VENDOR_PIN) == 40


def test_field_energy_uses_arithmetic_face_to_cell_charge_and_curved_weight():
    dr = 0.5
    centers = np.array([0.25, 0.75])
    geometry = {
        "X": np.array([1.0, 2.0]),
        "a": np.array([4.0, 1.0]),
        "b": np.array([1.0, 4.0]),
    }
    charge_faces = np.array([0.0, 0.5, 1.0])
    q_cell = np.array([0.25, 0.75])
    weight = geometry["X"] * np.sqrt(geometry["a"]) / (
        centers**2 * geometry["b"]
    )
    expected = 0.5 * dr * np.sum(weight * q_cell**2)
    observed = field_energy_functional(dr, centers, geometry, charge_faces)
    assert observed == pytest.approx(expected, rel=0.0, abs=0.0)
    assert observed > 0.0


def test_fixed_areal_radius_interpolation_is_bounded_and_monotone():
    areal_radii = np.array([1.0, 3.0, 5.0])
    profile = np.array([2.0, 6.0, 10.0])
    assert interpolate_at_areal_radius(areal_radii, profile, 4.0) == pytest.approx(8.0)
    with pytest.raises(ValueError, match="strictly increasing"):
        interpolate_at_areal_radius([1.0, 4.0, 3.0], profile, 3.5)
    with pytest.raises(ValueError, match="outside"):
        interpolate_at_areal_radius(areal_radii, profile, 7.0)


def test_cellwise_rhs_reconstruction_passes_roundoff_scale_and_reports_ratio():
    terms = (
        ("flux", np.array([1.0e8, 1.0e-10])),
        ("source", np.array([-1.0e8, -2.0e-10])),
    )
    observed = np.array([np.nextafter(0.0, 1.0), -1.0e-10])
    reconstructed, residual, ratio_max = _reconstruct_energy_rhs_with_roundoff_contract(
        observed, terms, label="test species start energy RHS"
    )
    assert reconstructed.shape == observed.shape
    assert np.all(np.isfinite(residual))
    assert ratio_max <= 1.0


def test_cellwise_rhs_reconstruction_fails_closed_above_roundoff_bound():
    with pytest.raises(RuntimeError, match="exceeded scale-aware roundoff bound"):
        _reconstruct_energy_rhs_with_roundoff_contract(
            np.array([1.0e-3]),
            (("flux", np.array([0.0])), ("source", np.array([0.0]))),
            label="charged_dust primary_predictor energy RHS",
        )


def test_cellwise_rhs_reconstruction_fails_closed_on_shape_mismatch():
    with pytest.raises(RuntimeError, match="does not match RHS shape"):
        _reconstruct_energy_rhs_with_roundoff_contract(
            np.array([0.0, 0.0]),
            (("flux", np.array([0.0])),),
            label="test RHS",
        )
