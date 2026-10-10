"""Baseline tripwires for the legacy reference initial-data diagnostic.

These tests preserve a known numerical mismatch as an OPEN diagnostic. A passing
test means the witness was reproduced, not that the physical model is validated.
"""
import pytest

from engine.reference_initial_data_consistency import audit_initial_slice


def test_unperturbed_reference_background_closes_to_roundoff():
    for n, r_max in ((40, 40.0), (80, 40.0), (160, 160.0)):
        case = audit_initial_slice(n, r_max, amplitude=0.0)
        assert case["finite"]
        assert case["metric_B_identity_max_abs_b3_minus_X6"] <= 1.0e-14
        assert case["source_balance_B_plus_r_dB_minus_source"]["max_abs_all"] <= 1.0e-11
        assert case["ricci_BSSN_minus_polar"]["max_abs_all"] <= 1.0e-11
        assert case["hamiltonian_source_algebraic_closure"]["max_abs_all"] <= 1.0e-11
        assert case["hamiltonian_BSSN"]["max_abs_all"] <= 1.0e-11
        assert case["hamiltonian_independent_polar"]["max_abs_all"] <= 1.0e-11
        assert case["hamiltonian_decomposition_closure"]["max_abs_all"] <= 1.0e-11


def test_perturbed_reference_slice_preserves_open_center_ricci_mismatch():
    coarse = audit_initial_slice(40, 40.0, amplitude=0.01)
    fine = audit_initial_slice(320, 40.0, amplitude=0.01)
    coarse_center = coarse["ricci_BSSN_minus_polar"]["first_cell"]
    fine_center = fine["ricci_BSSN_minus_polar"]["first_cell"]
    assert abs(coarse_center) > 3.0e-4
    assert abs(fine_center) > 3.0e-4
    assert abs(coarse_center - fine_center) < 6.0e-5
    assert abs(fine["hamiltonian_BSSN"]["first_cell"]) > 5.0e-4


def test_fixed_domain_refinement_improves_off_center_witness_but_not_center_one():
    n40 = audit_initial_slice(40, 40.0, amplitude=0.01)
    n80 = audit_initial_slice(80, 40.0, amplitude=0.01)
    n320 = audit_initial_slice(320, 40.0, amplitude=0.01)
    r40 = abs(n40["hamiltonian_BSSN"]["at_r_near_2p5"])
    r80 = abs(n80["hamiltonian_BSSN"]["at_r_near_2p5"])
    assert r40 > 10.0 * r80
    assert abs(n320["ricci_BSSN_minus_polar"]["first_cell"]) > 3.0e-4


def test_fixed_spacing_domain_growth_leaves_local_r25_residual_unchanged():
    cases = [
        audit_initial_slice(n, float(n), amplitude=0.01)
        for n in (40, 80, 160, 320)
    ]
    values = [case["hamiltonian_BSSN"]["at_r_near_2p5"] for case in cases]
    for value in values[1:]:
        assert value == pytest.approx(values[0], rel=0.0, abs=1.0e-12)
    assert values[0] == pytest.approx(-3.698065997959077e-5, abs=2.0e-10)


def test_hamiltonian_residual_is_resolved_into_source_and_ricci_terms():
    case = audit_initial_slice(40, 40.0, amplitude=0.01)
    assert case["hamiltonian_source_algebraic_closure"]["max_abs_all"] < 1.0e-11
    assert case["hamiltonian_decomposition_closure"]["max_abs_all"] < 1.0e-10
    rb_p = case["ricci_BSSN_minus_polar"]["first_cell"]
    rp_s = case["ricci_polar_minus_source"]["first_cell"]
    rb_s = case["ricci_BSSN_minus_source"]["first_cell"]
    assert rb_s == pytest.approx(rb_p + rp_s, abs=1.0e-12)
    rho_mismatch = case["matter_density_constraint_minus_source"]["first_cell"]
    h_bssn = case["hamiltonian_BSSN"]["first_cell"]
    assert h_bssn == pytest.approx(rb_s - 16.0 * 3.141592653589793 * rho_mismatch, abs=1.0e-10)
    assert abs(rho_mismatch) < 1.0e-8
