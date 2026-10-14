"""Regression tests for the independent fixed-R artifact recomputation."""
import copy

import numpy as np
import pytest

from tools.recompute_fixed_r_balance_artifact import (
    recompute_arm,
    recompute_em_stage_exchange_witness,
)


def _raw_stage():
    indices = list(range(6))
    dr = 1.0
    r = np.arange(0.5, 6.5, 1.0)
    a = 1.1 + 0.015 * r
    b = np.ones_like(r)
    X = np.ones_like(r)
    alpha = 1.0 + 0.03 * r
    a_dot = 0.001 + 0.0002 * r
    b_dot = 0.02 + 0.001 * r
    X_dot = -0.0005 + 0.0001 * r
    R = r * np.sqrt(b) / X
    Rdot = R * (0.5 * b_dot / b - X_dot / X)

    def d4(values, i):
        return (
            values[i - 2] - 8.0 * values[i - 1]
            + 8.0 * values[i + 1] - values[i + 2]
        ) / (12.0 * dr)

    Rr = np.ones_like(R)
    Rdot_r = np.array([d4(Rdot, i) if 2 <= i <= 3 else 0.0 for i in indices])

    K = 0.012 + 0.0004 * r
    Aa = 0.0003 - 0.00002 * r
    l2_K = 0.002 + 0.0001 * r
    l3_K = -0.0007 + 0.00004 * r
    l2_Aa = 0.00005 - 0.00001 * r
    l3_Aa = -0.00002 + 0.000005 * r
    Kdot = l2_K + l3_K
    Aadot = l2_Aa + l3_Aa
    Ktheta = K / 3.0 - Aa / 2.0
    Ktheta_dot = Kdot / 3.0 - Aadot / 2.0
    gamma_inv = X**2 / a
    gamma_inv_dot = gamma_inv * (2.0 * X_dot / X - a_dot / a)
    chi = gamma_inv * Rr**2 - (R * Ktheta)**2
    mass = 0.5 * R * (1.0 - chi)

    terms = (
        0.5 * Rdot * (1.0 - chi),
        -0.5 * R * gamma_inv_dot * Rr**2,
        -R * gamma_inv * Rr * Rdot_r,
        R**2 * Rdot * Ktheta**2,
        R**3 * Ktheta * Ktheta_dot,
    )
    mass_dot = terms[0] + terms[1] + terms[2] + terms[3] + terms[4]

    sectors = {
        "scalar": {
            "rho": 0.003 + 0.0001 * r,
            "pr": 0.0002 + 0.00003 * r,
            "pt": 0.00015 + 0.00002 * r,
            "j": 0.00001 * r,
        },
        "fluid": {
            "rho": 0.0004 + 0.00001 * r,
            "pr": 0.00006 + 0.000004 * r,
            "pt": 0.00005 + 0.000003 * r,
            "j": -0.000002 * r,
        },
        "electromagnetic": {
            "rho": 0.00003 + 0.000001 * r,
            "pr": -(0.00003 + 0.000001 * r),
            "pt": 0.00003 + 0.000001 * r,
            "j": np.zeros_like(r),
        },
    }
    aggregate = {
        key: sum((sector[key] for sector in sectors.values()), np.zeros_like(r))
        for key in ("rho", "pr", "pt", "j")
    }
    raw = {
        "schema_version": 1,
        "target_areal_radius": 3.0,
        "grid_dr": dr,
        "global_cell_indices": indices,
        "active_bracket": {"left_index": 2, "right_index": 3, "theta": 0.5},
        "profiles": {
            "r": r, "a": a, "b": b, "X": X, "alpha": alpha,
            "a_dot_explicit": a_dot, "b_dot_explicit": b_dot, "X_dot_explicit": X_dot,
            "R": R, "Rr_parity_aware": Rr, "Rdot": Rdot,
            "Rdot_r_parity_aware": Rdot_r, "gamma_inv": gamma_inv,
            "gamma_inv_dot": gamma_inv_dot, "K": K, "Aa": Aa,
            "Ktheta": Ktheta, "K_dot_stage_sum": Kdot, "Aa_dot_stage_sum": Aadot,
            "Ktheta_dot": Ktheta_dot, "chi": chi,
            "mass_native": mass, "mass_chain_rule": mass,
            "mass_dot_coordinate_chain_rule": mass_dot,
            **aggregate,
        },
        "rhs_split": {
            "explicit": {"a": a_dot, "b": b_dot, "X": X_dot},
            "primary_l2": {"K": l2_K, "Aa": l2_Aa},
            "primary_l3": {"K": l3_K, "Aa": l3_Aa},
            "vacuum_primary_l3": {"K": np.zeros_like(r), "Aa": np.zeros_like(r)},
        },
        "stress_sectors": sectors,
    }
    return {"t": 0.25, "tube_fixed_R_balance_raw_witness": raw}


def test_independent_fixed_r_reconstruction_keeps_operator_conventions_separate():
    result = recompute_arm(_raw_stage())
    values = result["interpolated"]
    assert result["active_bracket"]["left_index"] == 2
    assert result["active_bracket"]["right_index"] == 3
    assert result["active_bracket"]["theta"] == pytest.approx(0.5)
    assert values["moving_interpolant_residual_unthresholded"] == pytest.approx(
        values["composite_fixed_R_residual_unthresholded"], abs=2e-17
    )
    # Factor- and composite-interpolated products are separate discretizations.
    assert values["composite_fixed_R_residual_unthresholded"] != pytest.approx(
        values["factor_fixed_R_residual_unthresholded"], abs=1e-14
    )
    assert set(result["cells"]) == {"2", "3"}
    assert set(result["stencil_roundoff_checks"]) == {"2", "3"}
    assert "mass_gradient_identity_residual_unthresholded" in result["cells"]["2"]
    assert "Ktheta_dot_matter_source_from_live_minus_vacuum_L3" in result["cells"]["3"]
    assert "stress_j_sector_additivity" in result["algebraic_roundoff_checks"]
    assert "without a pass/fail threshold" in result["interpretation"]


def test_independent_fixed_r_reconstruction_fails_on_bad_saved_derivative():
    stage = _raw_stage()
    broken = copy.deepcopy(stage)
    broken["tube_fixed_R_balance_raw_witness"]["profiles"][
        "Rdot_r_parity_aware"
    ][2] += 1.0e-4
    with pytest.raises(RuntimeError, match=r"D_r\(Rdot\)"):
        recompute_arm(broken)


def test_independent_fixed_r_reconstruction_fails_on_rhs_stage_mismatch():
    stage = _raw_stage()
    broken = copy.deepcopy(stage)
    broken["tube_fixed_R_balance_raw_witness"]["rhs_split"]["primary_l3"]["K"][2] += 1.0e-4
    with pytest.raises(RuntimeError, match=r"K dot from primary L2\+L3"):
        recompute_arm(broken)



def _raw_em_stage_record():
    centers = np.array([0.5, 1.5, 2.5], dtype=np.float64)
    dr = 1.0
    eta = 0.2
    charge_faces = np.array([0.0, 0.08, 0.19, 0.27], dtype=np.float64)
    rest_flux = np.array([0.0, 0.01, -0.02, 0.03], dtype=np.float64)
    current = 4.0 * np.pi * eta * rest_flux
    qdot = -current
    qcell = 0.5 * (charge_faces[:-1] + charge_faces[1:])
    qcelldot = 0.5 * (qdot[:-1] + qdot[1:])
    X = np.array([1.1, 1.2, 1.3])
    a = np.array([1.2, 1.1, 1.4])
    b = np.array([1.0, 1.1, 1.3])
    Xdot = np.array([0.01, -0.01, 0.02])
    adot = np.array([0.005, 0.002, -0.003])
    bdot = np.array([0.01, -0.004, 0.005])
    weight = X * np.sqrt(a) / (centers**2 * b)
    log_weight_dot = Xdot / X + 0.5 * adot / a - bdot / b
    weight_dot = weight * log_weight_dot
    charge_terms = dr * weight * qcell * qcelldot
    metric_terms = 0.5 * dr * weight * qcell**2 * log_weight_dot
    charge_power = float(np.sum(charge_terms))
    metric_power = float(np.sum(metric_terms))
    lorentz = np.array([1.0e-5, -2.0e-5, 1.5e-5])
    lorentz_power = 4.0 * np.pi * dr * float(np.sum(lorentz))
    outer_power = -0.5 * dr * weight[-1] * qcell[-1] * current[-1]
    witness = {
        "grid_dr": dr,
        "charge_to_rest_mass": eta,
        "centers": centers.tolist(),
        "geometry": {"X": X.tolist(), "a": a.tolist(), "b": b.tolist()},
        "explicit_geometry_rhs": {
            "X": Xdot.tolist(), "a": adot.tolist(), "b": bdot.tolist()
        },
        "charge_faces": charge_faces.tolist(),
        "rest_mass_flux_faces": rest_flux.tolist(),
        "charge_current_faces": current.tolist(),
        "charge_face_dot": qdot.tolist(),
        "charge_cell_average": qcell.tolist(),
        "charge_cell_dot_average": qcelldot.tolist(),
        "field_energy_weight": weight.tolist(),
        "log_energy_weight_rate": log_weight_dot.tolist(),
        "energy_weight_dot": weight_dot.tolist(),
        "charge_current_tangent_cell_terms": charge_terms.tolist(),
        "metric_tangent_cell_terms": metric_terms.tolist(),
        "field_energy_at_stage": 0.5 * dr * float(np.sum(weight * qcell**2)),
        "charge_current_tangent": charge_power,
        "metric_tangent": metric_power,
        "total_field_energy_tangent": charge_power + metric_power,
        "lorentz_energy_source_profile": lorentz.tolist(),
        "matter_lorentz_power": lorentz_power,
        "outer_current_power": outer_power,
    }
    return {
        "stage": "start",
        "t": 0.01,
        "charge_power": charge_power,
        "metric_power": metric_power,
        "matter_lorentz_power": lorentz_power,
        "outer_current_power": outer_power,
        "stage_witness": witness,
    }


def test_independent_em_stage_recomputation_rebuilds_face_charge_and_energy_tangent():
    result = recompute_em_stage_exchange_witness(_raw_em_stage_record())
    assert result["field_energy_at_stage"] > 0.0
    assert result["physical_residual_threshold_applied"] is False
    assert result["total_field_energy_tangent"] == pytest.approx(
        result["charge_current_tangent"] + result["metric_tangent"], abs=1e-18
    )
    assert result["exchange_only_remainder_unthresholded"] != 0.0
    assert all(check["bound_ratio"] <= 1.0 for check in result["algebraic_roundoff_checks"].values())


def test_independent_em_stage_recomputation_fails_closed_on_wrong_charge_tangent():
    broken = _raw_em_stage_record()
    broken["stage_witness"]["charge_face_dot"][1] += 1e-4
    with pytest.raises(RuntimeError, match="charge-face tangent reconstruction"):
        recompute_em_stage_exchange_witness(broken)
