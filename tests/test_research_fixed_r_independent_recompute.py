"""Regression tests for the independent fixed-R artifact recomputation."""
import copy

import numpy as np
import pytest

from tools.recompute_fixed_r_balance_artifact import recompute_arm


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
