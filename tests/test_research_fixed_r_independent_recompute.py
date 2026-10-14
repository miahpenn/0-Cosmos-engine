"""Regression tests for the independent fixed-R artifact recomputation."""
import copy

import numpy as np
import pytest

from tools.recompute_fixed_r_balance_artifact import (
    _geometric_radial_einstein_rr_at_cell,
    _summarize_em_stage_exchange,
    recompute_arm,
    recompute_charged_dust_boundary_witness,
    recompute_combined_matter_em_accounting_candidates,
    recompute_em_stage_exchange_campaign,
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
    momentum_constraint_residual = np.zeros_like(r)
    for i in range(2, len(r) - 2):
        aa_r = d4(Aa, i)
        k_r = d4(K, i)
        x_r = d4(X, i)
        b_r = d4(b, i)
        chi_p = -0.5 * x_r / X[i]
        ab = -0.5 * Aa[i]
        vacuum_momentum = (
            aa_r - (2.0 / 3.0) * k_r + 6.0 * Aa[i] * chi_p
            + (Aa[i] - ab) * (2.0 / r[i] + b_r / b[i])
        )
        momentum_constraint_residual[i] = (
            vacuum_momentum - 8.0 * np.pi * aggregate["j"][i]
        )
    hamiltonian_constraint_residual = 1.0e-5 + 2.0e-6 * r
    raw = {
        "schema_version": 1,
        "target_areal_radius": 3.0,
        "grid_dr": dr,
        "global_cell_indices": indices,
        "active_bracket": {"left_index": 2, "right_index": 3, "theta": 0.5},
        "profiles": {
            "r": r, "a": a, "b": b, "X": X, "alpha": alpha,
            "beta": np.zeros_like(r),
            "a_dot_explicit": a_dot, "b_dot_explicit": b_dot, "X_dot_explicit": X_dot,
            "R": R, "Rr_parity_aware": Rr, "Rdot": Rdot,
            "Rdot_r_parity_aware": Rdot_r, "gamma_inv": gamma_inv,
            "gamma_inv_dot": gamma_inv_dot, "K": K, "Aa": Aa,
            "Ktheta": Ktheta, "K_dot_stage_sum": Kdot, "Aa_dot_stage_sum": Aadot,
            "Ktheta_dot": Ktheta_dot, "chi": chi,
            "mass_native": mass, "mass_chain_rule": mass,
            "mass_dot_coordinate_chain_rule": mass_dot,
            "hamiltonian_constraint_residual": hamiltonian_constraint_residual,
            "momentum_constraint_residual": momentum_constraint_residual,
            **aggregate,
        },
        "constraint_conventions": {
            "hamiltonian_residual": (
                "vacuum.constraints(grid, geometry)['hamiltonian'] - 16*pi*total.rho"
            ),
            "momentum_residual": (
                "vacuum.constraints(grid, geometry)['momentum'] - 8*pi*total.j"
            ),
            "mass_gradient_hamiltonian_term": "R**2 * D_r(R) * H_residual / 4",
            "mass_gradient_momentum_term": "R**2 * (Rdot/alpha) * M_residual / 2",
            "radial_mass_rate_identity": (
                "epsilon_M = -R**2*alpha*gamma_rr_inv*D_r(R)*M_residual/2 "
                "- R**2*Rdot*C^r_r/2"
            ),
            "residual_sign_policy": "signed values retained; no physical residual threshold",
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
    assert set(result["constraint_profile_roundoff_checks"]) == {"2", "3"}
    assert all(
        item["bound_ratio"] <= 1.0
        for item in result["constraint_profile_roundoff_checks"].values()
    )
    assert "C_radial_mixed_unthresholded" in result["cells"]["2"]["independent_geometric_radial_einstein"]
    assert "hamiltonian_constraint_residual_signed" in result["cells"]["2"]
    assert "constraint_identity_operator_remainder_unthresholded" in result["cells"]["2"]


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



def _raw_em_stage_record(stage_name="start", charge_faces=None, rest_flux=None):
    centers = np.array([0.5, 1.5, 2.5], dtype=np.float64)
    dr = 1.0
    eta = 0.2
    charge_faces = np.asarray(
        [0.0, 0.08, 0.19, 0.27] if charge_faces is None else charge_faces,
        dtype=np.float64,
    )
    rest_flux = np.asarray(
        [0.0, 0.01, -0.02, 0.03] if rest_flux is None else rest_flux,
        dtype=np.float64,
    )
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
    lapse = np.array([0.95, 0.97, 0.99])
    sqrt_gamma = np.sqrt(a) * b * centers**2 / X**3
    lorentz_current_r = lorentz / (lapse * sqrt_gamma * weight * qcell)
    witness = {
        "grid_dr": dr,
        "charge_to_rest_mass": eta,
        "centers": centers.tolist(),
        "lapse_profile": lapse.tolist(),
        "lorentz_current_r_profile": lorentz_current_r.tolist(),
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
        "stage": stage_name,
        "t": 0.01 if stage_name == "start" else 0.02,
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
    assert result["full_functional_charge_tangent_plus_lorentz_remainder_unthresholded"] == pytest.approx(
        result["charge_current_tangent"] + result["matter_lorentz_power"], abs=1e-18
    )
    charge_cell_terms = np.asarray(
        result["field_charge_tangent_cell_terms_unthresholded"], dtype=np.float64
    )
    cell_current_terms = np.asarray(
        result["cell_current_tangent_terms_unthresholded"], dtype=np.float64
    )
    lorentz_cell_terms = np.asarray(
        result["lorentz_power_cell_terms_unthresholded"], dtype=np.float64
    )
    observed_cell_defect = np.asarray(
        result["face_cell_current_defect_by_cell_unthresholded"], dtype=np.float64
    )
    full_pairing_terms = np.asarray(
        result["full_functional_cellwise_face_plus_lorentz_terms_unthresholded"],
        dtype=np.float64,
    )
    assert observed_cell_defect == pytest.approx(
        charge_cell_terms - cell_current_terms, abs=1e-30
    )
    assert full_pairing_terms == pytest.approx(
        charge_cell_terms + lorentz_cell_terms, abs=1e-30
    )
    assert result["full_functional_cellwise_face_plus_lorentz_sum_unthresholded"] == pytest.approx(
        result["full_functional_charge_tangent_plus_lorentz_remainder_unthresholded"],
        abs=1e-25,
    )
    assert result["full_functional_cellwise_face_lorentz_absmax_cell"] == int(
        np.argmax(np.abs(full_pairing_terms))
    )
    assert result["algebraic_roundoff_checks"][
        "per_cell_defect_sums_to_full_exchange_remainder"
    ]["bound_ratio"] <= 1.0
    assert result["boundary_separated_exchange_remainder_unthresholded"] == pytest.approx(
        result["exchange_only_remainder_unthresholded"], abs=0.0
    )
    assert result["full_minus_boundary_separated_remainder_unthresholded"] == pytest.approx(
        result["outer_current_power"], abs=1e-18
    )
    assert result["exchange_only_remainder_unthresholded"] != 0.0
    assert all(check["bound_ratio"] <= 1.0 for check in result["algebraic_roundoff_checks"].values())



def test_em_charge_tangent_is_reconstructed_face_by_face_including_outermost_face():
    result = recompute_em_stage_exchange_witness(_raw_em_stage_record("primary_predictor"))
    face_terms = np.asarray(
        result["field_charge_tangent_by_face_unthresholded"], dtype=np.float64
    )
    cell_terms = np.asarray(
        result["field_charge_tangent_cell_terms_unthresholded"], dtype=np.float64
    )
    assert face_terms.size == cell_terms.size + 1
    assert face_terms[-1] == pytest.approx(
        result["outer_current_power"], rel=0.0, abs=1e-28
    )
    assert face_terms[1:-1].sum() == pytest.approx(
        result["field_charge_tangent_interior_faces_unthresholded"] and
        sum(result["field_charge_tangent_interior_faces_unthresholded"]),
        rel=1e-14, abs=1e-28,
    )
    assert face_terms.sum() == pytest.approx(
        result["charge_current_tangent"], rel=1e-14, abs=1e-28
    )
    assert cell_terms.sum() == pytest.approx(
        result["charge_current_tangent"], rel=1e-14, abs=1e-28
    )
    assert max(
        check["bound_ratio"]
        for check in result["algebraic_roundoff_checks"].values()
    ) <= 1.0
    assert result["physical_residual_threshold_applied"] is False


def test_same_stage_em_witness_joins_actual_charged_dust_energy_boundary_flux():
    em_stage = _raw_em_stage_record("primary_predictor")
    em = recompute_em_stage_exchange_witness(em_stage)
    raw = em_stage["stage_witness"]
    rest = np.asarray(raw["rest_mass_flux_faces"], dtype=np.float64)
    energy = np.array([0.0, 1.0e-4, -2.0e-4, 3.0e-4])
    boundary_rate = 4.0 * np.pi * float(energy[0] - energy[-1])
    source_stage = {
        "stage": em_stage["stage"],
        "t": em_stage["t"],
        "resolution": len(raw["centers"]),
        "species": {
            "charged_dust": {
                "face_flux_rest_mass": rest.tolist(),
                "face_flux_energy": energy.tolist(),
                "face_flux_radial_momentum": [0.0] * len(rest),
                "boundary_energy_rate": boundary_rate,
            }
        },
    }
    witness = recompute_charged_dust_boundary_witness(
        em, source_stage, {"boundary_energy_rate": boundary_rate}
    )
    assert witness["outer_charge_current_from_face_flux"] == pytest.approx(
        witness["outer_charge_current_used_by_maxwell"], rel=0.0, abs=1e-16
    )
    assert witness["charged_dust_boundary_energy_rate_from_faces"] == pytest.approx(
        -4.0 * np.pi * energy[-1], rel=0.0, abs=1e-18
    )
    assert witness["charged_dust_boundary_energy_rate_from_faces"] == pytest.approx(
        witness["charged_dust_boundary_energy_rate_reported"], rel=0.0, abs=1e-18
    )
    # Keep the two boundary conventions separately measured, not assumed equal.
    assert "em_outer_current_power_unthresholded" in witness
    assert "outer_energy_flux_rate_stored_convention" in witness
    assert witness["physical_residual_threshold_applied"] is False


def test_same_stage_boundary_join_rejects_rest_flux_or_time_mismatch():
    em_stage = _raw_em_stage_record("start")
    em = recompute_em_stage_exchange_witness(em_stage)
    raw = em_stage["stage_witness"]
    rest = np.asarray(raw["rest_mass_flux_faces"], dtype=np.float64)
    energy = np.array([0.0, 0.0, 0.0, 2.0e-4])
    boundary_rate = 4.0 * np.pi * float(energy[0] - energy[-1])
    base = {
        "stage": em_stage["stage"],
        "t": em_stage["t"],
        "resolution": len(raw["centers"]),
        "species": {
            "charged_dust": {
                "face_flux_rest_mass": rest.tolist(),
                "face_flux_energy": energy.tolist(),
                "face_flux_radial_momentum": [0.0] * len(rest),
                "boundary_energy_rate": boundary_rate,
            }
        },
    }
    wrong_flux = copy.deepcopy(base)
    wrong_flux["species"]["charged_dust"]["face_flux_rest_mass"][2] += 1.0e-5
    with pytest.raises(RuntimeError, match="same-stage charged-dust / Maxwell rest-mass face flux"):
        recompute_charged_dust_boundary_witness(em, wrong_flux, {"boundary_energy_rate": boundary_rate})
    wrong_time = copy.deepcopy(base)
    wrong_time["t"] += 1.0e-6
    with pytest.raises(RuntimeError, match="stage times differ"):
        recompute_charged_dust_boundary_witness(em, wrong_time, {"boundary_energy_rate": boundary_rate})


def test_same_stage_boundary_join_rejects_mismatched_outer_energy_flux():
    em_stage = _raw_em_stage_record("primary_predictor")
    em = recompute_em_stage_exchange_witness(em_stage)
    raw = em_stage["stage_witness"]
    rest = np.asarray(raw["rest_mass_flux_faces"], dtype=np.float64)
    energy = np.array([0.0, 0.0, 0.0, -1.0e-4])
    boundary_rate = 4.0 * np.pi * float(energy[0] - energy[-1])
    source_stage = {
        "stage": em_stage["stage"],
        "t": em_stage["t"],
        "resolution": len(raw["centers"]),
        "species": {
            "charged_dust": {
                "face_flux_rest_mass": rest.tolist(),
                "face_flux_energy": energy.tolist(),
                "face_flux_radial_momentum": [0.0] * len(rest),
                "boundary_energy_rate": boundary_rate + 1.0e-8,
            }
        },
    }
    with pytest.raises(RuntimeError, match="charged-dust boundary energy rate reconstruction"):
        recompute_charged_dust_boundary_witness(
            em, source_stage, {"boundary_energy_rate": boundary_rate + 1.0e-8}
        )


def test_independent_em_stage_recomputation_fails_closed_on_wrong_charge_tangent():
    broken = _raw_em_stage_record()
    broken["stage_witness"]["charge_face_dot"][1] += 1e-4
    with pytest.raises(RuntimeError, match="charge-face tangent reconstruction"):
        recompute_em_stage_exchange_witness(broken)



def _valid_em_charge_transition():
    dt = 0.02
    q_before = np.array([0.0, 0.08, 0.19, 0.27])
    flux_start = np.array([0.0, 0.01, -0.02, 0.03])
    start = _raw_em_stage_record(
        "start", charge_faces=q_before, rest_flux=flux_start
    )
    i_start = np.asarray(start["stage_witness"]["charge_current_faces"])
    q_predictor = q_before - dt * i_start

    flux_predictor = np.array([0.0, -0.004, 0.012, 0.025])
    predictor = _raw_em_stage_record(
        "primary_predictor", charge_faces=q_predictor, rest_flux=flux_predictor
    )
    i_predictor = np.asarray(predictor["stage_witness"]["charge_current_faces"])
    q_after = q_before - 0.5 * dt * (i_start + i_predictor)

    # Preserve a same-stage finite-volume charged-dust boundary witness next
    # to the Maxwell face-current witness. These synthetic energy fluxes are
    # intentionally manufactured; they test identity pairing, not physics.
    source_stage_rows = []
    charged_dust_stage_summaries = []
    for idx, em_stage in enumerate((start, predictor)):
        em_raw = em_stage["stage_witness"]
        outer_energy_flux = (0.003 if idx == 0 else -0.002)
        energy_faces = np.array([0.0, 0.0, 0.0, outer_energy_flux])
        rest_faces = np.asarray(em_raw["rest_mass_flux_faces"], dtype=np.float64)
        boundary_rate = 4.0 * np.pi * float(energy_faces[0] - energy_faces[-1])
        source_stage_rows.append({
            "stage": em_stage["stage"],
            "t": em_stage["t"],
            "resolution": len(em_raw["centers"]),
            "species": {
                "charged_dust": {
                    "face_flux_rest_mass": rest_faces.tolist(),
                    "face_flux_energy": energy_faces.tolist(),
                    "face_flux_radial_momentum": np.zeros_like(rest_faces).tolist(),
                    "boundary_energy_rate": boundary_rate,
                }
            },
        })
        charged_dust_stage_summaries.append({"boundary_energy_rate": boundary_rate})

    step = {
        "dt": dt,
        "t_before": 0.01,
        "t_after": 0.03,
        "em_charge_faces_before": q_before.tolist(),
        "em_charge_faces_after": q_after.tolist(),
        "em_field_energy_before": 0.01,
        "em_field_energy_after": 0.0101,
        "em_field_energy_rate_measured": 0.005,
        "em_stage_exchange_witnesses": [start, predictor],
        "source_stage_rows": source_stage_rows,
        "matter_species": {
            "charged_dust": {"stage_terms": charged_dust_stage_summaries}
        },
    }
    run = {"resolution": 3, "q_over_m": 0.2, "step_rows": [step]}
    return {"runs": [run], "fixed_N_temporal_refinement_runs": [dict(run, cfl=0.0025)]}


def test_independent_em_recompute_verifies_actual_euler_and_accepted_trapezoid_charge_updates():
    report = recompute_em_stage_exchange_campaign(_valid_em_charge_transition())
    assert report["interval_count"] == 2
    for interval in report["intervals"]:
        update = interval["charge_face_update"]
        assert update["algebraic_update_identity_passed"] is True
        assert update["scheme"].startswith("Q_after = Q_before")
        assert all(
            check["bound_ratio"] <= 1.0
            for check in update["roundoff_checks"].values()
        )
        witnesses = interval["same_stage_charged_dust_boundary_witnesses"]
        assert len(witnesses) == 2
        for witness in witnesses:
            assert witness["physical_residual_threshold_applied"] is False
            assert witness["charged_dust_boundary_energy_rate_from_faces"] == pytest.approx(
                witness["charged_dust_boundary_energy_rate_reported"], abs=1e-18
            )
            assert all(
                check["bound_ratio"] <= 1.0
                for check in witness["algebraic_roundoff_checks"].values()
            )


def test_independent_em_recompute_rejects_accepted_charge_update_mismatch():
    data = _valid_em_charge_transition()
    data["runs"][0]["step_rows"][0]["em_charge_faces_after"][2] += 1.0e-8
    with pytest.raises(RuntimeError, match="EM accepted trapezoidal charge update"):
        recompute_em_stage_exchange_campaign(data)


def test_independent_em_recompute_rejects_predictor_charge_update_mismatch():
    data = _valid_em_charge_transition()
    step = data["runs"][0]["step_rows"][0]
    predictor = step["em_stage_exchange_witnesses"][1]
    altered_q = np.asarray(predictor["stage_witness"]["charge_faces"], dtype=float)
    altered_q[2] += 1.0e-8
    # Regenerate a self-consistent predictor witness so the temporal-state
    # contract (not a stale local algebraic field) is the failing check.
    flux = predictor["stage_witness"]["rest_mass_flux_faces"]
    altered_predictor = _raw_em_stage_record(
        "primary_predictor", charge_faces=altered_q, rest_flux=flux
    )
    step["em_stage_exchange_witnesses"][1] = altered_predictor
    with pytest.raises(RuntimeError, match="EM Euler predictor charge update"):
        recompute_em_stage_exchange_campaign(data)



def test_independent_metric_radial_einstein_residual_converges_on_static_schwarzschild():
    def compute(n):
        r_max = 20.0
        dr = r_max / n
        indices = list(range(n))
        r = (np.arange(n, dtype=float) + 0.5) * dr
        psi = 1.0 + 0.5 / r
        X = psi**-2
        alpha = (1.0 - 0.5 / r) / (1.0 + 0.5 / r)
        a = np.ones_like(r)
        b = np.ones_like(r)
        areal = r * np.sqrt(b) / X
        # Exact time-independent vacuum Schwarzschild slice in isotropic
        # coordinates: K=K^r_r=0 and the radial Einstein residual vanishes.
        profiles = {
            "r": r, "a": a, "b": b, "X": X, "alpha": alpha,
            "beta": np.zeros_like(r), "K": np.zeros_like(r),
            "Aa": np.zeros_like(r), "pr": np.zeros_like(r),
        }
        target_r = 8.0
        index = int(np.argmin(np.abs(r - target_r)))
        result = _geometric_radial_einstein_rr_at_cell(
            profiles, indices, index, dr,
            np.zeros_like(r), np.zeros_like(r),
        )
        return abs(result["C_radial_mixed_unthresholded"])

    coarse = compute(80)
    fine = compute(160)
    assert coarse > 0.0
    assert fine < coarse / 8.0




def test_combined_matter_em_candidates_keep_outer_term_choice_explicit():
    dt = 0.5
    measured_matter = (2.1 - 2.0) / dt
    measured_field = (3.2 - 3.0) / dt
    measured_combined = measured_matter + measured_field
    boundary, geometry, dm_exchange, metric, outer = 1.5, -0.2, 0.1, 0.3, -0.05
    known_with_outer = boundary + geometry + dm_exchange + metric + outer
    row = {
        "dt": dt,
        "matter_species": {
            "charged_dust": {"energy_before": 2.0, "energy_after": 2.1}
        },
        "em_field_energy_before": 3.0,
        "em_field_energy_after": 3.2,
        "combined_matter_plus_em_energy_rate_measured": measured_combined,
        "combined_rate_from_named_boundary_and_source_terms": known_with_outer,
        "unassigned_matter_plus_em_rate_remainder_not_global_conservation_residual": (
            measured_combined - known_with_outer
        ),
    }
    source_rows = [
        {
            "species": {
                "charged_dust": {
                    "boundary_energy_rate": 1.4,
                    "geometric_energy_source_rate": -0.1,
                    "dm_scalar_exchange_energy_rate": 0.08,
                }
            }
        },
        {
            "species": {
                "charged_dust": {
                    "boundary_energy_rate": 1.6,
                    "geometric_energy_source_rate": -0.3,
                    "dm_scalar_exchange_energy_rate": 0.12,
                }
            }
        },
    ]
    em_stage_rows = [
        {"metric_tangent": 0.2, "outer_current_power": -0.04},
        {"metric_tangent": 0.4, "outer_current_power": -0.06},
    ]

    result = recompute_combined_matter_em_accounting_candidates(
        row, source_rows, em_stage_rows
    )
    no_outer = result["remainder_without_separate_em_outer_candidate_unthresholded"]
    with_outer = result["remainder_with_separate_em_outer_candidate_unthresholded"]
    assert result["known_terms_without_separate_em_outer_unthresholded"] == pytest.approx(1.7)
    assert result["known_terms_with_separate_em_outer_unthresholded"] == pytest.approx(1.65)
    assert no_outer == pytest.approx(measured_combined - 1.7)
    assert with_outer == pytest.approx(measured_combined - 1.65)
    assert no_outer - with_outer == pytest.approx(outer)
    assert result["remainder_candidate_difference_minus_em_outer_unthresholded"] == pytest.approx(0.0)
    assert result["physical_residual_threshold_applied"] is False



def test_em_exchange_cli_summary_keeps_physical_remainders_unthresholded():
    report = {
        "stages": [
            {
                "run_group": "runs", "run_index": 0, "resolution": 24,
                "q_over_m": 0.25, "cfl": None,
                "stage": "predictor", "t": 0.04, "interval_index": 1,
                "charge_current_tangent": -1.0e-7,
                "cell_centered_current_tangent": -8.0e-8,
                "matter_lorentz_power": 8.0e-8,
                "outer_current_power": -3.0e-9,
                "full_functional_charge_tangent_plus_lorentz_remainder_unthresholded": -2.0e-8,
                "boundary_separated_exchange_remainder_unthresholded": -1.7e-8,
                "boundary_separation_adjustment_unthresholded": 3.0e-9,
                "full_minus_boundary_separated_remainder_unthresholded": -3.0e-9,
                "full_functional_exchange_relative_remainder_unthresholded": 2.0e-8 / 1.8e-7,
                "full_functional_cellwise_face_lorentz_absmax_unthresholded": 4.0e-9,
                "full_functional_cellwise_face_lorentz_absmax_cell": 7,
                "full_functional_cellwise_face_plus_lorentz_sum_unthresholded": -2.0e-8,
                "full_functional_cellwise_face_plus_lorentz_sum_minus_scalar_unthresholded": 2.0e-25,
                "face_minus_cell_current_tangent_mismatch_unthresholded": -2.0e-8,
                "exchange_from_face_cell_current_mismatch_unthresholded": -1.7e-8,
                "exchange_identity_gap_unthresholded": 1.0e-24,
                "exchange_only_remainder_unthresholded": -1.7e-8,
                "exchange_only_relative_remainder_unthresholded": 0.017,
                "algebraic_roundoff_checks": {"identity": {"bound_ratio": 0.5}},
            },
            {
                "run_group": "runs", "run_index": 0, "resolution": 24,
                "q_over_m": 0.25, "cfl": None,
                "stage": "start", "t": 0.03, "interval_index": 0,
                "charge_current_tangent": -5.0e-8,
                "cell_centered_current_tangent": -4.0e-8,
                "matter_lorentz_power": 4.0e-8,
                "outer_current_power": -2.0e-9,
                "full_functional_charge_tangent_plus_lorentz_remainder_unthresholded": -1.0e-8,
                "boundary_separated_exchange_remainder_unthresholded": -8.0e-9,
                "boundary_separation_adjustment_unthresholded": 2.0e-9,
                "full_minus_boundary_separated_remainder_unthresholded": -2.0e-9,
                "full_functional_exchange_relative_remainder_unthresholded": 1.0e-8 / 9.0e-8,
                "full_functional_cellwise_face_lorentz_absmax_unthresholded": 2.0e-9,
                "full_functional_cellwise_face_lorentz_absmax_cell": 3,
                "full_functional_cellwise_face_plus_lorentz_sum_unthresholded": -1.0e-8,
                "full_functional_cellwise_face_plus_lorentz_sum_minus_scalar_unthresholded": 0.0,
                "face_minus_cell_current_tangent_mismatch_unthresholded": -1.0e-8,
                "exchange_from_face_cell_current_mismatch_unthresholded": -8.0e-9,
                "exchange_identity_gap_unthresholded": 0.0,
                "exchange_only_remainder_unthresholded": -8.0e-9,
                "exchange_only_relative_remainder_unthresholded": 0.008,
                "algebraic_roundoff_checks": {"identity": {"bound_ratio": 0.75}},
            },
        ],
        "intervals": [
            {
                "run_group": "runs", "run_index": 0,
                "field_energy_time_update_remainder_unthresholded": -3.0e-8,
                "full_functional_exchange_start_predictor_average_remainder_unthresholded": -5.0e-9,
                "boundary_separation_adjustment_start_predictor_average_unthresholded": 1.0e-8,
                "full_minus_boundary_separated_interval_remainder_unthresholded": -1.0e-8,
                "exchange_only_start_predictor_average_remainder_unthresholded": 5.0e-9,
                "combined_matter_em_accounting_candidates": {
                    "remainder_without_separate_em_outer_candidate_unthresholded": 4.0e-9,
                    "remainder_with_separate_em_outer_candidate_unthresholded": 5.0e-9,
                    "em_outer_current_power_start_predictor_average_reconstructed_unthresholded": -1.0e-9,
                },
            }
        ],
    }
    summary = _summarize_em_stage_exchange(report)
    assert len(summary["groups"]) == 1
    row = summary["groups"][0]
    assert row["exchange_only_remainder_absmax_unthresholded"] == pytest.approx(1.7e-8)
    assert row["full_functional_charge_tangent_plus_lorentz_remainder_absmax_unthresholded"] == pytest.approx(2.0e-8)
    assert row["full_functional_exchange_relative_remainder_absmax_unthresholded"] == pytest.approx(1.0e-8 / 9.0e-8)
    assert row["full_functional_cellwise_face_lorentz_absmax_unthresholded"] == pytest.approx(4.0e-9)
    assert row["worst_cellwise_face_lorentz_stage_terms"]["stage"] == "predictor"
    assert row["worst_cellwise_face_lorentz_stage_terms"]["full_functional_cellwise_face_lorentz_absmax_cell"] == 7
    assert row["boundary_separation_adjustment_absmax_unthresholded"] == pytest.approx(3.0e-9)
    assert row["exchange_only_relative_remainder_absmax_unthresholded"] == pytest.approx(0.017)
    assert row["field_energy_time_update_remainder_absmax_unthresholded"] == pytest.approx(3.0e-8)
    assert row["interval_exchange_average_remainder_absmax_unthresholded"] == pytest.approx(5.0e-9)
    assert row["full_functional_interval_exchange_average_remainder_absmax_unthresholded"] == pytest.approx(5.0e-9)
    assert row["interval_boundary_separation_adjustment_absmax_unthresholded"] == pytest.approx(1.0e-8)
    assert row["combined_matter_em_remainder_without_separate_em_outer_candidate_absmax_unthresholded"] == pytest.approx(4.0e-9)
    assert row["combined_matter_em_remainder_with_separate_em_outer_candidate_absmax_unthresholded"] == pytest.approx(5.0e-9)
    assert row["combined_matter_em_outer_current_power_absmax_unthresholded"] == pytest.approx(1.0e-9)
    assert row["combined_matter_em_candidates_physical_residual_threshold_applied"] is False
    assert row["algebraic_roundoff_bound_ratio_absmax"] == pytest.approx(0.75)
    assert row["face_minus_cell_current_tangent_mismatch_absmax_unthresholded"] == pytest.approx(2.0e-8)
    assert row["exchange_identity_gap_absmax_unthresholded"] == pytest.approx(1.0e-24)
    assert row["worst_exchange_stage_terms"]["stage"] == "predictor"
    assert row["worst_exchange_stage_terms"]["exchange_identity_gap_unthresholded"] == pytest.approx(1.0e-24)
    assert row["worst_exchange_stage_terms"]["full_minus_boundary_separated_remainder_unthresholded"] == pytest.approx(-3.0e-9)
    assert row["physical_residual_threshold_applied"] is False



def test_independent_fixed_r_reconstruction_rejects_changed_archived_momentum_constraint_sign():
    stage = _raw_stage()
    broken = copy.deepcopy(stage)
    broken["tube_fixed_R_balance_raw_witness"]["profiles"][
        "momentum_constraint_residual"
    ][2] *= -1.0
    with pytest.raises(RuntimeError, match="archived signed momentum constraint residual at cell 2"):
        recompute_arm(broken)


def test_independent_fixed_r_reconstruction_requires_signed_constraint_profiles_and_conventions():
    stage = _raw_stage()
    del stage["tube_fixed_R_balance_raw_witness"]["profiles"]["hamiltonian_constraint_residual"]
    with pytest.raises(RuntimeError, match="lacks signed native-cell constraint profiles"):
        recompute_arm(stage)

    stage = _raw_stage()
    stage["tube_fixed_R_balance_raw_witness"]["constraint_conventions"]["momentum_residual"] = (
        "sign-ambiguous momentum residual"
    )
    with pytest.raises(RuntimeError, match="raw constraint convention 'momentum_residual'"):
        recompute_arm(stage)


def _manufactured_balanced_em_exchange_record():
    """Build a known-balanced donor/recipient exchange without engine imports.

    The face charge/current and geometry determine the field's charge tangent
    and the declared outer-face power. The manufactured Lorentz source is then
    chosen so charge-tangent + matter power - outer power is exactly zero.
    Its cell-current projection is built from the same source so that the
    independent checker has a fully self-consistent raw-operand witness.
    This fixture validates the ledger convention only; it is not production
    evidence and must never be used to threshold a production residual.
    """
    record = _raw_em_stage_record("manufactured_balanced")
    raw = record["stage_witness"]
    centers = np.asarray(raw["centers"], dtype=np.float64)
    dr = float(raw["grid_dr"])
    q_face = np.asarray(raw["charge_faces"], dtype=np.float64)
    current = np.asarray(raw["charge_current_faces"], dtype=np.float64)
    q_cell = 0.5 * (q_face[:-1] + q_face[1:])
    q_cell_dot = -0.5 * (current[:-1] + current[1:])
    geometry = raw["geometry"]
    X = np.asarray(geometry["X"], dtype=np.float64)
    a = np.asarray(geometry["a"], dtype=np.float64)
    b = np.asarray(geometry["b"], dtype=np.float64)
    lapse = np.asarray(raw["lapse_profile"], dtype=np.float64)
    weight = X * np.sqrt(a) / (centers**2 * b)
    charge_tangent = dr * float(np.sum(weight * q_cell * q_cell_dot))
    outer_power = -0.5 * dr * weight[-1] * q_cell[-1] * current[-1]

    # Set recipient gain to the exact opposite of the field charge tangent,
    # with the independently declared outer-current term retained explicitly.
    lorentz_power = -charge_tangent + outer_power
    lorentz_energy = np.zeros_like(q_cell)
    lorentz_energy[-1] = lorentz_power / (4.0 * np.pi * dr)
    sqrt_gamma = np.sqrt(a) * b * centers**2 / X**3
    lorentz_current_r = lorentz_energy / (
        lapse * sqrt_gamma * weight * q_cell
    )

    raw["lorentz_energy_source_profile"] = lorentz_energy.tolist()
    raw["lorentz_current_r_profile"] = lorentz_current_r.tolist()
    raw["matter_lorentz_power"] = float(lorentz_power)
    raw["outer_current_power"] = float(outer_power)
    record["matter_lorentz_power"] = float(lorentz_power)
    record["outer_current_power"] = float(outer_power)
    return record


def test_manufactured_em_donor_recipient_exchange_closes_at_roundoff():
    result = recompute_em_stage_exchange_witness(
        _manufactured_balanced_em_exchange_record()
    )
    remainder = result["exchange_only_remainder_unthresholded"]
    scale = (
        abs(result["charge_current_tangent"])
        + abs(result["matter_lorentz_power"])
        + abs(result["outer_current_power"])
    )
    # This roundoff comparison applies only to this manufactured exact-balance
    # fixture. No such criterion is applied to the measured production ledger.
    assert abs(remainder) <= 64.0 * np.finfo(np.float64).eps * scale
    # This fixture is balanced under the boundary-separated convention by
    # construction. The full-functional pairing differs by the explicit
    # outer-face summand; the test preserves that distinction rather than
    # declaring either convention physically correct.
    assert result["full_functional_charge_tangent_plus_lorentz_remainder_unthresholded"] == pytest.approx(
        result["outer_current_power"], abs=1e-18
    )
    assert result["boundary_separated_exchange_remainder_unthresholded"] == pytest.approx(
        remainder, abs=0.0
    )
    assert result["physical_residual_threshold_applied"] is False
    assert all(
        check["bound_ratio"] <= 1.0
        for check in result["algebraic_roundoff_checks"].values()
    )
    cell_terms = np.asarray(
        result["full_functional_cellwise_face_plus_lorentz_terms_unthresholded"],
        dtype=np.float64,
    )
    assert float(np.sum(cell_terms)) == pytest.approx(
        result["full_functional_cellwise_face_plus_lorentz_remainder_unthresholded"],
        abs=1e-25,
    )
    assert result["full_functional_cellwise_face_lorentz_vs_defect_gap_unthresholded"] == pytest.approx(
        np.zeros_like(cell_terms).tolist(), abs=1e-30
    )


def test_manufactured_nonzero_em_exchange_remains_an_unthresholded_measurement():
    record = _manufactured_balanced_em_exchange_record()
    raw = record["stage_witness"]
    injected_power_density = 1.0e-6
    raw["lorentz_energy_source_profile"][0] += injected_power_density
    centers = np.asarray(raw["centers"], dtype=np.float64)
    geometry = raw["geometry"]
    X = np.asarray(geometry["X"], dtype=np.float64)
    a = np.asarray(geometry["a"], dtype=np.float64)
    b = np.asarray(geometry["b"], dtype=np.float64)
    lapse = np.asarray(raw["lapse_profile"], dtype=np.float64)
    q_face = np.asarray(raw["charge_faces"], dtype=np.float64)
    q_cell = 0.5 * (q_face[:-1] + q_face[1:])
    weight = X * np.sqrt(a) / (centers**2 * b)
    sqrt_gamma = np.sqrt(a) * b * centers**2 / X**3
    lorentz_energy = np.asarray(raw["lorentz_energy_source_profile"], dtype=np.float64)
    lorentz_current_r = lorentz_energy / (
        lapse * sqrt_gamma * weight * q_cell
    )
    raw["lorentz_current_r_profile"] = lorentz_current_r.tolist()
    power_increment = 4.0 * np.pi * float(raw["grid_dr"]) * injected_power_density
    raw["matter_lorentz_power"] += power_increment
    record["matter_lorentz_power"] += power_increment

    result = recompute_em_stage_exchange_witness(record)
    assert result["exchange_only_remainder_unthresholded"] == pytest.approx(
        power_increment, rel=1.0e-10, abs=1.0e-18
    )
    assert result["physical_residual_threshold_applied"] is False
    assert all(
        check["bound_ratio"] <= 1.0
        for check in result["algebraic_roundoff_checks"].values()
    )
