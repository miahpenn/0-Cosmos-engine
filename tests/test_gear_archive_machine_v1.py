"""Focused contracts for the versioned archive-machine runner."""
import math

import pytest

from engine.cosmos import CosmosParams
from engine.local import LocalParams, local_rhs
from research.gear_archive_machine_v1 import (
    ARCHIVE_OPEN_GATES,
    LOCAL_ENERGY_CATEGORY_CONTRACT,
    archive_docx_local_rhs,
    run_machine,
    rk4_step,
)


def test_rk4_step_integrates_a_scalar_rhs_without_hidden_terms():
    result = rk4_step(lambda y: (y[0],), (1.0,), 0.01)
    assert len(result) == 1
    assert result[0] == pytest.approx(math.exp(0.01), rel=1.0e-9)


def test_archive_docx_local_rhs_matches_current_rhs_when_eta_is_zero():
    params = LocalParams(eta_G=0.0)
    state = (0.02, -0.003, 0.3, 0.4)
    archive = archive_docx_local_rhs(state, params)
    current = local_rhs(state, params, H=0.0)
    assert archive == pytest.approx(current, rel=0.0, abs=1.0e-15)


def test_archive_docx_local_rhs_contains_declared_signed_damping_term():
    params = LocalParams(eta_G=0.05)
    state = (0.02, -0.003, 0.3, 0.4)
    archive = archive_docx_local_rhs(state, params)
    current = local_rhs(state, params, H=0.0)
    assert archive[:3] == pytest.approx(current[:3])
    assert archive[3] == pytest.approx(
        current[3] - params.eta_G * state[2] ** 2 * state[3],
        rel=0.0,
        abs=1.0e-15,
    )
    assert archive[3] < current[3]


def test_runner_defaults_to_archive_local_rhs_and_allows_current_source_comparison():
    initial = (0.02, 0.0, 0.3, 0.4)
    archived = run_machine(
        local_tau_final=0.05,
        local_dt=0.01,
        cosmos_t_final=0.0,
        local_initial_state=initial,
    )
    current = run_machine(
        local_tau_final=0.05,
        local_dt=0.01,
        cosmos_t_final=0.0,
        local_initial_state=initial,
        local_equation_set="current_source",
    )
    assert archived["local"]["equation_set"] == "archive_docx"
    assert "preserved GEAR archive DOCX" in archived["local"]["equation_source"]
    assert current["local"]["equation_set"] == "current_source"
    assert current["local"]["equation_source"].startswith("engine.local.local_rhs")
    assert archived["local"]["state"][-1]["Dd"] != current["local"]["state"][-1]["Dd"]


def test_runner_preserves_two_independent_time_axes_and_ledger_shapes():
    result = run_machine(
        local_tau_final=0.02,
        local_dt=0.01,
        cosmos_t_final=0.015,
        cosmos_dt=0.01,
        sample_every=1,
    )

    assert result["status"] == "research_only_independent_lanes"
    assert result["scientific_admission"] == "not_claimed"
    assert result["lane_contract"]["cross_lane_source_terms"] is False

    local = result["local"]
    cosmos = result["cosmos"]
    assert local["time_name"] == "tau_local"
    assert cosmos["time_name"] == "t_cosmos"
    assert local["time"] == [0.0, 0.01, 0.02]
    assert cosmos["time"] == [0.0, 0.01, 0.015]
    assert len(local["state"]) == len(local["energy_and_activity_ledgers"])
    assert len(local["state"]) == len(local["oscillator_phase_theta"])
    assert len(cosmos["state"]) == len(cosmos["component_and_constraint_ledgers"])
    assert set(ARCHIVE_OPEN_GATES).issubset(result["open_gates"])


def test_expansion_and_contraction_select_only_initial_h_sign():
    expanding = run_machine(local_tau_final=0.0, cosmos_t_final=0.0, cosmic_branch=1)
    contracting = run_machine(local_tau_final=0.0, cosmos_t_final=0.0, cosmic_branch=-1)

    ex_state = expanding["cosmos"]["state"][0]
    co_state = contracting["cosmos"]["state"][0]
    assert ex_state["H"] > 0.0
    assert co_state["H"] == pytest.approx(-ex_state["H"])
    for key in ("a", "phi", "pi_phi", "rho_dm", "rho_b", "rho_r", "rho_shear"):
        assert co_state[key] == pytest.approx(ex_state[key])
    assert expanding["cosmos"]["initial_branch_sign"] == 1
    assert contracting["cosmos"]["initial_branch_sign"] == -1


def test_runner_rejects_a_local_parameter_set_that_breaks_gear03_lock():
    with pytest.raises(ValueError, match=r"alpha_D\*\*2 == g-kappa"):
        run_machine(
            local_tau_final=0.0,
            cosmos_t_final=0.0,
            local_params=LocalParams(kappa=0.1, g=1.0, alpha_D=1.0),
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"local_tau_final": -1.0},
        {"local_dt": 0.0},
        {"cosmos_t_final": math.inf},
        {"cosmic_branch": 0},
        {"sample_every": 0},
        {"local_equation_set": "unknown"},
    ],
)
def test_runner_fails_closed_on_invalid_controls(kwargs):
    options = {"local_tau_final": 0.0, "cosmos_t_final": 0.0}
    options.update(kwargs)
    with pytest.raises(ValueError):
        run_machine(**options)


def test_zero_length_run_still_reports_starting_state_and_open_gates():
    result = run_machine(local_tau_final=0.0, cosmos_t_final=0.0)
    assert result["local"]["time"] == [0.0]
    assert result["cosmos"]["time"] == [0.0]
    assert result["open_gates"]["D_mode_to_localized_matter"].startswith("OPEN")
    assert result["open_gates"]["local_source_to_global_cosmos_Q"].startswith("OPEN")


def test_cosmos_parameters_can_be_injected_without_changing_lane_contract():
    params = CosmosParams(beta=-0.04)
    result = run_machine(local_tau_final=0.0, cosmos_t_final=0.0, cosmos_params=params)
    assert result["cosmos"]["state"][0]["phi"] == pytest.approx(0.179055)
    assert result["lane_contract"]["local_time"] != result["lane_contract"]["cosmos_time"]


def test_growing_D_mode_separates_signed_energy_from_activity_proxy():
    """Pin the conditional signed-potential identity without inferring matter."""
    params = LocalParams(kappa=0.0, g=1.0, alpha_D=1.0)
    D = 0.3
    result = run_machine(
        local_tau_final=0.0,
        cosmos_t_final=0.0,
        local_params=params,
        local_initial_state=(0.0, 0.0, D, params.alpha_D * D),
    )

    ledger = result["local"]["energy_and_activity_ledgers"][0]
    # Conditional canonical scalar energy with the signed inverted potential:
    # 1/2 * Dd**2 - 1/2 * alpha_D**2 * D**2 = 0 on Dd=alpha_D*D.
    assert ledger["E_D_raw"] == pytest.approx(0.0, abs=1.0e-16)
    assert ledger["generated_activity_proxy"] == pytest.approx(
        params.alpha_D**2 * D**2, rel=0.0, abs=1.0e-15
    )
    assert ledger["generated_activity_proxy"] > abs(ledger["E_D_raw"])

    contract = result["interpretation_contract"]
    assert contract["D_energy_categories"] == LOCAL_ENERGY_CATEGORY_CONTRACT
    assert contract["D_mode_to_localized_matter_applied"] is False
    assert contract["local_source_to_global_cosmos_Q_applied"] is False
    assert contract["stress_energy_normalization"].startswith("OPEN")
    assert "not canonical" in contract["D_energy_categories"]["generated_activity_proxy"]
