"""Campaign-level contract for the opt-in, research-only radial Maxwell lane."""
import json

import numpy as np
import pytest

from engine.campaign import CampaignConfig
from engine.run_production import run_campaign


def _campaign(tmp_path, **overrides):
    values = dict(
        resolutions=(24,),
        r_max=24.0,
        cfl=0.01,
        final_time=0.01,
        checkpoint_interval=0.01,
        output_dir=str(tmp_path),
        radial_maxwell_charge_to_mass=0.25,
        # Caller-specified density knots are serialized and interpolated at
        # each resolution's cell centres; this is a test fixture, not a GEAR
        # carrier prescription.
        charged_dust_density_profile=((0.0, 2.0e-5), (24.0, 2.0e-5)),
    )
    values.update(overrides)
    return CampaignConfig(**values)


def test_radial_maxwell_campaign_inputs_fail_closed():
    with pytest.raises(ValueError, match="both"):
        CampaignConfig(radial_maxwell_charge_to_mass=0.25)

    with pytest.raises(ValueError, match="both"):
        CampaignConfig(
            charged_dust_density_profile=((0.0, 1.0e-6), (80.0, 1.0e-6))
        )

    with pytest.raises(ValueError, match="include radius 0 and radius r_max"):
        CampaignConfig(
            resolutions=(24,),
            r_max=24.0,
            charged_dust_density_profile=None,
            radial_maxwell_charge_to_mass=0.25,
        ) if False else _campaign(
            pytest.TempPathFactory(),  # replaced below without filesystem use
            charged_dust_density_profile=((1.0, 1.0e-6), (24.0, 1.0e-6)),
        )


def test_radial_maxwell_campaign_persists_accepted_state_and_stage_witnesses(tmp_path):
    config = _campaign(tmp_path)
    results = run_campaign(config=config)

    assert len(results) == 1
    assert results[0].status == "completed"
    assert results[0].t_final == pytest.approx(config.final_time)
    assert results[0].checkpoint_count == 1

    run_dir = tmp_path / "N24"
    summary = json.loads((tmp_path / "campaign_summary.json").read_text())
    em_summary = summary["physics_extensions"]["radial_maxwell"]
    assert em_summary["enabled"] is True
    assert em_summary["status"] == "research_only_not_admitted"
    assert "curved-slice donor-plus-field balance" in em_summary["open_gates"]

    ledger = json.loads((run_dir / "ledgers.json").read_text())
    history = ledger["history"]
    assert len(history) >= 2
    final = history[-1]
    for key in (
        "em_start_stage_current_absmax",
        "em_predictor_stage_current_absmax",
        "em_start_stage_lorentz_momentum_source_absmax",
        "em_predictor_stage_lorentz_momentum_source_absmax",
        "em_accepted_charge_change_absmax",
        "em_accepted_gauss_shell_residual_absmax",
        "em_stress_energy_density_absmax",
    ):
        assert key in final
        assert np.isfinite(final[key])
    assert final["em_predictor_stage_current_absmax"] > 0.0
    assert final["em_predictor_stage_lorentz_momentum_source_absmax"] > 0.0
    assert final["em_accepted_charge_change_absmax"] > 0.0

    checkpoints = list(run_dir.glob("checkpoint_t*.npz"))
    assert len(checkpoints) == 1
    with np.load(checkpoints[0], allow_pickle=False) as checkpoint:
        expected = {
            "charged_dust_rest",
            "charged_dust_energy_t",
            "charged_dust_momentum_r",
            "em_charge_to_rest_mass",
            "em_enclosed_charge_faces",
            "em_electric_field_hat",
            "em_stress_rho",
            "em_stress_pr",
            "em_stress_pt",
            "em_stress_j",
            "em_gauss_shell_residual",
            "em_charge_density",
            "em_current_r",
            "em_lorentz_energy_source",
            "em_lorentz_momentum_source_r",
        }
        assert expected.issubset(set(checkpoint.files))
        assert float(checkpoint["em_charge_to_rest_mass"]) == pytest.approx(0.25)
        assert checkpoint["em_enclosed_charge_faces"][0] == 0.0
        assert np.any(np.abs(checkpoint["em_electric_field_hat"]) > 0.0)
        assert np.any(np.abs(checkpoint["em_current_r"]) > 0.0)
        assert np.all(np.isfinite(checkpoint["em_gauss_shell_residual"]))


def test_default_campaign_config_keeps_the_em_extension_disabled():
    config = CampaignConfig()
    assert config.radial_maxwell_charge_to_mass is None
    assert config.charged_dust_density_profile is None
