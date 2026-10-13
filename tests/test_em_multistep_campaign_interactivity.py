"""Three-step full-sector regression for the opt-in radial-Maxwell campaign.

This checks repeated accepted evolution and checkpoint consistency on the
research-only candidate, not physical admission or long-horizon stability.
"""
import json

import numpy as np
import pytest

from engine.campaign import CampaignConfig
from engine.production_kernel import V55ProductionKernel
from engine.run_production import run_campaign


def _three_step_campaign(tmp_path):
    return CampaignConfig(
        resolutions=(24,),
        r_max=24.0,
        cfl=0.01,  # dr=1, so the intended bounded step is dt=0.01.
        final_time=0.03,
        checkpoint_interval=0.01,
        output_dir=str(tmp_path),
        radial_maxwell_charge_to_mass=0.25,
        # Explicit test-carrier profile; this is not inferred from GEAR.
        charged_dust_density_profile=((0.0, 2.0e-5), (24.0, 2.0e-5)),
    )


def _max_abs(*arrays):
    return max(float(np.max(np.abs(np.asarray(values, dtype=float)))) for values in arrays)


def test_opt_in_em_campaign_remains_interactive_over_three_accepted_steps(
    tmp_path, monkeypatch
):
    kernel = V55ProductionKernel()
    rhs_trace = []
    geometry_exchange = []

    rhs_original = kernel._rhs
    geometry_original = kernel._geometry_stage_terms

    def record_rhs(stage):
        scalar_rhs, matter_rhs, metric_derivatives = rhs_original(stage)
        scalar_names = ("S", "PS", "D", "PD", "phi", "Pi")
        matter_names = ("dark_matter", "baryons", "radiation", "charged_dust")
        row = {
            "scalar": _max_abs(*(getattr(scalar_rhs, name) for name in scalar_names)),
        }
        for name in matter_names:
            conserved_rhs = matter_rhs[name]
            row[name] = _max_abs(
                conserved_rhs.rest,
                conserved_rhs.energy_t,
                conserved_rhs.momentum_r,
            )
        for component in (metric_derivatives.time, metric_derivatives.radial):
            for values in component.values():
                assert np.all(np.isfinite(np.asarray(values, dtype=float)))
        assert all(np.isfinite(value) for value in row.values())
        rhs_trace.append(row)
        return scalar_rhs, matter_rhs, metric_derivatives

    def record_geometry(grid, geometry, scalars, matter, electromagnetic):
        terms = geometry_original(
            grid, geometry, scalars, matter, electromagnetic
        )
        if electromagnetic is not None:
            # Same exact stage and same matter state, with only the EM
            # projection removed, so this measures whether field stress reaches
            # the real geometry consumer on successive steps.
            baseline = geometry_original(
                grid, geometry, scalars, matter, None
            )
            projection = electromagnetic.stress_projection(grid, geometry)
            geometry_exchange.append({
                "rho_max": float(np.max(np.abs(projection["rho"]))),
                "Aa_l3_delta": _max_abs(
                    terms["primary_l3"]["Aa"] - baseline["primary_l3"]["Aa"]
                ),
                "K_l3_delta": _max_abs(
                    terms["primary_l3"]["K"] - baseline["primary_l3"]["K"]
                ),
            })
        return terms

    # Instance-level wrappers observe real internal call sites and leave the
    # candidate implementation untouched.
    monkeypatch.setattr(kernel, "_rhs", record_rhs)
    monkeypatch.setattr(kernel, "_geometry_stage_terms", record_geometry)

    config = _three_step_campaign(tmp_path)
    results = run_campaign(kernel=kernel, config=config)

    assert len(results) == 1
    result = results[0]
    assert result.status == "completed"
    assert result.steps == 3
    assert result.t_final == pytest.approx(0.03, abs=1.0e-12)
    assert result.checkpoint_count == 3

    # Two real matter/scalar RHS stages per accepted PIRK step. Requiring
    # nonzero RHS norms ensures these active sectors were actually evaluated,
    # while finite checks make silent non-finite propagation fail closed.
    assert len(rhs_trace) == 6
    for sector in ("scalar", "dark_matter", "baryons", "radiation", "charged_dust"):
        values = [row[sector] for row in rhs_trace]
        assert all(np.isfinite(value) for value in values), sector
        assert max(values) > 0.0, sector

    # Three geometry calls per step in the PIRK path (accepted start,
    # explicit predictor, primary predictor); all must include a nontrivial EM
    # stress contribution in the actual primary-L3 consumer.
    assert len(geometry_exchange) == 9
    for row in geometry_exchange:
        assert all(np.isfinite(value) for value in row.values())
        assert row["rho_max"] > 0.0
        assert row["Aa_l3_delta"] > 0.0
        assert row["K_l3_delta"] > 0.0

    run_dir = tmp_path / "N24"
    summary = json.loads((tmp_path / "campaign_summary.json").read_text())
    extension = summary["physics_extensions"]["radial_maxwell"]
    assert extension["enabled"] is True
    assert extension["status"] == "research_only_not_admitted"

    ledger = json.loads((run_dir / "ledgers.json").read_text())
    history = ledger["history"]
    accepted_em_rows = [
        row for row in history if "em_predictor_stage_current_absmax" in row
    ]
    assert len(history) == 4  # initial diagnostics plus three accepted steps
    assert len(accepted_em_rows) == 3
    for row in accepted_em_rows:
        for key in (
            "em_predictor_stage_current_absmax",
            "em_start_stage_lorentz_momentum_source_absmax",
            "em_predictor_stage_lorentz_momentum_source_absmax",
            "em_accepted_charge_change_absmax",
            "em_accepted_gauss_shell_residual_absmax",
            "em_stress_energy_density_absmax",
        ):
            assert key in row
            assert np.isfinite(row[key]), key
        assert row["em_predictor_stage_current_absmax"] > 0.0
        assert row["em_predictor_stage_lorentz_momentum_source_absmax"] > 0.0
        assert row["em_accepted_charge_change_absmax"] > 0.0
        assert row["em_stress_energy_density_absmax"] > 0.0

    checkpoints = sorted(run_dir.glob("checkpoint_t*.npz"))
    assert len(checkpoints) == 3
    snapshots = []
    for path in checkpoints:
        with np.load(path, allow_pickle=False) as saved:
            snap = {name: np.array(saved[name], copy=True) for name in saved.files}
        snapshots.append(snap)

    expected_times = (0.01, 0.02, 0.03)
    required_arrays = (
        "a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda",
        "S", "PS", "D", "PD", "phi", "Pi",
        "dm_rest", "dm_energy_t", "dm_momentum_r",
        "baryon_rest", "baryon_energy_t", "baryon_momentum_r",
        "radiation_rest", "radiation_energy_t", "radiation_momentum_r",
        "charged_dust_rest", "charged_dust_energy_t", "charged_dust_momentum_r",
        "em_enclosed_charge_faces", "em_electric_field_hat",
        "em_stress_rho", "em_stress_pr", "em_stress_pt", "em_stress_j",
        "em_gauss_shell_residual", "em_charge_density", "em_current_r",
        "em_lorentz_energy_source", "em_lorentz_momentum_source_r",
    )
    for expected_t, snap in zip(expected_times, snapshots):
        assert float(snap["t"]) == pytest.approx(expected_t, abs=1.0e-12)
        for name in required_arrays:
            assert name in snap, name
            assert np.all(np.isfinite(snap[name])), name
        assert float(snap["em_charge_to_rest_mass"]) == pytest.approx(0.25)
        assert snap["em_enclosed_charge_faces"][0] == 0.0
        assert np.any(np.abs(snap["em_electric_field_hat"]) > 0.0)
        assert np.any(np.abs(snap["em_current_r"]) > 0.0)

        # Independently recompute Gauss's shell invariant only from the saved
        # checkpoint variables; a finite stored residual alone is insufficient.
        q_faces = snap["em_enclosed_charge_faces"]
        charge_to_mass = float(snap["em_charge_to_rest_mass"])
        dr = float(snap["r"][1] - snap["r"][0])
        recomputed = np.diff(q_faces) - (
            4.0 * np.pi * charge_to_mass * dr * snap["charged_dust_rest"]
        )
        np.testing.assert_array_equal(snap["em_gauss_shell_residual"], recomputed)
        shell_charge = (
            4.0 * np.pi * charge_to_mass * dr * snap["charged_dust_rest"]
        )
        # The shell difference subtracts adjacent cumulative Q values, so the
        # roundoff scale must include those operands, not a unit-sized floor.
        # A conservative 64-epsilon bound covers the cumulative-charge,
        # stage-update, shell-difference and product roundoff without allowing
        # a small normalized charge state an unrelated absolute error budget.
        operand_scale = np.maximum.reduce((
            np.abs(q_faces[:-1]),
            np.abs(q_faces[1:]),
            np.abs(np.diff(q_faces)),
            np.abs(shell_charge),
        ))
        roundoff_bound = (
            64.0 * np.finfo(np.float64).eps
            * np.maximum(operand_scale, np.finfo(np.float64).tiny)
        )
        assert np.all(np.abs(recomputed) <= roundoff_bound)

    # This is specifically a repeated interaction test: the field charge and
    # carrier momentum must not be static copies across all accepted checkpoints.
    for previous, current in zip(snapshots, snapshots[1:]):
        assert np.any(
            previous["em_enclosed_charge_faces"]
            != current["em_enclosed_charge_faces"]
        )
        assert np.any(
            previous["charged_dust_momentum_r"]
            != current["charged_dust_momentum_r"]
        )
