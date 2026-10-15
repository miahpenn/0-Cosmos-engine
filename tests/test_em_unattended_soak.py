"""Unattended multi-step cross-engine soak for opt-in radial Maxwell.

This extends the existing three-step interaction regression to t=0.30 on two
resolutions. Inputs are fixed and explicit; there are no manual interventions,
parameter fits, production-source edits, or inferred GEAR carrier charges.
The test is a bounded research diagnostic, not scientific admission.
"""
import json
import math

import numpy as np

from engine.campaign import CampaignConfig
from engine.production_kernel import V55ProductionKernel
from engine.run_production import run_campaign


RESOLUTIONS = (24, 48)
R_MAX = 24.0
CFL = 0.01
FINAL_TIME = 0.30
CHECKPOINT_INTERVAL = 0.05
CHARGE_TO_MASS = 0.25
# A separate, explicit test-carrier prescription. Never inferred from GEAR.
DENSITY_PROFILE = ((0.0, 2.0e-5), (R_MAX, 2.0e-5))


def _max_abs(*arrays):
    return max(
        float(np.max(np.abs(np.asarray(values, dtype=float))))
        for values in arrays
    )


def _finite_or_none(value):
    if value is None:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _required_checkpoint_arrays():
    return (
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


def test_opt_in_em_cross_engine_soak_runs_unattended_to_t_030(
    tmp_path, monkeypatch, capsys
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
            "scalar": _max_abs(
                *(getattr(scalar_rhs, name) for name in scalar_names)
            ),
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
        terms = geometry_original(grid, geometry, scalars, matter, electromagnetic)
        if electromagnetic is not None:
            # Compare the real stage call with the same state absent the EM
            # projection; this observes actual geometry consumers, not a proxy.
            baseline = geometry_original(grid, geometry, scalars, matter, None)
            projection = electromagnetic.stress_projection(grid, geometry)
            geometry_exchange.append({
                "rho_max": _max_abs(projection["rho"]),
                "Aa_l3_delta": _max_abs(
                    terms["primary_l3"]["Aa"] - baseline["primary_l3"]["Aa"]
                ),
                "K_l3_delta": _max_abs(
                    terms["primary_l3"]["K"] - baseline["primary_l3"]["K"]
                ),
            })
        return terms

    monkeypatch.setattr(kernel, "_rhs", record_rhs)
    monkeypatch.setattr(kernel, "_geometry_stage_terms", record_geometry)

    config = CampaignConfig(
        resolutions=RESOLUTIONS,
        r_max=R_MAX,
        cfl=CFL,
        final_time=FINAL_TIME,
        checkpoint_interval=CHECKPOINT_INTERVAL,
        output_dir=str(tmp_path),
        radial_maxwell_charge_to_mass=CHARGE_TO_MASS,
        charged_dust_density_profile=DENSITY_PROFILE,
    )
    results = run_campaign(kernel=kernel, config=config)

    assert len(results) == len(RESOLUTIONS)
    per_resolution = {}
    all_checkpoint_times = {}
    required_arrays = _required_checkpoint_arrays()
    constraint_keys = (
        "hamiltonian_full_absmax",
        "momentum_full_absmax",
        "connection_full_absmax",
        "determinant_constraint_max",
        "M_MS",
        "lapse_min",
    )

    for resolution, result in zip(RESOLUTIONS, results):
        assert result.resolution == resolution
        assert result.status == "completed", (resolution, result.status)
        assert math.isclose(
            result.t_final, FINAL_TIME, rel_tol=0.0, abs_tol=1.0e-11
        )
        nominal_steps = math.ceil(
            FINAL_TIME / (CFL * R_MAX / resolution) - 1.0e-10
        )
        assert result.steps >= nominal_steps
        assert result.checkpoint_count >= 5

        run_dir = tmp_path / f"N{resolution}"
        ledger = json.loads((run_dir / "ledgers.json").read_text())
        history = ledger["history"]
        assert len(history) >= result.steps + 1
        accepted_em_rows = [
            row for row in history
            if "em_predictor_stage_current_absmax" in row
        ]
        assert len(accepted_em_rows) >= result.steps

        for row in accepted_em_rows:
            for key in (
                "em_predictor_stage_current_absmax",
                "em_start_stage_lorentz_momentum_source_absmax",
                "em_predictor_stage_lorentz_momentum_source_absmax",
                "em_accepted_gauss_shell_residual_absmax",
                "em_stress_energy_density_absmax",
            ):
                assert key in row
                assert row[key] is not None and math.isfinite(float(row[key])), (
                    resolution, key, row.get(key)
                )
            assert row["em_predictor_stage_current_absmax"] > 0.0
            assert row["em_predictor_stage_lorentz_momentum_source_absmax"] > 0.0
            assert row["em_stress_energy_density_absmax"] > 0.0

        checkpoint_paths = sorted(run_dir.glob("checkpoint_t*.npz"))
        assert len(checkpoint_paths) == result.checkpoint_count
        assert len(checkpoint_paths) >= 5
        snapshots = []
        gauss_absmax_values = []
        for path in checkpoint_paths:
            with np.load(path, allow_pickle=False) as saved:
                snap = {name: np.array(saved[name], copy=True) for name in saved.files}
            snapshots.append(snap)
            assert float(snap["em_charge_to_rest_mass"]) == CHARGE_TO_MASS
            assert snap["em_enclosed_charge_faces"][0] == 0.0
            for name in required_arrays:
                assert name in snap, (resolution, name)
                assert np.all(np.isfinite(snap[name])), (resolution, name)
            assert np.any(np.abs(snap["em_electric_field_hat"]) > 0.0)
            assert np.any(np.abs(snap["em_current_r"]) > 0.0)
            assert np.any(np.abs(snap["em_stress_rho"]) > 0.0)
            assert np.any(np.abs(snap["em_lorentz_momentum_source_r"]) > 0.0)

            # Recompute the Gauss shell identity from the saved Q faces,
            # explicit q/m, dr and carrier rest state; do not trust a stored norm.
            q_faces = snap["em_enclosed_charge_faces"]
            dr = float(snap["r"][1] - snap["r"][0])
            shell_charge = (
                4.0 * math.pi * CHARGE_TO_MASS * dr * snap["charged_dust_rest"]
            )
            recomputed = np.diff(q_faces) - shell_charge
            np.testing.assert_array_equal(
                snap["em_gauss_shell_residual"], recomputed
            )
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
            gauss_absmax_values.append(float(np.max(np.abs(recomputed))))

        times = [float(snap["t"]) for snap in snapshots]
        assert times == sorted(times)
        assert all(b > a for a, b in zip(times, times[1:]))
        assert math.isclose(times[-1], FINAL_TIME, rel_tol=0.0, abs_tol=1.0e-11)
        all_checkpoint_times[resolution] = times

        # Require evolving charge faces and carrier momentum at every saved
        # interval, so passing cannot be achieved by replaying static snapshots.
        for previous, current in zip(snapshots, snapshots[1:]):
            assert np.any(
                previous["em_enclosed_charge_faces"]
                != current["em_enclosed_charge_faces"]
            )
            assert np.any(
                previous["charged_dust_momentum_r"]
                != current["charged_dust_momentum_r"]
            )

        final_row = history[-1]
        per_resolution[resolution] = {
            "status": result.status,
            "accepted_steps": int(result.steps),
            "checkpoint_count": int(result.checkpoint_count),
            "t_final": float(result.t_final),
            "tau_final": _finite_or_none(result.tau_final),
            "max_recomputed_gauss_shell_residual": max(gauss_absmax_values),
            "final_constraints_and_mass": {
                key: _finite_or_none(final_row.get(key)) for key in constraint_keys
            },
        }

    total_accepted_steps = sum(result.steps for result in results)
    # Every accepted PIRK step evaluates the real RHS at start/predictor and
    # contributes three observed geometry-stage consumers. Extra calls are
    # allowed if the existing automatic admissibility retry path is exercised.
    assert len(rhs_trace) >= 2 * total_accepted_steps
    assert len(geometry_exchange) >= 3 * total_accepted_steps
    for sector in ("scalar", "dark_matter", "baryons", "radiation", "charged_dust"):
        values = [row[sector] for row in rhs_trace]
        assert all(math.isfinite(value) for value in values), sector
        assert max(values) > 0.0, sector
    for row in geometry_exchange:
        assert all(math.isfinite(value) for value in row.values())
        assert row["rho_max"] > 0.0
        assert row["Aa_l3_delta"] > 0.0
        assert row["K_l3_delta"] > 0.0

    summary = json.loads((tmp_path / "campaign_summary.json").read_text())
    extension = summary["physics_extensions"]["radial_maxwell"]
    assert extension["enabled"] is True
    assert extension["status"] == "research_only_not_admitted"

    report = {
        "classification": "bounded unattended research soak; not scientific admission",
        "engine_source_baseline": "44828c98c3752f3a44cac0167bfce667df0d83a4",
        "settings": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "cfl": CFL,
            "final_time": FINAL_TIME,
            "checkpoint_interval": CHECKPOINT_INTERVAL,
            "charge_to_mass": CHARGE_TO_MASS,
            "charged_dust_density_profile": [list(knot) for knot in DENSITY_PROFILE],
            "human_interventions_after_start": 0,
            "automatic_retry_path": "existing production behavior; no manual tuning",
        },
        "results": {str(key): value for key, value in per_resolution.items()},
        "checkpoint_times": {
            str(key): value for key, value in all_checkpoint_times.items()
        },
        "trace_counts": {
            "accepted_steps": total_accepted_steps,
            "scalar_and_matter_rhs_evaluations": len(rhs_trace),
            "geometry_stage_em_feedback_observations": len(geometry_exchange),
        },
        "limits": [
            "explicit test carrier is not a charge assignment to a GEAR species",
            "not the six-hour production campaign",
            "does not establish full field+matter+gravity conservation",
            "does not establish long-horizon spatial/temporal convergence or scientific admission",
        ],
    }
    payload = json.dumps(report, sort_keys=True, allow_nan=False)
    with capsys.disabled():
        print("EM_UNATTENDED_SOAK_REPORT_JSON=" + payload, flush=True)
