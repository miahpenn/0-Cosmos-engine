"""Quantitative, short-horizon radial-Maxwell behavior probe.

Paired control: the same separately supplied charged-dust profile at q/m=0
and q/m=0.25, at N=24/48/96 through t=0.03. This reports the observed
semi-discrete behavior; it is not a physics admission or convergence proof.
"""
import json
import math
import os

import numpy as np

from engine.campaign import CampaignConfig
from engine.run_production import run_campaign


RESOLUTIONS = (24, 48, 96)
R_MAX = 24.0
FINAL_TIME = 0.03
CFL = 0.01
DENSITY_PROFILE = ((0.0, 2.0e-5), (R_MAX, 2.0e-5))
CHARGE_TO_MASS_ARMS = (0.0, 0.25)


def _finite_or_none(value):
    number = float(value)
    return number if math.isfinite(number) else None


def _run_arm(root, eta):
    output_dir = root / f"q_over_m_{eta:.2f}"
    config = CampaignConfig(
        resolutions=RESOLUTIONS,
        r_max=R_MAX,
        cfl=CFL,
        final_time=FINAL_TIME,
        checkpoint_interval=0.01,
        output_dir=str(output_dir),
        radial_maxwell_charge_to_mass=eta,
        charged_dust_density_profile=DENSITY_PROFILE,
    )
    results = run_campaign(config=config)
    assert len(results) == len(RESOLUTIONS)
    rows = []

    for resolution, result in zip(RESOLUTIONS, results):
        assert result.status == "completed", (eta, resolution, result.status)
        assert math.isclose(result.t_final, FINAL_TIME, rel_tol=0.0, abs_tol=1.0e-11)

        run_dir = output_dir / f"N{resolution}"
        ledger = json.loads((run_dir / "ledgers.json").read_text())
        history = ledger["history"]
        assert len(history) >= 2
        initial, final = history[0], history[-1]

        checkpoint_paths = sorted(run_dir.glob("checkpoint_t*.npz"))
        assert checkpoint_paths
        snapshots = []
        for path in checkpoint_paths:
            with np.load(path, allow_pickle=False) as saved:
                snapshots.append({name: np.array(saved[name], copy=True) for name in saved.files})
        final_saved = max(snapshots, key=lambda state: float(state["t"]))

        q_faces = final_saved["em_enclosed_charge_faces"]
        saved_eta = float(final_saved["em_charge_to_rest_mass"])
        dr = float(final_saved["r"][1] - final_saved["r"][0])
        recomputed_gauss = np.diff(q_faces) - (
            4.0 * np.pi * saved_eta * dr * final_saved["charged_dust_rest"]
        )
        gauss_absmax = float(np.max(np.abs(recomputed_gauss)))
        shell_charge = (
            4.0 * np.pi * saved_eta * dr * final_saved["charged_dust_rest"]
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
        assert np.all(np.abs(recomputed_gauss) <= roundoff_bound)

        row = {
            "q_over_m": eta,
            "resolution": resolution,
            "steps": int(result.steps),
            "t_final": _finite_or_none(result.t_final),
            "tau_final": _finite_or_none(result.tau_final),
            "outer_enclosed_charge": _finite_or_none(q_faces[-1]),
            "electric_field_hat_absmax": float(np.max(np.abs(final_saved["em_electric_field_hat"]))),
            "em_energy_density_absmax": float(np.max(np.abs(final_saved["em_stress_rho"]))),
            "charge_current_absmax": float(np.max(np.abs(final_saved["em_current_r"]))),
            "lorentz_energy_source_absmax": float(np.max(np.abs(final_saved["em_lorentz_energy_source"]))),
            "lorentz_momentum_source_absmax": float(np.max(np.abs(final_saved["em_lorentz_momentum_source_r"]))),
            "gauss_shell_residual_absmax_recomputed": gauss_absmax,
            "gauss_shell_residual_absmax_reported": _finite_or_none(final.get("em_accepted_gauss_shell_residual_absmax")),
            "hamiltonian_full_absmax_initial": _finite_or_none(initial.get("hamiltonian_full_absmax")),
            "hamiltonian_full_absmax_final": _finite_or_none(final.get("hamiltonian_full_absmax")),
            "momentum_full_absmax_initial": _finite_or_none(initial.get("momentum_full_absmax")),
            "momentum_full_absmax_final": _finite_or_none(final.get("momentum_full_absmax")),
            "connection_full_absmax_initial": _finite_or_none(initial.get("connection_full_absmax")),
            "connection_full_absmax_final": _finite_or_none(final.get("connection_full_absmax")),
            "determinant_constraint_max_initial": _finite_or_none(initial.get("determinant_constraint_max")),
            "determinant_constraint_max_final": _finite_or_none(final.get("determinant_constraint_max")),
            "misner_sharp_mass_final": _finite_or_none(final.get("M_MS")),
            "lapse_min_final": _finite_or_none(final.get("lapse_min")),
            "checkpoint_times": [float(state["t"]) for state in snapshots],
        }
        rows.append(row)

        # The neutral carrier is a paired control for the same added dust mass;
        # its electric field, stress, current, and Lorentz sources should vanish.
        if eta == 0.0:
            for key in (
                "outer_enclosed_charge",
                "electric_field_hat_absmax",
                "em_energy_density_absmax",
                "charge_current_absmax",
                "lorentz_energy_source_absmax",
                "lorentz_momentum_source_absmax",
            ):
                assert abs(row[key]) <= 1.0e-30, (resolution, key, row[key])
        else:
            for key in (
                "outer_enclosed_charge",
                "electric_field_hat_absmax",
                "em_energy_density_absmax",
                "charge_current_absmax",
                "lorentz_momentum_source_absmax",
            ):
                assert row[key] > 0.0, (resolution, key, row[key])

    return rows


def test_report_short_horizon_opt_in_em_behavior(tmp_path, capsys):
    rows = []
    for eta in CHARGE_TO_MASS_ARMS:
        rows.extend(_run_arm(tmp_path, eta))

    neutral = {row["resolution"]: row for row in rows if row["q_over_m"] == 0.0}
    charged = {row["resolution"]: row for row in rows if row["q_over_m"] == 0.25}
    paired = []
    for resolution in RESOLUTIONS:
        control, active = neutral[resolution], charged[resolution]
        paired.append({
            "resolution": resolution,
            "outer_charge_delta": active["outer_enclosed_charge"] - control["outer_enclosed_charge"],
            "electric_field_absmax_delta": active["electric_field_hat_absmax"] - control["electric_field_hat_absmax"],
            "em_energy_density_absmax_delta": active["em_energy_density_absmax"] - control["em_energy_density_absmax"],
            "lorentz_momentum_source_absmax_delta": active["lorentz_momentum_source_absmax"] - control["lorentz_momentum_source_absmax"],
            "hamiltonian_full_absmax_delta": (
                None if active["hamiltonian_full_absmax_final"] is None or control["hamiltonian_full_absmax_final"] is None
                else active["hamiltonian_full_absmax_final"] - control["hamiltonian_full_absmax_final"]
            ),
            "momentum_full_absmax_delta": (
                None if active["momentum_full_absmax_final"] is None or control["momentum_full_absmax_final"] is None
                else active["momentum_full_absmax_final"] - control["momentum_full_absmax_final"]
            ),
        })

    report = {
        "classification": "bounded research-only short-horizon behavior probe; not scientific admission",
        "source_line": "main radial-Maxwell opt-in merged by PR #176",
        "settings": {
            "resolutions": list(RESOLUTIONS),
            "r_max": R_MAX,
            "final_time": FINAL_TIME,
            "cfl": CFL,
            "charged_dust_density_profile": [list(knot) for knot in DENSITY_PROFILE],
            "charge_to_mass_arms": list(CHARGE_TO_MASS_ARMS),
            "interpretation": "q/m=0 holds the extra dust mass fixed while removing the electric field; q/m=0.25 is the explicitly charged comparison",
        },
        "rows": rows,
        "charged_minus_neutral": paired,
        "open_gates": [
            "carrier/action derivation from GEAR",
            "physical field-to-engine normalization",
            "nonzero-shift coupling",
            "full donor+field+gravity conservation",
            "initial/evolved constraint convergence",
            "checkpoint restart equivalence",
            "longer-horizon and spatial/time convergence",
            "radiation-continuation gate",
        ],
    }
    payload = json.dumps(report, sort_keys=True, allow_nan=False)

    # Print outside pytest's output capture so the numbers survive in the
    # GitHub Actions log even when all assertions pass.
    with capsys.disabled():
        print("EM_SHORT_HORIZON_REPORT_JSON=" + payload, flush=True)

    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_path:
        def fmt(value):
            return "—" if value is None else f"{value:.6e}" if isinstance(value, float) else str(value)

        lines = [
            "## Radial Maxwell short-horizon behavior probe",
            "",
            "Research-only paired test, same separate dust profile, q/m=0 versus 0.25. Gaussian units, c=1. N=24/48/96, t=0.03. This is an early-time behavior readout, not a conservation/convergence or scientific-admission result.",
            "",
            "| q/m | N | steps | outer Q | max |Ê| | max EM ρ | max |j| | max Lorentz |S| source | max Gauss residual | max H | max M | max connection |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for row in rows:
            lines.append(
                f"| {row['q_over_m']:.2f} | {row['resolution']} | {row['steps']} | "
                f"{fmt(row['outer_enclosed_charge'])} | {fmt(row['electric_field_hat_absmax'])} | "
                f"{fmt(row['em_energy_density_absmax'])} | {fmt(row['charge_current_absmax'])} | "
                f"{fmt(row['lorentz_momentum_source_absmax'])} | "
                f"{fmt(row['gauss_shell_residual_absmax_recomputed'])} | "
                f"{fmt(row['hamiltonian_full_absmax_final'])} | "
                f"{fmt(row['momentum_full_absmax_final'])} | "
                f"{fmt(row['connection_full_absmax_final'])} |"
            )
        lines.extend(["", "### Charged minus neutral control", "", "| N | Δ outer Q | Δ max |Ê| | Δ max EM ρ | Δ Lorentz momentum source | Δ max H | Δ max M |", "|---:|---:|---:|---:|---:|---:|---:|"])
        for row in paired:
            lines.append(
                f"| {row['resolution']} | {fmt(row['outer_charge_delta'])} | "
                f"{fmt(row['electric_field_absmax_delta'])} | {fmt(row['em_energy_density_absmax_delta'])} | "
                f"{fmt(row['lorentz_momentum_source_absmax_delta'])} | "
                f"{fmt(row['hamiltonian_full_absmax_delta'])} | {fmt(row['momentum_full_absmax_delta'])} |"
            )
        with open(summary_path, "a", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")
