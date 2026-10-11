"""Full production/COSMOS bidirectional consistency audit.

Runs the canonical production spacetime and a synchronized homogeneous COSMOS
mirror driven only by the production spacetime's solved H_eff.

The production kernel already contains the reciprocal local geometry <-> COSMOS
field/matter interaction through the solved stress-energy and K/lapse terms.
This diagnostic does not add a new physical source, efficiency, closure, or
time-scale parameter. The homogeneous lane is a mirror/ledger channel only.
"""
from __future__ import annotations

import json
import math

import numpy as np
from pathlib import Path

from .cosmos import CosmosParams, construct_present_day, geometry_driven_step, state_rho_p
from .production_kernel import V55ProductionKernel
from .stress_energy import assemble_total_stress_energy


def _signed_friedmann_H(rho: float, H_ref: float) -> float:
    return math.copysign(math.sqrt(max(rho / 3.0, 0.0)), H_ref if H_ref != 0.0 else 1.0)


def run_case(
    D_amplitude: float,
    *,
    resolution: int = 80,
    r_max: float = 80.0,
    final_time: float = 24.0,
    cfl: float = 0.03,
    include_radiation: bool = True,
):
    if D_amplitude < 0.0:
        raise ValueError("D_amplitude must be non-negative")

    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        D_amplitude=D_amplitude,
        include_radiation=include_radiation,
    )

    params = CosmosParams()
    cosmos = construct_present_day(params)
    dt_nominal = cfl * state.grid.dr
    previous_H = float(kernel.diagnostics(state)["H_eff"])

    rows = []
    for _ in range(math.ceil(final_time / dt_nominal) + 1):
        if state.t >= final_time - 1.0e-15:
            break
        dt = min(dt_nominal, final_time - state.t)
        state = kernel.step(state, dt)
        obs = state.history[-1]
        H_eff = float(obs["H_eff"])

        cosmos = geometry_driven_step(
            cosmos, params, previous_H, H_eff, dt
        )
        rho_c, p_c = state_rho_p(cosmos, params)
        H_f = _signed_friedmann_H(rho_c, H_eff)
        total = assemble_total_stress_energy(
            state.grid, state.geometry, state.scalars, state.matter
        )
        w = state.grid.volumes
        rho_prod_volume = 8.0 * math.pi * float(np.sum(w * total.rho) / np.sum(w))
        rho_prod_outer_physical = 8.0 * math.pi * float(obs["rho_outer"])
        H_f_prod_volume = _signed_friedmann_H(rho_prod_volume, H_eff)
        a_from_geometry = math.exp(float(state.e_folds)) * 33.8983

        phi_gap = float(obs["phi_outer"]) - float(cosmos.phi)
        pi_gap = float(obs["Pi_outer"]) - float(cosmos.pi_phi)
        rho_gap = float(obs["rho_outer"]) - float(rho_c)

        rows.append({
            "hamiltonian_max": float(obs["hamiltonian_max"]),
            "hamiltonian_normalized_max": float(obs["hamiltonian_normalized_max"]),
            "hamiltonian_l2_inner": float(obs["hamiltonian_l2_inner"]),
            "hamiltonian_l2_outer": float(obs["hamiltonian_l2_outer"]),
            "momentum_max": float(obs["momentum_max"]),
            "momentum_normalized_max": float(obs["momentum_normalized_max"]),
            "momentum_l2_inner": float(obs["momentum_l2_inner"]),
            "momentum_l2_outer": float(obs["momentum_l2_outer"]),
            "cmc_residual_outer_max": float(obs["cmc_residual_outer_max"]),
            "cmc_residual_last_interior": float(obs["cmc_residual_last_interior"]),
            "determinant_min": float(obs["determinant_min"]),
            "determinant_constraint_max": float(obs["determinant_constraint_max"]),
            "trapping_min": float(obs["trapping_min"]),
            "trapped_roots_count": len(obs["trapped_roots"]),
            "alpha_r_max": float(obs["alpha_r_max"]),
            "alpha_sigma": float(obs["alpha_sigma"]),
            "t": float(state.t),
            "tau": float(state.tau),
            "H_eff": H_eff,
            "H_friedmann_cosmos": H_f,
            "friedmann_residual": float(3.0 * H_eff * H_eff - rho_c),
            "production_volume_rho_physical": rho_prod_volume,
            "production_volume_friedmann_residual": float(3.0 * H_eff * H_eff - rho_prod_volume),
            "production_volume_H_friedmann": H_f_prod_volume,
            "production_a_from_e_folds": a_from_geometry,
            "cosmos_a": float(cosmos.a),
            "a_log_gap": float(math.log(max(cosmos.a, 1.0e-300) / max(a_from_geometry, 1.0e-300))),
            "production_phi_outer": float(obs["phi_outer"]),
            "cosmos_phi": float(cosmos.phi),
            "phi_gap": phi_gap,
            "production_Pi_outer": float(obs["Pi_outer"]),
            "cosmos_pi_phi": float(cosmos.pi_phi),
            "pi_gap": pi_gap,
            "production_rho_outer": float(obs["rho_outer"]),
            "production_rho_outer_physical": rho_prod_outer_physical,
            "cosmos_rho": float(rho_c),
            "rho_gap": rho_gap,
            "rho_gap_rel": abs(rho_gap) / max(abs(rho_c), 1.0e-300),
            "rho_outer_physical_gap": float(rho_prod_outer_physical - rho_c),
            "rho_outer_physical_gap_rel": abs(float(rho_prod_outer_physical - rho_c)) / max(abs(rho_c), 1.0e-300),
            "lapse_min": float(obs["lapse_min"]),
            "lapse_min_r": float(obs["lapse_min_r"]),
            "D_center": float(state.scalars.D[0]),
            "cycle_events": len(state.cycle.events),
        })
        previous_H = H_eff

    final = rows[-1]
    avg_abs_rho_gap_rel = sum(r["rho_gap_rel"] for r in rows) / len(rows)
    avg_abs_phi_gap = sum(abs(r["phi_gap"]) for r in rows) / len(rows)
    avg_abs_pi_gap = sum(abs(r["pi_gap"]) for r in rows) / len(rows)
    max_abs_friedmann = max(abs(r["friedmann_residual"]) for r in rows)
    exact_checks = {}
    for key in ("hamiltonian_max", "hamiltonian_normalized_max", "momentum_max", "momentum_normalized_max", "cmc_residual_outer_max", "cmc_residual_last_interior", "determinant_constraint_max"):
        exact_checks[f"max_abs_{key}"] = float(max(abs(r[key]) for r in rows))
    max_abs_production_volume_friedmann = max(abs(r["production_volume_friedmann_residual"]) for r in rows)
    avg_abs_rho_outer_physical_gap_rel = sum(r["rho_outer_physical_gap_rel"] for r in rows) / len(rows)

    return {
        "D_amplitude": D_amplitude,
        "resolution": resolution,
        "r_max": r_max,
        "final_time": float(state.t),
        "steps": len(rows),
        "tau_final": float(state.tau),
        "H_eff_initial": float(rows[0]["H_eff"]),
        "H_eff_final": float(final["H_eff"]),
        "lapse_min_final": float(final["lapse_min"]),
        "lapse_min_r_final": float(final["lapse_min_r"]),
        "cosmos_a_final": float(cosmos.a),
        "production_a_from_e_folds_final": float(final["production_a_from_e_folds"]),
        "a_log_gap_final": float(final["a_log_gap"]),
        "cosmos_phi_final": float(cosmos.phi),
        "production_phi_outer_final": float(final["production_phi_outer"]),
        "phi_gap_final": float(final["phi_gap"]),
        "cosmos_pi_final": float(cosmos.pi_phi),
        "production_Pi_outer_final": float(final["production_Pi_outer"]),
        "pi_gap_final": float(final["pi_gap"]),
        "cosmos_rho_final": float(final["cosmos_rho"]),
        "production_rho_outer_final": float(final["production_rho_outer"]),
        "rho_gap_final": float(final["rho_gap"]),
        "avg_abs_rho_gap_rel": float(avg_abs_rho_gap_rel),
        "avg_abs_phi_gap": float(avg_abs_phi_gap),
        "avg_abs_pi_gap": float(avg_abs_pi_gap),
        "max_abs_friedmann_residual": float(max_abs_friedmann),
        "exact_constraint_witnesses": exact_checks,
        "final_exact_constraints": {k: final[k] for k in ("hamiltonian_max", "hamiltonian_normalized_max", "hamiltonian_l2_inner", "hamiltonian_l2_outer", "momentum_max", "momentum_normalized_max", "momentum_l2_inner", "momentum_l2_outer", "cmc_residual_outer_max", "cmc_residual_last_interior", "determinant_min", "determinant_constraint_max", "trapping_min", "trapped_roots_count", "alpha_r_max", "alpha_sigma")},
        "max_abs_production_volume_friedmann_residual": float(max_abs_production_volume_friedmann),
        "avg_abs_rho_outer_physical_gap_rel": float(avg_abs_rho_outer_physical_gap_rel),
        "cycle_events": len(state.cycle.events),
        "turnaround_events": sum(e.kind == "turnaround" for e in state.cycle.events),
        "reexpansion_events": sum(e.kind == "re-expansion_crossing" for e in state.cycle.events),
        "handoffs": len(state.handoffs),
        "no_phenomenological_Q": True,
    }, rows


def run(
    *,
    resolution: int = 80,
    r_max: float = 80.0,
    final_time: float = 24.0,
    cfl: float = 0.03,
):
    strong, strong_rows = run_case(
        1.0e-4,
        resolution=resolution,
        r_max=r_max,
        final_time=final_time,
        cfl=cfl,
        include_radiation=True,
    )
    control, control_rows = run_case(
        0.0,
        resolution=resolution,
        r_max=r_max,
        final_time=final_time,
        cfl=cfl,
        include_radiation=True,
    )

    return {
        "status": "completed",
        "architecture": {
            "production_cosmos_interaction": "shared production spacetime stress-energy and K/lapse coupling",
            "homogeneous_channel": "derived H_eff mirror only",
            "reciprocal_physical_source_added": False,
            "phenomenological_Q": False,
        },
        "strong_D": strong,
        "D0_control": control,
        "strong_minus_control_tau": float(strong["tau_final"] - control["tau_final"]),
        "strong_tau_ratio": float(strong["tau_final"] / control["tau_final"]),
        "strong_minus_control_H_final": float(
            strong["H_eff_final"] - control["H_eff_final"]
        ),
        "production_cosmos_response": {
            "phi_outer_change_strong_minus_control": float(
                strong["production_phi_outer_final"] - control["production_phi_outer_final"]
            ),
            "Pi_outer_change_strong_minus_control": float(
                strong["production_Pi_outer_final"] - control["production_Pi_outer_final"]
            ),
            "homogeneous_phi_change_strong_minus_control": float(
                strong["cosmos_phi_final"] - control["cosmos_phi_final"]
            ),
            "homogeneous_Pi_change_strong_minus_control": float(
                strong["cosmos_pi_final"] - control["cosmos_pi_final"]
            ),
            "extra_local_phi_shift_beyond_H_mirror": float(
                (strong["production_phi_outer_final"] - control["production_phi_outer_final"])
                - (strong["cosmos_phi_final"] - control["cosmos_phi_final"])
            ),
            "extra_local_Pi_shift_beyond_H_mirror": float(
                (strong["production_Pi_outer_final"] - control["production_Pi_outer_final"])
                - (strong["cosmos_pi_final"] - control["cosmos_pi_final"])
            ),
        },
        "sample_rows": {
            "strong_D_first": strong_rows[0],
            "strong_D_mid": strong_rows[len(strong_rows) // 2],
            "strong_D_last": strong_rows[-1],
        },
    }


def main(output="runs/production-cosmos-bidirectional/report.json"):
    result = run()
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
