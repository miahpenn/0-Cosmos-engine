"""Parameter-free production->COSMOS geometric handoff diagnostic.

This does not alter the production equations. It evolves the canonical V5.5
production kernel, then continuously drives the archive homogeneous COSMOS
equations with the production spacetime's solved H_eff. The result is a
geometric handoff/consistency test, not a phenomenological energy source.
"""
import json
from pathlib import Path

from .cosmos import (
    CosmosParams,
    construct_present_day,
    geometry_driven_step,
    state_rho_p,
)
from .diagnostics import friedmann_constraint
from .production_kernel import V55ProductionKernel


def run(
    resolution: int = 80,
    r_max: float = 80.0,
    final_time: float = 24.0,
    cfl: float = 0.03,
):
    if resolution <= 0 or r_max <= 0.0 or final_time <= 0.0 or cfl <= 0.0:
        raise ValueError("invalid production geometry bridge configuration")

    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=r_max,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    dt_nominal = cfl * state.grid.dr
    cosmos_params = CosmosParams()
    cosmos = construct_present_day()
    initial_H = float(kernel.diagnostics(state)["H_eff"])
    previous_H = initial_H

    while state.t < final_time:
        dt = min(dt_nominal, final_time - state.t)
        state = kernel.step(state, dt)
        row = state.history[-1]
        current_H = float(row["H_eff"])
        cosmos = geometry_driven_step(
            cosmos, cosmos_params, previous_H, current_H, dt
        )
        previous_H = current_H

    rho_cosmos, p_cosmos = state_rho_p(cosmos, cosmos_params)
    local_H = previous_H
    friedmann_residual = friedmann_constraint(local_H, rho_cosmos)
    return {
        "status": "completed",
        "resolution": resolution,
        "r_max": r_max,
        "final_time": float(state.t),
        "steps": len(state.history),
        "tau_final": float(state.tau),
        "H_eff_initial": initial_H,
        "H_eff_final": local_H,
        "cosmos_a_final": float(cosmos.a),
        "cosmos_phi_final": float(cosmos.phi),
        "cosmos_rho_final": float(rho_cosmos),
        "cosmos_pressure_final": float(p_cosmos),
        "geometry_to_cosmos_friedmann_residual": float(friedmann_residual),
        "cycle_events": len(state.cycle.events),
        "turnaround_events": sum(
            e.kind == "turnaround" for e in state.cycle.events
        ),
        "reexpansion_events": sum(
            e.kind == "re-expansion_crossing" for e in state.cycle.events
        ),
        "handoffs": len(state.handoffs),
        "no_phenomenological_Q": True,
    }


def main(output="runs/production-geometry-bridge/report.json"):
    result = run()
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
