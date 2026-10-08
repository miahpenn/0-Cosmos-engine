"""Common CFL=0.0075 four-way center-repair clock-profile diagnostic.

The four cases use this same runner/schema:
  D=0 or 1e-4, repair OFF or ON
  N=160, r_max=80, CFL=0.0075, final t=22.5
  radiation ON
  target times [0.01, 4, 8, 12, 16, 20, 22.5]

Diagnostic-only: production equations are unchanged.
Select the case with environment variables CASE_D, CASE_REPAIR, CASE_NAME.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from engine.production_kernel import V55ProductionKernel


TARGET_TIMES = [0.01, 4.0, 8.0, 12.0, 16.0, 20.0, 22.5]
N = 160
R_MAX = 80.0
CFL = 0.0075
FINAL_TIME = 22.5


class NoRepairKernel(V55ProductionKernel):
    @staticmethod
    def _enforce_center_regularity(grid, geometry):
        # Diagnostic A/B only: bypass the existing algebraic/numerical repair.
        return geometry


def run_case(D: float, repair: bool) -> dict:
    Kernel = V55ProductionKernel if repair else NoRepairKernel
    k = Kernel()
    state = k.initialize(
        resolution=N,
        r_max=R_MAX,
        D_amplitude=D,
        include_radiation=True,
    )

    dt_nom = CFL * state.grid.dr
    r = np.asarray(state.grid.centers, dtype=float).copy()
    tau_profile = np.zeros_like(r)
    alpha_prev = np.asarray(state.geometry.alpha, dtype=float).copy()
    samples = []
    target_idx = 0
    failed = None

    def append_sample() -> None:
        last = state.history[-1] if state.history else {}
        alpha_now = np.asarray(state.geometry.alpha, dtype=float).copy()
        samples.append(
            {
                "t": float(state.t),
                "tau_central": float(state.tau),
                "r": r.tolist(),
                "tau_r": tau_profile.tolist(),
                "alpha_r": alpha_now.tolist(),
                "hamiltonian_max": float(last.get("hamiltonian_max", np.nan)),
                "hamiltonian_max_r": float(last.get("hamiltonian_max_r", np.nan)),
                "hamiltonian_center": float(last.get("hamiltonian_at_max", np.nan)),
                "lapse_min": float(last.get("lapse_min", np.nan)),
                "lapse_min_r": float(last.get("lapse_min_r", np.nan)),
                "momentum_max": float(last.get("momentum_max", np.nan)),
                "cmc_residual_outer_max": float(last.get("cmc_residual_outer_max", np.nan)),
            }
        )

    while state.t < FINAL_TIME - 1e-14:
        dt = min(dt_nom, FINAL_TIME - state.t)
        try:
            state = k.step(state, dt)
            state.geometry.assert_finite_positive()
        except (FloatingPointError, ValueError) as exc:
            failed = {"t": float(state.t), "error": str(exc)}
            break

        alpha_now = np.asarray(state.geometry.alpha, dtype=float).copy()
        tau_profile += 0.5 * (alpha_prev + alpha_now) * dt
        alpha_prev = alpha_now

        while target_idx < len(TARGET_TIMES) and state.t + 0.5 * dt >= TARGET_TIMES[target_idx]:
            append_sample()
            target_idx += 1

    last = state.history[-1] if state.history else {}
    return {
        "status": "numerical_failure" if failed else "completed",
        "diagnostic_only": True,
        "purpose": "Common-schema four-way center-repair radial proper-time and lapse profile.",
        "parameters": {
            "D_amplitude": float(D),
            "repair": bool(repair),
            "resolution": N,
            "r_max": R_MAX,
            "cfl": CFL,
            "final_time_requested": FINAL_TIME,
            "include_radiation": True,
            "target_times": TARGET_TIMES,
        },
        "final_time": float(state.t),
        "tau_final": float(state.tau),
        "r": r.tolist(),
        "samples": samples,
        "failure": failed,
        "final_diagnostics": {
            "hamiltonian_max": float(last.get("hamiltonian_max", np.nan)),
            "hamiltonian_max_r": float(last.get("hamiltonian_max_r", np.nan)),
            "hamiltonian_at_max": float(last.get("hamiltonian_at_max", np.nan)),
            "lapse_min": float(last.get("lapse_min", np.nan)),
            "lapse_min_r": float(last.get("lapse_min_r", np.nan)),
            "momentum_max": float(last.get("momentum_max", np.nan)),
            "cmc_residual_outer_max": float(last.get("cmc_residual_outer_max", np.nan)),
        },
    }


def main() -> int:
    try:
        D = float(os.environ.get("CASE_D", "0"))
        repair = os.environ.get("CASE_REPAIR", "0") == "1"
        case_name = os.environ.get("CASE_NAME", "unnamed")
    except ValueError as exc:
        raise SystemExit(f"Invalid case configuration: {exc}")

    result = run_case(D, repair)
    out = Path(f"GEAR_group/four_way_{case_name}_N160.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
