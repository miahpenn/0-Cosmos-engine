"""Quarantined strong-D center-repair-OFF spatial clock-profile diagnostic.

Diagnostic-only. This runs exactly one case:
  D=1e-4, repair OFF, N=160, r_max=80, CFL=0.0075, final t=22.5.

The production equations are not modified. The output path is intentionally
fixed to GEAR_group/strong_D1e-4_repair_off_N160.json so the artifact can be
validated independently before any physics interpretation.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from engine.center_repair_ab_clock_profile import NoRepairKernel


OUTPUT = Path("GEAR_group/strong_D1e-4_repair_off_N160.json")
TARGET_TIMES = [0.01, 4.0, 8.0, 12.0, 16.0, 20.0, 22.5]


def run_case() -> dict:
    D = 1e-4
    N = 160
    r_max = 80.0
    cfl = 0.0075
    final_time = 22.5

    k = NoRepairKernel()
    state = k.initialize(
        resolution=N,
        r_max=r_max,
        D_amplitude=D,
        include_radiation=True,
    )

    dt_nom = cfl * state.grid.dr
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

    while state.t < final_time - 1e-14:
        dt = min(dt_nom, final_time - state.t)
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
    result = {
        "status": "numerical_failure" if failed else "completed",
        "diagnostic_only": True,
        "purpose": "Strong-D center-repair-OFF radial proper-time and lapse profile.",
        "parameters": {
            "D_amplitude": D,
            "repair": False,
            "resolution": N,
            "r_max": r_max,
            "cfl": cfl,
            "final_time_requested": final_time,
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
    return result


def main() -> int:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result = run_case()
    OUTPUT.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "completed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
