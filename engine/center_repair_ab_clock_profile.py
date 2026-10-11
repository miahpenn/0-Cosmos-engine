"""Fixed-domain center-repair A/B plus radial clock-profile diagnostic.

Diagnostic-only: the production equations are untouched. The only A/B switch is
whether the existing algebraic center-regularity projection is applied.
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from engine.production_kernel import V55ProductionKernel


class NoRepairKernel(V55ProductionKernel):
    @staticmethod
    def _enforce_center_regularity(grid, geometry):
        # Diagnostic A/B only: bypass the existing algebraic/numerical repair.
        return geometry


def run_case(D, repair, N=160, r_max=80.0, cfl=0.0075, final_time=24.0):
    Kernel = V55ProductionKernel if repair else NoRepairKernel
    k = Kernel()
    state = k.initialize(
        resolution=N,
        r_max=r_max,
        D_amplitude=D,
        include_radiation=True,
    )
    dt_nom = cfl * state.grid.dr
    target_times = [0.0, 6.0, 12.0, 18.0, 24.0]
    sample_idx = 1
    r = np.asarray(state.grid.centers, dtype=float).copy()
    tau_profile = np.zeros_like(r)
    alpha_prev = np.asarray(state.geometry.alpha, dtype=float).copy()

    samples = [{
        "t": 0.0,
        "tau_central": 0.0,
        "r": r.tolist(),
        "tau_r": tau_profile.tolist(),
        "alpha_r": alpha_prev.tolist(),
    }]

    failed = None
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

        if sample_idx < len(target_times) and state.t + 0.5 * dt >= target_times[sample_idx]:
            last = state.history[-1]
            samples.append({
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
            })
            sample_idx += 1

    last = state.history[-1] if state.history else {}
    return {
        "D_amplitude": float(D),
        "repair": bool(repair),
        "resolution": N,
        "r_max": r_max,
        "cfl": cfl,
        "final_time": float(state.t),
        "status": "numerical_failure" if failed else "completed",
        "failure": failed,
        "tau_final": float(state.tau),
        "r": r.tolist(),
        "samples": samples,
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


def main(output="runs/0star-center-repair-ab/report.json"):
    rows = []
    for D in (0.0, 1e-4):
        for repair in (True, False):
            rows.append(run_case(D=D, repair=repair))
    out = {
        "status": "completed",
        "diagnostic_only": True,
        "purpose": "A center-regularity repair A/B plus B radial tau(r), alpha(r)",
        "cases": rows,
    }
    p = Path(output)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
