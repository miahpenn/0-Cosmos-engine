#!/usr/bin/env python3
"""Diagnostic-only radiation four-force counterfactual.

Null: beta_rad = 0.0 (current production)
Candidate: beta_rad = -0.04 (same structural form as DM)

Frozen-state + one ordinary production step. No long campaign, coefficient
fitting, or merge to main. This runner is intentionally not launched by CI.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from engine.production_kernel import V55ProductionKernel

RESOLUTIONS = (80, 160)
BETA_NULL = 0.0
BETA_CAND = -0.04
CFL = 0.03


def _snapshot(kernel, state, label: str) -> dict:
    obs = kernel.diagnostics(state, profiles=True)
    return {
        "label": label,
        "t": float(state.t),
        "tau": float(state.tau),
        "H_eff": float(obs.get("H_eff", np.nan)),
        "max_abs_H_admitted": float(
            obs.get("hamiltonian_max_admitted", obs.get("hamiltonian_max", np.nan))
        ),
        "max_abs_M_admitted": float(
            obs.get("momentum_max_admitted", obs.get("momentum_max", np.nan))
        ),
        "max_abs_C_Lambda": float(obs.get("connection_max", np.nan)),
        # The current production diagnostics do not yet define a
        # Misner-Sharp current; retain explicit NaNs rather than inventing one.
        "MS_current_max": float(obs.get("misner_sharp_current_max", np.nan)),
        "MS_current_rms": float(obs.get("misner_sharp_current_rms", np.nan)),
        "rad_margin_i0": float(obs.get("radiation_margin_i0", np.nan)),
        "rad_margin_admitted": float(obs.get("radiation_margin_admitted", np.nan)),
    }


def run_one(resolution: int, beta_rad: float) -> dict:
    kernel = V55ProductionKernel()
    kernel.set_diagnostic_beta_rad(beta_rad)
    state = kernel.initialize(
        resolution=resolution, r_max=40.0, include_radiation=True
    )
    dr = state.grid.centers[1] - state.grid.centers[0]
    dt = CFL * dr

    pre = _snapshot(kernel, state, "t0")
    state = kernel.step(state, dt)
    post = _snapshot(kernel, state, "t0+dt")

    return {
        "resolution": resolution,
        "beta_rad": beta_rad,
        "dt": float(dt),
        "pre": pre,
        "post": post,
    }


def main() -> None:
    results = []
    for resolution in RESOLUTIONS:
        for beta in (BETA_NULL, BETA_CAND):
            print(f"Running N={resolution} beta_rad={beta}")
            results.append(run_one(resolution, beta))

    out = Path("artifacts/radiation_four_force_counterfactual.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(results, indent=2, sort_keys=True)
    out.write_text(payload + "\\n", encoding="utf-8")
    digest = hashlib.sha256((payload + "\\n").encode("utf-8")).hexdigest()
    print(f"Wrote {out} SHA-256={digest}")


if __name__ == "__main__":
    main()
