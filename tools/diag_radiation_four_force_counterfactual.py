#!/usr/bin/env python3
"""Diagnostic-only radiation four-force counterfactual.

Null: beta_rad = 0.0 (current production)
Candidate: beta_rad = -0.04 (same structural form as DM)

Frozen-state + one ordinary production step. No long campaign, coefficient
fitting, or merge to main. This runner is executed only by the diagnostic
branch's gated CI workflow after the focused regression and existing suite
pass. It records each case's exception, writes and hashes the complete
artifact, then exits nonzero if any case failed.
"""
from __future__ import annotations

import hashlib
import json
import traceback
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
    result = {
        "resolution": resolution,
        "beta_rad": beta_rad,
        "status": "failed",
        "dt": None,
        "pre": None,
        "post": None,
        "error_type": None,
        "error": None,
        "traceback": None,
    }
    try:
        kernel = V55ProductionKernel()
        kernel.set_diagnostic_beta_rad(beta_rad)
        state = kernel.initialize(
            resolution=resolution, r_max=40.0, include_radiation=True
        )
        dr = state.grid.centers[1] - state.grid.centers[0]
        dt = CFL * dr
        result["dt"] = float(dt)
        result["pre"] = _snapshot(kernel, state, "t0")
        state = kernel.step(state, dt)
        result["post"] = _snapshot(kernel, state, "t0+dt")
        result["status"] = "passed"
    except Exception as exc:
        result["error_type"] = type(exc).__name__
        result["error"] = str(exc)
        result["traceback"] = traceback.format_exc()
        print(
            f"CASE FAILED N={resolution} beta_rad={beta_rad}: "
            f"{type(exc).__name__}: {exc}"
        )
    return result


def main() -> None:
    results = []
    for resolution in RESOLUTIONS:
        for beta in (BETA_NULL, BETA_CAND):
            print(f"Running N={resolution} beta_rad={beta}")
            results.append(run_one(resolution, beta))

    out = Path("artifacts/radiation_four_force_counterfactual.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(results, indent=2, sort_keys=True, allow_nan=True)
    out.write_text(payload + "\n", encoding="utf-8")
    digest = hashlib.sha256((payload + "\n").encode("utf-8")).hexdigest()
    print(f"Wrote {out} SHA-256={digest}")

    failures = [r for r in results if r["status"] != "passed"]
    print(f"Counterfactual cases: {len(results) - len(failures)} passed, {len(failures)} failed")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
