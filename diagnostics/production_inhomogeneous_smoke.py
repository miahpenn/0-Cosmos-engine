"""Short inhomogeneous production-kernel convergence test.

A small localized S-field pulse is evolved with the actual V55ProductionKernel.
Three timestep resolutions are compared against a finer-step production run at
the same final coordinate time. This is a numerical convergence diagnostic,
not an independent physical solution or a campaign validation.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel

OUT = ROOT / "runs" / "production_inhomogeneous_smoke"
T_FINAL = 2.0e-4
STEPS = (1.0e-4, 5.0e-5, 2.5e-5)
REFERENCE_DT = 6.25e-6


def evolve(dt: float):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=64,
        r_max=32.0,
        amplitude=1.0e-3,
        width=7.0,
        D_amplitude=0.0,
        include_radiation=False,
    )
    initial_diag = kernel.diagnostics(state, profiles=False)
    nsteps = round(T_FINAL / dt)
    if not np.isclose(nsteps * dt, T_FINAL, rtol=0.0, atol=1.0e-15):
        raise ValueError(f"dt={dt} does not divide final time")
    for _ in range(nsteps):
        state = kernel.step(state, dt)
    final_diag = kernel.diagnostics(state, profiles=False)
    arrays = {}
    for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        arrays[f"geometry.{name}"] = np.asarray(getattr(state.geometry, name), dtype=float)
    for name in ("S", "PS", "D", "PD", "phi", "Pi"):
        arrays[f"scalars.{name}"] = np.asarray(getattr(state.scalars, name), dtype=float)
    for species in ("dark_matter", "baryons", "radiation"):
        item = getattr(state.matter, species)
        for name in ("rest", "energy_t", "momentum_r"):
            arrays[f"matter.{species}.{name}"] = np.asarray(getattr(item, name), dtype=float)
    finite = all(np.all(np.isfinite(v)) for v in arrays.values())
    finite = finite and all(np.isfinite(float(final_diag[k])) for k in (
        "H_eff", "hamiltonian_max", "momentum_max", "connection_max",
        "lapse_min", "determinant_min",
    ))
    return arrays, {
        "dt": dt,
        "steps": nsteps,
        "t_final": float(state.t),
        "finite": bool(finite),
        "H_eff": float(final_diag["H_eff"]),
        "initial_hamiltonian_max": float(initial_diag["hamiltonian_max"]),
        "final_hamiltonian_max": float(final_diag["hamiltonian_max"]),
        "initial_momentum_max": float(initial_diag["momentum_max"]),
        "final_momentum_max": float(final_diag["momentum_max"]),
        "final_connection_max": float(final_diag["connection_max"]),
        "final_lapse_min": float(final_diag["lapse_min"]),
        "final_determinant_min": float(final_diag["determinant_min"]),
        "S_center_initial": float(1.0e-3 * np.exp(-(state.grid.centers[0] / 7.0) ** 2)),
        "S_center_final": float(state.scalars.S[0]),
    }


def normalized_state_error(candidate, reference):
    per_field = {}
    for key in sorted(reference):
        scale = max(float(np.sqrt(np.mean(reference[key] ** 2))), 1.0e-12)
        delta = candidate[key] - reference[key]
        per_field[key] = float(np.sqrt(np.mean(delta ** 2)) / scale)
    aggregate = float(np.sqrt(np.mean(np.square(list(per_field.values())))))
    return aggregate, per_field


def main():
    runs = {}
    meta = {}
    for dt in (*STEPS, REFERENCE_DT):
        arrays, info = evolve(dt)
        runs[dt] = arrays
        meta[dt] = info

    ref = runs[REFERENCE_DT]
    errors = []
    for dt in STEPS:
        error, per_field = normalized_state_error(runs[dt], ref)
        errors.append(error)
        meta[dt]["normalized_state_error_vs_fine_production_run"] = error
        meta[dt]["largest_field_error"] = max(per_field.items(), key=lambda item: item[1])
    ratios = [errors[i] / errors[i + 1] for i in range(len(errors) - 1)]
    finite = all(item["finite"] for item in meta.values())
    convergence = all(2.0 < ratio < 6.5 for ratio in ratios)
    passed = finite and convergence
    payload = {
        "test": "inhomogeneous S-pulse full production-kernel timestep convergence",
        "scope": "actual V55ProductionKernel; localized S pulse amplitude 1e-3; D amplitude zero; radiation off; same final coordinate time; fine production run used as numerical reference, not an independent analytic solution",
        "grid": {"resolution": 64, "r_max": 32.0, "dr": 0.5},
        "final_time": T_FINAL,
        "reference_dt": REFERENCE_DT,
        "runs": [meta[dt] for dt in (*STEPS, REFERENCE_DT)],
        "normalized_state_error_ratios_when_dt_halved": ratios,
        "checks": {
            "all_fields_and_diagnostics_finite": finite,
            "observed_second_order_convergence": convergence,
        },
        "tolerances": {"error_ratio_for_convergence": [2.0, 6.5]},
        "pass": bool(passed),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "production_inhomogeneous_smoke.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    if not passed:
        raise SystemExit("FAIL: inhomogeneous production pulse failed finiteness or convergence gate")
    print("PASS: inhomogeneous production pulse is finite and shows timestep convergence")


if __name__ == "__main__":
    main()
