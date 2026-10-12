"""Production-kernel cosmology check against an independent homogeneous reference.

Runs the archive-normalized homogeneous operating point through the actual
V55ProductionKernel stepper at three decreasing timesteps. Compares H_eff
against an independent RK4 integration of engine.cosmos.cosmos_rhs and checks
the expected second-order convergence of the production stepper.

Local S/D amplitudes are zero and radiation is disabled. The native
exponential COSMOS potential remains dynamic; this is not constant-Lambda
LCDM validation.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import cosmos_potential
from engine.cosmos import CosmosParams, cosmos_rhs

OUT = ROOT / "runs" / "production_cosmology_smoke"
STEPS = (1.0e-4, 5.0e-5, 2.5e-5)


def add(y, k, scale):
    return tuple(v + scale * kv for v, kv in zip(y, k))


def rk4_step(y, dt, params):
    k1 = cosmos_rhs(y, params)
    k2 = cosmos_rhs(add(y, k1, 0.5 * dt), params)
    k3 = cosmos_rhs(add(y, k2, 0.5 * dt), params)
    k4 = cosmos_rhs(add(y, k3, dt), params)
    return tuple(
        v + dt * (k1i + 2*k2i + 2*k3i + k4i) / 6.0
        for v, k1i, k2i, k3i, k4i in zip(y, k1, k2, k3, k4)
    )


def main():
    resolution, r_max = 64, 32.0
    phi0, pi0 = 0.179055, 0.00341913
    rho_dm0, rho_b0 = 2.5857e-5, 4.0306e-6
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += rho_dm0 + rho_b0
    expected_h0 = float(np.sqrt(rho0 / 3.0))

    params = CosmosParams(
        lam=1.0, D=0.10, beta=-0.04, V0=8.242522415500654,
        omega_de_today=0.69, omega_r_today=0.0,
        omega_shear_today=0.0, G=1.0,
    )
    y0 = (33.8983, expected_h0, phi0, pi0, rho_dm0, rho_b0, 0.0, 0.0)

    results = []
    for dt in STEPS:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=resolution,
            r_max=r_max,
            amplitude=0.0,
            width=7.0,
            D_amplitude=0.0,
            include_radiation=False,
        )
        initial = kernel.diagnostics(state, profiles=False)
        h0_error = abs(float(initial["H_eff"]) - expected_h0)
        y_reference = rk4_step(y0, dt, params)

        final = kernel.step(state, dt)
        after = kernel.diagnostics(final, profiles=False)
        h_error = abs(float(after["H_eff"]) - float(y_reference[1]))
        names = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
        finite = (
            all(np.all(np.isfinite(getattr(final.geometry, n))) for n in names)
            and all(np.all(np.isfinite(getattr(final.scalars, n)))
                    for n in ("S", "PS", "D", "PD", "phi", "Pi"))
            and all(np.all(np.isfinite(getattr(getattr(final.matter, n), q)))
                    for n in ("dark_matter", "baryons", "radiation")
                    for q in ("rest", "energy_t", "momentum_r"))
            and all(np.isfinite(v) for v in (
                float(after["H_eff"]), h_error,
                float(after["hamiltonian_max"]),
                float(after["momentum_max"]), float(after["connection_max"]),
            ))
        )
        results.append({
            "dt": dt,
            "H_eff_initial": float(initial["H_eff"]),
            "H0_from_initializer_inputs": expected_h0,
            "H0_abs_error": h0_error,
            "H_eff_production": float(after["H_eff"]),
            "H_eff_homogeneous_RK4_reference": float(y_reference[1]),
            "H_eff_abs_error_vs_homogeneous_RK4": h_error,
            "hamiltonian_max": float(after["hamiltonian_max"]),
            "momentum_max": float(after["momentum_max"]),
            "connection_max": float(after["connection_max"]),
            "finite": bool(finite),
        })

    errors = [row["H_eff_abs_error_vs_homogeneous_RK4"] for row in results]
    ratios = [errors[i] / errors[i + 1] for i in range(len(errors) - 1)]
    init_matches = all(row["H0_abs_error"] <= 1.0e-12 for row in results)
    finite = all(row["finite"] for row in results)
    reference_close = all(error < 1.0e-9 for error in errors)
    # For a second-order production method, halving dt should reduce the
    # leading temporal error by about four in the asymptotic regime.
    second_order = all(2.5 < ratio < 5.5 for ratio in ratios)
    passed = finite and init_matches and reference_close and second_order

    payload = {
        "test": "full production homogeneous archive cosmology versus independent RK4",
        "scope": "actual V55ProductionKernel.step at three timesteps; local S/D amplitudes zero; radiation off; native exponential COSMOS potential remains dynamic; not constant-Lambda LCDM validation",
        "grid": {"resolution": resolution, "r_max": r_max, "dr": r_max / resolution},
        "steps": results,
        "H_error_ratios_when_dt_halved": ratios,
        "checks": {
            "all_fields_finite": finite,
            "initial_Friedmann_normalization_matches": init_matches,
            "all_H_errors_below_1e-9": reference_close,
            "second_order_convergence": second_order,
        },
        "tolerances": {
            "initial_H_abs": 1.0e-12,
            "reference_H_abs": 1.0e-9,
            "error_ratio_for_second_order": [2.5, 5.5],
        },
        "pass": bool(passed),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "production_cosmology_smoke.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    if not passed:
        raise SystemExit("FAIL: production cosmology did not match the homogeneous reference with expected convergence")
    print("PASS: production kernel matches the homogeneous reference with second-order convergence")


if __name__ == "__main__":
    main()
