"""Short homogeneous cosmology smoke test through the full production kernel.

Uses the engine's archive-normalized homogeneous operating point with local
S and D perturbations disabled. Checks the initialized Hubble rate against
the initializer's Friedmann normalization, then executes one actual
V55ProductionKernel.step and records finite-state/constraint diagnostics.

This is not a claim that the production scalar potential is constant-Lambda
LCDM; the production phi potential is exponential and rolls.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.production_kernel import V55ProductionKernel
from engine.v55_initial import build_initial_data
from engine import v55_pirk_adapter as adapter
from engine.scalar_system import cosmos_potential
from engine.cosmos import CosmosParams, cosmos_rhs

OUT = ROOT / "runs" / "production_cosmology_smoke"


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
    dt = 1.0e-4
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

    # Match the archive initializer's own H0 normalization independently
    # from the geometry K field, using its published operating-point inputs.
    phi0 = 0.179055
    pi0 = 0.00341913
    rho_dm0 = 2.5857e-5
    rho_b0 = 4.0306e-6
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += rho_dm0 + rho_b0
    expected_h0 = float(np.sqrt(rho0 / 3.0))
    h0_error = abs(float(initial["H_eff"]) - expected_h0)

    # Independent homogeneous reference integration from the same archive
    # operating point. This checks the full production step's H response,
    # not merely that it returned finite arrays.
    params = CosmosParams(
        lam=1.0, D=0.10, beta=-0.04, V0=8.242522415500654,
        omega_de_today=0.69, omega_r_today=0.0,
        omega_shear_today=0.0, G=1.0,
    )
    y0 = (
        33.8983, expected_h0, phi0, pi0, rho_dm0, rho_b0, 0.0, 0.0,
    )
    y_reference = rk4_step(y0, dt, params)

    final = kernel.step(state, dt)
    after = kernel.diagnostics(final, profiles=False)
    fields = {
        "a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"
    }
    max_geometry_abs = {
        name: float(np.max(np.abs(getattr(final.geometry, name))))
        for name in sorted(fields)
    }
    all_finite = (
        all(np.all(np.isfinite(getattr(final.geometry, n))) for n in fields)
        and all(np.all(np.isfinite(getattr(final.scalars, n)))
                for n in ("S", "PS", "D", "PD", "phi", "Pi"))
        and all(np.all(np.isfinite(getattr(getattr(final.matter, n), q)))
                for n in ("dark_matter", "baryons", "radiation")
                for q in ("rest", "energy_t", "momentum_r"))
        and all(np.isfinite(v) for v in (
            float(initial["H_eff"]), expected_h0, h0_error,
            float(after["H_eff"]), float(after["hamiltonian_max"]),
            float(after["momentum_max"]), float(after["connection_max"]),
        ))
    )
    reference_h_error = abs(float(after["H_eff"]) - float(y_reference[1]))
    # The production integrator is second-order; this short-step comparison
    # uses a conservative 1e-9 absolute H tolerance against independent RK4.
    reference_tolerance = 1.0e-9
    reference_matches = reference_h_error <= reference_tolerance
    init_tolerance = 1.0e-12
    init_matches = h0_error <= init_tolerance
    payload = {
        "test": "full production-kernel homogeneous archive cosmology one-step smoke",
        "scope": "V55ProductionKernel.initialize + actual V55ProductionKernel.step; local S/D amplitudes zero; radiation off; native exponential COSMOS potential remains dynamic; not constant-Lambda LCDM validation",
        "grid": {"resolution": resolution, "r_max": r_max, "dr": float(state.grid.dr)},
        "dt": dt,
        "initial": {
            "H_eff": float(initial["H_eff"]),
            "H0_from_initializer_inputs": expected_h0,
            "H0_abs_error": h0_error,
            "hamiltonian_max": float(initial["hamiltonian_max"]),
            "momentum_max": float(initial["momentum_max"]),
            "connection_max": float(initial["connection_max"]),
        },
        "after_one_step": {
            "t": float(final.t),
            "H_eff": float(after["H_eff"]),
            "H_eff_homogeneous_RK4_reference": float(y_reference[1]),
            "H_eff_abs_error_vs_homogeneous_RK4": reference_h_error,
            "hamiltonian_max": float(after["hamiltonian_max"]),
            "momentum_max": float(after["momentum_max"]),
            "connection_max": float(after["connection_max"]),
            "max_abs_geometry_fields": max_geometry_abs,
        },
        "finite": all_finite,
        "initial_friedmann_normalization_matches": init_matches,
        "initialization_tolerance": init_tolerance,
        "homogeneous_RK4_reference_tolerance": reference_tolerance,
        "homogeneous_reference_matches": reference_matches,
        "pass": bool(all_finite and init_matches and reference_matches),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "production_cosmology_smoke.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    if not payload["pass"]:
        raise SystemExit("FAIL: production cosmology disagrees with the homogeneous reference or produced nonfinite state")
    print("PASS: production kernel matches the homogeneous reference over one short step")


if __name__ == "__main__":
    main()
