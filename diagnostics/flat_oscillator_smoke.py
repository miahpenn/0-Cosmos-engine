"""Cheap known-physics smoke test: [flat-oscillator-smoke] homogeneous massive scalar in flat spacetime.

For the S field with flat geometry, K=0, beta=0 and spatially constant S,
the implemented equations reduce exactly to dS/dt=PS, dPS/dt=-S.
Thus S=A cos(t), PS=-A sin(t). This deliberately isolates the scalar RHS;
it is not a full cosmology or a validation of coupled geometry evolution.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine import v55_pirk_adapter as adapter
from engine.scalar_system import ScalarFields, scalar_rhs_arrays

OUT = ROOT / "runs" / "flat_oscillator_smoke"
AMPLITUDE = 1.0e-3
T_FINAL = 2.0 * np.pi
STEPS = (400, 800, 1600)


def rhs(grid, geom, S, PS):
    zeros = np.zeros_like(S)
    # Other fields are zero; only the S,PS subsystem is measured.
    fields = ScalarFields(S, PS, zeros, zeros, zeros, zeros)
    out = scalar_rhs_arrays(grid, geom, fields)
    return out.S, out.PS


def run(nsteps):
    grid_ops, vacuum, _ = adapter.vendor_modules()
    grid = grid_ops.SphericalCellGrid(40, 40.0)
    geom = vacuum.flat_state(grid)
    geom.alpha[:] = 1.0
    geom.beta[:] = 0.0
    geom.K[:] = 0.0
    geom.a[:] = 1.0
    geom.b[:] = 1.0
    geom.X[:] = 1.0

    S = np.full(len(grid.centers), AMPLITUDE)
    PS = np.zeros_like(S)
    dt = T_FINAL / nsteps

    for _ in range(nsteps):
        k1s, k1p = rhs(grid, geom, S, PS)
        k2s, k2p = rhs(grid, geom, S + 0.5 * dt * k1s, PS + 0.5 * dt * k1p)
        k3s, k3p = rhs(grid, geom, S + 0.5 * dt * k2s, PS + 0.5 * dt * k2p)
        k4s, k4p = rhs(grid, geom, S + dt * k3s, PS + dt * k3p)
        S = S + (dt / 6.0) * (k1s + 2*k2s + 2*k3s + k4s)
        PS = PS + (dt / 6.0) * (k1p + 2*k2p + 2*k3p + k4p)

    exact_s = AMPLITUDE * np.cos(T_FINAL)
    exact_ps = -AMPLITUDE * np.sin(T_FINAL)
    err_s = float(np.max(np.abs(S - exact_s)))
    err_ps = float(np.max(np.abs(PS - exact_ps)))
    return {
        "steps": nsteps,
        "dt": dt,
        "t_final": T_FINAL,
        "S_exact": float(exact_s),
        "PS_exact": float(exact_ps),
        "S_numeric_center": float(S[0]),
        "PS_numeric_center": float(PS[0]),
        "max_abs_S_error": err_s,
        "max_abs_PS_error": err_ps,
        "max_error": max(err_s, err_ps),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    results = [run(n) for n in STEPS]
    ratios = [
        results[i]["max_error"] / results[i + 1]["max_error"]
        for i in range(len(results) - 1)
    ]
    # RK4 should approach a factor of 16 error reduction when dt halves.
    passed = all(np.isfinite(r["max_error"]) for r in results) and all(
        ratio > 10.0 for ratio in ratios
    )
    payload = {
        "test": "homogeneous massive scalar oscillator in flat spacetime",
        "scope": "scalar RHS plus independent RK4 harness; not full kernel validation",
        "expected": "S=A*cos(t), PS=-A*sin(t), period=2*pi",
        "results": results,
        "error_reduction_ratios_when_dt_halved": ratios,
        "pass": passed,
    }
    (OUT / "flat_oscillator_smoke.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    if not passed:
        raise SystemExit("FAIL: oscillator accuracy/convergence check did not pass")
    print("PASS: flat-space oscillator follows the analytic solution with fourth-order convergence")


if __name__ == "__main__":
    main()
