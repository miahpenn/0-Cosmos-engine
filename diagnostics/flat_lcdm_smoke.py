"""Known-physics cosmology check: flat dust + cosmological constant.

Runs engine.cosmos.cosmos_rhs with beta=0, radiation/shear/baryons off,
constant scalar potential, and Omega_m=0.3, Omega_Lambda=0.7 at H0=1.
Compares a(t) and H(t) to the exact flat LCDM solution and checks Friedmann
closure and dust conservation. This exercises the homogeneous COSMOS lane,
not the full spatial production kernel.
"""
from __future__ import annotations

import json
from math import asinh, sinh, sqrt
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from engine.cosmos import CosmosParams, cosmos_rhs

OMEGA_M = 0.30
OMEGA_L = 0.70
T_FINAL = 1.0
STEPS = (100, 200, 400)


def exact_a(t):
    # Normalize a(0)=1; exact flat matter+Lambda solution.
    phase0 = asinh(sqrt(OMEGA_L / OMEGA_M))
    return (OMEGA_M / OMEGA_L) ** (1.0 / 3.0) * sinh(
        phase0 + 1.5 * sqrt(OMEGA_L) * t
    ) ** (2.0 / 3.0)


def exact_h(a):
    return sqrt(OMEGA_M / a**3 + OMEGA_L)


def rhs(y, p):
    return cosmos_rhs(y, p)


def add(y, k, scale):
    return tuple(v + scale * kv for v, kv in zip(y, k))


def rk4_step(y, dt, p):
    k1 = rhs(y, p)
    k2 = rhs(add(y, k1, 0.5 * dt), p)
    k3 = rhs(add(y, k2, 0.5 * dt), p)
    k4 = rhs(add(y, k3, dt), p)
    return tuple(
        v + dt * (k1i + 2*k2i + 2*k3i + k4i) / 6.0
        for v, k1i, k2i, k3i, k4i in zip(y, k1, k2, k3, k4)
    )


def run(nsteps):
    # rho_m(today)=3*Omega_m; V0=rho_Lambda=3*Omega_Lambda because H0=1,
    # and the code's flat Friedmann normalization is H^2=rho/3.
    p = CosmosParams(
        lam=0.0, D=0.0, beta=0.0, V0=3.0 * OMEGA_L,
        omega_de_today=OMEGA_L, omega_r_today=0.0,
        omega_shear_today=0.0, G=1.0,
    )
    y = (1.0, 1.0, 0.0, 0.0, 3.0 * OMEGA_M, 0.0, 0.0, 0.0)
    dt = T_FINAL / nsteps
    max_friedmann = 0.0
    max_dust_conservation = 0.0
    for _ in range(nsteps):
        y = rk4_step(y, dt, p)
        a, H, phi, pi, rho_dm, rho_b, rho_r, rho_s = y
        rho, _ = __import__("engine.cosmos", fromlist=["cosmos_rho_p"]).cosmos_rho_p(y, p)
        max_friedmann = max(max_friedmann, abs(H * H - rho / 3.0))
        max_dust_conservation = max(
            max_dust_conservation, abs(rho_dm * a**3 - 3.0 * OMEGA_M)
        )
    a, H = y[0], y[1]
    return {
        "steps": nsteps,
        "dt": dt,
        "a_numeric": a,
        "a_exact": exact_a(T_FINAL),
        "a_abs_error": abs(a - exact_a(T_FINAL)),
        "H_numeric": H,
        "H_exact": exact_h(exact_a(T_FINAL)),
        "H_abs_error": abs(H - exact_h(exact_a(T_FINAL))),
        "max_friedmann_abs_residual": max_friedmann,
        "max_dust_conservation_abs_error": max_dust_conservation,
    }


def main():
    results = [run(n) for n in STEPS]
    a_ratios = [
        results[i]["a_abs_error"] / results[i + 1]["a_abs_error"]
        for i in range(len(results) - 1)
    ]
    h_ratios = [
        results[i]["H_abs_error"] / results[i + 1]["H_abs_error"]
        for i in range(len(results) - 1)
    ]
    passed = (
        all(12.0 < x < 20.0 for x in a_ratios)
        and all(12.0 < x < 20.0 for x in h_ratios)
        and results[-1]["max_friedmann_abs_residual"] < 1e-9
        and results[-1]["max_dust_conservation_abs_error"] < 1e-10
    )
    payload = {
        "test": "flat LCDM homogeneous background",
        "scope": "engine.cosmos.cosmos_rhs + independent RK4 driver; not full spatial production-kernel validation",
        "parameters": {
            "H0": 1.0, "Omega_m": OMEGA_M, "Omega_Lambda": OMEGA_L,
            "Omega_r": 0.0, "Omega_shear": 0.0, "beta": 0.0,
            "potential_is_constant": True,
        },
        "expected": "exact flat dust+Lambda solution with a(0)=1 and H(a)^2=0.3/a^3+0.7",
        "results": results,
        "a_error_ratios_when_dt_halved": a_ratios,
        "H_error_ratios_when_dt_halved": h_ratios,
        "pass": passed,
    }
    out = ROOT / "runs" / "flat_lcdm_smoke"
    out.mkdir(parents=True, exist_ok=True)
    (out / "flat_lcdm_smoke.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
    if not passed:
        raise SystemExit("FAIL: flat LCDM analytic comparison/convergence/closure failed")
    print("PASS: homogeneous COSMOS lane matches the exact flat LCDM solution")


if __name__ == "__main__":
    main()
