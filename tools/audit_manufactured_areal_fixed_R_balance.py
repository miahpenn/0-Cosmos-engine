#!/usr/bin/env python3
"""Manufactured exact fixed-R mass/current interpolation control.

The analytic continuum identity is exact by construction:
    m(t,R) = (M + A*t) f(R)
    chi(t,R) = 1 - 2*m(t,R)/R
    F(t,R) = d_t m|_R = A*f(R)
so -R/2 * d_t chi|_R == F at every radius. Profiles are sampled on a
smoothly moving areal-radius map R(t,r), then queried on common fixed-R
samples using the same four-point Lagrange interpolation pattern as the
accepted-state diagnostics. This isolates interpolation/query error; it
does not exercise the BSSN tensor transformation or admit a physical law.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
from pathlib import Path
from typing import Any

import numpy as np

SOURCE = "4aec88d0e3ff028e0a65a8ac36b30948625df792"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RMAX = 80.0
T0 = 0.60
CFL = 0.03
M0 = 0.30
M_RATE = 0.12
PROFILE_SCALE = 7.0
MAP_AMPLITUDE = 0.08
MAP_SCALE = 7.0
ANNULUS = (2.0, 10.0)
RESOLUTIONS = (60, 120, 240)


def finite(name: str, value: Any) -> np.ndarray:
    arr = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(arr)):
        raise RuntimeError(f"{name} contains non-finite values")
    return arr


def shape_function(R: np.ndarray | float) -> np.ndarray | float:
    rr = np.asarray(R, dtype=np.float64)
    L3 = PROFILE_SCALE**3
    out = rr**3 / (L3 + rr**3)
    return float(out) if out.ndim == 0 else out


def areal_map(r: np.ndarray, t: float) -> np.ndarray:
    x = r / MAP_SCALE
    return finite("mapped areal radius", r * (1.0 + MAP_AMPLITUDE * t * np.exp(-(x * x))))


def analytic_fields(R: np.ndarray, t: float) -> dict[str, np.ndarray]:
    f = finite("shape function", shape_function(R))
    mass = finite("analytic mass", (M0 + M_RATE * t) * f)
    chi = finite("analytic chi", 1.0 - 2.0 * mass / R)
    current = finite("analytic fixed-R current", M_RATE * f)
    return {"R": finite("R", R), "mass": mass, "chi": chi, "current": current}


def interp4(x: float, xs: np.ndarray, ys: np.ndarray) -> float:
    """Four-point Lagrange interpolation, matching the production diagnostic helper."""
    if not (xs[0] <= x <= xs[-1]) or xs.size < 4:
        raise RuntimeError("cubic query is outside the in-range four-point domain")
    idx = int(np.searchsorted(xs, x))
    start = min(max(idx - 2, 0), xs.size - 4)
    nodes = xs[start : start + 4]
    weights = np.ones(4, dtype=np.float64)
    for k in range(4):
        for ell in range(4):
            if ell != k:
                weights[k] *= (x - nodes[ell]) / (nodes[k] - nodes[ell])
    return float(weights @ ys[start : start + 4])


def sample(values: np.ndarray, profile_R: np.ndarray, radii: np.ndarray, method: str) -> np.ndarray:
    if method == "linear":
        return finite(f"linear interpolation {method}", np.interp(radii, profile_R, values))
    if method == "cubic":
        return finite(
            f"four-point interpolation {method}",
            np.array([interp4(float(x), profile_R, values) for x in radii], dtype=np.float64),
        )
    raise ValueError(f"unsupported interpolation method: {method}")


def stats(values: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    selected = finite("selected residual", values)[mask]
    if selected.size == 0:
        raise RuntimeError("empty fixed-radius analysis region")
    return {
        "rms": float(np.sqrt(np.mean(selected * selected))),
        "max_abs": float(np.max(np.abs(selected))),
    }


def observed_order(coarse: float, fine: float) -> float:
    if coarse <= 0 or fine <= 0:
        raise RuntimeError("orders require positive norms")
    return float(math.log(coarse / fine, 2.0))


def run_arm(n: int, step: float) -> dict[str, Any]:
    dr = RMAX / n
    h = step
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr
    times = (T0, T0 + h, T0 + 2.0 * h)
    profiles = []
    for t in times:
        R = areal_map(r, t)
        if np.any(np.diff(R) <= 0.0):
            raise RuntimeError("manufactured areal map is not monotone")
        profiles.append(analytic_fields(R, t))

    lo = max(float(p["R"][0]) for p in profiles)
    hi = min(float(p["R"][-1]) for p in profiles)
    lo = max(lo, ANNULUS[0])
    hi = min(hi, ANNULUS[1])
    if not hi > lo:
        raise RuntimeError("no common fixed-R annulus")
    radii = np.linspace(ANNULUS[0], ANNULUS[1], max(300, 2 * n))
    mask = (radii >= ANNULUS[0]) & (radii <= ANNULUS[1])
    truth = finite("exact current on query radii", M_RATE * shape_function(radii))

    out: dict[str, Any] = {
        "resolution": n,
        "dr": dr,
        "dt": h,
        "time_samples": list(times),
        "continuum_identity_max_abs": 0.0,
        "methods": {},
    }
    for method in ("linear", "cubic"):
        chi = [sample(p["chi"], p["R"], radii, method) for p in profiles]
        current = [sample(p["current"], p["R"], radii, method) for p in profiles]
        forward_rate = finite("forward mass-rate", -0.5 * radii * (chi[1] - chi[0]) / h)
        centered_rate = finite("centered mass-rate", -0.5 * radii * (chi[2] - chi[0]) / (2.0 * h))
        endpoint_current = finite("endpoint-average current", 0.5 * (current[0] + current[1]))
        midpoint_current = current[1]
        forward_balance = finite("forward balance residual", forward_rate - endpoint_current)
        centered_balance = finite("centered balance residual", centered_rate - midpoint_current)
        forward_vs_truth = finite("forward minus exact current", forward_rate - truth)
        centered_vs_truth = finite("centered minus exact current", centered_rate - truth)
        current_interp_error = finite("interpolated-current error", midpoint_current - truth)
        out["methods"][method] = {
            "forward_minus_endpoint_average_current": stats(forward_balance, mask),
            "centered_minus_midpoint_current": stats(centered_balance, mask),
            "forward_mass_rate_minus_exact_current": stats(forward_vs_truth, mask),
            "centered_mass_rate_minus_exact_current": stats(centered_vs_truth, mask),
            "interpolated_current_minus_exact_current": stats(current_interp_error, mask),
        }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/manufactured-areal-fixed-R-balance-2026-10-25/summary.json"))
    args = parser.parse_args()
    results = []
    for n in RESOLUTIONS:
        dr = RMAX / n
        for label, factor in (("dt", 1.0), ("dt_half", 0.5)):
            results.append({"step_arm": label, **run_arm(n, factor * CFL * dr)})

    # Observed spatial orders for the matched dt/dt-half controls.
    order_summary: dict[str, Any] = {}
    for step_arm in ("dt", "dt_half"):
        for method in ("linear", "cubic"):
            for witness in ("forward_minus_endpoint_average_current", "centered_minus_midpoint_current"):
                values = [
                    next(item for item in results if item["resolution"] == n and item["step_arm"] == step_arm)["methods"][method][witness]["rms"]
                    for n in RESOLUTIONS
                ]
                order_summary[f"{step_arm}.{method}.{witness}"] = {
                    "rms_by_resolution": values,
                    "observed_orders": [observed_order(values[0], values[1]), observed_order(values[1], values[2])],
                }

    document = {
        "diagnostic": "manufactured exact fixed-areal-radius mass/current interpolation control",
        "status": "research-only; interpolation/query calibration, not BSSN or conservation admission",
        "source_pin": SOURCE,
        "vendor_pin": VENDOR,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__},
        "analytic_profile": {
            "mass": "(M0 + M_RATE*t) * R^3/(L^3 + R^3)",
            "chi": "1 - 2*mass/R",
            "fixed_R_current": "M_RATE * R^3/(L^3 + R^3)",
            "exact_identity": "-R/2 * d_t(chi)|_R == fixed_R_current",
            "M0": M0, "M_RATE": M_RATE, "L": PROFILE_SCALE,
            "moving_areal_map": "R(t,r)=r*(1+MAP_AMPLITUDE*t*exp(-(r/MAP_SCALE)^2))",
            "MAP_AMPLITUDE": MAP_AMPLITUDE, "MAP_SCALE": MAP_SCALE,
            "rmax": RMAX, "annulus": list(ANNULUS),
            "dt_rule": "dt=CFL*dr, compared with dt/2", "CFL": CFL, "t0": T0,
        },
        "results": results,
        "observed_order_summary": order_summary,
        "interpretation": "A linear-vs-cubic difference quantifies fixed-R query sensitivity for a continuum balance that is exactly zero by construction. It does not test the BSSN tensor transformation, the physical source law, or the parity-aware FV face reconstruction.",
        "production_changed": False,
        "scientific_admission": False,
        "full_matter_coupled_propagation_identity_closed": False,
        "total_stress_energy_sector_conservation_closed": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "n_resolutions": len(RESOLUTIONS), "results": results, "orders": order_summary}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
