#!/usr/bin/env python3
"""Manufactured-metric control for the metric-Lambda Bianchi commutator.

Research-only. No evolution, accepted-state artifact, production equation, or
vendor code is imported or changed. Reproduces the pinned parity-aware
cell-centred fourth-order D1/D2 stencils and decomposes the induced Bianchi
commutator when D1(Lambda_metric) is compared with the primitive chain rule.
"""
from __future__ import annotations

import argparse
import json
import math
from typing import Any

import numpy as np

R_MAX = 80.0
DEFAULT_RESOLUTIONS = (60, 120, 240)
DEFAULT_CONTROLS = (
    ("reference", 0.1, 7.0),
    ("half_amplitude", 0.05, 7.0),
    ("shorter_scale", 0.1, 5.0),
)


def finite(name: str, values: Any) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(array)):
        raise FloatingPointError(f"non-finite values in {name}")
    return array


def d1(values: np.ndarray, dr: float, parity: int) -> np.ndarray:
    """Pinned SphericalCellGrid fourth-order D1 with parity ghosts."""
    u = finite("D1 input", values)
    if u.ndim != 1 or u.size < 5 or dr <= 0 or parity not in (-1, 1):
        raise ValueError("invalid D1 inputs")
    n = u.size
    ext = np.empty(n + 4, dtype=np.float64)
    ext[2:-2] = u
    ext[1], ext[0] = parity * u[0], parity * u[1]
    out = np.empty_like(u)
    for i in range(n - 2):
        stencil = ext[i:i + 5]
        out[i] = (
            stencil[0] - 8 * stencil[1] + 8 * stencil[3] - stencil[4]
        ) / (12 * dr)
    indices = np.arange(n - 5, n)
    for i in (n - 2, n - 1):
        offsets = indices.astype(np.float64) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[1] = 1.0
        out[i] = float(np.linalg.solve(matrix, target) @ u[indices]) / dr
    return finite("D1 output", out)


def d2(values: np.ndarray, dr: float, parity: int) -> np.ndarray:
    """Pinned SphericalCellGrid fourth-order D2 with parity ghosts."""
    u = finite("D2 input", values)
    if u.ndim != 1 or u.size < 5 or dr <= 0 or parity not in (-1, 1):
        raise ValueError("invalid D2 inputs")
    n = u.size
    ext = np.empty(n + 4, dtype=np.float64)
    ext[2:-2] = u
    ext[1], ext[0] = parity * u[0], parity * u[1]
    out = np.empty_like(u)
    for i in range(n - 2):
        stencil = ext[i:i + 5]
        out[i] = (
            -stencil[0] + 16 * stencil[1] - 30 * stencil[2]
            + 16 * stencil[3] - stencil[4]
        ) / (12 * dr**2)
    indices = np.arange(n - 5, n)
    for i in (n - 2, n - 1):
        offsets = indices.astype(np.float64) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[2] = 2.0
        out[i] = float(np.linalg.solve(matrix, target) @ u[indices]) / dr**2
    return finite("D2 output", out)


def one_resolution(n: int, amplitude: float, length: float) -> dict[str, Any]:
    dr = R_MAX / n
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr
    q = 1.0 + amplitude * (1.0 - np.exp(-((r / length) ** 2)))
    a, b, X = q ** (-4 / 3), q ** (2 / 3), q ** (-2 / 3)
    ap, bp = d1(a, dr, +1), d1(b, dr, +1)
    app, bpp = d2(a, dr, +1), d2(b, dr, +1)
    delta = 1.0 / b - 1.0 / a

    # Three pieces of Lambda_metric, assembled from the pinned D1(a), D1(b).
    terms = {
        "a_derivative": (
            ap / (2 * a * a),
            app / (2 * a * a) - ap * ap / (a * a * a),
        ),
        "b_derivative": (
            -bp / (a * b),
            -bpp / (a * b)
            + bp * (ap * b + a * bp) / (a * a * b * b),
        ),
        "radial_connection": (
            2 * delta / r,
            -2 * delta / (r * r)
            + 2 * (-bp / (b * b) + ap / (a * a)) / r,
        ),
    }
    chi_r = -0.5 * d1(X, dr, +1) / X
    P = 6 * chi_r + 1.5 * (2 / r + bp / b)
    channels: dict[str, np.ndarray] = {}
    for key, (assembled, chain) in terms.items():
        delta_lambda_prime = d1(assembled, dr, -1) - chain
        # Exact linear response in B_h to this Delta(Lambda_metric') channel.
        channels[key] = (
            0.5 * d1(X * X * delta_lambda_prime, dr, +1)
            + (2 / 3) * P * X * X * delta_lambda_prime
        )
    total = sum(channels.values())

    # Independently reconstruct the unsplit nested-vs-chain commutator.
    metric_lambda = sum(assembled for assembled, _chain in terms.values())
    full_chain = sum(chain for _assembled, chain in terms.values())
    full_delta_lambda_prime = d1(metric_lambda, dr, -1) - full_chain
    full_bianchi_commutator = (
        0.5 * d1(X * X * full_delta_lambda_prime, dr, +1)
        + (2 / 3) * P * X * X * full_delta_lambda_prime
    )

    areal_r = r * q
    mask = (areal_r >= 2.0) & (areal_r < 10.0)
    if int(mask.sum()) < 4:
        raise RuntimeError(
            f"too few cells in annulus for N={n}: {int(mask.sum())}"
        )
    rms = {
        key: float(np.sqrt(np.mean(values[mask] ** 2)))
        for key, values in channels.items()
    }
    rms["signed_sum"] = float(np.sqrt(np.mean(total[mask] ** 2)))
    scaled = {key: value * n**4 for key, value in rms.items()}
    closure = float(np.max(np.abs(total[mask] - full_bianchi_commutator[mask])))
    if closure > 5.0e-13:
        raise RuntimeError(
            f"signed channel decomposition failed identity closure: {closure:.6e}"
        )
    ratios = {
        "radial_over_a": rms["radial_connection"] / max(rms["a_derivative"], 1e-300),
        "radial_over_b": rms["radial_connection"] / max(rms["b_derivative"], 1e-300),
    }

    all_channels = {**channels, "signed_sum": total}
    origin_native_cells = {
        str(i): {key: float(values[i]) for key, values in all_channels.items()}
        for i in range(min(3, n))
    }
    region_rms: dict[str, Any] = {}
    for region, region_mask in (
        ("origin_R_lt_1", areal_r < 1.0),
        ("inner_shell_1_to_2", (areal_r >= 1.0) & (areal_r < 2.0)),
    ):
        region_rms[region] = {"cell_count": int(region_mask.sum()), "rms": {}}
        for key, values in all_channels.items():
            region_rms[region]["rms"][key] = (
                float(np.sqrt(np.mean(values[region_mask] ** 2)))
                if np.any(region_mask) else None
            )
    fixed_R = {}
    for radius in (0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0, 5.0, 8.0):
        fixed_R[str(radius)] = {
            key: float(np.interp(radius, areal_r, values))
            for key, values in channels.items()
        }
        fixed_R[str(radius)]["signed_sum"] = float(np.interp(radius, areal_r, total))
    return {
        "N": n,
        "dr": dr,
        "r_max": R_MAX,
        "amplitude": amplitude,
        "length_scale": length,
        "annulus_areal_R": [2.0, 10.0],
        "annulus_native_cell_count": int(mask.sum()),
        "rms": rms,
        "N4_scaled_rms": scaled,
        "radial_to_other_channel_ratios": ratios,
        "origin_native_cells_first_three": origin_native_cells,
        "region_rms": region_rms,
        "identity_closure_max_abs": closure,
        "fixed_R_signed_channels": fixed_R,
    }


def observed_order(e0: float, e1: float, n0: int, n1: int) -> float:
    return math.log(abs(e0) / abs(e1)) / math.log(n1 / n0)


def run_control(
    name: str, amplitude: float, length: float, resolutions: tuple[int, ...]
) -> dict[str, Any]:
    rows = [one_resolution(n, amplitude, length) for n in resolutions]
    keys = ("a_derivative", "b_derivative", "radial_connection", "signed_sum")
    orders = []
    for coarse, fine in zip(rows, rows[1:]):
        orders.append({
            "from_N": coarse["N"],
            "to_N": fine["N"],
            "observed_order": {
                key: observed_order(
                    coarse["rms"][key], fine["rms"][key],
                    coarse["N"], fine["N"],
                )
                for key in keys
            },
        })
    regional_orders = []
    for coarse, fine in zip(rows, rows[1:]):
        pair: dict[str, Any] = {
            "from_N": coarse["N"],
            "to_N": fine["N"],
            "regions": {},
        }
        for region in ("origin_R_lt_1", "inner_shell_1_to_2"):
            pair["regions"][region] = {}
            for key in keys:
                e0 = coarse["region_rms"][region]["rms"][key]
                e1 = fine["region_rms"][region]["rms"][key]
                pair["regions"][region][key] = (
                    observed_order(float(e0), float(e1), coarse["N"], fine["N"])
                    if e0 is not None and e1 is not None and e0 > 0 and e1 > 0
                    else None
                )
        regional_orders.append(pair)
    return {
        "control": name,
        "profile": {
            "q": "1 + A*(1-exp(-(r/L)^2))",
            "amplitude": amplitude,
            "length_scale": length,
            "a": "q^(-4/3)",
            "b": "q^(2/3)",
            "X": "q^(-2/3)",
            "exact_areal_R": "r*q",
        },
        "resolutions": rows,
        "observed_order_between_resolutions": orders,
        "observed_order_origin_and_inner_shell_between_resolutions": regional_orders,
        "interpretation_limits": [
            "Analytic prescribed metric only; no target accepted-state profiles and no evolution.",
            "This is the Bianchi response to replacing nested D1(Lambda_metric) with its D1/D2 chain-rule form; not a matter-coupled propagation identity.",
            "Channel decomposition is exact signed bookkeeping under the specified grouping, not independent physical causation.",
            "No derivative or Ricci representation is selected for production based on residual size.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolutions", nargs="+", type=int, default=list(DEFAULT_RESOLUTIONS))
    parser.add_argument("--output", type=str, default="")
    args = parser.parse_args()
    controls = [
        run_control(*control, tuple(args.resolutions))
        for control in DEFAULT_CONTROLS
    ]
    result = {
        "schema": "manufactured_lambda_metric_bianchi_commutator_v1",
        "status": "diagnostic_only",
        "source_anchor": "4aec88d0e3ff028e0a65a8ac36b30948625df792",
        "stencil_reference": (
            "miahpenn/0-Cosmos-engine tools/analyze_centre_ricci_derivative_channels.py "
            "(parity-aware D1/D2 contract, locally reproduced)"
        ),
        "vendor_reference": "d6052d605673ce9d82cdc99b0df79855fa2ea215",
        "method": (
            "For Lambda_metric = D1(a)/(2a^2)-D1(b)/(ab)+2/r*(1/b-1/a), "
            "compare D1 of each assembled term with its primitive D1/D2 chain derivative; "
            "map each signed delta Lambda-prime term into Delta B = "
            "0.5 D1(X^2 deltaLambdaPrime)+(2/3)P X^2 deltaLambdaPrime, "
            "P=6chi_r+3/2(2/r+b_r/b)."
        ),
        "controls": controls,
        "production_changed": False,
        "scientific_admission": False,
        "full_matter_coupled_propagation_identity_closed": False,
        "total_stress_energy_sector_conservation_closed": False,
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        with open(args.output, "w", encoding="utf-8") as stream:
            stream.write(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
