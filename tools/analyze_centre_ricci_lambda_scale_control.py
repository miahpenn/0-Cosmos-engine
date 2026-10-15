#!/usr/bin/env python3
"""Manufactured length-scale control for the regular-centre Lambda_r channel.

This is a read-only analytic/numerical experiment. It uses the fourth-order
cell operators already transcribed from the pinned vendor implementation,
but its smooth metric profile and continuum derivative are independently
prescribed. It does not initialize or evolve the engine.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from analyze_centre_ricci_derivative_channels import (
    cell_derivative_fourth,
    metric_lambda_term_responses,
)

SCHEMA = "centre_ricci_lambda_radial_scale_control_v1"
RESOLUTIONS = (80, 160, 320)
LENGTH_SCALES = (4.0, 8.0, 16.0)
R_MAX = 32.0
AMPLITUDE = 0.1
ANNULUS = (1.0, 4.0)


def _finite(name: str, value: np.ndarray | float) -> None:
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} contains non-finite values")


def analyze_profile(
    n: int,
    length_scale: float,
    *,
    r_max: float = R_MAX,
    amplitude: float = AMPLITUDE,
    annulus: tuple[float, float] = ANNULUS,
) -> dict[str, Any]:
    """Compare two discrete derivative representations with a closed form."""
    if n < 16 or r_max <= 0 or length_scale <= 0 or amplitude <= 0:
        raise ValueError("n, r_max, length_scale and amplitude must be positive and n >= 16")
    if not (0 < annulus[0] < annulus[1] < r_max):
        raise ValueError("annulus must be positive and lie strictly inside the domain")

    dr = r_max / n
    r = (np.arange(n, dtype=float) + 0.5) * dr

    # Smooth, even, determinant-consistent manufactured metric family.
    # q(0)=1, q'(0)=0; the length scale changes the same smooth shape without
    # altering the fixed physical domain or the analytic formula.
    exponential = np.exp(-((r / length_scale) ** 2))
    q = 1.0 + amplitude * (1.0 - exponential)
    q_prime = 2.0 * amplitude * r / (length_scale**2) * exponential
    a = q ** (-4.0 / 3.0)
    b = q ** (2.0 / 3.0)
    X = q ** (1.0 / 3.0)
    delta = 1.0 / b - 1.0 / a

    # Third term in the metric-defined conformal connection.
    lambda_radial = 2.0 * delta / r

    # Independent continuum derivative, obtained by differentiating the
    # prescribed q(r), not by applying either discrete derivative route.
    delta_prime = (
        -(2.0 / 3.0) * q ** (-5.0 / 3.0)
        -(4.0 / 3.0) * q ** (1.0 / 3.0)
    ) * q_prime
    exact_lambda_prime = 2.0 * delta_prime / r - 2.0 * delta / (r * r)

    sampled_d4 = cell_derivative_fourth(lambda_radial, dr, parity=-1)
    d1_a = cell_derivative_fourth(a, dr, parity=+1)
    d1_b = cell_derivative_fourth(b, dr, parity=+1)
    chain_rule = (
        -2.0 * delta / (r * r)
        + 2.0 * (-d1_b / (b * b) + d1_a / (a * a)) / r
    )

    err_d4 = sampled_d4 - exact_lambda_prime
    err_chain = chain_rule - exact_lambda_prime
    representation_gap = sampled_d4 - chain_rule

    annulus_mask = (r >= annulus[0]) & (r <= annulus[1])
    if int(annulus_mask.sum()) < 4:
        raise ValueError("fixed annulus has fewer than four cell centres")
    _finite("q", q)
    _finite("exact derivative", exact_lambda_prime)
    _finite("D4 error", err_d4)
    _finite("chain-rule error", err_chain)

    def channel_summary(error: np.ndarray) -> dict[str, Any]:
        return {
            "cell0": {
                "r": float(r[0]),
                "signed_error": float(error[0]),
                "error_over_r2": float(error[0] / r[0]**2),
            },
            "cell1": {
                "r": float(r[1]),
                "signed_error": float(error[1]),
                "error_over_r2": float(error[1] / r[1]**2),
            },
            "fixed_annulus_max_abs": float(np.max(np.abs(error[annulus_mask]))),
            "fixed_annulus_signed_min": float(np.min(error[annulus_mask])),
            "fixed_annulus_signed_max": float(np.max(error[annulus_mask])),
            "fixed_annulus_rms": float(np.sqrt(np.mean(error[annulus_mask]**2))),
        }

    term_responses, full_channel_closure = metric_lambda_term_responses(
        a, b, X, r, dr
    )
    radial_channel = (2.0 / 3.0) * X * X * representation_gap
    radial_term_response = term_responses["2/r*(1/b-1/a)"]
    radial_identity_residual = float(
        np.max(np.abs(radial_channel - radial_term_response))
    )

    return {
        "resolution": int(n),
        "length_scale": float(length_scale),
        "amplitude": float(amplitude),
        "r_max": float(r_max),
        "dr": float(dr),
        "fixed_annulus": [float(annulus[0]), float(annulus[1])],
        "grid_points_in_annulus": int(annulus_mask.sum()),
        "all_profiles_finite": True,
        "D4_sampled_odd_Lambda_minus_exact": channel_summary(err_d4),
        "D1_chain_rule_minus_exact": channel_summary(err_chain),
        "D4_minus_chain_rule": channel_summary(representation_gap),
        "channel_response": {
            "cell0": float(radial_channel[0]),
            "cell1": float(radial_channel[1]),
            "fixed_annulus_max_abs": float(np.max(np.abs(radial_channel[annulus_mask]))),
            "full_three_term_recombination_max_abs": float(full_channel_closure),
            "radial_term_identity_max_abs": radial_identity_residual,
        },
    }


def _observed_orders(errors: list[float], resolutions: tuple[int, ...]) -> list[dict[str, float]]:
    rows = []
    for i in range(len(errors) - 1):
        coarse, fine = float(errors[i]), float(errors[i + 1])
        ratio = coarse / fine if fine > 0 else math.inf
        order = math.log(ratio, resolutions[i + 1] / resolutions[i]) if fine > 0 else math.inf
        rows.append({
            "from_N": int(resolutions[i]),
            "to_N": int(resolutions[i + 1]),
            "error_ratio": float(ratio),
            "observed_order": float(order),
        })
    return rows


def analyze_scales() -> dict[str, Any]:
    """Run only the deterministic manufactured-profile matrix."""
    rows_by_scale = []
    for length in LENGTH_SCALES:
        rows = [
            analyze_profile(n, length)
            for n in RESOLUTIONS
        ]
        order_rows = {}
        for channel, field in (
            ("D4_sampled_odd_Lambda_minus_exact", "D4_sampled_odd_Lambda_minus_exact"),
            ("D1_chain_rule_minus_exact", "D1_chain_rule_minus_exact"),
            ("D4_minus_chain_rule", "D4_minus_chain_rule"),
        ):
            errors = [row[field]["fixed_annulus_max_abs"] for row in rows]
            order_rows[channel] = _observed_orders(errors, RESOLUTIONS)
        rows_by_scale.append({
            "length_scale": float(length),
            "resolutions": rows,
            "fixed_annulus_observed_orders": order_rows,
        })
    return {
        "schema": SCHEMA,
        "experiment": "manufactured profile only; no engine initialization, state capture, or evolution",
        "profile": "q(r)=1+A*(1-exp(-(r/L)^2)); a=q^(-4/3), b=q^(2/3), X=q^(1/3)",
        "amplitude": AMPLITUDE,
        "r_max": R_MAX,
        "resolutions": list(RESOLUTIONS),
        "length_scales": list(LENGTH_SCALES),
        "fixed_annulus": list(ANNULUS),
        "method": {
            "sampled_route": "fourth-order D4 of the sampled odd Lambda_r with odd parity",
            "chain_route": "D1(a), D1(b) differentiated chain expression for d/dr [2*(1/b-1/a)/r]",
            "independent_reference": "closed-form derivative of the analytically prescribed q(r)",
            "centre_witness": "signed cell-0 and cell-1 errors and each error divided by r_i^2",
            "formal_order_rule": "orders are reported from fixed-physical-annulus max errors only; centre cells are not used as a same-coordinate convergence fit",
        },
        "scales": rows_by_scale,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="optional JSON output path")
    args = parser.parse_args()
    result = analyze_scales()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
