"""Manufactured-source consistency audit for origin half-cell quadrature.

This is a research-only numerical discriminator. It does not import or modify
production engine code. The goal is to check whether a polynomial fit in r^2
recovers the moment over [0, dr/2] at the expected order under grid refinement.
"""
from __future__ import annotations

import json
import math
from typing import Callable

import numpy as np

ArraySource = Callable[[np.ndarray], np.ndarray]


def polynomial_even_source(r: np.ndarray) -> np.ndarray:
    """An even source exactly cubic in r^2."""
    r = np.asarray(r, dtype=float)
    return 1.0 + 0.3 * r**2 - 0.04 * r**4 + 0.004 * r**6


def smooth_even_source(r: np.ndarray) -> np.ndarray:
    """A smooth even source that is not cubic in r^2 on a fixed window."""
    r = np.asarray(r, dtype=float)
    return (
        1.0
        + 0.3 * r**2
        - 0.04 * r**4
        + 0.004 * r**6
        + 0.001 * r**8
        + 0.0002 * np.exp(-(r**2))
    )


def exact_polynomial_moment(r0: float) -> float:
    """Exact integral of polynomial_even_source from zero to r0."""
    return (
        r0
        + 0.3 * r0**3 / 3.0
        - 0.04 * r0**5 / 5.0
        + 0.004 * r0**7 / 7.0
    )


def exact_smooth_moment(r0: float) -> float:
    """Exact integral of smooth_even_source from zero to r0."""
    return (
        r0
        + 0.3 * r0**3 / 3.0
        - 0.04 * r0**5 / 5.0
        + 0.004 * r0**7 / 7.0
        + 0.001 * r0**9 / 9.0
        + 0.0002 * math.sqrt(math.pi) * math.erf(r0) / 2.0
    )


def fit_half_cell_moment(
    dr: float,
    fit_radius: float,
    *,
    degree_in_r_squared: int = 3,
    source: ArraySource = smooth_even_source,
) -> dict:
    """Fit an even polynomial on cell centres and integrate to the origin.

    The scaled coordinate z=(r/fit_radius)^2 spans the same polynomial space
    as raw powers of r^2 but avoids unnecessary conditioning problems when
    the fit window shrinks with dr.
    """
    if not math.isfinite(dr) or dr <= 0.0:
        raise ValueError("dr must be finite and positive")
    if not math.isfinite(fit_radius) or fit_radius <= 0.0:
        raise ValueError("fit_radius must be finite and positive")
    if degree_in_r_squared < 0:
        raise ValueError("degree_in_r_squared must be nonnegative")

    count = max(0, int(math.floor(fit_radius / dr + 0.5)))
    r = (np.arange(count, dtype=float) + 0.5) * dr
    r = r[r <= fit_radius + 1.0e-12 * max(1.0, fit_radius)]
    if len(r) < degree_in_r_squared + 2:
        raise ValueError("fit requires at least degree+2 cell-centre samples")

    z = (r / fit_radius) ** 2
    design = np.column_stack([z**k for k in range(degree_in_r_squared + 1)])
    values = np.asarray(source(r), dtype=float)
    if values.shape != r.shape or not np.all(np.isfinite(values)):
        raise ValueError("source must return finite values matching the sample shape")

    coeff, *_ = np.linalg.lstsq(design, values, rcond=None)
    r0 = 0.5 * dr
    moment = sum(
        float(coeff[k]) * r0 * (r0 / fit_radius) ** (2 * k) / (2 * k + 1)
        for k in range(degree_in_r_squared + 1)
    )
    return {
        "dr": float(dr),
        "fit_radius": float(fit_radius),
        "degree_in_r_squared": int(degree_in_r_squared),
        "sample_count": int(len(r)),
        "moment": float(moment),
        "fitted_intercept": float(coeff[0]),
        "source_at_origin": float(np.asarray(source(np.array([0.0])))[0]),
        "max_fit_residual": float(np.max(np.abs(values - design @ coeff))),
    }


def observed_order(coarse_error: float, fine_error: float) -> float:
    """Observed order for a factor-two refinement."""
    if coarse_error == 0.0 or fine_error == 0.0:
        return float("inf")
    return math.log(abs(coarse_error / fine_error), 2.0)


def run_audit() -> dict:
    """Return exact-polynomial controls and fixed/local-window error sweeps."""
    dr_values = (0.125, 0.0625, 0.03125)
    rows = {"fixed_physical_window": [], "local_six_cell_window": []}

    for mode in rows:
        errors = []
        for dr in dr_values:
            fit_radius = 2.5 if mode == "fixed_physical_window" else 6.0 * dr
            result = fit_half_cell_moment(dr, fit_radius, source=smooth_even_source)
            exact = exact_smooth_moment(0.5 * dr)
            error = result["moment"] - exact
            result["exact_moment"] = exact
            result["signed_moment_error"] = error
            errors.append(error)
            rows[mode].append(result)
        for i in range(1, len(rows[mode])):
            rows[mode][i]["observed_order_vs_previous"] = observed_order(
                errors[i - 1], errors[i]
            )

    exact_controls = []
    for dr in (0.5, 0.25, 0.125, 0.0625):
        for mode, fit_radius in (
            ("fixed_physical_window", 2.5),
            ("local_six_cell_window", 6.0 * dr),
        ):
            result = fit_half_cell_moment(
                dr, fit_radius, source=polynomial_even_source
            )
            error = result["moment"] - exact_polynomial_moment(0.5 * dr)
            exact_controls.append({
                "window": mode,
                "dr": dr,
                "signed_moment_error": float(error),
                "sample_count": result["sample_count"],
            })

    return {
        "method": "degree-3 least squares in r^2, with a numerically scaled basis",
        "production_source_changed": False,
        "production_evolution_steps": 0,
        "smooth_even_source": "1 + 0.3*r^2 - 0.04*r^4 + 0.004*r^6 + 0.001*r^8 + 0.0002*exp(-r^2)",
        "fixed_vs_local_window": rows,
        "exact_degree3_polynomial_controls": exact_controls,
        "interpretation": [
            "The fixed physical fit window tests the same support policy as the current copied-initializer counterfactual.",
            "The local window keeps six cell-centre samples as dr decreases.",
            "The manufactured source is diagnostic only and is not claimed to reproduce the engine's exact source profile.",
            "A result for this source is a method discriminator, not a unique production quadrature prescription.",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(run_audit(), indent=2, sort_keys=True))
