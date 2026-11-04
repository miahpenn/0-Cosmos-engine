"""Reference witnesses for the accepted-state full-Bianchi regression."""
from __future__ import annotations

import math
import numpy as np

REFERENCE = {
    60: {
        "pointwise": (1.8406e-5, 2.0498e-6, (1.8406e-5, 1.4017e-5, 7.7777e-6)),
        "finite_volume": (9.6404e-6, 1.8089e-6, (-1.2950e-6, 9.6404e-6, 6.9830e-6)),
        "tangent_epsilon_delta": 1.38e-13,
        "pointwise_predictor_epsilon_delta": 7.11e-15,
    },
    120: {
        "pointwise": (1.5272e-6, 1.6403e-7, (1.5272e-6, 1.4480e-6, 9.8044e-7)),
        "finite_volume": (1.2827e-6, 1.3547e-7, (-1.2827e-6, 5.9546e-7, 5.6046e-7)),
        "tangent_epsilon_delta": 1.93e-13,
        "pointwise_predictor_epsilon_delta": 1.93e-14,
    },
    240: {
        "pointwise": (3.9897e-7, 1.0667e-8, (3.9897e-7, 1.2524e-7, 9.3559e-8)),
        "finite_volume": (4.0771e-8, 8.9879e-9, (3.5862e-8, 7.1513e-9, 2.5872e-8)),
        "tangent_epsilon_delta": 5.60e-13,
        "pointwise_predictor_epsilon_delta": 4.36e-14,
    },
}


def assert_reference_regression(result: dict) -> dict:
    """Check both signed channels and the archived directional epsilon sweep."""
    n = int(result["resolution"])
    if n not in REFERENCE:
        raise ValueError(f"no reference regression for N={n}")
    expected = REFERENCE[n]
    radii = np.asarray(result["grid_centers"], dtype=np.float64)
    mask = radii <= 10.0
    if not np.any(mask):
        raise ValueError("r<=10 regression mask is empty")
    checks = []
    for channel in ("pointwise", "finite_volume"):
        values = np.asarray(result["channels"]["eps2e-05"][channel]["residual_tangent_minus_predictor"])
        exp_max, exp_annulus, exp_cells = expected[channel]
        annulus = (radii >= 5.0) & (radii <= 10.0)
        witnesses = [
            ("max_abs_r_le_10", float(np.max(np.abs(values[mask]))), exp_max),
            ("max_abs_annulus_5_to_10", float(np.max(np.abs(values[annulus]))), exp_annulus),
        ]
        witnesses.extend((f"signed_cell_{i}", float(values[i]), exp_cells[i]) for i in range(3))
        for name, actual, target in witnesses:
            passed = math.isclose(actual, target, rel_tol=5.0e-3, abs_tol=1.0e-12)
            checks.append({"channel": channel, "witness": name, "actual": actual, "expected": target, "passed": passed})
    tangent_delta = float(np.max(np.abs(
        np.asarray(result["channels"]["eps2e-05"]["measured_tangent"])
        - np.asarray(result["channels"]["eps1e-05"]["measured_tangent"])
    )[mask]))
    predictor_delta = float(np.max(np.abs(
        np.asarray(result["channels"]["eps2e-05"]["pointwise"]["predictor"])
        - np.asarray(result["channels"]["eps1e-05"]["pointwise"]["predictor"])
    )[mask]))
    for name, actual, target in (
        ("tangent_epsilon_delta", tangent_delta, expected["tangent_epsilon_delta"]),
        ("pointwise_predictor_epsilon_delta", predictor_delta, expected["pointwise_predictor_epsilon_delta"]),
    ):
        passed = math.isclose(actual, target, rel_tol=0.10, abs_tol=1.0e-15)
        checks.append({"channel": "epsilon_sweep", "witness": name, "actual": actual, "expected": target, "passed": passed})
    failures = [item for item in checks if not item["passed"]]
    if failures:
        raise RuntimeError(f"accepted-state Bianchi regression mismatch: {failures}")
    return {"passed": True, "residual_sign": "measured_tangent_minus_predictor", "checks": checks}
