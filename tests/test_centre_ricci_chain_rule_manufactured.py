"""Manufactured regular-profile control for Ricci derivative channels."""
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from analyze_centre_ricci_derivative_channels import (
    cell_derivative_fourth,
    cell_second_derivative_fourth,
)


def manufactured_chain_channel(n: int) -> tuple[float, float, float]:
    """Return centre and fixed-annulus response for a smooth regular metric.

    This is a no-evolution control. q(0)=1 makes a(0)=b(0)=X(0)=1;
    a, b and X are even regular fields and the reconstructed Lambda is odd.
    """
    rmax = 80.0
    dr = rmax / n
    r = (np.arange(n, dtype=float) + 0.5) * dr
    q = 1.0 + 0.1 * (1.0 - np.exp(-(r / 10.0) ** 2))
    a = q ** (-4.0 / 3.0)
    b = q ** (2.0 / 3.0)
    X = q ** (1.0 / 3.0)

    ap = cell_derivative_fourth(a, dr, parity=+1)
    bp = cell_derivative_fourth(b, dr, parity=+1)
    app = cell_second_derivative_fourth(a, dr, parity=+1)
    bpp = cell_second_derivative_fourth(b, dr, parity=+1)
    lambda_metric = (
        ap / (2.0 * a * a)
        - bp / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )

    # Compare two representations of the derivative of the same sampled
    # metric reconstruction: D4(Lambda_metric) versus its analytic
    # product/chain-rule expression using the pinned D1/D2 operands.
    d4_lambda = cell_derivative_fourth(lambda_metric, dr, parity=-1)
    d_chain = (
        app / (2.0 * a * a)
        - ap * ap / (a ** 3)
        - bpp / (a * b)
        + bp * (ap * b + a * bp) / (a * a * b * b)
        - 2.0 * (1.0 / b - 1.0 / a) / (r * r)
        + 2.0 * (-bp / (b * b) + ap / (a * a)) / r
    )
    # Pinned vendor TF-Ricci derivative response coefficient.
    channel = (2.0 / 3.0) * X * X * (d4_lambda - d_chain)
    if not np.all(np.isfinite(channel)):
        raise AssertionError("manufactured derivative channel is not finite")

    annulus = (r >= 5.0) & (r <= 10.0)
    if not np.any(annulus):
        raise AssertionError("manufactured grid did not sample r in [5, 10]")
    return (
        float(channel[0]),
        float(np.max(np.abs(channel[annulus]))),
        dr,
    )


def test_regular_manufactured_chain_channel_decays_under_refinement():
    results = [manufactured_chain_channel(n) for n in (40, 80, 160, 320)]
    centre_abs = [abs(row[0]) for row in results]
    annulus_max = [row[1] for row in results]

    # The first cell moves toward the origin as dr shrinks, so this is an
    # origin-regularity control rather than a same-coordinate pointwise fit.
    assert all(new < old / 8.0 for old, new in zip(centre_abs, centre_abs[1:]))
    # The annulus is held fixed in physical coordinates; a ratio near 16
    # is consistent with fourth-order truncation on this smooth profile.
    assert all(new < old / 8.0 for old, new in zip(annulus_max, annulus_max[1:]))
    for old, new in zip(annulus_max, annulus_max[1:]):
        ratio = old / new
        assert 8.0 < ratio < 32.0
