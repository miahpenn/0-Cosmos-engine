"""Central-proper-time CMC gauge.

This is an archive-derived time reparameterization, not a 0* prescription.
The archived clock relation is d tau = alpha_c dt on the chosen central
worldline.  Starting from the existing CMC lapse, divide the entire lapse by
alpha_c so the evolution coordinate is central proper time.

A spatially constant lapse rescaling preserves the CMC spatial shape:
if L[alpha] = -Kdot/inv, then L[alpha/alpha_c]
= -(Kdot/alpha_c)/inv.  No spatial boundary condition, source, coefficient,
or physical term is changed.
"""
from __future__ import annotations

import numpy as np

from .cmc_gauge import solve_cmc_lapse


def solve_central_proper_time_lapse(
    grid,
    geometry,
    scalars,
    matter,
    outer_frac: float = 0.20,
) -> tuple[np.ndarray, float]:
    """Return the existing CMC lapse in central proper-time coordinates."""
    alpha, kdot = solve_cmc_lapse(
        grid, geometry, scalars, matter, outer_frac=outer_frac
    )
    alpha_c = float(alpha[0])
    if not np.isfinite(alpha_c) or alpha_c <= 0.0:
        raise ValueError(
            "central proper-time reparameterization requires positive central lapse"
        )

    alpha_tau = np.asarray(alpha, dtype=float) / alpha_c
    if not np.all(np.isfinite(alpha_tau)) or float(np.min(alpha_tau)) <= 0.0:
        raise FloatingPointError(
            "central proper-time lapse reparameterization is non-finite/non-positive"
        )

    return alpha_tau, float(kdot / alpha_c)
