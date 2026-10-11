"""Central-proper-time CMC gauge.

This is an archive-derived time reparameterization, not a 0* prescription.
The archived clock relation is d tau = alpha_c dt on the chosen central
worldline.  Starting from the existing CMC lapse, divide the entire lapse by
alpha_c so the evolution coordinate is central proper time.

The stored production lapse is the central-proper-time lapse. Before solving
the CMC equation again, recover the canonically normalized CMC lapse from the
outer value: the CMC normalization is alpha(R)=1, while the transformed lapse
has alpha_tau(R)=1/alpha_c. This keeps the CMC solve in its original clock
normalization and applies the coordinate transformation only after that solve.

No spatial boundary condition, source, coefficient, or physical term is
changed.
"""
from __future__ import annotations

import numpy as np

from .cmc_gauge import solve_cmc_lapse, target_kdot


def recover_cmc_clock(
    geometry,
) -> tuple[object, float]:
    """Recover the outer-normalized CMC lapse and old-clock scale.

    The central-time lapse satisfies
        alpha_tau = alpha_cmc / alpha_cmc(0)
    and therefore
        alpha_tau(R) = 1 / alpha_cmc(0)
    because the underlying CMC gauge uses alpha_cmc(R)=1.

    Returns a geometry copy carrying the canonically normalized CMC lapse and
    the corresponding central proper-time scale alpha_cmc(0).
    """
    alpha_tau = np.asarray(geometry.alpha, dtype=float)
    alpha_boundary = float(alpha_tau[-1])
    if not np.isfinite(alpha_boundary) or alpha_boundary <= 0.0:
        raise ValueError(
            "central proper-time gauge requires a positive outer lapse"
        )

    alpha_cmc = alpha_tau / alpha_boundary
    alpha_c = 1.0 / alpha_boundary

    recovered = geometry.copy()
    recovered.alpha = alpha_cmc

    if not np.all(np.isfinite(alpha_cmc)) or float(np.min(alpha_cmc)) <= 0.0:
        raise FloatingPointError(
            "recovered CMC lapse is non-finite/non-positive"
        )
    return recovered, alpha_c


def central_proper_time_target_kdot(
    grid,
    geometry,
    scalars,
    matter,
) -> float:
    """Return the CMC Kdot target expressed in central proper time."""
    recovered, alpha_c = recover_cmc_clock(geometry)
    kdot_cmc = target_kdot(grid, recovered, scalars, matter)
    return float(kdot_cmc / alpha_c)


def solve_central_proper_time_lapse(
    grid,
    geometry,
    scalars,
    matter,
    outer_frac: float = 0.20,
) -> tuple[np.ndarray, float]:
    """Return the existing CMC lapse in central proper-time coordinates."""
    recovered, _ = recover_cmc_clock(geometry)

    alpha_cmc, kdot_cmc = solve_cmc_lapse(
        grid, recovered, scalars, matter, outer_frac=outer_frac
    )
    alpha_c = float(alpha_cmc[0])
    if not np.isfinite(alpha_c) or alpha_c <= 0.0:
        raise ValueError(
            "central proper-time reparameterization requires positive central lapse"
        )

    alpha_tau = np.asarray(alpha_cmc, dtype=float) / alpha_c
    if not np.all(np.isfinite(alpha_tau)) or float(np.min(alpha_tau)) <= 0.0:
        raise FloatingPointError(
            "central proper-time lapse reparameterization is non-finite/non-positive"
        )

    return alpha_tau, float(kdot_cmc / alpha_c)
