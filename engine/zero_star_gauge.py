"""0* temporal-gauge candidate (diagnostic only).

This is not yet a production gauge.  It removes the worldtube lapse
normalization entirely and uses the already-derived central proper-time
meaning of t: with beta=0, alpha=1 gives dt=d tau on the central worldline.

No source, coefficient, boundary condition, bounce rule, or physical term is
added.  The existing CMC Kdot projection is retained only as a diagnostic so
the 0* criterion can be checked against the same geometric slice.
"""
import numpy as np

from .cmc_gauge import target_kdot


def solve_zero_star_lapse(grid, geometry, scalars, matter, outer_frac=0.20):
    """Return the synchronous proper-time candidate lapse.

    The candidate is deliberately alpha=1 everywhere and beta=0 in the
    production kernel.  This is a diagnostic temporal gauge, not a promoted
    physical assumption.
    """
    del outer_frac
    alpha = np.ones(grid.n, dtype=float)
    if not np.all(np.isfinite(alpha)):
        raise FloatingPointError("0* candidate lapse is non-finite")
    return alpha, target_kdot(grid, geometry, scalars, matter)
