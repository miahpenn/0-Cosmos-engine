"""First-party geometry-only primary-L2 RHS block.

The mixed radial trace-free curvature component and trace evolve in this block as

    Aa_L2 = -(D_r D_r alpha - (D^2 alpha)/3)
            + alpha * (R^r_r - R/3)
    K_L2  = -D^2 alpha

Here R^r_r is the native reference-metric BSSN radial mixed component returned
as Rr_r; R is the independently returned BSSN scalar witness. They are not
reconstructed from the metric-only Ricci fields, because that would erase the
evolved-Lambda principal-structure distinction off the connection-constraint
surface.

This is deliberately one RHS block, not a complete evolution stage. It does
not include primary-L3 terms (shift advection, nonlinear curvature terms,
matter projection or dissipation), explicit/gauge evolution, Lambda L2/L3, or
PIRK2 stage composition. It is a pure candidate and is not wired into the
production adapter.
"""

from __future__ import annotations

import numpy as np

from .geometry import native_bssn_geometry_terms


def native_primary_l2_rhs(grid, state) -> dict[str, np.ndarray]:
    """Return the native Aa and K primary-L2 RHS profiles.

    The function reads a caller-owned state without modifying it and returns
    new finite cell profiles. The state must expose the native geometry fields
    a, b, X, alpha and evolved Lambda with shape (n,). The radial trace-free
    component Aa and trace K are not used as inputs by this L2 block.
    """
    geometry = native_bssn_geometry_terms(grid, state)
    centers = np.asarray(grid.centers, dtype=float)
    alpha = np.asarray(state.alpha, dtype=float)
    if alpha.shape != centers.shape:
        raise ValueError(f"alpha must have shape {centers.shape}")
    if not np.all(np.isfinite(alpha)):
        raise ValueError("alpha must be finite")

    with np.errstate(over="ignore", invalid="ignore"):
        aa_rhs = -(
            geometry["DrDr_alpha"] - geometry["lap_alpha"] / 3.0
        ) + alpha * (
            geometry["Rr_r"] - geometry["R"] / 3.0
        )
        k_rhs = -geometry["lap_alpha"]

    out = {
        "Aa": np.asarray(aa_rhs, dtype=float),
        "K": np.asarray(k_rhs, dtype=float),
    }
    for name, values in out.items():
        if values.shape != centers.shape:
            raise ValueError(f"{name} RHS must have shape {centers.shape}")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} RHS produced non-finite values")
    return out
