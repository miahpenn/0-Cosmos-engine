"""First-party native spherical geometry and isolated evolution routines.

The production CMC kernel consumes native_explicit_rhs for the explicit metric
fields a/b/X. Its alpha output is not integrated by the CMC stage path; lapse
remains CMC-solved, and beta/B remain pinned by the production stage contract.
Native primary/Lambda, constraints, gauge, and full PIRK2 composition remain
isolated candidates until each has its own bounded integration contract.
The native package itself has no dependency on the vendor implementation.
"""

from .operators import NativeSphericalCellGrid
from .geometry import metric_connection_lambda, native_bssn_geometry_terms
from .constraints import native_bssn_constraints
from .evolution import (
    native_cmc_explicit_metric_rhs,
    native_explicit_rhs,
    native_gauge_b_rhs,
    native_lambda_l2_rhs,
    native_lambda_l3_rhs,
    native_primary_l2_rhs,
    native_primary_l3_rhs,
    native_pirk2_step,
)

__all__ = [
    "NativeSphericalCellGrid",
    "metric_connection_lambda",
    "native_bssn_geometry_terms",
    "native_bssn_constraints",
    "native_cmc_explicit_metric_rhs",
    "native_explicit_rhs",
    "native_gauge_b_rhs",
    "native_lambda_l2_rhs",
    "native_lambda_l3_rhs",
    "native_primary_l2_rhs",
    "native_primary_l3_rhs",
    "native_pirk2_step",
]
