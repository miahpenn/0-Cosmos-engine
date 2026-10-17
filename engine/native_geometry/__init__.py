"""First-party native spherical geometry and isolated evolution candidates.

This package contains independent grid/operator, metric/connection/Ricci,
constraint, and narrowly scoped evolution RHS routines. It has no vendor
dependency and is not selected by production evolution; vendor comparisons
remain in the separate diagnostic/test layer.
"""

from .operators import NativeSphericalCellGrid
from .geometry import metric_connection_lambda, native_bssn_geometry_terms
from .constraints import native_bssn_constraints
from .evolution import native_lambda_l2_rhs, native_lambda_l3_rhs, native_primary_l2_rhs, native_primary_l3_rhs

__all__ = [
    "NativeSphericalCellGrid",
    "metric_connection_lambda",
    "native_bssn_geometry_terms",
    "native_bssn_constraints",
    "native_lambda_l2_rhs",
    "native_lambda_l3_rhs",
    "native_primary_l2_rhs",
    "native_primary_l3_rhs",
]
