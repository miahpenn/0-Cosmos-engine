"""First-party native spherical geometry core.

This package contains independent grid/operator and metric/connection/Ricci
routines. It has no vendor dependency and is not yet selected by production
evolution; vendor comparisons remain in the separate diagnostic/test layer.
"""

from .operators import NativeSphericalCellGrid
from .geometry import metric_connection_lambda, native_bssn_geometry_terms
from .constraints import native_bssn_constraints

__all__ = [
    "NativeSphericalCellGrid",
    "metric_connection_lambda",
    "native_bssn_geometry_terms",
    "native_bssn_constraints",
]
