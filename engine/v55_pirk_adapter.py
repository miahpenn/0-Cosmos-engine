"""Production V5.5 adapter around the pinned reference-metric BSSN/PIRK kernel.

Only the reference numerical operators and PIRK ordering are imported from the
vendor repository. V5.5 stress-energy and scalar equations remain local to
this repository.

This module is intentionally a stage/RHS adapter rather than a claimed
completed production integrator. The final campaign must use a kernel that
also advances the conservative matter states and records invariant witnesses.
"""
from dataclasses import dataclass
from pathlib import Path
import importlib
import sys
import numpy as np

from . import v55_matter as vm


VENDOR = Path(__file__).resolve().parents[1] / "vendor" / "bb-palatini-unified-r0"


def vendor_modules():
    if not VENDOR.exists():
        raise RuntimeError(
            "Pinned reference kernel is absent. Initialize git submodules first."
        )
    path = str(VENDOR)
    if path not in sys.path:
        sys.path.insert(0, path)
    return (
        importlib.import_module("r0_operators"),
        importlib.import_module("vacuum_bssn"),
        importlib.import_module("moving_puncture_pirk"),
    )


@dataclass(frozen=True)
class V55KernelSlice:
    geometry: object
    matter: vm.V55MatterState


def primary_l3_with_matter(grid, geometry, fields, matter, radiation_recovery_metric=None):
    _, vacuum, _ = vendor_modules()
    out = vacuum.primary_l3_rhs(grid, geometry)
    if radiation_recovery_metric is None:
        total = vm.total_matter_projection(grid, geometry, fields, matter)
    else:
        total = vm.total_matter_projection(
            grid, geometry, fields, matter, radiation_recovery_metric=radiation_recovery_metric
        )
    out["Aa"] = out["Aa"] - (16.0 * np.pi / 3.0) * geometry.alpha * (total["pr"] - total["pt"])
    out["K"] = out["K"] + 4.0 * np.pi * geometry.alpha * (total["rho"] + total["pr"] + 2.0 * total["pt"])
    return out
def lambda_l3_with_matter(
    grid, geometry, fields, matter, lambda_m=2.0, radiation_recovery_metric=None
):
    _, vacuum, _ = vendor_modules()
    if radiation_recovery_metric is None:
        total = vm.total_matter_projection(grid, geometry, fields, matter)
    else:
        total = vm.total_matter_projection(
            grid, geometry, fields, matter, radiation_recovery_metric=radiation_recovery_metric
        )
    return vacuum.lambda_l3_rhs(grid, geometry) - 8.0 * np.pi * lambda_m * geometry.alpha * total["j"] / geometry.a
def geometry_stage_terms(
    grid, geometry, fields, matter, lambda_m=2.0, radiation_recovery_metric=None
):
    """Return geometry blocks with stage metrics and optional radiation inversion metric."""
    _, vacuum, moving = vendor_modules()
    explicit = moving.moving_puncture_explicit_rhs(grid, geometry)
    l2 = vacuum.primary_l2_rhs(grid, geometry)
    if radiation_recovery_metric is None:
        l3 = primary_l3_with_matter(grid, geometry, fields, matter)
        ll3 = lambda_l3_with_matter(grid, geometry, fields, matter, lambda_m=lambda_m)
    else:
        l3 = primary_l3_with_matter(
            grid, geometry, fields, matter, radiation_recovery_metric=radiation_recovery_metric
        )
        ll3 = lambda_l3_with_matter(
            grid, geometry, fields, matter, lambda_m=lambda_m,
            radiation_recovery_metric=radiation_recovery_metric
        )
    ll2 = vacuum.lambda_l2_rhs(grid, geometry, lambda_m=lambda_m)
    return {"explicit": explicit, "primary_l2": l2, "primary_l3": l3, "lambda_l2": ll2, "lambda_l3": ll3}
def geometry_constraints(grid, geometry):
    _, vacuum, _ = vendor_modules()
    return vacuum.constraints(grid, geometry)
