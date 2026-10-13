"""Single total spherical stress-energy assembly.

All matter sources presented to Einstein are assembled here from the solved
local scalars, COSMOS scalar, conservative fluids, and their metric
projection. There is no separate interface reservoir.
"""
from dataclasses import dataclass
import numpy as np

from .scalar_system import ScalarFields
from . import v55_matter
from .v55_matter import V55MatterState


@dataclass(frozen=True)
class TotalStressEnergy:
    rho: np.ndarray
    pr: np.ndarray
    pt: np.ndarray
    j: np.ndarray
    scalar_rho: np.ndarray
    scalar_pr: np.ndarray
    scalar_pt: np.ndarray
    scalar_j: np.ndarray
    fluid_rho: np.ndarray
    fluid_pr: np.ndarray
    fluid_pt: np.ndarray
    fluid_j: np.ndarray

    def trace_spatial(self) -> np.ndarray:
        return self.pr + 2.0 * self.pt

    def finite(self) -> bool:
        return all(
            np.all(np.isfinite(x))
            for x in (
                self.rho, self.pr, self.pt, self.j,
                self.scalar_rho, self.scalar_pr, self.scalar_pt, self.scalar_j,
                self.fluid_rho, self.fluid_pr, self.fluid_pt, self.fluid_j,
            )
        )


def assemble_total_stress_energy(grid, geometry, scalars: ScalarFields,
                                 matter: V55MatterState) -> TotalStressEnergy:
    # Keep one projection assembly for geometry/CMC and constraint,
    # boundary, and diagnostic consumers. Future sector additions (including
    # any approved EM projection) must enter this shared path once.
    projection = v55_matter.total_matter_projection(
        grid, geometry, scalars, matter
    )
    scalar = projection["scalar"]
    fluid = projection["fluid"]
    return TotalStressEnergy(
        rho=projection["rho"],
        pr=projection["pr"],
        pt=projection["pt"],
        j=projection["j"],
        scalar_rho=scalar["rho"],
        scalar_pr=scalar["pr"],
        scalar_pt=scalar["pt"],
        scalar_j=scalar["j"],
        fluid_rho=fluid["rho"],
        fluid_pr=fluid["pr"],
        fluid_pt=fluid["pt"],
        fluid_j=fluid["j"],
    )
