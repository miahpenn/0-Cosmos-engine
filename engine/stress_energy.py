"""Single total spherical stress-energy assembly.

All Einstein/diagnostic/boundary consumers share the canonical scalar+fluid
projection in :mod:`engine.v55_matter`. Keeping one projection entry point is
important for compatibility-preserving sector extensions: a future approved
EM contribution must enter the shared projection once, not be independently
re-created in a second consumer path. There is no separate interface reservoir.
"""
from dataclasses import dataclass
import numpy as np

from . import v55_matter
from .scalar_system import ScalarFields
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
    electromagnetic_rho: np.ndarray
    electromagnetic_pr: np.ndarray
    electromagnetic_pt: np.ndarray
    electromagnetic_j: np.ndarray

    def trace_spatial(self) -> np.ndarray:
        return self.pr + 2.0 * self.pt

    def finite(self) -> bool:
        return all(
            np.all(np.isfinite(x))
            for x in (
                self.rho, self.pr, self.pt, self.j,
                self.scalar_rho, self.scalar_pr, self.scalar_pt, self.scalar_j,
                self.fluid_rho, self.fluid_pr, self.fluid_pt, self.fluid_j,
                self.electromagnetic_rho, self.electromagnetic_pr,
                self.electromagnetic_pt, self.electromagnetic_j,
            )
        )


def assemble_total_stress_energy(
    grid,
    geometry,
    scalars: ScalarFields,
    matter: V55MatterState,
    *,
    electromagnetic_projection: dict[str, np.ndarray] | None = None,
) -> TotalStressEnergy:
    """Adapt the canonical projection to the diagnostic and boundary contract.

    The optional EM arrays must already use consumer-normalized units. None
    preserves the pre-EM projection exactly.
    """
    projection = v55_matter.total_matter_projection(
        grid, geometry, scalars, matter,
        electromagnetic_projection=electromagnetic_projection,
    )
    scalar = projection["scalar"]
    fluid = projection["fluid"]
    em = projection["electromagnetic"]
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
        electromagnetic_rho=em["rho"],
        electromagnetic_pr=em["pr"],
        electromagnetic_pt=em["pt"],
        electromagnetic_j=em["j"],
    )
