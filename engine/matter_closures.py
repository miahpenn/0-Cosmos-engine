"""Covariant matter-sector interfaces for the spherical production engine.

These formulas are kinematic stress-energy projections only.  Evolution laws
must be supplied by the corresponding conservative/primitive integrator;
nothing here invents an interface source or matter-identification rule.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class PerfectFluid:
    rho: float
    pressure: float
    radial_velocity: float = 0.0

    def stress_normal(self, lorentz_gamma=1.0):
        v = self.radial_velocity
        W = lorentz_gamma
        h = self.rho + self.pressure
        return {
            "rho": h * W * W - self.pressure,
            "j_r": h * W * W * v,
            "S_r": h * W * W * v * v + self.pressure,
            "S_t": self.pressure,
        }


@dataclass(frozen=True)
class RadiationState:
    rho: float

    @property
    def pressure(self):
        return self.rho / 3.0

    def stress_normal(self):
        return PerfectFluid(self.rho, self.pressure).stress_normal()


@dataclass(frozen=True)
class DustState:
    rest_density: float
    radial_velocity: float = 0.0

    def stress_normal(self, lorentz_gamma=1.0):
        return PerfectFluid(
            self.rest_density, 0.0, self.radial_velocity
        ).stress_normal(lorentz_gamma)


def dark_matter_exchange(beta, rho_dm, phi_gradient):
    """Archive-locked DM four-force density in the chosen radial frame."""
    q = beta * rho_dm * phi_gradient
    return q


def total_exchange_pair(beta, rho_dm, phi_gradient):
    q = dark_matter_exchange(beta, rho_dm, phi_gradient)
    return q, -q
