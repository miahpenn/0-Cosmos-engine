"""Covariant spherical matter evolution interfaces.

This module defines the conservative/primitive identities required by the
V5.5 production matter sector. It deliberately does not choose a numerical
reconstruction, Riemann solver, or artificial interface source.

Locked V5.5 exchange:
    ∇_μ T_DM^{μν}    = + beta rho_DM ∇^ν phi
    ∇_μ T_phi^{μν}   = - beta rho_DM ∇^ν phi

Baryons are pressureless dust. Radiation is a perfect fluid with p=rho/3.
Shear is outside the exact spherical branch.
"""
from dataclasses import dataclass
from math import sqrt
from typing import Mapping


@dataclass(frozen=True)
class Primitive:
    """Eulerian primitive variables for a radial perfect fluid."""
    rho: float
    pressure: float
    velocity: float

    def lorentz_factor(self) -> float:
        v2 = self.velocity * self.velocity
        if v2 < 0.0 or v2 >= 1.0:
            raise ValueError("radial Eulerian velocity must satisfy |v| < 1")
        return 1.0 / sqrt(1.0 - v2)


@dataclass(frozen=True)
class Conservative:
    """Valencia conservative variables (D, S_r, tau)."""
    D: float
    S_r: float
    tau: float


def primitive_to_conservative(q: Primitive) -> Conservative:
    W = q.lorentz_factor()
    h = q.rho + q.pressure
    D = q.rho * W
    S_r = h * W * W * q.velocity
    tau = h * W * W - q.pressure - D
    return Conservative(D, S_r, tau)


def conservative_to_primitive(U: Conservative, pressure: float,
                              velocity: float) -> Primitive:
    """Construct primitives after an external conservative inversion.

    The inversion itself is intentionally supplied by the production
    integrator; no root finder or EOS closure is hidden here.
    """
    if U.D <= 0.0:
        raise ValueError("conserved rest-mass density must be positive")
    W = 1.0 / sqrt(1.0 - velocity * velocity)
    rho = U.D / W
    return Primitive(rho=rho, pressure=pressure, velocity=velocity)


def radial_flux(q: Primitive, U: Conservative, alpha: float,
                beta: float) -> Conservative:
    """Valencia radial coordinate flux for lapse alpha and shift beta."""
    transport = alpha * q.velocity - beta
    return Conservative(
        D=U.D * transport,
        S_r=U.S_r * transport + alpha * q.pressure,
        tau=U.tau * transport + alpha * q.pressure * q.velocity,
    )


@dataclass(frozen=True)
class DarkMatterFourForce:
    """Archive-locked DM/scalar exchange before 3+1 projection."""
    beta: float
    rho_dm: float
    dphi_contravariant: float

    @property
    def Q_dm(self) -> float:
        return self.beta * self.rho_dm * self.dphi_contravariant

    @property
    def Q_phi(self) -> float:
        return -self.Q_dm


def dark_matter_four_force(beta: float, rho_dm: float,
                           dphi_contravariant: float) -> DarkMatterFourForce:
    return DarkMatterFourForce(beta, rho_dm, dphi_contravariant)


def perfect_fluid_source_projection(
    four_force: Mapping[str, float] | None,
) -> Mapping[str, float]:
    """Pass a covariant four-force projection to the conservative integrator.

    Expected keys are supplied by the metric/source layer (for example normal
    energy and radial momentum projections). No projection is guessed here.
    """
    return {} if four_force is None else dict(four_force)


def radiation_primitive(rho: float, velocity: float = 0.0) -> Primitive:
    if rho < 0.0:
        raise ValueError("radiation density must be non-negative")
    return Primitive(rho=rho, pressure=rho / 3.0, velocity=velocity)


def baryon_primitive(rest_density: float, velocity: float = 0.0) -> Primitive:
    if rest_density < 0.0:
        raise ValueError("baryon rest density must be non-negative")
    return Primitive(rho=rest_density, pressure=0.0, velocity=velocity)


def validate_spherical_scope(*, include_shear: bool = False) -> None:
    """Guard the exact spherical branch against silently isotropizing shear."""
    if include_shear:
        raise ValueError(
            "Bianchi-I shear is not an isotropic spherical perfect-fluid "
            "component; lift the geometry before enabling the shear sector."
        )
