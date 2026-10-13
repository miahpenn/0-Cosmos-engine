"""Archive-derived V5.5 unified spherical diagnostics.

This module contains only relations explicitly established by the V5.5
covariant global gate. It is a reference layer, not yet a time integrator.
"""

from dataclasses import dataclass
from math import pi


@dataclass(frozen=True)
class UnifiedUnits:
    """Bridge corrected COSMOS archive normalization to G=1 Einstein units."""
    G: float = 1.0

    @property
    def kappa(self) -> float:
        return 8.0 * pi * self.G

    def rho_phys(self, rho_norm: float) -> float:
        return rho_norm / self.kappa

    def phi_phys(self, phi_norm: float) -> float:
        return phi_norm / self.kappa**0.5

    def V_phys(self, V_norm: float) -> float:
        return V_norm / self.kappa


@dataclass(frozen=True)
class ScalarStress:
    rho: float
    p_r: float
    p_t: float
    j_r: float


def scalar_stress(phi_r: float, Pi: float, potential: float, a: float) -> ScalarStress:
    """Normal-frame stress tensor for one scalar in the archive convention."""
    grad2 = (phi_r * phi_r) / (a * a)
    kin2 = Pi * Pi
    return ScalarStress(
        rho=0.5 * (kin2 + grad2) + potential,
        p_r=0.5 * (kin2 + grad2) - potential,
        p_t=0.5 * (kin2 - grad2) - potential,
        j_r=-Pi * phi_r / (a * a),
    )


def sum_stress(*parts: ScalarStress) -> ScalarStress:
    return ScalarStress(
        rho=sum(x.rho for x in parts),
        p_r=sum(x.p_r for x in parts),
        p_t=sum(x.p_t for x in parts),
        j_r=sum(x.j_r for x in parts),
    )


def momentum_curvature(K_theta: float, dK_theta_dr: float, r: float, j_r: float,
                       G: float = 1.0) -> float:
    """Hamiltonian/momentum-gate relation for K_r."""
    return K_theta + r * dK_theta_dr + 4.0 * pi * G * r * j_r


def radial_B_rhs(r: float, rho: float, K_r: float, K_theta: float,
                 G: float = 1.0) -> float:
    """Right side of (r B)' from the V5.5 Hamiltonian constraint."""
    return (
        1.0
        - 8.0 * pi * G * r * r * rho
        + 2.0 * r * r * K_r * K_theta
        + r * r * K_theta * K_theta
    )


def misner_sharp_mass(r: float, B: float, beta: float) -> float:
    return 0.5 * r * (1.0 - B + beta * beta)


def trapping_indicator(B: float, beta: float) -> float:
    """chi = g^{ab}(dr)_a(dr)_b in the areal-coordinate gate."""
    return B - beta * beta


def misner_sharp_current_residual(dm_dt: float, r: float, T_t_r: float,
                                  G: float = 1.0) -> float:
    return dm_dt - 4.0 * pi * G * r * r * T_t_r
