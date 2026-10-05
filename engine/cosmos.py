"""Archive-locked homogeneous COSMOS background lane.

The homogeneous control uses the same beta dark-sector exchange later used
by the unified spacetime:

    phi_ddot + 3 H phi_dot + V_,phi = beta rho_DM
    rho_DM_dot + 3 H rho_DM = -beta rho_DM phi_dot

Radiation has p=rho/3. The archived shear contribution has effective
background pressure p_shear=rho_shear and evolves as a^-6. It is NOT inserted
into the exact spherical local stress tensor.
"""
from dataclasses import dataclass
from math import exp, sqrt


@dataclass(frozen=True)
class CosmosParams:
    lam: float = 1.0
    D: float = 0.10
    beta: float = -0.04
    omega_de_today: float = 0.69
    omega_r_today: float = 9.2e-5
    omega_shear_today: float = 1.745e-8
    f_b: float = 0.157
    n_contract: float = 47.203
    V0: float = 8.242522415500654e-5
    G: float = 1.0


@dataclass
class CosmosState:
    a: float
    H: float
    phi: float
    pi_phi: float
    rho_dm: float
    rho_b: float
    rho_r: float
    rho_shear: float


def potential(phi, p: CosmosParams):
    return p.V0 * (exp(-p.lam * phi) - p.D)


def dV_dphi(phi, p: CosmosParams):
    return -p.V0 * p.lam * exp(-p.lam * phi)


def cosmos_rho_p(y, p: CosmosParams):
    a, H, phi, pi, rho_dm, rho_b, rho_r, rho_s = y
    rho_phi = 0.5 * pi * pi + potential(phi, p)
    p_phi = 0.5 * pi * pi - potential(phi, p)
    rho = rho_phi + rho_dm + rho_b + rho_r + rho_s
    pressure = p_phi + rho_r / 3.0 + rho_s
    return rho, pressure


def construct_present_day(p: CosmosParams = CosmosParams()) -> CosmosState:
    """Return the archive-reconstructed present-day homogeneous state."""
    a = 33.8983
    phi = 0.179055
    phidot = 0.00341913
    rho_dm = 2.5857e-5
    rho_b = 4.0306e-6
    rho_total = 9.64116e-5
    rho_r = p.omega_r_today * rho_total
    rho_shear = p.omega_shear_today * rho_total
    H = sqrt(max(rho_total / 3.0, 0.0))
    return CosmosState(
        a=a,
        H=H,
        phi=phi,
        pi_phi=phidot,
        rho_dm=rho_dm,
        rho_b=rho_b,
        rho_r=rho_r,
        rho_shear=rho_shear,
    )


def state_rho_p(state: CosmosState, p: CosmosParams):
    y = (
        state.a, state.H, state.phi, state.pi_phi,
        state.rho_dm, state.rho_b, state.rho_r, state.rho_shear
    )
    return cosmos_rho_p(y, p)


def normalized_budgets(state: CosmosState, p: CosmosParams = CosmosParams()):
    rho, _ = state_rho_p(state, p)
    if rho == 0.0:
        raise ZeroDivisionError("cannot normalize a zero total density")
    rho_phi = 0.5 * state.pi_phi**2 + potential(state.phi, p)
    return {
        "Omega_phi": rho_phi / rho,
        "Omega_dm": state.rho_dm / rho,
        "Omega_b": state.rho_b / rho,
        "Omega_r": state.rho_r / rho,
        "Omega_shear": state.rho_shear / rho,
    }


def cosmos_rhs(y, p: CosmosParams, exchange_Q=0.0):
    """Homogeneous RHS in one signed-H dynamical system.

    H is an evolved state variable, so H can cross zero from the equations.
    exchange_Q is retained only as a compatibility argument and is ignored;
    no phenomenological reciprocal source enters this lane.
    """
    a, H, phi, pi, rho_dm, rho_b, rho_r, rho_s = y
    rho, pressure = cosmos_rho_p(y, p)

    da = H * a
    dH = -0.5 * (rho + pressure)
    dphi = pi
    dpi = -3.0 * H * pi - dV_dphi(phi, p) + p.beta * rho_dm
    drho_dm = -3.0 * H * rho_dm - p.beta * rho_dm * pi
    drho_b = -3.0 * H * rho_b
    drho_r = -4.0 * H * rho_r
    drho_s = -6.0 * H * rho_s

    return (
        da, dH, dphi, dpi,
        drho_dm, drho_b, drho_r, drho_s
    )
