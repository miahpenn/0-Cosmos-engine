"""Archive-locked homogeneous COSMOS background lane.

Densities in this lane are Einstein-normalized, rho_E = 8*pi*G*rho_physical,
with G fixed to 1 and c=1. Consequently the background equation is
H**2 = rho_E/3. In the local production lane the matching archive-normalized
fluid density is divided by KAPPA=8*pi at the stress-projection boundary,
and the BSSN constraints use the corresponding 16*pi source coefficient.
See docs/UNITS_CONVENTION.md and its production-path Friedmann regression.

The homogeneous control uses the same beta dark-sector exchange later used
by the unified spacetime:

    phi_ddot + 3 H phi_dot + V_,phi = beta rho_DM
    rho_DM_dot + 3 H rho_DM = -beta rho_DM phi_dot

Radiation has p=rho/3. The archived shear contribution has effective
background pressure p_shear=rho_shear and evolves as a^-6. It is NOT inserted
into the exact spherical local stress tensor.

Floor (hand-floor removal step 2, Nov 16 2026; labels in engine/model_constants.py).
The hand floor D = 0.10 of V0 is off. The floor is the S-D sea's vacuum, D at
rest, felt through phi by the GEAR-91 coupling, carried adiabatically (no D
field in this lane; D stays at its phi-dependent rest, research PRs #500-#501):

    rho_floor(phi) = -KAPPA e^(-4 beta_D phi) / (2 lambda),   p_floor = -rho_floor,
    phi's equation gains -d rho_floor/dphi.

It enters as a term of phi's potential, so it is counted once in rho, p and
phi's force. Floor: mechanism derived, size set through lambda (a ruler matched
today through the borrowed V0 and phi_today; C6 gives lambda* ~ 10-1000, the
relocated hierarchy). beta_D: form derived, value not. All times read from this
lane are chart-dependent through phi_today.

The old path (hand floor, no coupled floor) is kept only as the labelled
reference HAND_FLOOR_REFERENCE_PARAMS, bit-identical to the lane before step 2.
"""
from dataclasses import dataclass
from math import exp, isfinite, log, sqrt

from . import model_constants as _mc


@dataclass(frozen=True)
class CosmosParams:
    lam: float = _mc.LAMBDA_PHI
    D: float = _mc.D_POT                     # hand floor: 0 in production; 0.10 on the reference only
    beta: float = _mc.BETA_DM
    # Adiabatic coupled floor of the S-D sea (step 2). coupled_floor=False is the reference path.
    coupled_floor: bool = True
    sd_lambda: float = _mc.SD_QUARTIC_LAMBDA   # ruler, matched today (not derived)
    beta_D: float = _mc.BETA_D                 # form derived, value not
    omega_de_today: float = 0.69
    omega_r_today: float = 9.2e-5
    omega_shear_today: float = 1.745e-8
    f_b: float = 0.157
    n_contract: float = 47.203
    V0: float = _mc.V0
    # Legacy metadata kept for archived parameter payload compatibility.
    # G is absorbed into rho_E=8*pi*G*rho_physical; this implementation fixes
    # G=1 rather than supporting an independent gravitational-coupling sweep.
    G: float = 1.0

    def __post_init__(self):
        if self.coupled_floor and (
            isinstance(self.sd_lambda, bool)
            or not isfinite(float(self.sd_lambda))
            or float(self.sd_lambda) <= 0.0
        ):
            raise ValueError(
                "the coupled floor needs lambda > 0 (at lambda = 0 D is unbounded and has no floor); "
                "use HAND_FLOOR_REFERENCE_PARAMS for the old path"
            )
        if (
            isinstance(self.G, bool)
            or not isinstance(self.G, (int, float))
            or not isfinite(float(self.G))
            or float(self.G) != 1.0
        ):
            raise ValueError(
                "COSMOS uses Einstein-normalized densities rho_E=8*pi*G*rho_physical "
                "with G fixed to 1; varying G requires a consistent density "
                "renormalization and is unsupported"
            )


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


# REFERENCE ONLY: the pre-step-2 lane (hand floor 0.10 of V0, no coupled floor), bit-identical to it.
HAND_FLOOR_REFERENCE_PARAMS = CosmosParams(D=_mc.HAND_FLOOR_REFERENCE, coupled_floor=False)


def floor_density(phi, p: CosmosParams):
    """rho_floor(phi) = -KAPPA e^(-4 beta_D phi) / (2 lambda): the sea's vacuum at D's rest (0 on the reference)."""
    if not p.coupled_floor:
        return 0.0
    return -_mc.KAPPA * exp(-4.0 * p.beta_D * phi) / (2.0 * p.sd_lambda)


def potential(phi, p: CosmosParams):
    if not p.coupled_floor:
        return p.V0 * (exp(-p.lam * phi) - p.D)
    return p.V0 * (exp(-p.lam * phi) - p.D) + floor_density(phi, p)


def dV_dphi(phi, p: CosmosParams):
    if not p.coupled_floor:
        return -p.V0 * p.lam * exp(-p.lam * phi)
    # d rho_floor / dphi = -4 beta_D rho_floor: phi's extra force (GEAR-91 at the sea's rest).
    return -p.V0 * p.lam * exp(-p.lam * phi) - 4.0 * p.beta_D * floor_density(phi, p)


def potential_node(p: CosmosParams):
    """phi where phi's potential (with the floor) vanishes: ln(1/D)/lam on the reference; with the coupled floor and
    D = 0, V0 e^(-lam phi) = KAPPA e^(-4 beta_D phi)/(2 lambda), i.e. phi = ln(2 lambda V0 / KAPPA)/(lam - 4 beta_D)."""
    if not p.coupled_floor:
        return log(1.0 / p.D) / p.lam
    if p.D != 0.0:
        raise ValueError("the node has no closed form with both the hand floor and the coupled floor on")
    return log(2.0 * p.sd_lambda * p.V0 / _mc.KAPPA) / (p.lam - 4.0 * p.beta_D)


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
