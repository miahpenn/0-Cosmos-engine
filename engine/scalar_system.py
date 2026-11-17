"""V5.5 scalar-field evolution and stress projections.

The scalar equations are the frozen spherical BSSN equations from the
archive. The helper accepts any grid exposing fourth-order cell derivatives
with the parity convention used by the pinned reference kernel.

Hand-floor removal, step 2 (Nov 16 2026; labels in engine/model_constants.py)
----------------------------------------------------------------------------
Production now runs the S-D pair with the GEAR-67 quartic ON (SD_QUARTIC_LAMBDA, a ruler matched today) and the
GEAR-91 coupling to phi ON (BETA_S = 0, BETA_D = BETA_DM), with the hand floor OFF (D_POT = 0):

    V_SD(S, D; phi) = e^(-2 BETA_S phi) S^2/2 - e^(-2 BETA_D phi) D^2/2 + (lambda/8)(S^4 + 6 S^2 D^2 + D^4)

S and D see the phi-dependent masses through dV_SD/dS, dV_SD/dD; phi's equation gains KAPPA dV_SD/dphi (phi at
8 pi G = 1, S and D at G = 1); the stress of S/D carries the full V_SD. The floor is D's vacuum at rest in the sea
(mechanism derived, size set through lambda). The physics is the research reference implementation
research/gear91_coupling_v1.py (potential, gradient, coupled()), moved into the kernel.

State variables: the engine keeps S and D (not the scaled s = sqrt(lambda) S, d = sqrt(lambda) D). lambda scales out
of the S-D equations (PR #486), and at lambda ~ 1.5e6 the unscaled path agrees with the scaled one to 1e-13 and
resolves D's rest value at N = 20 (PR #501, results/floor-swap-2026-11-15/). Keeping S, D leaves every consumer
(initializer, projections, diagnostics, stored runs) reading the same fields; see the PR for the choice.

The old path (hand floor, lambda = 0, no coupling) is kept only as a labelled reference: hand_floor_reference().

Dark matter as S-D matter whose split follows phi (Nov 17 2026; labels in engine/model_constants.py)
------------------------------------------------------------------------------------------------
The kernel's dust is S-D matter: ONE transported density whose coupling to phi is beta_eff(phi) = BETA_DM f_D(phi), with
f_D(phi) DERIVED from GEAR-91 (model_constants.sd_matter_coupled_fraction; D_r's share today READ). Its two parts are the
coupled part f_D rho (D_r plus D's share in the stones) and the uncoupled part (1 - f_D) rho (S's own share):
dark_matter_parts(). They are parts of the same stones, so they are transported together (one velocity); carrying them as
two independently transported fluids would let phi's fifth force pull D's share out of the stones. The COSMOS lane carries
the same physics as two slots with the exchange Gamma (engine/cosmos.py); the sums are identical.
The old imported, fully coupled dust (beta = BETA_DM for all of it) is kept only as the labelled reference
imported_dust_reference(); hand_floor_reference() implies it.
"""
import contextlib
from dataclasses import dataclass
import numpy as np

from . import model_constants as _mc
from .model_constants import BETA_DM, D_POT, KAPPA, V0  # noqa: F401  (one copy: engine/model_constants.py)
from .model_constants import LAMBDA_PHI as LAMBDA

# GEAR-67 on-site quartic of the S-D pair.
#
# GEAR-67 models the pair as two identical elements with an on-site quartic lambda x^4 / 4. In pair
# coordinates (engine labels: S the surviving oscillator, D the ejected mode) the potential at phi's reference is
#     V_SD = 1/2 S^2 - 1/2 D^2 + (lambda/8)(S^4 + 6 S^2 D^2 + D^4)      (g = 1, kappa = 0).
# lambda has no derived value (GEAR-67 3(e)). Production's lambda is a RULER, matched today to the old hand floor
# through the borrowed "today" (engine/model_constants.py); C6 gives lambda* ~ 10-1000: the relocated hierarchy.
SD_QUARTIC_LAMBDA = _mc.SD_QUARTIC_LAMBDA

# GEAR-91 couplings to phi (switches; BETA_D's form derived, value borrowed; BETA_S = 0, S is tangible matter).
BETA_S = _mc.BETA_S
BETA_D = _mc.BETA_D

# Initial structure resolved by engine.v55_initial when the caller does not choose one: production's default is the
# background alone ("no structure placed"); under hand_floor_reference() it is the old reference pocket.
DEFAULT_INITIAL_STRUCTURE = "none"

# Dark matter as S-D matter whose coupled share follows phi (production). False: the imported, fully coupled dust
# (labelled reference only; set by imported_dust_reference() and hand_floor_reference()).
SD_MATTER_SPLIT_FOLLOWS_PHI = True


@contextlib.contextmanager
def hand_floor_reference():
    """REFERENCE ONLY (labelled old path). Restore the pre-step-2 production machine for the block.

    Hand floor on (D_POT = 0.10 of V0, put in by hand), lambda = 0 (S and D interact only through gravity; D has no
    floor of its own and is unbounded), no GEAR-91 coupling, and the old default initial data (the S pocket of
    amplitude 0.01, width 7, D seed 1e-10). Kept so the old path stays reproducible bit for bit; not live physics.
    """
    g = globals()
    keys = ("SD_QUARTIC_LAMBDA", "BETA_S", "BETA_D", "D_POT", "DEFAULT_INITIAL_STRUCTURE", "SD_MATTER_SPLIT_FOLLOWS_PHI")
    saved = {k: g[k] for k in keys}
    g.update(SD_QUARTIC_LAMBDA=0.0, BETA_S=0.0, BETA_D=0.0, D_POT=_mc.HAND_FLOOR_REFERENCE,
             DEFAULT_INITIAL_STRUCTURE="hand_floor_reference", SD_MATTER_SPLIT_FOLLOWS_PHI=False)
    try:
        yield
    finally:
        g.update(saved)


@contextlib.contextmanager
def imported_dust_reference():
    """REFERENCE ONLY (labelled old path). The kernel's dust as the imported, fully coupled dust (beta = BETA_DM for all
    of it), bit-identical to the kernel before the split followed phi. Everything else stays production."""
    g = globals()
    saved = g["SD_MATTER_SPLIT_FOLLOWS_PHI"]
    g["SD_MATTER_SPLIT_FOLLOWS_PHI"] = False
    try:
        yield
    finally:
        g["SD_MATTER_SPLIT_FOLLOWS_PHI"] = saved


def dark_matter_coupling(phi, beta_dm=BETA_DM):
    """The dust's coupling to phi, per cell. Production: beta_eff(phi) = beta_dm f_D(phi) (S-D matter, the derived
    phi-dependent split; an array). Reference (imported_dust_reference / hand_floor_reference): beta_dm itself, the same
    float as before, so the old arithmetic is unchanged."""
    if not SD_MATTER_SPLIT_FOLLOWS_PHI:
        return beta_dm
    return beta_dm * _mc.sd_matter_coupled_fraction(np.asarray(phi, dtype=float))


def dark_matter_parts(phi, rho_dm):
    """(coupled, uncoupled) parts of the kernel's dust: (f_D rho, (1 - f_D) rho) in production; (rho, 0) on the
    reference."""
    rho = np.asarray(rho_dm, dtype=float)
    if not SD_MATTER_SPLIT_FOLLOWS_PHI:
        return rho, np.zeros_like(rho)
    f = _mc.sd_matter_coupled_fraction(np.asarray(phi, dtype=float))
    return f * rho, (1.0 - f) * rho


def coupling_active():
    """True when the GEAR-91 coupling is on (production); False on the reference path or research beta = 0."""
    return BETA_S != 0.0 or BETA_D != 0.0


def sd_rest_value(phi, lam=None, beta_D=None):
    """D's rest value in the sea (S = 0): D_s(phi) = sqrt(2 e^(-2 beta_D phi) / lambda). Held-step readings with the
    coupling on are made against this phi-dependent value. Undefined (inf) at lambda = 0: D is unbounded there."""
    lam = SD_QUARTIC_LAMBDA if lam is None else lam
    bD = BETA_D if beta_D is None else beta_D
    with np.errstate(divide="ignore"):
        return np.sqrt(2.0 * np.exp(-2.0 * bD * np.asarray(phi, dtype=float)) / lam)


def sea_floor_density(phi, lam=None, beta_D=None):
    """The sea's vacuum energy density at D's rest, Einstein-normalised: KAPPA V_SD(0, D_s; phi)
    = -KAPPA e^(-4 beta_D phi) / (2 lambda). Mechanism derived; size set through lambda."""
    lam = SD_QUARTIC_LAMBDA if lam is None else lam
    bD = BETA_D if beta_D is None else beta_D
    return -KAPPA * np.exp(-4.0 * bD * np.asarray(phi, dtype=float)) / (2.0 * lam)


@dataclass(frozen=True)
class ScalarFields:
    S: np.ndarray
    PS: np.ndarray
    D: np.ndarray
    PD: np.ndarray
    phi: np.ndarray
    Pi: np.ndarray

    def copy(self) -> "ScalarFields":
        return ScalarFields(
            self.S.copy(), self.PS.copy(), self.D.copy(),
            self.PD.copy(), self.phi.copy(), self.Pi.copy()
        )


def _d1(grid, x, parity):
    if hasattr(grid, "cell_derivative_fourth"):
        return grid.cell_derivative_fourth(x, parity=parity)
    return grid.d1(x, parity)


def _d2(grid, x, parity):
    if hasattr(grid, "cell_second_derivative_fourth"):
        return grid.cell_second_derivative_fourth(x, parity=parity)
    return grid.d2(x, parity)


def sd_potential(S, D, lam=None, phi=None, beta_S=None, beta_D=None):
    """V_SD(S, D; phi) of GEAR-67 with the GEAR-91 coupling, in engine labels (see SD_QUARTIC_LAMBDA).

    phi=None gives the uncoupled form (phi's reference, e^0 = 1), as before step 2."""
    lam = SD_QUARTIC_LAMBDA if lam is None else lam
    if phi is None:
        return 0.5 * S**2 - 0.5 * D**2 + (lam / 8.0) * (S**4 + 6.0 * S**2 * D**2 + D**4)
    bS = BETA_S if beta_S is None else beta_S
    bD = BETA_D if beta_D is None else beta_D
    return (np.exp(-2.0 * bS * phi) * 0.5 * S**2 - np.exp(-2.0 * bD * phi) * 0.5 * D**2
            + (lam / 8.0) * (S**4 + 6.0 * S**2 * D**2 + D**4))


def sd_potential_gradient(S, D, lam=None):
    """(dV_SD/dS, dV_SD/dD), uncoupled form. For lambda = 0 these are exactly (S, -D)."""
    lam = SD_QUARTIC_LAMBDA if lam is None else lam
    if lam == 0.0:
        return S, -D
    dS = S + (lam / 2.0) * S**3 + (1.5 * lam) * S * D**2
    dD = -D + (lam / 2.0) * D**3 + (1.5 * lam) * S**2 * D
    return dS, dD


def sd_coupled_gradient(S, D, phi, lam=None, beta_S=None, beta_D=None):
    """(dV_SD/dS, dV_SD/dD, dV_SD/dphi) with the GEAR-91 coupling (research/gear91_coupling_v1.gradient)."""
    lam = SD_QUARTIC_LAMBDA if lam is None else lam
    bS = BETA_S if beta_S is None else beta_S
    bD = BETA_D if beta_D is None else beta_D
    eS, eD = np.exp(-2.0 * bS * phi), np.exp(-2.0 * bD * phi)
    dS = eS * S + (lam / 2.0) * S**3 + 1.5 * lam * S * D**2
    dD = -eD * D + (lam / 2.0) * D**3 + 1.5 * lam * S**2 * D
    dphi = -bS * eS * S**2 + bD * eD * D**2
    return dS, dD, dphi


def cosmos_potential(phi):
    """phi's own potential. D_POT = 0 in production (the floor is the sea's, carried by D); 0.10 on the reference."""
    return V0 * (np.exp(-LAMBDA * phi) - D_POT)


def cosmos_potential_prime(phi):
    return -V0 * LAMBDA * np.exp(-LAMBDA * phi)


def scalar_rhs_arrays(grid, geometry, fields, beta_dm=BETA_DM,
                      rho_dm=None):
    """Return d_t(S,PS,D,PD,phi,Pi) for the frozen V5.5 scalar system."""
    r = np.asarray(grid.centers)
    invr = geometry.X**2 / geometry.a
    ap = _d1(grid, geometry.a, 1)
    bp = _d1(grid, geometry.b, 1)
    Xp = _d1(grid, geometry.X, 1)
    alphap = _d1(grid, geometry.alpha, 1)

    def one(f, p, Vp):
        fp = _d1(grid, f, 1)
        fpp = _d2(grid, f, 1)
        lap = invr * (
            fpp - fp * (
                0.5 * ap / geometry.a
                - bp / geometry.b
                + Xp / geometry.X
                - 2.0 / r
            )
        )
        field_t = geometry.alpha * p + geometry.beta * fp
        p_t = (
            geometry.beta * _d1(grid, p, 1)
            + geometry.alpha * lap
            + geometry.alpha * geometry.K * p
            + invr * alphap * fp
            - geometry.alpha * Vp
        )
        return field_t, p_t

    if coupling_active():
        # GEAR-91 (step 2): S, D feel the phi-dependent masses; phi gains KAPPA dV_SD/dphi.
        VpS, VpD, Vphi_SD = sd_coupled_gradient(fields.S, fields.D, fields.phi)
        Vp_phi = cosmos_potential_prime(fields.phi) + KAPPA * Vphi_SD
    else:
        # Reference path / uncoupled research setting: arithmetic unchanged from before step 2.
        VpS, VpD = sd_potential_gradient(fields.S, fields.D)
        Vp_phi = cosmos_potential_prime(fields.phi)
    St, PSt = one(fields.S, fields.PS, VpS)
    Dt, Pdt = one(fields.D, fields.PD, VpD)
    phit, Pit = one(
        fields.phi, fields.Pi, Vp_phi
    )
    if rho_dm is not None:
        # Pi is a normal derivative, so coordinate-time evolution carries alpha.
        Pit = Pit + geometry.alpha * beta_dm * rho_dm
    return ScalarFields(St, PSt, Dt, Pdt, phit, Pit)


def scalar_projection(grid, geometry, fields):
    """Return the archive-normalized scalar Einstein projections."""
    invr = geometry.X**2 / geometry.a
    e = np.zeros_like(fields.S)
    pr = np.zeros_like(fields.S)
    pt = np.zeros_like(fields.S)
    j = np.zeros_like(fields.S)

    # The S-D potential is joint (GEAR-67): its quadratic parts are attributed to S and D as before
    # and the quartic, including the S-D cross term, is added once with D.
    if coupling_active():
        # GEAR-91 (step 2): the full V_SD(S, D; phi), counted once; the floor is D's vacuum in the sea.
        V_S = np.exp(-2.0 * BETA_S * fields.phi) * 0.5 * fields.S**2
        V_D = -np.exp(-2.0 * BETA_D * fields.phi) * 0.5 * fields.D**2 + (SD_QUARTIC_LAMBDA / 8.0) * (
            fields.S**4 + 6.0 * fields.S**2 * fields.D**2 + fields.D**4)
    else:
        # Reference path / uncoupled research setting: arithmetic unchanged from before step 2.
        V_quartic = sd_potential(fields.S, fields.D) - (0.5 * fields.S**2 - 0.5 * fields.D**2)
        if SD_QUARTIC_LAMBDA == 0.0:
            V_quartic = 0.0
        V_S = 0.5 * fields.S**2
        V_D = -0.5 * fields.D**2 + V_quartic
    for f, p, V in (
        (fields.S, fields.PS, V_S),
        (fields.D, fields.PD, V_D),
        (fields.phi, fields.Pi, cosmos_potential(fields.phi)),
    ):
        fp = _d1(grid, f, 1)
        gr2 = invr * fp * fp
        scale = KAPPA if f is fields.phi else 1.0
        e += (0.5 * p * p + 0.5 * gr2 + V) / scale
        pr += (0.5 * p * p + 0.5 * gr2 - V) / scale
        pt += (0.5 * p * p - 0.5 * gr2 - V) / scale
        j += (-p * fp) / scale

    return e, pr, pt, j
