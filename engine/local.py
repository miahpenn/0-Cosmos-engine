"""Local GEAR-03 S/D sector with the GEAR-63/65 geometric barrier.

The cubic barrier -(2 pi/3) G alpha_D^2 D^3 is derived from a closure and is
conditional: conformal force -(R/6) D with dust-remnant curvature
R = 8 pi G (1/2 alpha_D^2 D^2) (GEAR-67E converted branch). It assumes D -> matter
conversion, which is open; with D's own stress-energy the same coupling is
anti-restoring, +(4 pi/3) G alpha_D^2 D^3. It is applied at all times, so the
post-conversion force also acts on the pre-conversion mode. See docs/PHYSICS_LOCK.md.
eta_G is declared for parameter compatibility and is not used here.

GEAR-67 S-D quartic: this reduced lane has no (lambda/8)(S^4 + 6 S^2 D^2 + D^4) term; it is the
lambda = 0 case, like the production engine (engine.scalar_system.SD_QUARTIC_LAMBDA). Its own cubic
barrier above is a different, conditional closure. local_rhs therefore refuses to run when
SD_QUARTIC_LAMBDA is nonzero, so the two lanes cannot silently disagree; a lambda > 0 reduced model
would need its own derivation.
"""
from dataclasses import dataclass

from . import scalar_system as _scalar_system

@dataclass(frozen=True)
class LocalParams:
    kappa: float = 0.0
    g: float = 1.0
    eta_G: float = 0.05
    G: float = 1.0
    alpha_D: float = 1.0

@dataclass
class LocalState:
    S: float
    Sd: float
    D: float
    Dd: float

def local_rhs(y, p: LocalParams, H: float = 0.0):
    if _scalar_system.SD_QUARTIC_LAMBDA != 0.0:
        raise NotImplementedError(
            "the reduced local lane is the GEAR-67 lambda = 0 case; "
            "SD_QUARTIC_LAMBDA is nonzero (research setting)"
        )
    S, Sd, D, Dd = y
    omega2 = p.kappa + p.g
    alpha2 = p.g - p.kappa
    geom = (2.0*3.141592653589793/3.0) * p.G * p.alpha_D*p.alpha_D * D**3
    Sdd = -omega2*S - 3.0*H*Sd
    Ddd = alpha2*D - 3.0*H*Dd - geom
    return (Sd, Sdd, Dd, Ddd)

def local_energy(y):
    S,Sd,D,Dd=y
    return 0.5*Sd*Sd + 0.5*S*S + 0.5*Dd*Dd - 0.5*D*D
