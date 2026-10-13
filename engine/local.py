"""Local GEAR-03 S/D sector with the GEAR-63/65 geometric barrier."""
from dataclasses import dataclass

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
