"""PROVENANCE ONLY. Reduced GEAR-03 S/D lane with the GEAR-63/65 cubic barrier (retired Nov 11 2026).

Retired from live physics by DECISION-map-row-5 (Jeremiah, Nov 11 2026; docs/GEAR_TO_CODE_MAP.md row 5): the GEAR-67
S-D quartic (engine.scalar_system, research settings at lambda > 0) is the one arrest closure. This module is kept,
unchanged in behaviour, so archive comparisons (research/gear_archive_machine_v1.py) and the reduced control runner
(engine/run_controls.py) still reproduce it; no production kernel path imports it. See docs/PHYSICS_LOCK.md.

Why retired: the barrier is conditional on D -> matter conversion, which is not shown (with D's own stress-energy the
same coupling is anti-restoring), and the archive's Term A damping eta_G D^2 Ddot has no derived origin; "always derive"
and "never force a component" exclude both as live physics. eta_G is declared for parameter compatibility only.

Kept on record: the cubic force -(2 pi/3) G alpha_D^2 D^3 is the force of the quartic (pi/6) G alpha_D^2 D^4, i.e. the
GEAR-67 D^4 term with lambda = (4 pi/3) G alpha_D^2 = 4.19 (D*^2 = 2/lambda = 3/(2 pi), D* = 0.690988). It is the only
place the archive fixed lambda from something deeper; it has no S-D cross term (6 S^2 D^2), so matter cannot hold or
pump D through it. The lambda guard below stays: this lane is the lambda = 0 case of the production quartic.
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
