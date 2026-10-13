"""Canonical dynamical handoff and cycle ledger.

All quantities are observations of the solved state.  No reset, bounce,
branch flip, fitted source, or clock conversion is imposed.
"""
from dataclasses import dataclass
from math import atan2, isfinite

@dataclass(frozen=True)
class HandoffPoint:
    t: float
    tau: float
    H_eff: float
    phi: float
    theta_D: float
    theta_S: float
    radius: float
    misner_sharp_mass: float
    chi: float
    flux_T: float
    work_pR: float

@dataclass
class CycleEvent:
    index: int
    t: float
    tau: float
    kind: str
    H_before: float
    H_after: float

@dataclass
class CycleLedger:
    events: list

    def __init__(self):
        self.events = []

    def observe(self, index, t, tau, H_before, H_after):
        if not (isfinite(H_before) and isfinite(H_after)):
            return None
        kind = None
        if H_before > 0.0 and H_after <= 0.0:
            kind = "turnaround"
        elif H_before < 0.0 and H_after >= 0.0:
            kind = "re-expansion_crossing"
        if kind is None:
            return None
        event = CycleEvent(index, t, tau, kind, H_before, H_after)
        self.events.append(event)
        return event

def local_phase(S, Sd):
    return atan2(Sd, S)

def ejected_phase(D, Dd):
    return atan2(Dd, D)

def handoff_from_ledger(t, tau, observables, S, Sd, D, Dd):
    return HandoffPoint(
        t=float(t),
        tau=float(tau),
        H_eff=float(observables["H_eff"]),
        phi=float(observables["phi_outer"]),
        theta_D=ejected_phase(D, Dd),
        theta_S=local_phase(S, Sd),
        radius=float(observables["R_sigma"]),
        misner_sharp_mass=float(observables["M_MS"]),
        chi=float(observables["chi_sigma"]),
        flux_T=float(observables["flux_T"]),
        work_pR=float(observables["work_pR"]),
    )
