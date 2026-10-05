"""Cycle coordinates and event observation; no imposed reset rule."""
from dataclasses import dataclass

@dataclass
class CycleObservation:
    theta: float=None
    phi: float=None
    lambda_cycle: float=None
    turnaround: bool=False
    handoff: bool=False

def observe(H,H_prev,phi=None,theta=None):
    turnaround=H_prev is not None and H_prev>0.0 and H<=0.0
    return CycleObservation(theta=theta,phi=phi,turnaround=turnaround)
