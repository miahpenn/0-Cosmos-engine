"""Explicit local/COSMOS interface quantities.

No fitted feedback law or artificial source closure is contained here.
"""
from dataclasses import dataclass

@dataclass
class InterfaceState:
    tau: float = 0.0
    radius: float = 30.0
    radius_dot: float = 0.0
    local_energy: float = 0.0
    outward_flux: float = 0.0
    stress_energy_flux: float = 0.0
    ledger_exchange: float = 0.0

def exchange_from_flux(flux: float, donor_energy_rate=None):
    recipient = flux
    donor = -flux if donor_energy_rate is None else donor_energy_rate
    return {"recipient": recipient, "donor": donor, "net": recipient + donor}
