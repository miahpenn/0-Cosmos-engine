"""COSMOS lane: corrected archive-derived homogeneous background skeleton."""
from dataclasses import dataclass
from math import exp

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
    V0: float = 8.2425e-5
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

def potential(phi,p):
    return p.V0*(exp(-p.lam*phi)-p.D)

def dV_dphi(phi,p):
    return -p.V0*p.lam*exp(-p.lam*phi)

def cosmos_rho_p(y,p):
    a,H,phi,pi,rho_dm,rho_b,rho_r,rho_s=y
    rho_phi=0.5*pi*pi+potential(phi,p)
    p_phi=0.5*pi*pi-potential(phi,p)
    return rho_phi+rho_dm+rho_b+rho_r+rho_s, p_phi+rho_r/3.0+rho_s/3.0

def cosmos_rhs(y,p,exchange_Q=0.0):
    a,H,phi,pi,rho_dm,rho_b,rho_r,rho_s=y
    rho,p_tot=cosmos_rho_p(y,p)
    return (a*H,-0.5*(rho+p_tot),pi,-3.0*H*pi-dV_dphi(phi,p),-3.0*H*rho_dm,-3.0*H*rho_b,-4.0*H*rho_r,-6.0*H*rho_s)
