"""RESEARCH ONLY. phi's links to S and D from the machine's own coupling (GEAR-91, research switch of PR #465), checked
against the unpaid items (Nov 15 2026: floor, seed, second reset, kappa*/lambda*, dust). No runs; algebra checked
numerically on the switch's own potential(). Note: results/phi-links-2026-11-15/README.md.

GEAR-91:  V_SD = e^(-2 bS phi) S^2/2 - e^(-2 bD phi) D^2/2 + (lam/8)(S^4 + 6 S^2 D^2 + D^4),  Box phi = V_phi' + 8 pi dV_SD/dphi
(phi at 8 pi G = 1, S and D at G = 1; bS = 0 default, bD = -0.04 locked with dark matter).
Usage: python research/phi_links_v1.py
"""
import math
import sys

sys.path.insert(0, ".")
import numpy as np  # noqa: E402

from engine.model_constants import D_POT, V0  # noqa: E402
from research.gear91_coupling_v1 import BETA, gradient, potential  # noqa: E402

bD = BETA


def sea(phi, lam):
    """D's rest value and vacuum energy in the sea (S = 0) at given phi: D^2 = 2 e^(-2 bD phi)/lam, V = -e^(-4 bD phi)/(2 lam)."""
    return 2 * np.exp(-2 * bD * phi) / lam, -np.exp(-4 * bD * phi) / (2 * lam)


if __name__ == "__main__":
    lam = 5000.0
    print("1. D's vacuum as a function of phi (checked against potential()):")
    for phi in (-1.0, 0.0, 1.0, 2.0, 3.0):
        D2, Vv = sea(phi, lam)
        Dg = np.linspace(0.5, 1.5, 200001) * math.sqrt(D2)
        Vn = potential(0.0, Dg, phi, beta_S=0.0, beta_D=bD, lam=lam)
        print(f"   phi {phi:+.1f}: D_s^2 = {D2:.6e} (numeric min at {Dg[np.argmin(Vn)]**2:.6e}),  V_vac = {Vv:.6e} (numeric {Vn.min():.6e})")
    print("\n2. The floor from D's step, felt by phi:  V_eff(phi) = V0 e^(-phi) + 8 pi V_vac(phi)  (hand floor off)")
    lam_f = 4 * math.pi / (D_POT * V0)
    print(f"   matching the hand floor at phi = 0 (8 pi / (2 lam) = D_POT V0): lam = 4 pi/(D_POT V0) = {lam_f:.4e}")
    node_hand = math.log(1 / D_POT); node_D = math.log(1 / D_POT) / (1 - 4 * bD)
    print(f"   phi's node (V_eff = 0): hand floor ln(1/D_POT) = {node_hand:.4f};  D's floor ln(1/D_POT)/(1 - 4 bD) = {node_D:.4f}")
    phis = np.linspace(-2, 6, 8001)
    Veff = V0 * np.exp(-phis) - 8 * math.pi * np.exp(-4 * bD * phis) / (2 * lam_f)
    dV = np.gradient(Veff, phis)
    print(f"   numeric node {phis[np.argmin(np.abs(Veff))]:.4f};  V_eff' < 0 everywhere on [-2, 6]: {bool(np.all(dV < 0))} (no minimum: phi does not park)")
    print(f"   floor slope / e^(-phi) slope at the node: {(-4*bD)*D_POT*math.exp(-4*bD*node_D)/math.exp(-node_D):.3f}")
    print("\n3. Seed: is D = 0 invariant whatever phi does?  dV/dD at D = 0, random S, phi:")
    rng = np.random.default_rng(0); S = rng.normal(size=1000); ph = rng.normal(size=1000) * 3
    print(f"   max |dV/dD(S, 0, phi)| = {np.max(np.abs(gradient(S, 0.0, ph, beta_S=0.0, beta_D=bD, lam=lam)[1])):.1e}  (D -> -D symmetry: phi cannot seed D)")
    print("\n4. Second reset: D's mass^2 at the saddle (S = 0, D = 0) is -e^(-2 bD phi):")
    print(f"   over phi in [-50, 50]: max = {np.max(-np.exp(-2*bD*np.linspace(-50,50,1001))):.3e} < 0 always. phi cannot make the saddle stable.")
    print("   The saddle is stable for D only if (3 lam/2) S^2 > e^(-2 bD phi): it needs S's energy. With bS = 0, phi has no")
    print("   channel into S; with bS != 0 phi changes S's mass only at phi's rate (~H << S's mass): adiabatic, no S made.")
    print("\n5. kappa: GEAR-03 split K_D - K_S over two, at bS = 0: kappa(phi) = (1 - e^(-2 bD phi))/2 ~ bD phi near phi = 0.")
    for phi in (0.0, 0.5, 1.0, node_D):
        print(f"   phi {phi:.3f}: kappa = {(1-np.exp(-2*bD*phi))/2:+.4f}")
