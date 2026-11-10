"""RESEARCH ONLY. Contraction depth of the production COSMOS lane for different physical endpoints (GEAR-93 §2).

The turnaround has H = 0, so in this flat lane rho_total = 3 H^2 = 0 there: a depth cannot be counted from the
turnaround density by w = 1 scaling. It is measured by integrating (research/cosmos_efold_v1.py) until
rho_total / rho_today reaches the endpoint ratio. |H| = mu (GEAR-93's matter endpoint) is rho_total = 3 mu^2, i.e.
ratio (mu / H_today)^2. Usage (repo root): python research/cosmos_endpoint_depths_v1.py
"""
import sys

sys.path.insert(0, ".")

from engine.cosmos import CosmosParams, construct_present_day  # noqa: E402
from research.cosmos_efold_v1 import PLANCK_TO_CRITICAL_DENSITY_RATIO, integrate_cycle  # noqa: E402

HBAR_H0_EV = 1.054571817e-34 * 2.20e-18 / 1.602176634e-19   # hbar H0 with H0 = 2.20e-18 s^-1 (CORRECTION-GEAR-50)


def main():
    s = construct_present_day(CosmosParams())
    cases = [
        ("|H| = mu, production units (mu = 1, H_today = %.4e)" % s.H, 1.0 / s.H**2),
        ("|H| = mu, mu = m_proton, observed H0", (938.272e6 / HBAR_H0_EV) ** 2),
        ("|H| = mu, mu = m_electron, observed H0", (0.51099895e6 / HBAR_H0_EV) ** 2),
        ("Planck density (full; CORRECTION-GEAR-50 ratio)", PLANCK_TO_CRITICAL_DENSITY_RATIO),
    ]
    print("endpoint | rho_end/rho_today | e-folds below turnaround | e-folds below today")
    for name, ratio in cases:
        tr = integrate_cycle(density_ratio=ratio)
        print(f"{name} | {ratio:.3e} | {-tr.N_turn[-1]:.4f} | {-tr.ln_a_today[-1]:.4f}   ({tr.status})")


if __name__ == "__main__":
    main()
