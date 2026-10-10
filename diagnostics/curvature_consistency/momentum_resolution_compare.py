"""Spatial-refinement check for the evolved, matter-sourced momentum residual.

Runs the same discrete-consistent S/D candidate at N=40 and N=80 to t=1.
At each resolution, measures signed one-step changes to the momentum and
Hamiltonian constraints at several timestep factors, separates projection
jumps, and linearly extrapolates rates against dt to estimate their dt->0
limits. This is diagnostic only; no production equations or parameters change.
"""
import copy
import pathlib
import sys
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

from engine.discrete_consistent_initial import discrete_consistent_state
from engine.production_kernel import V55ProductionKernel
from engine import v55_pirk_adapter as adapter
from evolved_curvature_compare import hamiltonian_residual, momentum_components

_, vacuum, _ = adapter.vendor_modules()

R_MAX = 40.0
TARGET_T = 1.0
AMPLITUDE = 0.01
D_AMPLITUDE = 1.0e-10
EVOLUTION_CFL = 0.015
PROBE_CFL_FACTORS = (0.015, 0.0075, 0.00375)


def maxabs(values):
    return float(np.max(np.abs(np.asarray(values))))


def run_resolution(n):
    kernel = V55ProductionKernel()
    state, _, _ = discrete_consistent_state(
        kernel,
        resolution=n,
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=7.0,
        D_amplitude=D_AMPLITUDE,
        include_radiation=True,
    )
    dr = state.grid.dr
    dt_evolve = EVOLUTION_CFL * dr
    while state.t < TARGET_T - 1.0e-12:
        state = kernel.step(state, min(dt_evolve, TARGET_T - state.t))

    # Record the natural evolved state, then identify any state-map offset
    # before taking the one-step probes.
    H_raw = hamiltonian_residual(state)
    M_raw = momentum_components(state)["M_total"]
    projected = copy.deepcopy(state)
    projected.geometry = vacuum.enforce_algebraic_regularity(
        projected.grid, projected.geometry.copy()
    )
    H_base = hamiltonian_residual(projected)
    M_base = momentum_components(projected)["M_total"]
    print(f"[RESOLUTION N={n}] dr={dr:.9g} t={state.t:.12g} evolution_dt/dr={EVOLUTION_CFL}")
    print(
        f"  pre-projection center H={H_raw[0]:+.6e}, M={M_raw[0]:+.6e}; "
        f"post-projection center H={H_base[0]:+.6e}, M={M_base[0]:+.6e}"
    )
    print(
        f"  projection jump: cell0 dH={H_base[0]-H_raw[0]:+.3e}, "
        f"dM={M_base[0]-M_raw[0]:+.3e}; "
        f"maxabs dM cells0-4={maxabs(M_base[:5]-M_raw[:5]):.3e}"
    )
    M_terms_base = momentum_components(projected)
    rates_m = []
    rates_h = []
    probe_dt_values = []

    for factor in PROBE_CFL_FACTORS:
        dt = factor * dr
        advanced = kernel.step(copy.deepcopy(projected), dt)
        H_after = hamiltonian_residual(advanced)
        M_terms_after = momentum_components(advanced)
        delta_M = M_terms_after["M_total"] - M_terms_base["M_total"]
        delta_H = H_after - H_base
        rate_M = delta_M / dt
        rate_H = delta_H / dt
        rates_m.append(float(rate_M[0]))
        rates_h.append(float(rate_H[0]))
        probe_dt_values.append(dt)
        print(
            f"  probe dt/dr={factor:g} dt={dt:.9g}: "
            f"dM0={delta_M[0]:+.6e} rateM0={rate_M[0]:+.6e}; "
            f"dH0={delta_H[0]:+.6e} rateH0={rate_H[0]:+.6e}"
        )
        for name in ("M_DAa", "M_DK", "M_chi", "M_radial_Aa",
                     "M_S_source", "M_D_source", "M_phi_source",
                     "M_fluid_source"):
            dterm = M_terms_after[name] - M_terms_base[name]
            print(
                f"    d{name}: cell0={dterm[0]:+.6e}; "
                f"rate0={dterm[0]/dt:+.6e}; "
                f"maxabs rate cells0-4={maxabs(dterm[:5])/dt:.3e}"
            )
        closure = (
            sum(
                M_terms_after[name] - M_terms_base[name]
                for name in ("M_DAa", "M_DK", "M_chi", "M_radial_Aa",
                             "M_S_source", "M_D_source", "M_phi_source",
                             "M_fluid_source")
            )
            - delta_M
        )
        print(f"    momentum-budget delta closure maxabs={maxabs(closure):.3e}")

    # A simple first-order-in-dt diagnostic extrapolation. Values are displayed
    # at all dt values too; this intercept is not a fitted physics parameter.
    m_intercept = float(np.polyfit(np.asarray(probe_dt_values), np.asarray(rates_m), 1)[1])
    h_intercept = float(np.polyfit(np.asarray(probe_dt_values), np.asarray(rates_h), 1)[1])
    print(
        f"  dt->0 linear rate estimate: Mdot_center={m_intercept:+.6e}; "
        f"Hdot_center={h_intercept:+.6e}"
    )
    print(
        f"  current constraint norms: maxabs M all={maxabs(M_base):.6e}, "
        f"cells0-4={maxabs(M_base[:5]):.6e}; "
        f"maxabs H all={maxabs(H_base):.6e}, cells0-4={maxabs(H_base[:5]):.6e}"
    )
    return {
        "n": n,
        "dr": dr,
        "Mdot0": m_intercept,
        "Hdot0": h_intercept,
        "M0": float(M_base[0]),
        "H0": float(H_base[0]),
        "Mnorm0_4": maxabs(M_base[:5]),
        "Hnorm0_4": maxabs(H_base[:5]),
    }


if __name__ == "__main__":
    results = [run_resolution(n) for n in (40, 80, 160)]
    print("=== SPATIAL REFINEMENT SUMMARY ===")
    for row in results:
        print(
            f"N={row['n']} dr={row['dr']:.9g} M0={row['M0']:+.6e} "
            f"H0={row['H0']:+.6e} Mdot0={row['Mdot0']:+.6e} "
            f"Hdot0={row['Hdot0']:+.6e} "
            f"maxM(cells0-4)={row['Mnorm0_4']:.6e} "
            f"maxH(cells0-4)={row['Hnorm0_4']:.6e}"
        )
    if abs(results[0]["Mdot0"]) > 0.0 and abs(results[1]["Mdot0"]) > 0.0:
        print(
            "  estimated |Mdot| refinement ratio N80/N40="
            f"{abs(results[1]['Mdot0']/results[0]['Mdot0']):.6g}"
        )
    if abs(results[0]["Hdot0"]) > 0.0 and abs(results[1]["Hdot0"]) > 0.0:
        print(
            "  estimated |Hdot| refinement ratio N80/N40="
            f"{abs(results[1]['Hdot0']/results[0]['Hdot0']):.6g}"
        )
