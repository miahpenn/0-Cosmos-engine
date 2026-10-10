"""Evolved-state curvature consistency audit.

Compares the pinned vendor BSSN curvature with the repository's independent
general-BSSN ricci_terms() on identical geometry arrays. Both use evolved
Lambda in the Ricci principal part. No polar-areal identity is used.
The initial-slice agreement is a gate; failure stops before evolution.
Measurement only: no production equations or trajectory settings are changed.
"""
import math
import pathlib
import sys
from types import SimpleNamespace

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import reference_pirk_unified as ref
from engine import v55_pirk_adapter as adapter
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

_, vacuum, _ = adapter.vendor_modules()

N = 40
R_MAX = 40.0
AMPLITUDE = 0.01
D_AMPLITUDE = 1.0e-10
CFL = 0.03
TARGET_TIMES = (1.0, 2.0, 3.0)
INITIAL_AGREEMENT_TOL = 1.0e-12


def reference_state(geometry):
    return SimpleNamespace(
        a=np.asarray(geometry.a, dtype=float),
        b=np.asarray(geometry.b, dtype=float),
        X=np.asarray(geometry.X, dtype=float),
        alpha=np.asarray(geometry.alpha, dtype=float),
        Lambda=np.asarray(geometry.Lambda, dtype=float),
    )


def compare_derivative_operators(state, reference_grid):
    """Compare every first/second derivative used by the curvature formula."""
    grid = state.grid
    geometry = state.geometry
    results = {}
    first_order = {
        "a": (geometry.a, +1),
        "b": (geometry.b, +1),
        "X": (geometry.X, +1),
        "Lambda": (geometry.Lambda, -1),
    }
    second_order = {
        "a": (geometry.a, +1),
        "b": (geometry.b, +1),
        "X": (geometry.X, +1),
    }
    # Alpha derivatives do not enter R itself, but are checked as a useful
    # companion because the same vendor geometry_terms routine also returns
    # lapse-Hessian terms used by the geometry RHS.
    first_order["alpha"] = (geometry.alpha, +1)
    second_order["alpha"] = (geometry.alpha, +1)

    for name, (values, parity) in first_order.items():
        vendor = np.asarray(grid.cell_derivative_fourth(values, parity=parity))
        reference = np.asarray(ref.D(reference_grid, values, parity))
        results[f"D1:{name}"] = vendor - reference

    for name, (values, parity) in second_order.items():
        vendor = np.asarray(grid.cell_second_derivative_fourth(values, parity=parity))
        reference = np.asarray(ref.D2(reference_grid, values, parity))
        results[f"D2:{name}"] = vendor - reference

    return results


def hamiltonian_residual(state):
    raw = vacuum.constraints(state.grid, state.geometry)
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    return np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(total.rho)


def report(state, reference_grid, label, initial_gate=False):
    geometry = state.geometry
    radius = np.asarray(state.grid.centers)
    vendor_R = np.asarray(vacuum.geometry_terms(state.grid, geometry)["R"])
    reference_R = np.asarray(ref.ricci_terms(reference_grid, reference_state(geometry))[0])
    curvature_diff = vendor_R - reference_R
    derivative_diffs = compare_derivative_operators(state, reference_grid)
    operator_maxima = {
        name: float(np.max(np.abs(values)))
        for name, values in derivative_diffs.items()
    }
    operator_max = max(operator_maxima.values(), default=0.0)
    center_curvature_max = float(np.max(np.abs(curvature_diff[:5])))
    all_curvature_max = float(np.max(np.abs(curvature_diff)))
    idx = int(np.argmax(np.abs(curvature_diff)))

    H = hamiltonian_residual(state)
    areal_departure = np.asarray(geometry.b) / np.asarray(geometry.X) ** 2 - 1.0

    print(f"[{label}] t={state.t:.12g}")
    print("  derivative operator max-absolute differences:")
    for name in sorted(operator_maxima):
        d = derivative_diffs[name]
        print(
            f"    {name}: all={operator_maxima[name]:.3e}; "
            f"cells0-4={np.max(np.abs(d[:5])):.3e}"
        )
    print(f"  curvature R vendor cells0-4 = {np.array2string(vendor_R[:5], precision=9)}")
    print(f"  curvature R ref    cells0-4 = {np.array2string(reference_R[:5], precision=9)}")
    print(f"  signed R difference cells0-4 = {np.array2string(curvature_diff[:5], precision=4)}")
    print(
        f"  |delta R| max: all={all_curvature_max:.3e} at cell {idx} "
        f"(r={radius[idx]:.6g}); cells0-4={center_curvature_max:.3e}; "
        f"cells>=5={np.max(np.abs(curvature_diff[5:])):.3e}"
    )
    print(
        f"  b/X^2 - 1: cells0-4={np.array2string(areal_departure[:5], precision=4)}; "
        f"max|.| all={np.max(np.abs(areal_departure)):.3e}"
    )
    print(
        f"  Hamiltonian residual: cells0-4={np.array2string(H[:5], precision=4)}; "
        f"max|H| all={np.max(np.abs(H)):.3e}"
    )

    if initial_gate:
        if operator_max > INITIAL_AGREEMENT_TOL:
            raise AssertionError(
                f"initial derivative-operator agreement failed: {operator_max:.3e} "
                f"> {INITIAL_AGREEMENT_TOL:.1e}"
            )
        if all_curvature_max > INITIAL_AGREEMENT_TOL:
            raise AssertionError(
                f"initial curvature agreement failed: {all_curvature_max:.3e} "
                f"> {INITIAL_AGREEMENT_TOL:.1e}"
            )
        print(
            f"  INITIAL_GATE=PASS (operator max {operator_max:.3e}, "
            f"curvature max {all_curvature_max:.3e}, tolerance "
            f"{INITIAL_AGREEMENT_TOL:.1e})"
        )
    return all_curvature_max, operator_max


if __name__ == "__main__":
    kernel = V55ProductionKernel()
    state, _, _ = discrete_consistent_state(
        kernel,
        resolution=N,
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=7.0,
        D_amplitude=D_AMPLITUDE,
        include_radiation=True,
    )
    reference_grid = ref.Grid(N, R_MAX)
    report(state, reference_grid, "initial slice", initial_gate=True)

    dt = CFL * state.grid.dr
    for target in TARGET_TIMES:
        while state.t < target - 1.0e-12:
            state = kernel.step(state, min(dt, target - state.t))
        report(state, reference_grid, "evolved slice")
