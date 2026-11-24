"""Compare production and direct areal curvature on accepted W1 slices."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from engine.stress_energy import assemble_total_stress_energy
from research.w1_fring_smoke.origin_trace import SAMPLES, run

WIDTHS = (7, 9)
TRACE_DIR = Path("results/research/w1_fring_origin_trace")


def derivative(grid, values, width, order):
    # Same shifted polynomial stencil as the frozen areal diagnostic.
    half = width // 2
    padded = np.r_[values[:half][::-1], values]
    cells = np.arange(len(values)) + half
    starts = np.minimum(cells - half, len(padded) - width)
    offsets = starts - cells
    result = np.empty_like(values)
    for offset in np.unique(offsets):
        positions = offset + np.arange(width)
        moments = np.zeros(width)
        moments[order] = 1 if order == 1 else 2
        weights = np.linalg.solve(np.vander(positions, increasing=True).T, moments)
        selected = offsets == offset
        samples = padded[starts[selected, None] + np.arange(width)]
        result[selected] = (samples - values[selected, None]) @ weights / grid.dr**order
    return result


def areal_scalar(grid, q, u, width):
    """Scalar curvature for inverse radial metric 1+q and areal radius r(1+u)."""
    r = grid.centers
    radius = r * (1 + u)
    u_prime = derivative(grid, u, width, 1)
    q_prime = derivative(grid, q, width, 1)
    v = u + r * u_prime  # R' - 1
    radius_second = 2 * u_prime + r * derivative(grid, u, width, 2)
    # 1 - (1+q)R'^2 without subtracting two nearly unit quantities.
    angular = -q - (1 + q) * v * (2 + v)
    return 2 * angular / radius**2 - (
        4 * (1 + q) * radius_second + 2 * q_prime * (1 + v)
    ) / radius


def manufacture(grid):
    r = grid.centers
    zero = np.zeros_like(r)
    k = 1e-5
    u = 1e-7 * r**2
    v = 3 * u
    q = -v * (2 + v) / (1 + v)**2
    for width in WIDTHS:
        np.testing.assert_array_equal(areal_scalar(grid, zero, zero, width), zero)
        np.testing.assert_allclose(areal_scalar(grid, -k*r**2, zero, width), 6*k, atol=2e-14, rtol=0)
        np.testing.assert_allclose(areal_scalar(grid, q, u, width), zero, atol=2e-14, rtol=0)


def observe(state, vacuum, reference):
    grid = state.grid
    metric = state.geometry.copy()
    if state.electromagnetic is not None:
        raise ValueError("expected the uncharged W1 replay")
    if state.t == 0:
        manufacture(grid)
    vacuum_constraints = vacuum.constraints(grid, metric)
    production = vacuum.geometry_terms(grid, metric)["R"]
    rho = assemble_total_stress_energy(grid, metric, state.scalars, state.matter).rho
    source = 16 * np.pi * rho
    kinematic = (-metric.Aa**2 - 2 * (-0.5 * metric.Aa)**2) + (2/3) * metric.K**2
    official = vacuum_constraints["hamiltonian"] - source
    np.testing.assert_array_equal(official[SAMPLES], reference["cells"]["H"])
    connection = metric.X**2 * grid.cell_derivative_fourth(vacuum_constraints["connection"], parity=-1)

    # Factor the near-identity differences before division.
    x_squared_offset = (metric.X - 1) * (metric.X + 1)
    q = (x_squared_offset - (metric.a - 1)) / metric.a
    u = ((metric.b - 1) - x_squared_offset) / (metric.X * (np.sqrt(metric.b) + metric.X))
    fields = {
        "H_production": official, "R_production": production, "R_connection": connection,
        "H_source": source, "H_kinematic": kinematic, "areal_radius_over_r_minus_one": u,
    }
    for width in WIDTHS:
        direct = areal_scalar(grid, q, u, width)
        adjusted = direct + connection
        fields[f"R_metric_{width}"] = direct
        fields[f"H_metric_{width}"] = direct + kinematic - source
        fields[f"H_adjusted_{width}"] = adjusted + kinematic - source
        fields[f"R_gap_{width}"] = production - adjusted
        np.testing.assert_allclose(
            official - fields[f"H_adjusted_{width}"], production - adjusted, atol=2e-19, rtol=0,
        )

    for field in vacuum.ALL_FIELDS:
        np.testing.assert_array_equal(getattr(metric, field), getattr(state.geometry, field))
    maxima = {}
    for name, values in fields.items():
        if not np.all(np.isfinite(values)):
            raise FloatingPointError(name)
        i = int(np.argmax(np.abs(values)))
        maxima[name] = {"cell": i, "signed": float(values[i])}
    return {
        "t": state.t, "cells": {name: values[SAMPLES].tolist() for name, values in fields.items()},
        "maxima": maxima,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dr", type=float, choices=(0.025, 0.05))
    parser.add_argument("arm", choices=("A", "B"))
    parser.add_argument("--half", action="store_true")
    args = parser.parse_args()
    suffix = "-half" if args.half else ""
    original = TRACE_DIR / f"{args.arm}-{args.dr:g}{suffix}.json"
    trace = json.loads(original.read_text())
    observations = []

    def record(state, vacuum, row):
        if row != trace["rows"][len(observations)]:
            raise ValueError(f"accepted slice disagrees with the pinned trace at t={state.t}")
        observations.append(observe(state, vacuum, row))

    replay = run(args.dr, args.arm, 2 if args.half else 1, accepted_observer=record)
    if replay["first_step_stages"] != trace["first_step_stages"] or len(observations) != len(trace["rows"]):
        raise ValueError("replay is missing rows or does not reproduce the first-step brackets")
    metadata = {key: value for key, value in replay.items() if key not in ("rows", "first_step_stages")}
    metadata["areal_script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    metadata["original_trace_sha256"] = hashlib.sha256(original.read_bytes()).hexdigest()
    print(json.dumps({"metadata": metadata, "rows": observations}, allow_nan=False))


if __name__ == "__main__":
    main()
