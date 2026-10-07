"""Space-time D->geometry lock diagnostic.

This is diagnostic-only: it records the production engine without modifying
the evolution equations or adding source terms.

For each case it stores full radial profiles at regular system-time samples.
The aggregate analysis then computes lag surfaces

    C_DX(dr, dt) = corr(|d_r D_active|(r,t),
                       |d_r X|(r+dr,t+dt))

and the same maps for D=0-subtracted geometry fields on the matched N=160
grid. A positive (dr, dt) peak is interpreted only as a directional
space-time lock consistent with D leading X; it is not by itself proof of
causality.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.v55_initial import build_initial_data
from engine.v55_matter import metric_slice_from_q
from engine.stress_energy import assemble_total_stress_energy
from engine.matter_system import Species, project_species
from engine import v55_pirk_adapter as adapter


N = int(os.environ["RESOLUTION"])
RMAX = float(os.environ["RMAX"])
CFL = float(os.environ["CFL"])
TARGET = float(os.environ.get("TARGET", "22.5"))
SAMPLE_DT = float(os.environ.get("SAMPLE_DT", "0.25"))
DAMP = float(os.environ["DAMP"])
INCLUDE_RADIATION = os.environ["RADIATION"] == "1"
CASE = os.environ["CASE_LABEL"]


class ProbeKernel(V55ProductionKernel):
    def initialize(self, resolution=N, r_max=RMAX, **kwargs):
        grid_ops, vacuum, _ = adapter.vendor_modules()
        grid = grid_ops.SphericalCellGrid(resolution, r_max)
        init = build_initial_data(
            grid,
            vacuum.flat_state,
            D_amplitude=DAMP,
            include_radiation=INCLUDE_RADIATION,
        )
        init.geometry.alpha = self._solve_lapse(
            grid, init.geometry, init.scalars, init.matter
        )[0]
        init.geometry.beta.fill(0.0)
        init.geometry.B.fill(0.0)
        return ProductionState(
            grid=grid,
            geometry=init.geometry,
            scalars=init.scalars,
            matter=init.matter,
            e_folds=0.0,
        )


def deriv(g, y):
    return np.asarray(
        g.cell_derivative_fourth(np.asarray(y), parity=1), dtype=float
    )


kernel = ProbeKernel()
state = kernel.initialize()

outdir = Path("runs/0star-spacetime-lock") / CASE
outdir.mkdir(parents=True, exist_ok=True)

sample_times = np.arange(0.0, TARGET + 0.5 * SAMPLE_DT, SAMPLE_DT)
times: list[float] = []
taus: list[float] = []
profiles: dict[str, list[np.ndarray]] = {
    "D_active_grad_abs": [],
    "alpha_grad_abs": [],
    "theta_grad_abs": [],
    "Rdot_grad_abs": [],
    "rho_total_grad_abs": [],
    "j_total_grad_abs": [],
    "rho_rad_grad_abs": [],
    "j_rad_grad_abs": [],
    "D_active": [],
    "alpha": [],
    "theta": [],
    "Rdot": [],
    "rho_total": [],
    "j_total": [],
}
r_grid = None


def capture(st):
    global r_grid

    g = st.grid
    geom = st.geometry
    fields = st.scalars
    metric = metric_slice_from_q(g, geom)
    total = assemble_total_stress_energy(g, geom, fields, st.matter)
    _, _, moving = adapter.vendor_modules()
    explicit = moving.moving_puncture_explicit_rhs(g, geom)

    r = np.asarray(g.centers, dtype=float)
    if r_grid is None:
        r_grid = r.copy()

    R = r * np.sqrt(geom.b) / geom.X
    Rdot = R * (
        0.5 * explicit["b"] / geom.b - explicit["X"] / geom.X
    )
    theta = -geom.alpha * geom.K
    D_active = 2.0 * fields.PD**2 + fields.D**2

    rad = project_species(
        metric, st.matter.radiation, Species.RADIATION
    )

    base = {
        "D_active": np.asarray(D_active, dtype=float),
        "alpha": np.asarray(geom.alpha, dtype=float),
        "theta": np.asarray(theta, dtype=float),
        "Rdot": np.asarray(Rdot, dtype=float),
        "rho_total": np.asarray(total.rho, dtype=float),
        "j_total": np.asarray(total.j, dtype=float),
    }
    for name, arr in base.items():
        profiles[name].append(arr)
        profiles[name + "_grad_abs"].append(np.abs(deriv(g, arr)))

    profiles["rho_rad_grad_abs"].append(
        np.abs(deriv(g, np.asarray(rad["rho"], dtype=float)))
    )
    profiles["j_rad_grad_abs"].append(
        np.abs(deriv(g, np.asarray(rad["j"], dtype=float)))
    )

    times.append(float(st.t))
    taus.append(float(st.tau))


capture(state)
next_sample = 1
dt = CFL * state.grid.dr
failure = None

while state.t < TARGET:
    try:
        state = kernel.step(state, min(dt, TARGET - state.t))
        state.geometry.assert_finite_positive()
        while (
            next_sample < len(sample_times)
            and state.t >= sample_times[next_sample] - 0.5 * dt
        ):
            capture(state)
            next_sample += 1
    except Exception as exc:
        failure = {
            "type": type(exc).__name__,
            "message": str(exc),
            "t": float(state.t),
            "tau": float(state.tau),
        }
        break

np.savez_compressed(
    outdir / "profiles.npz",
    r=np.asarray(r_grid, dtype=float),
    t=np.asarray(times, dtype=float),
    tau=np.asarray(taus, dtype=float),
    **{k: np.asarray(v, dtype=float) for k, v in profiles.items()},
)

meta = {
    "case": CASE,
    "D_amplitude": DAMP,
    "include_radiation": INCLUDE_RADIATION,
    "resolution": N,
    "r_max": RMAX,
    "cfl": CFL,
    "target_time": TARGET,
    "sample_dt_requested": SAMPLE_DT,
    "sample_count": len(times),
    "t_final": float(state.t),
    "tau_final": float(state.tau),
    "failure": failure,
    "turnaround_events": sum(
        e.kind == "turnaround" for e in state.cycle.events
    ),
    "reexpansion_events": sum(
        e.kind == "re-expansion_crossing" for e in state.cycle.events
    ),
    "handoff_events": len(state.handoffs),
}
(outdir / "meta.json").write_text(json.dumps(meta, indent=2, default=float))
print(json.dumps(meta, indent=2, default=float))
