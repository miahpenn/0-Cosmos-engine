"""Envelope field-lock diagnostic for the 0-star coupled production engine.

Diagnostic only: does not alter evolution equations or add source terms.
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


def grad(g, y, parity=1):
    return np.asarray(g.cell_derivative_fourth(np.asarray(y), parity=parity), dtype=float)


def radius_at_fraction(r, y, frac):
    y = np.asarray(y, dtype=float)
    m = float(np.max(y))
    if not np.isfinite(m) or m <= 0.0:
        return float("nan")
    ids = np.where(y >= frac * m)[0]
    return float(r[ids[-1]]) if len(ids) else float("nan")


def peak_info(r, y, deriv):
    a = np.abs(deriv)
    i = int(np.argmax(a))
    return float(a[i]), float(r[i])


def normalized_shape(r, y, edge, n=201):
    if not np.isfinite(edge) or edge <= 0:
        return None
    x = np.asarray(r, dtype=float) / edge
    mask = (x >= 0.0) & (x <= 1.0)
    if mask.sum() < 4:
        return None
    grid = np.linspace(0.0, 1.0, n)
    vals = np.interp(grid, x[mask], np.asarray(y, dtype=float)[mask])
    lo, hi = float(vals.min()), float(vals.max())
    if hi > lo:
        vals = (vals - lo) / (hi - lo)
    return vals


def shape_metrics(a, b):
    if a is None or b is None:
        return {"corr": float("nan"), "rms": float("nan")}
    aa = np.asarray(a) - np.mean(a)
    bb = np.asarray(b) - np.mean(b)
    den = float(np.linalg.norm(aa) * np.linalg.norm(bb))
    corr = float(np.dot(aa, bb) / den) if den > 0 else float("nan")
    rms = float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))
    return {"corr": corr, "rms": rms}


kernel = ProbeKernel()
state = kernel.initialize()

outdir = Path("runs/0star-envelope-lock") / CASE
outdir.mkdir(parents=True, exist_ok=True)

sample_times = [float(x) for x in np.arange(0.0, TARGET + 0.5, 1.0)]
rows = []
prev_shapes = {}

def capture(st):
    g = st.grid
    geom = st.geometry
    fields = st.scalars
    metric = metric_slice_from_q(g, geom)
    total = assemble_total_stress_energy(g, geom, fields, st.matter)
    _, _, moving = adapter.vendor_modules()
    explicit = moving.moving_puncture_explicit_rhs(g, geom)

    r = np.asarray(g.centers, dtype=float)
    R = r * np.sqrt(geom.b) / geom.X
    Rdot = R * (0.5 * explicit["b"] / geom.b - explicit["X"] / geom.X)
    theta = -geom.alpha * geom.K

    D_active = 2.0 * fields.PD**2 + fields.D**2
    alpha_r = grad(g, geom.alpha)
    D_active_r = grad(g, D_active)
    theta_r = grad(g, theta)
    Rdot_r = grad(g, Rdot)
    rho_r = grad(g, total.rho)
    j_r = grad(g, total.j)

    rad = project_species(metric, st.matter.radiation, Species.RADIATION)
    rad_r = grad(g, rad["rho"])
    rad_j_r = grad(g, rad["j"])

    i_D = int(np.argmax(D_active))
    Dmax = float(D_active[i_D])

    peaks = {
        "D": peak_info(r, D_active, D_active_r),
        "alpha": peak_info(r, geom.alpha, alpha_r),
        "theta": peak_info(r, theta, theta_r),
        "Rdot": peak_info(r, Rdot, Rdot_r),
        "rho_total": peak_info(r, total.rho, rho_r),
        "j_total": peak_info(r, total.j, j_r),
        "rho_rad": peak_info(r, rad["rho"], rad_r),
        "j_rad": peak_info(r, rad["j"], rad_j_r),
    }

    thresholds = {
        f"r_D_{int(frac*1000)/10:g}pct": radius_at_fraction(r, D_active, frac)
        for frac in (0.5, 0.1, 0.05, 0.01, 0.001)
    }
    r_D1 = thresholds["r_D_1pct"]

    shapes = {
        "D": normalized_shape(r, D_active / max(Dmax, 1e-300), r_D1),
        "alpha": normalized_shape(r, geom.alpha, r_D1),
    }
    sh = {}
    for name, cur in shapes.items():
        sh[name] = shape_metrics(prev_shapes.get(name), cur)
        prev_shapes[name] = cur

    row = {
        "t": float(st.t),
        "tau": float(st.tau),
        "D_active_max": Dmax,
        "r_D_peak": float(r[i_D]),
        "alpha_min": float(np.min(geom.alpha)),
        "alpha_min_r": float(r[np.argmin(geom.alpha)]),
        **thresholds,
    }

    for name, (amp, rr) in peaks.items():
        row[f"grad_{name}_max"] = amp
        row[f"r_grad_{name}"] = rr

    rDg = peaks["D"][1]
    for name in ("alpha", "theta", "Rdot", "rho_total", "j_total", "rho_rad", "j_rad"):
        row[f"offset_{name}_minus_Dgrad"] = peaks[name][1] - rDg

    row.update({
        "shape_D_corr_prev": sh["D"]["corr"],
        "shape_D_rms_prev": sh["D"]["rms"],
        "shape_alpha_corr_prev": sh["alpha"]["corr"],
        "shape_alpha_rms_prev": sh["alpha"]["rms"],
    })

    # Direct field values at the D-gradient peak and D 1% support edge.
    def interp(y, rr):
        return float(np.interp(rr, r, np.asarray(y, dtype=float)))

    row.update({
        "alpha_at_Dgrad": interp(geom.alpha, rDg),
        "theta_at_Dgrad": interp(theta, rDg),
        "Rdot_at_Dgrad": interp(Rdot, rDg),
        "rho_total_at_Dgrad": interp(total.rho, rDg),
        "j_total_at_Dgrad": interp(total.j, rDg),
        "alpha_at_D1pct": interp(geom.alpha, r_D1) if np.isfinite(r_D1) else float("nan"),
        "theta_at_D1pct": interp(theta, r_D1) if np.isfinite(r_D1) else float("nan"),
        "Rdot_at_D1pct": interp(Rdot, r_D1) if np.isfinite(r_D1) else float("nan"),
        "rho_total_at_D1pct": interp(total.rho, r_D1) if np.isfinite(r_D1) else float("nan"),
        "j_total_at_D1pct": interp(total.j, r_D1) if np.isfinite(r_D1) else float("nan"),
        "rho_rad_max": float(np.max(rad["rho"])),
        "j_rad_max_abs": float(np.max(np.abs(rad["j"]))),
    })
    rows.append(row)

capture(state)
next_sample = 1
failure = None

while state.t < TARGET:
    try:
        state = kernel.step(state, min(CFL * state.grid.dr, TARGET - state.t))
        state.geometry.assert_finite_positive()
        while next_sample < len(sample_times) and state.t >= sample_times[next_sample] - 0.5 * CFL * state.grid.dr:
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

def fit_speed(key):
    t = np.array([x["t"] for x in rows], dtype=float)
    y = np.array([x[key] for x in rows], dtype=float)
    mask = (t >= min(7.5, TARGET / 2.0)) & np.isfinite(y)
    if mask.sum() < 2:
        return float("nan")
    return float(np.polyfit(t[mask], y[mask], 1)[0])

summary = {
    "Dgrad_speed_late": fit_speed("r_grad_D"),
    "alpha_grad_speed_late": fit_speed("r_grad_alpha"),
    "theta_grad_speed_late": fit_speed("r_grad_theta"),
    "Rdot_grad_speed_late": fit_speed("r_grad_Rdot"),
    "D1pct_speed_late": fit_speed("r_D_1pct"),
    "final_offsets": {k: rows[-1][k] for k in rows[-1] if k.startswith("offset_")},
    "final_radii": {k: rows[-1][k] for k in rows[-1] if k.startswith("r_grad_") or k.startswith("r_D_")},
    "final_lapse": {k: rows[-1][k] for k in ("alpha_min", "alpha_min_r")},
}

result = {
    "case": CASE,
    "D_amplitude": DAMP,
    "include_radiation": INCLUDE_RADIATION,
    "resolution": N,
    "r_max": RMAX,
    "cfl": CFL,
    "target_time": TARGET,
    "t_final": float(state.t),
    "tau_final": float(state.tau),
    "failure": failure,
    "turnaround_events": sum(e.kind == "turnaround" for e in state.cycle.events),
    "reexpansion_events": sum(e.kind == "re-expansion_crossing" for e in state.cycle.events),
    "handoff_events": len(state.handoffs),
    "summary": summary,
    "rows": rows,
}

(outdir / "envelope_lock.json").write_text(json.dumps(result, indent=2, default=float))
print(json.dumps({
    "case": CASE,
    "t_final": result["t_final"],
    "tau_final": result["tau_final"],
    "failure": failure,
    "summary": summary,
}, indent=2, default=float))
