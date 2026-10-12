#!/usr/bin/env python3
"""Matched short A/B test for the opt-in discrete-consistent initial metric.

Both arms use the same pinned production equations and integrator, N=40/60/80/120,
rmax=80, CFL=.03 and target t=2.16. The candidate arm solves only the initial
metric B so the exact discrete Hamiltonian constraint closes at t=0. Matter,
K, Aa, source terms, gauge equations and evolution code are unchanged.

The script runs each resolution independently so an initialization failure is
recorded without deleting completed results from other resolutions.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys

import numpy as np

from engine.campaign import CampaignConfig
from engine.run_production import run_campaign
from engine.production_kernel import V55ProductionKernel
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.scalar_system import ScalarFields
from engine.matter_system import ConservedSpecies
from engine.v55_matter import V55MatterState
from engine import v55_pirk_adapter as adapter
from engine.stress_energy import assemble_total_stress_energy

SOURCE_COMMIT = os.environ.get("GITHUB_SHA", "unavailable")
BASELINE_SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
RESOLUTIONS = (40, 60, 80, 120)
R_MAX = 80.0
CFL = 0.03
FINAL_TIME = 2.16
CHECKPOINT_INTERVAL = 0.24
FIXED_RADII = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)
ROOT = Path("runs/discrete-consistent-initial-production-test")
BASELINE_ROOT = ROOT / "baseline"
CANDIDATE_ROOT = ROOT / "candidate"
ARCHIVED_BASELINE_ROOT = Path("baseline-spatial-artifact")
OUT = Path("results/discrete-consistent-initial-production-test")

grid_ops, vacuum, _ = adapter.vendor_modules()
GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = {"dark_matter": "dm", "baryons": "baryon", "radiation": "radiation"}
CONSERVED = ("rest", "energy_t", "momentum_r")


def output_path(root: Path, resolution: int, t: float) -> Path:
    return root / f"N{resolution}" / f"checkpoint_t{t:012.6f}.npz"


class RecordingKernel(V55ProductionKernel):
    """Record t=0 alongside the strict production runner's scheduled checkpoints."""
    def __init__(self, output_root: Path):
        self.output_root = Path(output_root)
        super().__init__()

    def initialize(self, resolution=80, r_max=40.0, *, amplitude=0.01, width=7.0,
                   D_amplitude=1.0e-10, include_radiation=True):
        state = super().initialize(
            resolution=resolution, r_max=r_max, amplitude=amplitude, width=width,
            D_amplitude=D_amplitude, include_radiation=include_radiation,
        )
        self.checkpoint(state, output_path(self.output_root, resolution, state.t))
        return state


class DiscreteConsistentRecordingKernel(RecordingKernel):
    """Opt-in B-solve candidate; recursive initializer call is guarded."""
    def __init__(self, output_root: Path):
        super().__init__(output_root)
        self._inside_discrete_solve = False
        self.solver_reports = {}

    def initialize(self, resolution=80, r_max=40.0, *, amplitude=0.01, width=7.0,
                   D_amplitude=1.0e-10, include_radiation=True):
        if self._inside_discrete_solve:
            # discrete_consistent_state obtains an otherwise standard production
            # state before solving B. RecordingKernel writes it first; the
            # completed candidate below overwrites the same t=0 checkpoint.
            return super().initialize(
                resolution=resolution, r_max=r_max, amplitude=amplitude, width=width,
                D_amplitude=D_amplitude, include_radiation=include_radiation,
            )

        self._inside_discrete_solve = True
        try:
            state, B, residual_history, info = discrete_consistent_state(
                self, resolution=resolution, r_max=r_max, amplitude=amplitude,
                width=width, D_amplitude=D_amplitude, include_radiation=include_radiation,
                max_iter=12, tol=1.0e-13, max_cond=1.0e12, return_info=True,
            )
        finally:
            self._inside_discrete_solve = False

        # The metric solve alters only a,b,X,Lambda. Retain the original K,
        # Aa, scalar fields and matter; re-solve the same CMC gauge on the new slice.
        state.geometry.alpha = np.asarray(
            self._solve_lapse(state.grid, state.geometry, state.scalars, state.matter)[0],
            dtype=float,
        ).copy()
        state.geometry.beta.fill(0.0)
        state.geometry.B.fill(0.0)
        self.solver_reports[str(resolution)] = {
            **info,
            "residual_history": [float(x) for x in residual_history],
            "B_min": float(np.min(B)),
            "B_max": float(np.max(B)),
            "initial_H_max": float(residual_history[-1]),
            "metric_identity_max": float(np.max(np.abs(
                state.geometry.a * state.geometry.b**2 - 1.0
            ))),
            "t0_checkpoint": str(output_path(self.output_root, resolution, 0.0)),
        }
        self.checkpoint(state, output_path(self.output_root, resolution, state.t))
        return state


def load_arrays(path: Path):
    if not path.is_file():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=False) as z:
        return {key: np.asarray(z[key]).copy() for key in z.files}


def state_from_checkpoint(path: Path):
    data = load_arrays(path)
    r = np.asarray(data["r"], dtype=float)
    n = int(r.size)
    dr = float(r[1] - r[0])
    grid = grid_ops.SphericalCellGrid(n, n * dr)
    if not np.allclose(grid.centers, r, rtol=0.0, atol=1.0e-12 * max(1.0, n * dr)):
        raise ValueError(f"{path}: reconstructed grid does not match saved radii")
    _, geom, _ = adapter.vendor_modules()
    geometry = geom.flat_state(grid)
    for name in GEOMETRY_FIELDS:
        setattr(geometry, name, np.asarray(data[name], dtype=float).copy())
    scalars = ScalarFields(*(np.asarray(data[name], dtype=float).copy() for name in SCALAR_FIELDS))
    parts = {}
    for name, prefix in SPECIES.items():
        parts[name] = ConservedSpecies(*(
            np.asarray(data[f"{prefix}_{field}"], dtype=float).copy()
            for field in CONSERVED
        ))
    matter = V55MatterState(
        dark_matter=parts["dark_matter"], baryons=parts["baryons"], radiation=parts["radiation"]
    )
    return grid, geometry, scalars, matter, data


def metric_connection(grid, geometry):
    r = np.asarray(grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    d_a = grid.cell_derivative_fourth(a, parity=+1)
    d_b = grid.cell_derivative_fourth(b, parity=+1)
    return d_a/(2.0*a*a) - d_b/(a*b) + 2.0/r*(1.0/b - 1.0/a)


def summarize(values, radius, mask):
    values = np.asarray(values, dtype=float)
    ids = np.flatnonzero(np.asarray(mask, dtype=bool) & np.isfinite(values))
    if not ids.size:
        return {"count": 0}
    index = int(ids[np.argmax(np.abs(values[ids]))])
    return {
        "count": int(ids.size),
        "max_abs": float(abs(values[index])),
        "peak_signed": float(values[index]),
        "peak_r": float(radius[index]),
        "peak_i": index,
        "rms": float(np.sqrt(np.mean(values[ids]**2))),
    }


def witness(path: Path):
    grid, g, scalars, matter, raw_data = state_from_checkpoint(path)
    radius = np.asarray(grid.centers, dtype=float)
    total = assemble_total_stress_energy(grid, g, scalars, matter)
    raw = vacuum.constraints(grid, g)
    fields = {
        "H": np.asarray(raw["hamiltonian"], dtype=float)
              - 16.0*math.pi*np.asarray(total.rho, dtype=float),
        "M": np.asarray(raw["momentum"], dtype=float)
              - 8.0*math.pi*np.asarray(total.j, dtype=float),
        "C_lambda": np.asarray(g.Lambda, dtype=float) - metric_connection(grid, g),
        "determinant": np.asarray(g.a, dtype=float)*np.asarray(g.b, dtype=float)**2 - 1.0,
    }
    masks = {
        "full_domain": np.ones(radius.shape, dtype=bool),
        "first_two_cells": np.arange(radius.size) < 2,
        "first_five_cells": np.arange(radius.size) < 5,
        "formal_stencil_i_ge_2": np.arange(radius.size) >= 2,
        "admitted_inner_i_ge_2_r_lt_5": (np.arange(radius.size) >= 2) & (radius < 5.0),
        "fixed_annulus_2_le_r_lt_5": (radius >= 2.0) & (radius < 5.0),
        "r_lt_5": radius < 5.0,
        "r_ge_5": radius >= 5.0,
    }
    profile_summary = {
        name: {region: summarize(values, radius, mask) for region, mask in masks.items()}
        for name, values in fields.items()
    }
    fixed_radius = {}
    for r_target in FIXED_RADII:
        if radius[0] <= r_target <= radius[-1]:
            fixed_radius[str(r_target)] = {
                name: float(np.interp(r_target, radius, values))
                for name, values in fields.items()
            }
    first_cells = {}
    for i in range(min(5, radius.size)):
        first_cells[str(i)] = {
            "r": float(radius[i]),
            **{name: float(values[i]) for name, values in fields.items()},
        }
    return {
        "checkpoint": str(path),
        "time": float(raw_data["t"]),
        "tau": float(raw_data["tau"]),
        "resolution": int(radius.size),
        "r_max": float(radius.size * grid.dr),
        "summary": profile_summary,
        "fixed_radius": fixed_radius,
        "first_five_cells": first_cells,
    }


def compare_witnesses(baseline, candidate):
    delta = {}
    for field in ("H", "M", "C_lambda", "determinant"):
        delta[field] = {}
        for region in (
            "full_domain", "first_two_cells", "formal_stencil_i_ge_2",
            "admitted_inner_i_ge_2_r_lt_5", "fixed_annulus_2_le_r_lt_5", "r_ge_5",
        ):
            b = baseline["summary"][field][region]
            c = candidate["summary"][field][region]
            if not b.get("count") or not c.get("count"):
                delta[field][region] = {
                    "baseline_count": b.get("count"), "candidate_count": c.get("count"),
                    "delta_max_abs": None, "delta_rms": None,
                }
            else:
                delta[field][region] = {
                    "baseline_max_abs": b["max_abs"],
                    "candidate_max_abs": c["max_abs"],
                    "candidate_minus_baseline_max_abs": c["max_abs"] - b["max_abs"],
                    "baseline_rms": b["rms"],
                    "candidate_rms": c["rms"],
                    "candidate_minus_baseline_rms": c["rms"] - b["rms"],
                }
    first_delta = {}
    for i in range(min(5, len(candidate["first_five_cells"]))):
        key = str(i)
        first_delta[key] = {
            field: {
                "baseline": baseline["first_five_cells"][key][field],
                "candidate": candidate["first_five_cells"][key][field],
                "candidate_minus_baseline": (
                    candidate["first_five_cells"][key][field] - baseline["first_five_cells"][key][field]
                ),
            }
            for field in ("H", "M", "C_lambda", "determinant")
        }
    return {
        "baseline": baseline,
        "candidate": candidate,
        "candidate_minus_baseline": delta,
        "first_five_cell_delta": first_delta,
    }


def run_arm(kernel, root: Path, resolution: int):
    config = CampaignConfig(
        resolutions=(resolution,),
        r_max=R_MAX,
        cfl=CFL,
        final_time=FINAL_TIME,
        checkpoint_interval=CHECKPOINT_INTERVAL,
        output_dir=str(root),
    )
    result = run_campaign(kernel=kernel, config=config)
    if len(result) != 1 or result[0].status != "completed" or abs(result[0].t_final - FINAL_TIME) > 1.0e-9:
        raise RuntimeError(f"strict bounded admission failed for N{resolution}: {result}")
    return result[0]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ROOT.mkdir(parents=True, exist_ok=True)
    cases = {}
    arm_results = {"baseline": {}, "candidate": {}}
    failures = []
    for resolution in RESOLUTIONS:
        try:
            kernel = RecordingKernel(BASELINE_ROOT)
            result = run_arm(kernel, BASELINE_ROOT, resolution)
            arm_results["baseline"][str(resolution)] = {
                "status": result.status, "t_final": result.t_final,
                "tau_final": result.tau_final, "steps": result.steps,
                "misner_sharp_current_max": result.misner_sharp_current_max,
                "misner_sharp_current_rms": result.misner_sharp_current_rms,
            }
        except Exception as exc:
            failures.append({
                "arm": "baseline", "resolution": resolution,
                "error_type": type(exc).__name__, "error": str(exc),
            })

        try:
            kernel = DiscreteConsistentRecordingKernel(CANDIDATE_ROOT)
            result = run_arm(kernel, CANDIDATE_ROOT, resolution)
            arm_results["candidate"][str(resolution)] = {
                "status": result.status, "t_final": result.t_final,
                "tau_final": result.tau_final, "steps": result.steps,
                "misner_sharp_current_max": result.misner_sharp_current_max,
                "misner_sharp_current_rms": result.misner_sharp_current_rms,
            }
            solver_path = CANDIDATE_ROOT / f"N{resolution}" / "initial_solver.json"
            solver_path.write_text(json.dumps(
                kernel.solver_reports[str(resolution)], indent=2, sort_keys=True, allow_nan=False
            ))
        except Exception as exc:
            failures.append({
                "arm": "candidate", "resolution": resolution,
                "error_type": type(exc).__name__, "error": str(exc),
            })

    comparisons = []
    for resolution in RESOLUTIONS:
        case = {"resolution": resolution, "states": {}, "comparison": {}}
        for epoch, t in (("t0", 0.0), ("t2p16", FINAL_TIME)):
            bp = output_path(BASELINE_ROOT, resolution, t)
            cp = output_path(CANDIDATE_ROOT, resolution, t)
            if bp.exists():
                case["states"].setdefault(epoch, {})["baseline"] = witness(bp)
            if cp.exists():
                case["states"].setdefault(epoch, {})["candidate"] = witness(cp)
            if bp.exists() and cp.exists():
                pair = compare_witnesses(
                    case["states"][epoch]["baseline"], case["states"][epoch]["candidate"]
                )
                case["comparison"][epoch] = pair
        # Cross-check the baseline trajectory at t=2.16 against the archived,
        # exact-source spatial refinement checkpoint if present.
        archived = list(ARCHIVED_BASELINE_ROOT.rglob(
            f"checkpoint_t{FINAL_TIME:012.6f}.npz"
        ))
        archived = [p for p in archived if p.parent.name == f"N{resolution}"]
        if len(archived) == 1:
            try:
                case["archived_exact_source_baseline_t2p16"] = witness(archived[0])
                if output_path(BASELINE_ROOT, resolution, FINAL_TIME).exists():
                    case["same_run_vs_archived_baseline"] = compare_witnesses(
                        case["archived_exact_source_baseline_t2p16"],
                        case["states"]["t2p16"]["baseline"],
                    )
            except Exception as exc:
                failures.append({
                    "arm": "archived_baseline_audit", "resolution": resolution,
                    "error_type": type(exc).__name__, "error": str(exc),
                })
        comparisons.append(case)

    baseline_complete = len(arm_results["baseline"]) == len(RESOLUTIONS)
    candidate_complete = len(arm_results["candidate"]) == len(RESOLUTIONS)
    report = {
        "schema": "discrete_consistent_initial_production_ab_t2p16_v1",
        "status": (
            "completed" if baseline_complete and candidate_complete and not failures
            else "incomplete_or_failed"
        ),
        "kind": "matched_short_production_A_B_diagnostic_only",
        "source_commit": SOURCE_COMMIT,
        "baseline_exact_source_commit_for_archived_crosscheck": BASELINE_SOURCE_COMMIT,
        "settings": {
            "resolutions": list(RESOLUTIONS), "r_max": R_MAX, "CFL": CFL,
            "final_time": FINAL_TIME, "checkpoint_interval": CHECKPOINT_INTERVAL,
            "radiation": True, "kernel_equations_changed": False,
            "matter_state_or_source_changed": False, "production_default_initializer_changed": False,
            "candidate_change": "opt-in solve for initial B so the discrete H residual closes at t=0; K/Aa/scalars/matter retained from the standard initializer",
            "t0_and_t2p16_saved_for_both_arms": True,
        },
        "arm_results": arm_results,
        "candidate_solver_reports": {
            str(n): json.loads((CANDIDATE_ROOT / f"N{n}" / "initial_solver.json").read_text())
            for n in RESOLUTIONS
            if (CANDIDATE_ROOT / f"N{n}" / "initial_solver.json").exists()
        },
        "comparisons": comparisons,
        "failures": failures,
        "interpretation_guardrail": (
            "The initial H residual is enforced by the candidate solve, so t=0 H is a solver check, not independent validation. "
            "The decision is whether H/M/C/determinant improve together at t=2.16, including admitted-inner and fixed-radius norms. "
            "Do not admit or promote an initial-data candidate on H alone."
        ),
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    print("DISCRETE_CONSISTENT_INITIAL_PRODUCTION_AB " + json.dumps({
        "status": report["status"], "settings": report["settings"],
        "arm_results": report["arm_results"],
        "candidate_solver_reports": report["candidate_solver_reports"],
        "compact_comparison": [{
            "N": c["resolution"],
            "t0": {
                field: c.get("comparison", {}).get("t0", {}).get("candidate_minus_baseline", {}).get(field, {})
                for field in ("H","M","C_lambda","determinant")
            },
            "t2p16": {
                field: c.get("comparison", {}).get("t2p16", {}).get("candidate_minus_baseline", {}).get(field, {})
                for field in ("H","M","C_lambda","determinant")
            },
            "candidate_t2p16_first_cells": c.get("states",{}).get("t2p16",{}).get("candidate",{}).get("first_five_cells",{}),
        } for c in comparisons],
        "failures": failures,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
