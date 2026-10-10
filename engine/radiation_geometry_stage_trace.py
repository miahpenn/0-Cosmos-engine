"""Diagnostic-only trace of the V5.5 radiation predictor / CMC interface.

This module does not modify production equations, states, admissibility checks,
or stage ordering. It wraps the existing stage calls and records the existing
RHS decomposition into conservative flux transport and geometric source.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from . import matter_rhs as matter_rhs_module
from . import matter_system as matter_system_module
from . import production_kernel as production_kernel_module
from .matter_system import Species


REPOSITORY = "miahpenn/0-Cosmos-engine"
BRANCH = "0star-central-clock"
RESOLUTION = 80
R_MAX = 80.0
CFL = 0.03
TARGET_TIME = 50.0
TAIL_WINDOW = 2.0
OUTER_CELLS = 5


def _git(command: list[str]) -> str:
    try:
        return subprocess.check_output(
            ["git", *command], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return "unavailable"


def _safe_float(value) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _array(values) -> np.ndarray:
    return np.asarray(values, dtype=float)


class RadiationStageTrace:
    def __init__(
        self, output_dir: Path,
        use_accepted_metric_for_predictor_radiation_recovery: bool = False,
    ):
        self.output_dir = output_dir
        self.use_accepted_metric_for_predictor_radiation_recovery = bool(
            use_accepted_metric_for_predictor_radiation_recovery
        )
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.started_wall = time.time()
        self.step_index = 0
        self.current_step = None
        self.active_radiation_evolve = None
        self.last_radiation_sources: dict[int, tuple[float, float]] = {}
        self.stage_snapshots: list[dict] = []
        self.predictor_budgets: list[dict] = []
        self.progress_samples: list[dict] = []
        self.global_geometry_minima: dict[str, dict] = {}
        self.failure: dict | None = None
        self.instrumentation_errors: list[dict] = []
        self.latest_state = None
        self.latest_dt = None
        self._patch_production_kernel()

    def geometry_and_metrics(self, grid, geometry):
        metric_slice = production_kernel_module.metric_slice_from_q(
            grid, geometry
        )
        metrics = matter_system_module._metric_arrays(metric_slice)
        sqrt_gamma = np.asarray([m.sqrt_gamma for m in metrics], dtype=float)
        gamma_rr_inv = np.asarray(
            [m.gamma_rr_inv for m in metrics], dtype=float
        )
        return metric_slice, metrics, sqrt_gamma, gamma_rr_inv

    def cone_values(self, grid, geometry, radiation):
        _, metrics, sqrt_gamma, gamma_rr_inv = self.geometry_and_metrics(
            grid, geometry
        )
        U_E = _array(radiation.energy_t)
        U_r = _array(radiation.momentum_r)
        C = U_E - np.sqrt(gamma_rr_inv) * np.abs(U_r)
        physical_E = U_E / sqrt_gamma
        physical_S_abs = (
            np.sqrt(gamma_rr_inv) * np.abs(U_r) / sqrt_gamma
        )
        ratio = np.full_like(U_E, np.nan)
        positive = U_E > 0.0
        ratio[positive] = (
            np.sqrt(gamma_rr_inv[positive])
            * np.abs(U_r[positive])
            / U_E[positive]
        )
        return {
            "metrics": metrics,
            "sqrt_gamma": sqrt_gamma,
            "gamma_rr_inv": gamma_rr_inv,
            "U_E": U_E,
            "U_r": U_r,
            "C": C,
            "physical_E": physical_E,
            "physical_S_abs": physical_S_abs,
            "ratio": ratio,
        }

    def _geometry_minima(self, grid, geometry, t, stage, step_index):
        r = _array(grid.centers)
        fields = ("a", "b", "X", "alpha")
        context_values = {
            name: _array(getattr(geometry, name)) for name in fields
        }
        for field, vals in context_values.items():
            if not vals.size or not np.any(np.isfinite(vals)):
                continue
            safe = np.where(np.isfinite(vals), vals, np.inf)
            i = int(np.argmin(safe))
            candidate = {
                "value": float(vals[i]),
                "t": float(t),
                "step_index": int(step_index),
                "stage": stage,
                "field": field,
                "index": i,
                "r": float(r[i]),
                "a": float(context_values["a"][i]),
                "b": float(context_values["b"][i]),
                "X": float(context_values["X"][i]),
                "alpha": float(context_values["alpha"][i]),
            }
            old = self.global_geometry_minima.get(field)
            if old is None or candidate["value"] < old["value"]:
                self.global_geometry_minima[field] = candidate

    def _should_expand(self, cone, t, force=False):
        finite_ratio = cone["ratio"][np.isfinite(cone["ratio"])]
        max_ratio = float(np.max(finite_ratio)) if finite_ratio.size else None
        # Keep full stage rows in a rolling two-time-unit buffer. The
        # failing time is not assumed in advance, so do not key this to TARGET_TIME.
        return True

    def record_stage(self, stage, grid, geometry, matter, t, tau,
                     step_index, dt, force=False, exception=None):
        try:
            cone = self.cone_values(grid, geometry, matter.radiation)
            self._geometry_minima(
                grid, geometry, t, stage, step_index
            )
        except Exception as exc:
            self.instrumentation_errors.append({
                "where": "record_stage",
                "stage": stage,
                "t": float(t),
                "error": f"{type(exc).__name__}: {exc}",
            })
            return None

        ratio = cone["ratio"]
        C = cone["C"]
        idx_c = int(np.nanargmin(C)) if C.size else 0
        finite_ratio = ratio[np.isfinite(ratio)]
        idx_r = (
            int(np.nanargmax(ratio)) if finite_ratio.size else idx_c
        )
        tmin = float(t)
        current = self.current_step or {}
        stage_row = {
            "t": tmin,
            "tau": _safe_float(tau),
            "step_index": int(step_index),
            "dt": _safe_float(dt),
            "stage": stage,
            "grid_index_min_C": idx_c,
            "r_min_C": float(grid.centers[idx_c]),
            "min_C": _safe_float(C[idx_c]),
            "grid_index_max_ratio": idx_r,
            "r_max_ratio": float(grid.centers[idx_r]),
            "max_ratio": _safe_float(ratio[idx_r]),
            "negative_energy_cells": np.where(cone["U_E"] < 0.0)[0].astype(int).tolist(),
            "negative_cone_cells": np.where(C < 0.0)[0].astype(int).tolist(),
            "exception": exception,
        }
        if self._should_expand(cone, tmin, force=force or exception is not None):
            n = len(cone["U_E"])
            lo = max(0, n - OUTER_CELLS)
            cells = []
            for i in range(lo, n):
                g = geometry
                m = cone["metrics"][i]
                cells.append({
                    "i": i,
                    "r": float(grid.centers[i]),
                    "U_E": float(cone["U_E"][i]),
                    "U_r": float(cone["U_r"][i]),
                    "sqrt_gamma": float(cone["sqrt_gamma"][i]),
                    "gamma_rr_inv": float(cone["gamma_rr_inv"][i]),
                    "C": float(C[i]),
                    "ratio_absS_over_E": _safe_float(ratio[i]),
                    "E": _safe_float(cone["physical_E"][i]),
                    "abs_S": _safe_float(cone["physical_S_abs"][i]),
                    "a": float(g.a[i]),
                    "b": float(g.b[i]),
                    "X": float(g.X[i]),
                    "alpha": float(g.alpha[i]),
                    "beta": float(g.beta[i]),
                })
            stage_row["outer_cell_snapshots"] = cells
            if current.get("rhs0_record") is not None and stage == "explicit_predictor_first_CMC_lapse":
                try:
                    budget = self.predictor_budget(
                        grid, geometry, matter, current["rhs0_record"],
                        current["dt"], tmin, step_index, stage
                    )
                    self.predictor_budgets.append(budget)
                    stage_row["predictor_budget_index"] = len(self.predictor_budgets) - 1
                except Exception as exc:
                    err = {
                        "where": "predictor_budget",
                        "stage": stage,
                        "t": tmin,
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                    self.instrumentation_errors.append(err)
                    stage_row["predictor_budget_error"] = err["error"]
        if self.current_step is not None:
            self.stage_snapshots.append(stage_row)
        return stage_row

    def predictor_budget(self, grid, predictor_geometry, predictor_matter,
                         rhs0, dt, t, step_index, stage):
        """Ordered, same-cell cone-margin budget; no state is changed."""
        g0 = rhs0["geometry"]
        m0 = rhs0["matter_radiation"]
        rhs_e = rhs0["rhs_energy"]
        rhs_s = rhs0["rhs_momentum"]
        sources = rhs0["sources"]
        cone0 = self.cone_values(grid, g0, m0)
        conep = self.cone_values(grid, predictor_geometry, predictor_matter.radiation)
        n = len(cone0["U_E"])
        lo = max(0, n - OUTER_CELLS)
        rows = []
        for i in range(lo, n):
            src = sources.get(i)
            if src is None:
                rows.append({
                    "i": i, "r": float(grid.centers[i]),
                    "source_capture_missing": True,
                    "C_accepted": float(cone0["C"][i]),
                    "C_predictor": float(conep["C"][i]),
                })
                continue
            src_e, src_s = float(src[0]), float(src[1])
            total_e, total_s = float(rhs_e[i]), float(rhs_s[i])
            flux_e, flux_s = total_e - src_e, total_s - src_s
            U_E0, U_r0 = float(cone0["U_E"][i]), float(cone0["U_r"][i])
            U_E_flux = U_E0 + dt * flux_e
            U_r_flux = U_r0 + dt * flux_s
            root_ginv0 = math.sqrt(float(cone0["gamma_rr_inv"][i]))
            C0 = U_E0 - root_ginv0 * abs(U_r0)
            C_flux = U_E_flux - root_ginv0 * abs(U_r_flux)
            U_E_source = U_E_flux + dt * src_e
            U_r_source = U_r_flux + dt * src_s
            C_source = U_E_source - root_ginv0 * abs(U_r_source)
            U_E_pred_actual = float(predictor_matter.radiation.energy_t[i])
            U_r_pred_actual = float(predictor_matter.radiation.momentum_r[i])
            root_ginv_p = math.sqrt(float(conep["gamma_rr_inv"][i]))
            C_pred = U_E_pred_actual - root_ginv_p * abs(U_r_pred_actual)
            dflux = C_flux - C0
            dsource = C_source - C_flux
            dmetric = C_pred - C_source
            dtotal = C_pred - C0
            closure = dflux + dsource + dmetric - dtotal
            rows.append({
                "i": i,
                "r": float(grid.centers[i]),
                "source_capture_missing": False,
                "C_accepted": C0,
                "U_E_accepted": U_E0,
                "U_r_accepted": U_r0,
                "delta_flux_U_E": float(dt * flux_e),
                "delta_flux_U_r": float(dt * flux_s),
                "delta_geometric_source_U_E": float(dt * src_e),
                "delta_geometric_source_U_r": float(dt * src_s),
                "total_rhs_U_E": total_e,
                "total_rhs_U_r": total_s,
                "source_U_E_rate": src_e,
                "source_U_r_rate": src_s,
                "flux_U_E_rate": flux_e,
                "flux_U_r_rate": flux_s,
                "C_after_flux": C_flux,
                "C_after_sources": C_source,
                "C_predictor": C_pred,
                "delta_C_flux": dflux,
                "delta_C_sources": dsource,
                "delta_C_metric": dmetric,
                "delta_C_total": dtotal,
                "closure_error": closure,
                "predictor_U_E_actual_minus_rhs": U_E_pred_actual - (U_E0 + dt * total_e),
                "predictor_U_r_actual_minus_rhs": U_r_pred_actual - (U_r0 + dt * total_s),
                "accepted_metric": {
                    "sqrt_gamma": float(cone0["sqrt_gamma"][i]),
                    "gamma_rr_inv": float(cone0["gamma_rr_inv"][i]),
                    "a": float(g0.a[i]), "b": float(g0.b[i]),
                    "X": float(g0.X[i]), "alpha": float(g0.alpha[i]),
                },
                "predictor_metric": {
                    "sqrt_gamma": float(conep["sqrt_gamma"][i]),
                    "gamma_rr_inv": float(conep["gamma_rr_inv"][i]),
                    "a": float(predictor_geometry.a[i]),
                    "b": float(predictor_geometry.b[i]),
                    "X": float(predictor_geometry.X[i]),
                    "alpha_before_solve": float(predictor_geometry.alpha[i]),
                },
                "predictor_U_E_actual": U_E_pred_actual,
                "predictor_U_r_actual": U_r_pred_actual,
            })
        return {
            "t": float(t),
            "step_index": int(step_index),
            "stage": stage,
            "dt": float(dt),
            "budget_order": [
                "accepted cone on accepted metric",
                "flux-only conservative update evaluated on accepted metric",
                "geometric-source update evaluated on accepted metric",
                "same predictor conservative state evaluated on predictor metric",
            ],
            "cells": rows,
            "max_abs_closure_error": max(
                (abs(r["closure_error"]) for r in rows if "closure_error" in r),
                default=None,
            ),
            "max_abs_predictor_U_E_reconstruction_error": max(
                (abs(r["predictor_U_E_actual_minus_rhs"]) for r in rows if "predictor_U_E_actual_minus_rhs" in r),
                default=None,
            ),
            "max_abs_predictor_U_r_reconstruction_error": max(
                (abs(r["predictor_U_r_actual_minus_rhs"]) for r in rows if "predictor_U_r_actual_minus_rhs" in r),
                default=None,
            ),
        }

    def _patch_production_kernel(self):
        kernel_cls = production_kernel_module.V55ProductionKernel
        original_solve = kernel_cls._solve_lapse
        original_rhs = kernel_cls._rhs
        original_step = kernel_cls.step
        original_evolve = matter_rhs_module.evolve_species
        original_radiation_source = matter_system_module._radiation_source
        trace = self

        def traced_radiation_source(*args, **kwargs):
            result = original_radiation_source(*args, **kwargs)
            active = trace.active_radiation_evolve
            if active is not None:
                i = active["source_index"]
                if i >= len(active["metric_r"]) - OUTER_CELLS:
                    active["sources"][i] = (float(result[0]), float(result[1]))
                active["source_index"] += 1
            return result

        def traced_evolve(*args, **kwargs):
            species = kwargs.get("species")
            if species is None and len(args) >= 4:
                species = args[3]
            if species is not Species.RADIATION:
                return original_evolve(*args, **kwargs)
            metric = args[0] if len(args) >= 1 else kwargs["metric"]
            self.active_radiation_evolve = {
                "source_index": 0,
                "sources": {},
                "metric_r": _array(metric.r),
            }
            try:
                out = original_evolve(*args, **kwargs)
                self.last_radiation_sources = dict(
                    self.active_radiation_evolve["sources"]
                )
                return out
            finally:
                self.active_radiation_evolve = None

        def traced_solve(grid, geometry, scalars, matter, radiation_recovery_metric=None):
            step = trace.current_step
            t = float(step["t0"]) if step is not None else 0.0
            tau = float(step["tau0"]) if step is not None else 0.0
            idx = int(step["solve_count"]) if step is not None else 0
            if step is not None:
                step["solve_count"] += 1
            labels = (
                "accepted_start_CMC_lapse",
                "explicit_predictor_first_CMC_lapse",
                "primary_predictor_CMC_lapse",
                "completed_state_CMC_lapse",
            )
            stage = labels[idx] if idx < len(labels) else f"CMC_lapse_call_{idx}"
            row = trace.record_stage(
                stage, grid, geometry, matter, t, tau,
                trace.step_index, step["dt"] if step is not None else None,
                force=(step is not None and idx == 1),
            )
            try:
                if radiation_recovery_metric is None:
                    result = original_solve(grid, geometry, scalars, matter)
                else:
                    result = original_solve(
                        grid, geometry, scalars, matter,
                        radiation_recovery_metric=radiation_recovery_metric,
                    )
                if row is not None and (
                    t >= TARGET_TIME - TAIL_WINDOW
                    or row.get("max_ratio") is not None and row["max_ratio"] >= 0.98
                ):
                    row["lapse_solve_status"] = "returned"
                return result
            except Exception as exc:
                detail = f"{type(exc).__name__}: {exc}"
                trace.record_stage(
                    stage + "_failure", grid, geometry, matter, t, tau,
                    trace.step_index, step["dt"] if step is not None else None,
                    force=True, exception=detail,
                )
                if step is not None:
                    step["failing_stage"] = stage
                raise

        def traced_rhs(kernel_self, state, radiation_recovery_metric=None):
            step = trace.current_step
            t = float(state.t)
            tau = float(state.tau)
            rhs_idx = int(step["rhs_count"]) if step is not None else 0
            if step is not None:
                step["rhs_count"] += 1
            stage = "rhs0_accepted_geometry" if rhs_idx == 0 else (
                "rhs1_primary_predictor" if rhs_idx == 1 else f"rhs_call_{rhs_idx}"
            )
            trace.record_stage(
                stage, state.grid, state.geometry, state.matter, t, tau,
                trace.step_index, step["dt"] if step is not None else None,
            )
            if radiation_recovery_metric is None:
                out = original_rhs(kernel_self, state)
            else:
                out = original_rhs(
                    kernel_self, state, radiation_recovery_metric=radiation_recovery_metric
                )
            srhs, mrhs, md = out
            if step is not None and rhs_idx == 0:
                step["rhs0_record"] = {
                    "geometry": state.geometry.copy(),
                    "matter_radiation": state.matter.radiation.copy(),
                    "rhs_energy": _array(mrhs["radiation"].energy_t).copy(),
                    "rhs_momentum": _array(mrhs["radiation"].momentum_r).copy(),
                    "sources": dict(trace.last_radiation_sources),
                }
            return out

        def traced_step(kernel_self, state, dt):
            trace.latest_state = state
            trace.latest_dt = float(dt)
            trace.current_step = {
                "step_index": trace.step_index,
                "t0": float(state.t),
                "tau0": float(state.tau),
                "dt": float(dt),
                "solve_count": 0,
                "rhs_count": 0,
                "rhs0_record": None,
                "failing_stage": None,
            }
            try:
                candidate = original_step(kernel_self, state, dt)
                trace.record_progress(candidate)
                trace.step_index += 1
                trace.prune_tail(float(candidate.t) - TAIL_WINDOW)
                return candidate
            except Exception as exc:
                trace.failure = {
                    "t_accepted": float(state.t),
                    "tau_accepted": float(state.tau),
                    "step_index": int(trace.step_index),
                    "dt": float(dt),
                    "t_attempted": float(state.t + dt),
                    "stage": trace.current_step.get("failing_stage"),
                    "exception_class": type(exc).__name__,
                    "exception_message": str(exc),
                    "wall_elapsed_seconds": float(time.time() - trace.started_wall),
                }
                raise
            finally:
                trace.current_step = None

        matter_rhs_module.evolve_species = traced_evolve
        matter_system_module._radiation_source = traced_radiation_source
        kernel_cls._solve_lapse = staticmethod(traced_solve)
        kernel_cls._rhs = traced_rhs
        kernel_cls.step = traced_step

        self._restore_hooks = (
            kernel_cls, original_solve, original_rhs, original_step,
            original_evolve, original_radiation_source,
        )

    def prune_tail(self, cutoff_t):
        self.stage_snapshots = [
            row for row in self.stage_snapshots
            if float(row.get("t", -math.inf)) >= cutoff_t
        ]
        self.predictor_budgets = [
            row for row in self.predictor_budgets
            if float(row.get("t", -math.inf)) >= cutoff_t
        ]

    def record_progress(self, state):
        try:
            cone = self.cone_values(
                state.grid, state.geometry, state.matter.radiation
            )
            n = len(cone["U_E"])
            lo = max(0, n - OUTER_CELLS)
            band = range(lo, n)
            ratios = [
                _safe_float(cone["ratio"][i]) for i in band
            ]
            Cs = [float(cone["C"][i]) for i in band]
            min_local = min(band, key=lambda i: cone["C"][i])
            max_local = max(
                (i for i in band if np.isfinite(cone["ratio"][i])),
                key=lambda i: cone["ratio"][i],
                default=min_local,
            )
            self.progress_samples.append({
                "t": float(state.t),
                "tau": float(state.tau),
                "accepted_step_index": int(self.step_index),
                "outer_band_i_range": [lo, n - 1],
                "outer_band_min_C": float(cone["C"][min_local]),
                "outer_band_min_C_i": int(min_local),
                "outer_band_max_ratio": _safe_float(cone["ratio"][max_local]),
                "outer_band_max_ratio_i": int(max_local),
                "outer_band_ratios": ratios,
                "outer_band_cone_margins": Cs,
                "alpha_min": float(np.min(_array(state.geometry.alpha))),
                "alpha_min_i": int(np.argmin(_array(state.geometry.alpha))),
            })
            self._geometry_minima(
                state.grid, state.geometry, state.t,
                "accepted_state_after_step", self.step_index
            )
        except Exception as exc:
            self.instrumentation_errors.append({
                "where": "record_progress",
                "t": float(state.t),
                "error": f"{type(exc).__name__}: {exc}",
            })

    def run(self):
        kernel = production_kernel_module.V55ProductionKernel()
        kernel.use_accepted_metric_for_predictor_radiation_recovery = (
            self.use_accepted_metric_for_predictor_radiation_recovery
        )
        state = kernel.initialize(
            resolution=RESOLUTION,
            r_max=R_MAX,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
        dt_nominal = CFL * float(state.grid.dr)
        start = time.time()
        last_report_step = -1
        status = "completed_without_reproducing_expected_failure"
        target_steps = int(math.ceil(TARGET_TIME / dt_nominal)) + 5
        while state.t < TARGET_TIME and self.step_index < target_steps:
            dt = min(dt_nominal, TARGET_TIME - state.t)
            try:
                state = kernel.step(state, dt)
                state.geometry.assert_finite_positive()
            except (FloatingPointError, ValueError) as exc:
                # The unmodified numerical failure is data for this diagnostic.
                # Do not repair or suppress it; the trace is the run product.
                status = "expected_radiation_failure_captured" if (
                    self.failure is not None
                    and (
                        "radiation" in self.failure["exception_message"].lower()
                        or self.failure.get("stage") == "explicit_predictor_first_CMC_lapse"
                    )
                ) else "numerical_failure_captured"
                if self.failure is None:
                    self.failure = {
                        "t_accepted": float(state.t),
                        "tau_accepted": float(state.tau),
                        "step_index": int(self.step_index),
                        "dt": float(dt),
                        "t_attempted": float(state.t + dt),
                        "stage": None,
                        "exception_class": type(exc).__name__,
                        "exception_message": str(exc),
                        "wall_elapsed_seconds": float(time.time() - self.started_wall),
                    }
                break
            if self.step_index - last_report_step >= 100:
                print(json.dumps({
                    "trace_progress": True,
                    "step": self.step_index,
                    "t": state.t,
                    "tau": state.tau,
                    "outer_band_max_ratio": self.progress_samples[-1].get("outer_band_max_ratio") if self.progress_samples else None,
                    "elapsed_seconds": time.time() - start,
                }), flush=True)
                last_report_step = self.step_index
        else:
            if state.t >= TARGET_TIME:
                status = "target_time_completed_without_failure"
            elif self.step_index >= target_steps:
                status = "step_cap_reached_without_failure"

        # Anchor the final trace window to the attempted failing step, not
        # to the planned target time.
        if self.failure is not None:
            attempted_t = float(self.failure.get("t_attempted", state.t))
            self.prune_tail(attempted_t - TAIL_WINDOW)
        return self.write_outputs(state, status, dt_nominal)

    def write_outputs(self, state, status, dt_nominal):
        hashes = {}
        for rel in (
            "engine/production_kernel.py",
            "engine/matter_system.py",
            "engine/matter_rhs.py",
            "engine/valencia.py",
            "engine/v55_matter.py",
            "engine/v55_pirk_adapter.py",
            "engine/cmc_gauge.py",
            "engine/stress_energy.py",
            "engine/true_cmc_pirk_kernel.py",
            "engine/radiation_geometry_stage_trace.py",
            ".github/workflows/zero_star_radiation_stage_trace.yml",
            ".github/workflows/zero_star_extended_b_recovery_metric.yml",
        ):
            path = Path(rel)
            hashes[rel] = (
                hashlib.sha256(path.read_bytes()).hexdigest()
                if path.is_file() else None
            )
        payload = {
            "schema_name": "0star_radiation_predictor_trace_v1",
            "purpose": "diagnostic-only trace; production defaults unchanged; optional predictor radiation recovery metric experiment",
            "status": status,
            "provenance": {
                "repository": REPOSITORY,
                "branch": os.environ.get("GITHUB_REF_NAME", _git(["branch", "--show-current"])),
                "instrumentation_commit": os.environ.get("GITHUB_SHA", _git(["rev-parse", "HEAD"])),
                "source_commit_before_trace_workflow": _git(["rev-parse", "HEAD~2"]),
                "workflow_commit": _git(["rev-parse", "HEAD^"]),
                "source_file_sha256": hashes,
                "run_id": os.environ.get("GITHUB_RUN_ID"),
                "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
                "python": sys.version,
                "platform": platform.platform(),
                "numpy": np.__version__,
                "configuration": {
                    "resolution": RESOLUTION,
                    "r_max": R_MAX,
                    "dr": float(state.grid.dr),
                    "cfl": CFL,
                    "dt_nominal": float(dt_nominal),
                    "target_time": TARGET_TIME,
                    "tail_window": TAIL_WINDOW,
                    "outer_cells_traced": OUTER_CELLS,
                    "amplitude": 0.01,
                    "width": 7.0,
                    "D_amplitude": 1.0e-10,
                    "include_radiation": True,
                    "use_accepted_metric_for_predictor_radiation_recovery":
                        self.use_accepted_metric_for_predictor_radiation_recovery,
                    "outer_boundary": "existing production light-constraint boundary and radiation outer closure",
                },
            },
            "failure": self.failure,
            "final_state": {
                "t": float(state.t),
                "tau": float(state.tau),
                "steps_completed": int(self.step_index),
                "steps_recorded_in_history": len(state.history),
            },
            "global_geometry_minima": self.global_geometry_minima,
            "progress_samples": self.progress_samples,
            "stage_snapshots": self.stage_snapshots,
            "predictor_budgets": self.predictor_budgets,
            "instrumentation_errors": self.instrumentation_errors,
            "wall_elapsed_seconds": float(time.time() - self.started_wall),
            "notes": [
                "The conservative radiation variables remain unmodified by the trace.",
                "The source is captured directly from the existing production radiation source routine.",
                "Flux transport rate is reconstructed as the existing total radiation RHS minus that captured source; evolve_species sums only those two components for radiation.",
                "The cone budget is ordered and uses the accepted spatial metric for the accepted, flux-only, and source-updated stages, then the predictor metric for the final predictor margin.",
                (
                    "Extended-B mode uses the accepted start-of-step metric only for predictor-stage radiation primitive inversion; stage metrics still enter the fluxes, geometric sources, stress projections, and CMC spatial operator."
                    if self.use_accepted_metric_for_predictor_radiation_recovery
                    else "Extended-B mode is disabled; the original predictor-stage radiation recovery metric is used."
                ),
                "No clipping, floor, damping, physical source insertion, or timestep adjustment is introduced.",
            ],
        }
        trace_path = self.output_dir / "trace.json"
        trace_path.write_text(json.dumps(payload, indent=2, allow_nan=False))
        summary = {
            "schema_name": payload["schema_name"],
            "status": status,
            "failure": self.failure,
            "final_state": payload["final_state"],
            "configuration": payload["provenance"]["configuration"],
            "instrumentation_commit": payload["provenance"]["instrumentation_commit"],
            "source_commit_before_trace_workflow": payload["provenance"]["source_commit_before_trace_workflow"],
            "source_file_sha256": hashes,
            "stage_snapshot_count": len(self.stage_snapshots),
            "predictor_budget_count": len(self.predictor_budgets),
            "progress_sample_count": len(self.progress_samples),
            "instrumentation_error_count": len(self.instrumentation_errors),
            "wall_elapsed_seconds": payload["wall_elapsed_seconds"],
            "trace_sha256": hashlib.sha256(trace_path.read_bytes()).hexdigest(),
        }
        summary_path = self.output_dir / "summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False))
        checksum_path = self.output_dir / "SHA256SUMS.txt"
        checksum_path.write_text(
            hashlib.sha256(trace_path.read_bytes()).hexdigest() + "  trace.json\n"
            + hashlib.sha256(summary_path.read_bytes()).hexdigest() + "  summary.json\n"
        )
        print(json.dumps(summary, indent=2), flush=True)
        return 0


def main():
    output = Path("runs/radiation-stage-trace")
    use_accepted_recovery = (
        os.environ.get("EXTENDED_B_ACCEPTED_METRIC_RECOVERY", "0") == "1"
    )
    trace = RadiationStageTrace(
        output,
        use_accepted_metric_for_predictor_radiation_recovery=use_accepted_recovery,
    )
    try:
        return trace.run()
    except Exception as exc:
        # Preserve an instrumentation failure separately; never call it a
        # physical failure or rewrite a production outcome.
        if trace.failure is None:
            trace.failure = {
                "stage": "trace_harness",
                "exception_class": type(exc).__name__,
                "exception_message": str(exc),
                "traceback": __import__("traceback").format_exc(),
            }
        try:
            trace.write_outputs(trace.latest_state, "trace_harness_failure", float("nan"))
        except Exception:
            print(__import__("traceback").format_exc(), file=sys.stderr)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
