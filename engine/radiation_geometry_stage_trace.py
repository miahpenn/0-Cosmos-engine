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
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from . import matter_rhs as matter_rhs_module
from . import matter_system as matter_system_module
from . import production_kernel as production_kernel_module
from .adaptive_step import (
    RadiationStepSizeUnderflow,
    advance_with_radiation_admissibility_retries,
)
from .matter_system import Species


REPOSITORY = "miahpenn/0-Cosmos-engine"
BRANCH = "0star-central-clock"
RESOLUTION = 80
R_MAX = 80.0
CFL = 0.03
TARGET_TIME = 50.0
# Optional central proper-time endpoint for comparisons across differently gauged domains.
TARGET_TAU: float | None = None
TAIL_WINDOW = 2.0
OUTER_CELLS = 5


def _completion_status(
    t: float, tau: float, target_time: float, target_tau: float | None
) -> str | None:
    if target_tau is not None:
        if tau >= target_tau:
            return "target_proper_time_completed_without_failure"
        if t >= target_time:
            return "coordinate_cap_reached_before_proper_time_target"
        return None
    return "target_time_completed_without_failure" if t >= target_time else None


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


def _pminus_boundary_terms(grid, geometry, energy_density, momentum_radial) -> dict:
    """Read the pinned outer P-minus characteristic and its existing target.

    This is a pure diagnostic: it reads the live state and never modifies it.
    """
    production_kernel_module.adapter.vendor_modules()
    from bssn_characteristic_boundary import (
        light_characteristic_minus,
        pminus_source,
    )

    omega = _array(light_characteristic_minus(grid, geometry))
    source = _array(pminus_source(
        grid, geometry, energy_density, momentum_radial
    ))
    dr = float(grid.dr)
    target = float(omega[-3] + 2.0 * dr * source[-2])
    return {
        "omega_outer": float(omega[-1]),
        "omega_third_last": float(omega[-3]),
        "source_penultimate": float(source[-2]),
        "target_omega_outer": target,
        "target_minus_omega_outer": float(target - omega[-1]),
        "pminus_reconstruction_residual": float(
            (omega[-1] - omega[-3]) / (2.0 * dr) - source[-2]
        ),
        "Aa_outer": float(geometry.Aa[-1]),
        "K_outer": float(geometry.K[-1]),
        "Lambda_outer": float(geometry.Lambda[-1]),
        "a_outer": float(geometry.a[-1]),
        "b_outer": float(geometry.b[-1]),
        "X_outer": float(geometry.X[-1]),
        "alpha_outer": float(geometry.alpha[-1]),
        "alpha_penultimate": float(geometry.alpha[-2]),
    }


class RadiationStageTrace:
    def __init__(
        self, output_dir: Path,
        use_accepted_metric_for_predictor_radiation_recovery: bool = False,
        use_radiation_admissibility_retry: bool = False,
    ):
        self.output_dir = output_dir
        self.use_accepted_metric_for_predictor_radiation_recovery = bool(
            use_accepted_metric_for_predictor_radiation_recovery
        )
        self.use_radiation_admissibility_retry = bool(use_radiation_admissibility_retry)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.started_wall = time.time()
        self.step_index = 0
        self.current_step = None
        self.active_radiation_evolve = None
        self.active_boundary_context = None
        self.last_radiation_sources: dict[int, tuple[float, float]] = {}
        self.boundary_update_snapshots: list[dict] = []
        self.stage_snapshots: list[dict] = []
        self.predictor_budgets: list[dict] = []
        self.completed_state_budgets: list[dict] = []
        self.geometry_predictor_budgets: list[dict] = []
        self.progress_samples: list[dict] = []
        self.global_geometry_minima: dict[str, dict] = {}
        self.failure: dict | None = None
        self.rejected_step_attempts: list[dict] = []
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


    @staticmethod
    def _budget_metric_provenance(geometry, cone, i):
        return {
            "sqrt_gamma": float(cone["sqrt_gamma"][i]),
            "gamma_rr_inv": float(cone["gamma_rr_inv"][i]),
            "a": float(geometry.a[i]),
            "b": float(geometry.b[i]),
            "X": float(geometry.X[i]),
            "alpha": float(geometry.alpha[i]),
            "beta": float(geometry.beta[i]),
        }


    def geometry_predictor_budget(
        self, grid, predictor_geometry, predictor_matter, step, t, stage, stage_row
    ):
        """Audit actual explicit a/b/X rates versus the post-regularity predictor.

        Uses the explicit RHS returned by the live first geometry_stage_terms call.
        The raw Euler candidate and the algebraically regularized predictor are
        recorded separately. This method is observational and changes no state.
        """
        accepted = step.get("geometry_rhs0_input")
        rhs = step.get("geometry_rhs0_explicit")
        rhs0 = step.get("rhs0_record")
        if accepted is None or rhs is None or rhs0 is None:
            raise ValueError("missing accepted geometry or actual explicit RHS provenance")

        dt = float(step["dt"])
        U0 = rhs0["matter_radiation"]
        U1 = predictor_matter.radiation
        U0e, U0r = _array(U0.energy_t), _array(U0.momentum_r)
        U1e, U1r = _array(U1.energy_t), _array(U1.momentum_r)
        n = len(grid.centers)
        lo = max(0, n - OUTER_CELLS)
        cells = []
        for i in range(lo, n):
            a0, b0, X0 = (float(_array(getattr(accepted, name))[i]) for name in ("a", "b", "X"))
            a_raw = a0 + dt * float(rhs["a"][i])
            b_raw = b0 + dt * float(rhs["b"][i])
            X_raw = X0 + dt * float(rhs["X"][i])
            ap = float(_array(predictor_geometry.a)[i])
            bp = float(_array(predictor_geometry.b)[i])
            Xp = float(_array(predictor_geometry.X)[i])

            gamma0 = X0 * X0 / a0
            gamma_a_only = X0 * X0 / a_raw
            gamma_raw = X_raw * X_raw / a_raw
            gamma_pred = Xp * Xp / ap
            dgamma_a = gamma_a_only - gamma0
            dgamma_X = gamma_raw - gamma_a_only
            dgamma_projection = gamma_pred - gamma_raw
            dgamma_total = gamma_pred - gamma0
            gamma_closure = dgamma_a + dgamma_X + dgamma_projection - dgamma_total

            Aa0 = float(_array(accepted.Aa)[i])
            K0 = float(_array(accepted.K)[i])
            alpha0 = float(_array(accepted.alpha)[i])
            beta0 = float(_array(accepted.beta)[i])
            explicit_a = float(rhs["a"][i])
            explicit_X = float(rhs["X"][i])
            # On this zero-shift, KO_EPSILON=0 lane, these are the pinned
            # explicit BSSN equations. Residuals reveal any changed assumption.
            rhs_a_formula = -2.0 * alpha0 * a0 * Aa0
            rhs_X_formula = alpha0 * X0 * K0 / 3.0

            metric0 = gamma0
            metric_pred = gamma_pred
            C0 = float(U0e[i] - math.sqrt(metric0) * abs(U0r[i]))
            C1_same_metric = float(U1e[i] - math.sqrt(metric0) * abs(U1r[i]))
            C1_predictor_metric = float(U1e[i] - math.sqrt(metric_pred) * abs(U1r[i]))
            snapshot = next(
                (row for row in (stage_row or {}).get("outer_cell_snapshots", [])
                 if int(row.get("i", -1)) == i),
                None,
            )
            cells.append({
                "i": int(i),
                "r": float(grid.centers[i]),
                "alpha_accepted": alpha0,
                "beta_accepted": beta0,
                "Aa_accepted": Aa0,
                "K_accepted": K0,
                "a_accepted": a0,
                "b_accepted": b0,
                "X_accepted": X0,
                "rhs_explicit_a_actual": explicit_a,
                "rhs_explicit_b_actual": float(rhs["b"][i]),
                "rhs_explicit_X_actual": explicit_X,
                "rhs_explicit_a_formula_zero_shift": rhs_a_formula,
                "rhs_explicit_X_formula_zero_shift": rhs_X_formula,
                "rhs_a_formula_residual": explicit_a - rhs_a_formula,
                "rhs_X_formula_residual": explicit_X - rhs_X_formula,
                "a_predictor_raw_euler": a_raw,
                "b_predictor_raw_euler": b_raw,
                "X_predictor_raw_euler": X_raw,
                "a_predictor_after_regularity": ap,
                "b_predictor_after_regularity": bp,
                "X_predictor_after_regularity": Xp,
                "projection_delta_a": ap - a_raw,
                "projection_delta_b": bp - b_raw,
                "projection_delta_X": Xp - X_raw,
                "gamma_rr_inv_accepted": gamma0,
                "gamma_rr_inv_after_a_update": gamma_a_only,
                "gamma_rr_inv_after_raw_X_update": gamma_raw,
                "gamma_rr_inv_predictor_after_regularity": gamma_pred,
                "delta_gamma_rr_inv_from_a_update": dgamma_a,
                "delta_gamma_rr_inv_from_X_update": dgamma_X,
                "delta_gamma_rr_inv_from_regularity_projection": dgamma_projection,
                "delta_gamma_rr_inv_total": dgamma_total,
                "gamma_rr_inv_decomposition_closure_error": gamma_closure,
                "C_accepted_on_accepted_metric": C0,
                "C_predictor_conserved_state_on_accepted_metric": C1_same_metric,
                "C_predictor_conserved_state_on_predictor_metric": C1_predictor_metric,
                "delta_C_metric_on_same_predictor_matter": C1_predictor_metric - C1_same_metric,
                "stage_snapshot_C": None if snapshot is None else float(snapshot["C"]),
                "stage_snapshot_gamma_rr_inv": None if snapshot is None else float(snapshot["gamma_rr_inv"]),
            })
        return {
            "t": float(t),
            "t_attempted": float(t + dt),
            "step_index": int(step["step_index"]),
            "dt": dt,
            "stage": stage,
            "cell_band": [lo, n - 1],
            "rhs_source": "actual first adapter.geometry_stage_terms explicit RHS",
            "metric_identity": "gamma_rr_inv = X**2 / a",
            "pinned_zero_shift_explicit_identity": "a_t = -2 alpha a Aa; X_t = alpha X K / 3 (KO_EPSILON=0)",
            "cells": cells,
            "max_abs_gamma_rr_inv_decomposition_closure_error": max(
                (abs(cell["gamma_rr_inv_decomposition_closure_error"]) for cell in cells),
                default=None,
            ),
            "max_abs_rhs_a_formula_residual": max(
                (abs(cell["rhs_a_formula_residual"]) for cell in cells), default=None
            ),
            "max_abs_rhs_X_formula_residual": max(
                (abs(cell["rhs_X_formula_residual"]) for cell in cells), default=None
            ),
            "max_abs_projection_delta_X": max(
                (abs(cell["projection_delta_X"]) for cell in cells), default=None
            ),
        }


    def completed_state_budget(
        self, grid, completed_geometry, completed_matter,
        rhs0, rhs1, dt, t, step_index, stage,
    ):
        """Decompose the accepted-to-completed conservative radiation update.

        This is observational instrumentation. It uses the actual RHS0/RHS1
        arrays and captured source values. It does not alter equations, state,
        recovery, boundaries, or accept/reject decisions.
        """
        dt = float(dt)
        g0, g1 = rhs0["geometry"], rhs1["geometry"]
        U0, U1, Uc = (
            rhs0["matter_radiation"],
            rhs1["matter_radiation"],
            completed_matter.radiation,
        )
        U0e, U0r = _array(U0.energy_t), _array(U0.momentum_r)
        U1e, U1r = _array(U1.energy_t), _array(U1.momentum_r)
        Uce, Ucr = _array(Uc.energy_t), _array(Uc.momentum_r)
        r0e, r0r = _array(rhs0["rhs_energy"]), _array(rhs0["rhs_momentum"])
        r1e, r1r = _array(rhs1["rhs_energy"]), _array(rhs1["rhs_momentum"])
        arrays = (U0e, U0r, U1e, U1r, Uce, Ucr, r0e, r0r, r1e, r1r)
        if any(arr.shape != U0e.shape for arr in arrays):
            raise ValueError("completed radiation budget arrays have inconsistent shapes")

        c0 = self.cone_values(grid, g0, U0)
        c1 = self.cone_values(grid, g1, U1)
        cc = self.cone_values(grid, completed_geometry, Uc)
        cc_g0 = self.cone_values(grid, g0, Uc)
        cc_g1 = self.cone_values(grid, g1, Uc)
        sources0, sources1 = rhs0.get("sources", {}), rhs1.get("sources", {})
        n = len(U0e)
        lo = max(0, n - OUTER_CELLS)
        rows = []

        for i in range(lo, n):
            pred_e = U0e[i] + dt * r0e[i]
            pred_r = U0r[i] + dt * r0r[i]
            comp_e = U0e[i] + 0.5 * dt * (r0e[i] + r1e[i])
            comp_r = U0r[i] + 0.5 * dt * (r0r[i] + r1r[i])
            row = {
                "i": int(i), "r": float(grid.centers[i]),
                "source_capture_missing": False,
                "U_E_accepted": float(U0e[i]), "U_r_accepted": float(U0r[i]),
                "U_E_predictor_actual": float(U1e[i]), "U_r_predictor_actual": float(U1r[i]),
                "U_E_predictor_reconstructed_from_rhs0": float(pred_e),
                "U_r_predictor_reconstructed_from_rhs0": float(pred_r),
                "predictor_U_E_actual_minus_rhs0": float(U1e[i] - pred_e),
                "predictor_U_r_actual_minus_rhs0": float(U1r[i] - pred_r),
                "U_E_completed_actual": float(Uce[i]), "U_r_completed_actual": float(Ucr[i]),
                "U_E_completed_reconstructed_from_trapezoid": float(comp_e),
                "U_r_completed_reconstructed_from_trapezoid": float(comp_r),
                "completed_U_E_actual_minus_trapezoid": float(Uce[i] - comp_e),
                "completed_U_r_actual_minus_trapezoid": float(Ucr[i] - comp_r),
                "rhs0_U_E_rate": float(r0e[i]), "rhs0_U_r_rate": float(r0r[i]),
                "rhs1_U_E_rate": float(r1e[i]), "rhs1_U_r_rate": float(r1r[i]),
                "C_accepted": float(c0["C"][i]),
                "C_predictor_on_rhs1_metric": float(c1["C"][i]),
                "C_completed_on_accepted_metric": float(cc_g0["C"][i]),
                "C_completed_on_rhs1_metric": float(cc_g1["C"][i]),
                "C_completed_final_metric": float(cc["C"][i]),
                "accepted_metric": self._budget_metric_provenance(g0, c0, i),
                "rhs1_stage_metric": self._budget_metric_provenance(g1, c1, i),
                "completed_metric": self._budget_metric_provenance(completed_geometry, cc, i),
            }
            s0, s1 = sources0.get(i), sources1.get(i)
            if s0 is None or s1 is None:
                row["source_capture_missing"] = True
                row["missing_source_stages"] = [
                    label for label, value in (("rhs0", s0), ("rhs1", s1))
                    if value is None
                ]
                row["closure_error"] = None
                rows.append(row)
                continue

            s0e, s0r = float(s0[0]), float(s0[1])
            s1e, s1r = float(s1[0]), float(s1[1])
            f0e, f0r = float(r0e[i] - s0e), float(r0r[i] - s0r)
            f1e, f1r = float(r1e[i] - s1e), float(r1r[i] - s1r)
            dfe, dfr = 0.5 * dt * (f0e + f1e), 0.5 * dt * (f0r + f1r)
            dse, dsr = 0.5 * dt * (s0e + s1e), 0.5 * dt * (s0r + s1r)
            Ufe, Ufr = float(U0e[i] + dfe), float(U0r[i] + dfr)
            Use, Usr = float(Ufe + dse), float(Ufr + dsr)
            root_ginv0 = math.sqrt(float(c0["gamma_rr_inv"][i]))
            C_flux = Ufe - root_ginv0 * abs(Ufr)
            C_sources = Use - root_ginv0 * abs(Usr)
            dC_flux = C_flux - float(c0["C"][i])
            dC_sources = C_sources - C_flux
            dC_metric_rhs1 = float(cc_g1["C"][i] - cc_g0["C"][i])
            dC_metric_completed = float(cc["C"][i] - cc_g1["C"][i])
            dC_total = float(cc["C"][i] - c0["C"][i])
            closure = dC_flux + dC_sources + dC_metric_rhs1 + dC_metric_completed - dC_total
            dC_matter = float(cc_g1["C"][i] - c1["C"][i])
            dC_pred_completed = float(cc["C"][i] - c1["C"][i])
            pred_closure = dC_matter + dC_metric_completed - dC_pred_completed
            row.update({
                "rhs0_source_U_E_rate": s0e, "rhs0_source_U_r_rate": s0r,
                "rhs1_source_U_E_rate": s1e, "rhs1_source_U_r_rate": s1r,
                "rhs0_flux_U_E_rate": f0e, "rhs0_flux_U_r_rate": f0r,
                "rhs1_flux_U_E_rate": f1e, "rhs1_flux_U_r_rate": f1r,
                "delta_flux_U_E_trapezoidal": float(dfe),
                "delta_flux_U_r_trapezoidal": float(dfr),
                "delta_source_U_E_trapezoidal": float(dse),
                "delta_source_U_r_trapezoidal": float(dsr),
                "U_E_after_flux_on_accepted_metric": Ufe,
                "U_r_after_flux_on_accepted_metric": Ufr,
                "U_E_after_sources_on_accepted_metric": Use,
                "U_r_after_sources_on_accepted_metric": Usr,
                "C_after_flux_on_accepted_metric": float(C_flux),
                "C_after_sources_on_accepted_metric": float(C_sources),
                "delta_C_flux_trapezoidal": float(dC_flux),
                "delta_C_sources_trapezoidal": float(dC_sources),
                "delta_C_metric_accepted_to_rhs1": dC_metric_rhs1,
                "delta_C_metric_rhs1_to_completed": dC_metric_completed,
                "delta_C_total_accepted_to_completed": dC_total,
                "closure_error": float(closure),
                "delta_C_trapezoid_matter_correction_on_rhs1_metric": dC_matter,
                "delta_C_predictor_to_completed": dC_pred_completed,
                "predictor_to_completed_closure_error": float(pred_closure),
                "completed_conservative_energy_split_error": float(Use - Uce[i]),
                "completed_conservative_momentum_split_error": float(Usr - Ucr[i]),
            })
            rows.append(row)

        valid = [row for row in rows if row.get("closure_error") is not None]
        missing = [int(row["i"]) for row in rows if row["source_capture_missing"]]
        def max_abs(key, rows_for_key=valid):
            return max((abs(row[key]) for row in rows_for_key), default=None)
        return {
            "t_accepted": float(t), "t_completed_candidate": float(t + dt),
            "step_index": int(step_index), "stage": stage, "dt": dt,
            "cell_band": [lo, n - 1],
            "budget_order": [
                "accepted conservative state on accepted spatial metric",
                "trapezoidal flux-only update on accepted spatial metric",
                "trapezoidal source update on accepted spatial metric",
                "same completed conservative state on RHS1 stage spatial metric",
                "same completed conservative state on completed candidate spatial metric",
            ],
            "rhs_metric_provenance": {
                "rhs0_stage": self._budget_metric_provenance(g0, c0, min(n - 1, lo)),
                "rhs1_stage": self._budget_metric_provenance(g1, c1, min(n - 1, lo)),
            },
            "source_capture_missing_cells": missing, "cells": rows,
            "max_abs_budget_closure_error": max_abs("closure_error"),
            "max_abs_predictor_U_E_reconstruction_error": max_abs(
                "predictor_U_E_actual_minus_rhs0", rows
            ),
            "max_abs_predictor_U_r_reconstruction_error": max_abs(
                "predictor_U_r_actual_minus_rhs0", rows
            ),
            "max_abs_completed_U_E_reconstruction_error": max_abs(
                "completed_U_E_actual_minus_trapezoid", rows
            ),
            "max_abs_completed_U_r_reconstruction_error": max_abs(
                "completed_U_r_actual_minus_trapezoid", rows
            ),
            "max_abs_source_flux_split_energy_error": max_abs("completed_conservative_energy_split_error"),
            "max_abs_source_flux_split_momentum_error": max_abs("completed_conservative_momentum_split_error"),
            "max_abs_predictor_to_completed_closure_error": max_abs("predictor_to_completed_closure_error"),
        }


    def _patch_production_kernel(self):
        kernel_cls = production_kernel_module.V55ProductionKernel
        original_solve = kernel_cls._solve_lapse
        original_rhs = kernel_cls._rhs
        original_step = kernel_cls.step
        original_evolve = matter_rhs_module.evolve_species
        original_radiation_source = matter_system_module._radiation_source
        original_boundary = kernel_cls._apply_outer_light_boundary
        original_assemble_total = production_kernel_module.assemble_total_stress_energy
        original_geometry_stage_terms = production_kernel_module.adapter.geometry_stage_terms
        trace = self

        def note_failure_identity(step, stage, caller, grid, exc):
            if step is None:
                return
            step["failing_stage"] = stage
            step["failing_caller"] = caller
            match = re.search(r"\bcell i=(\d+)", str(exc))
            if match is None:
                return
            index = int(match.group(1))
            step["failing_grid_index"] = index
            centers = _array(grid.centers)
            if 0 <= index < len(centers):
                step["failing_radius"] = float(centers[index])

        def traced_geometry_stage_terms(*args, **kwargs):
            result = original_geometry_stage_terms(*args, **kwargs)
            step = trace.current_step
            if step is not None:
                idx = int(step.get("geometry_terms_count", 0))
                step["geometry_terms_count"] = idx + 1
                if idx == 0:
                    geometry = args[1]
                    step["geometry_rhs0_input"] = geometry.copy()
                    step["geometry_rhs0_explicit"] = {
                        name: _array(result["explicit"][name]).copy()
                        for name in ("a", "b", "X")
                    }
            return result

        def traced_assemble_total(*args, **kwargs):
            total = original_assemble_total(*args, **kwargs)
            context = trace.active_boundary_context
            if context is not None and context.get("pre") is None:
                try:
                    rho = _array(total.rho).copy()
                    momentum = _array(total.j).copy()
                    context["rho"] = rho
                    context["momentum"] = momentum
                    context["pre"] = _pminus_boundary_terms(
                        context["grid"], context["geometry"], rho, momentum
                    )
                except Exception as exc:
                    detail = f"{type(exc).__name__}: {exc}"
                    context["capture_error"] = detail
                    trace.instrumentation_errors.append({
                        "where": "pminus_boundary_pre_capture",
                        "stage": context.get("stage"),
                        "t": context.get("t"),
                        "error": detail,
                    })
            return total

        def traced_radiation_source(*args, **kwargs):
            result = original_radiation_source(*args, **kwargs)
            active = trace.active_radiation_evolve
            if active is not None:
                i = active["source_index"]
                if i >= len(active["metric_r"]) - OUTER_CELLS:
                    active["sources"][i] = (float(result[0]), float(result[1]))
                active["source_index"] += 1
            return result

        def traced_boundary(grid, geometry, scalars, matter,
                            radiation_recovery_metric=None):
            step = trace.current_step
            if step is None:
                if radiation_recovery_metric is None:
                    return original_boundary(grid, geometry, scalars, matter)
                return original_boundary(
                    grid, geometry, scalars, matter,
                    radiation_recovery_metric=radiation_recovery_metric,
                )

            t = float(step["t0"])
            tau = float(step["tau0"])
            idx = int(step.get("boundary_count", 0))
            step["boundary_count"] = idx + 1
            labels = (
                "primary_predictor_outer_light_boundary",
                "completed_state_outer_light_boundary",
            )
            stage = labels[idx] if idx < len(labels) else f"outer_light_boundary_call_{idx}"
            stage_row = trace.record_stage(
                stage, grid, geometry, matter, t, tau,
                trace.step_index, step["dt"], force=(idx == 1),
            )
            # Observe the completed candidate before normal boundary recovery.
            # This block does not modify matter/geometry or affect admissibility.
            if idx == 1 and trace.latest_state is not None and stage_row is not None:
                if step.get("rhs0_record") is None or step.get("rhs1_record") is None:
                    error = "completed-state budget lacks rhs0/rhs1 provenance"
                    stage_row["completed_state_budget_error"] = error
                    trace.instrumentation_errors.append({
                        "where": "completed_state_budget",
                        "t": t, "step_index": trace.step_index, "error": error,
                    })
                else:
                    try:
                        budget = trace.completed_state_budget(
                            grid, geometry, matter, step["rhs0_record"],
                            step["rhs1_record"], step["dt"], t,
                            trace.step_index, stage,
                        )
                        trace.completed_state_budgets.append(budget)
                        stage_row["completed_state_budget_index"] = len(trace.completed_state_budgets) - 1
                        if budget["source_capture_missing_cells"]:
                            error = "completed-state source capture missing at cells " + ",".join(
                                map(str, budget["source_capture_missing_cells"])
                            )
                            stage_row["completed_state_budget_error"] = error
                            trace.instrumentation_errors.append({
                                "where": "completed_state_budget",
                                "t": t, "step_index": trace.step_index, "error": error,
                            })
                    except Exception as exc:
                        error = f"{type(exc).__name__}: {exc}"
                        stage_row["completed_state_budget_error"] = error
                        trace.instrumentation_errors.append({
                            "where": "completed_state_budget",
                            "t": t, "step_index": trace.step_index, "error": error,
                        })
            try:
                if radiation_recovery_metric is None:
                    return original_boundary(grid, geometry, scalars, matter)
                return original_boundary(
                    grid, geometry, scalars, matter,
                    radiation_recovery_metric=radiation_recovery_metric,
                )
            except Exception as exc:
                detail = f"{type(exc).__name__}: {exc}"
                trace.record_stage(
                    stage + "_failure", grid, geometry, matter, t, tau,
                    trace.step_index, step["dt"], force=True, exception=detail,
                )
                note_failure_identity(
                    step, stage,
                    "engine.production_kernel.V55ProductionKernel._apply_outer_light_boundary",
                    grid, exc,
                )
                raise

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
            if step is not None and idx == 1:
                if step.get("geometry_rhs0_explicit") is None:
                    error = "missing actual RHS0 geometry explicit terms"
                    if row is not None:
                        row["geometry_predictor_budget_error"] = error
                    trace.instrumentation_errors.append({
                        "where": "geometry_predictor_budget",
                        "t": t, "step_index": trace.step_index, "error": error,
                    })
                else:
                    try:
                        budget = trace.geometry_predictor_budget(
                            grid, geometry, matter, step, t, stage, row
                        )
                        trace.geometry_predictor_budgets.append(budget)
                        if row is not None:
                            row["geometry_predictor_budget_index"] = len(trace.geometry_predictor_budgets) - 1
                    except Exception as exc:
                        error = f"{type(exc).__name__}: {exc}"
                        if row is not None:
                            row["geometry_predictor_budget_error"] = error
                        trace.instrumentation_errors.append({
                            "where": "geometry_predictor_budget",
                            "t": t, "step_index": trace.step_index, "error": error,
                        })
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
                note_failure_identity(
                    step, stage,
                    "engine.production_kernel.V55ProductionKernel._solve_lapse",
                    grid, exc,
                )
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
            try:
                if radiation_recovery_metric is None:
                    out = original_rhs(kernel_self, state)
                else:
                    out = original_rhs(
                        kernel_self, state,
                        radiation_recovery_metric=radiation_recovery_metric,
                    )
            except Exception as exc:
                detail = f"{type(exc).__name__}: {exc}"
                trace.record_stage(
                    stage + "_failure",
                    state.grid, state.geometry, state.matter, t, tau,
                    trace.step_index, step["dt"] if step is not None else None,
                    force=True, exception=detail,
                )
                note_failure_identity(
                    step, stage,
                    "engine.production_kernel.V55ProductionKernel._rhs",
                    state.grid, exc,
                )
                raise
            srhs, mrhs, md = out
            if step is not None and rhs_idx in (0, 1):
                rhs_record = {
                    "geometry": state.geometry.copy(),
                    "matter_radiation": state.matter.radiation.copy(),
                    "rhs_energy": _array(mrhs["radiation"].energy_t).copy(),
                    "rhs_momentum": _array(mrhs["radiation"].momentum_r).copy(),
                    "sources": dict(trace.last_radiation_sources),
                }
                step["rhs0_record" if rhs_idx == 0 else "rhs1_record"] = rhs_record
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
                "boundary_count": 0,
                "rhs0_record": None,
                "rhs1_record": None,
                "geometry_terms_count": 0,
                "geometry_rhs0_input": None,
                "geometry_rhs0_explicit": None,
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
                    "stage_name": trace.current_step.get("failing_stage"),
                    "caller": trace.current_step.get("failing_caller"),
                    "grid_index": trace.current_step.get("failing_grid_index"),
                    "radius": trace.current_step.get("failing_radius"),
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
        kernel_cls._apply_outer_light_boundary = staticmethod(traced_boundary)
        production_kernel_module.adapter.geometry_stage_terms = traced_geometry_stage_terms
        kernel_cls.step = traced_step

        self._restore_hooks = (
            kernel_cls, original_solve, original_rhs, original_step,
            original_evolve, original_radiation_source, original_boundary,
        )
        self._restore_geometry_stage_terms = original_geometry_stage_terms

    def prune_tail(self, cutoff_t):
        self.stage_snapshots = [
            row for row in self.stage_snapshots
            if float(row.get("t", -math.inf)) >= cutoff_t
        ]
        self.predictor_budgets = [
            row for row in self.predictor_budgets
            if float(row.get("t", -math.inf)) >= cutoff_t
        ]
        self.completed_state_budgets = [
            row for row in self.completed_state_budgets
            if float(row.get("t_completed_candidate", -math.inf)) >= cutoff_t
        ]
        self.geometry_predictor_budgets = [
            row for row in self.geometry_predictor_budgets
            if float(row.get("t_attempted", -math.inf)) >= cutoff_t
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
        dt_proposal = dt_nominal
        start = time.time()
        last_report_step = -1
        status = "completed_without_reproducing_expected_failure"
        target_steps = int(math.ceil(TARGET_TIME / dt_nominal)) + 5

        def record_rejected_attempt(record):
            # traced_step stores stage/cell identity in trace.failure. The
            # rejected candidate is discarded; every retry starts from state.
            if self.failure is not None:
                record.update(self.failure)
            record["rejected_for_radiation_admissibility"] = True
            self.rejected_step_attempts.append(dict(record))
            self.failure = None

        if TARGET_TAU is not None and (
            not math.isfinite(float(TARGET_TAU)) or float(TARGET_TAU) <= 0.0
        ):
            raise ValueError("TARGET_TAU must be finite and positive when set")

        while state.t < TARGET_TIME and (
            self.use_radiation_admissibility_retry or self.step_index < target_steps
        ) and (TARGET_TAU is None or state.tau < TARGET_TAU):
            dt = min(dt_proposal, TARGET_TIME - state.t)
            try:
                if self.use_radiation_admissibility_retry:
                    state, dt_used, _ = advance_with_radiation_admissibility_retries(
                        kernel, state, dt, on_reject=record_rejected_attempt
                    )
                    dt_proposal = min(dt_nominal, 2.0 * dt_used)
                else:
                    state = kernel.step(state, dt)
                    dt_used = dt
                    dt_proposal = dt_nominal
                state.geometry.assert_finite_positive()
            except RadiationStepSizeUnderflow as exc:
                status = "radiation_admissibility_step_underflow"
                if self.rejected_step_attempts:
                    self.failure = dict(self.rejected_step_attempts[-1])
                    self.failure["underflow_message"] = str(exc)
                else:
                    self.failure = {
                        "t_accepted": float(state.t),
                        "tau_accepted": float(state.tau),
                        "step_index": int(self.step_index),
                        "dt": float(dt),
                        "t_attempted": float(state.t + dt),
                        "stage": "radiation_admissibility_step_underflow",
                        "exception_class": type(exc).__name__,
                        "exception_message": str(exc),
                        "wall_elapsed_seconds": float(time.time() - self.started_wall),
                    }
                break
            except (FloatingPointError, ValueError) as exc:
                # Only typed radiation admissibility failures are retried.
                # Other numerical or geometry failures remain visible here.
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
                    "dt_accepted": dt_used,
                    "rejected_attempt_count": len(self.rejected_step_attempts),
                    "elapsed_seconds": time.time() - start,
                }), flush=True)
                last_report_step = self.step_index
        else:
            completed_status = _completion_status(
                float(state.t), float(state.tau), float(TARGET_TIME), TARGET_TAU
            )
            if completed_status is not None:
                status = completed_status
            elif not self.use_radiation_admissibility_retry and self.step_index >= target_steps:
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
            "engine/adaptive_step.py",
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
                    "target_proper_time": TARGET_TAU,
                    "tail_window": TAIL_WINDOW,
                    "outer_cells_traced": OUTER_CELLS,
                    "amplitude": 0.01,
                    "width": 7.0,
                    "D_amplitude": 1.0e-10,
                    "include_radiation": True,
                    "use_accepted_metric_for_predictor_radiation_recovery":
                        self.use_accepted_metric_for_predictor_radiation_recovery,
                    "radiation_admissibility_retry_enabled":
                        self.use_radiation_admissibility_retry,
                    "outer_boundary": "existing production light-constraint boundary and radiation outer closure",
                },
            },
            "failure": self.failure,
            "rejected_step_attempts": self.rejected_step_attempts,
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
            "completed_state_budgets": self.completed_state_budgets,
            "geometry_predictor_budgets": self.geometry_predictor_budgets,
            "instrumentation_errors": self.instrumentation_errors,
            "wall_elapsed_seconds": float(time.time() - self.started_wall),
            "notes": [
                "The conservative radiation variables remain unmodified by the trace.",
                (
                    "The trace stops at the first accepted state reaching the requested central proper time."
                    if TARGET_TAU is not None
                    else "The trace stops at the configured coordinate-time endpoint."
                ),
                "The source is captured directly from the existing production radiation source routine.",
                "Flux transport rate is reconstructed as the existing total radiation RHS minus that captured source; evolve_species sums only those two components for radiation.",
                "The cone budget is ordered and uses the accepted spatial metric for the accepted, flux-only, and source-updated stages, then the predictor metric for the final predictor margin.",
                "Completed-state budgets reconstruct the trapezoidal radiation update from actual RHS0/RHS1 evaluations and close the accepted-to-completed cone-margin change across stage-specific metrics.",
                "Geometry predictor budgets record the actual first geometry-stage explicit RHS and compare its raw Euler candidate with the predictor geometry after algebraic regularity projection.",
                (
                    "Extended-B mode uses the accepted start-of-step metric only for predictor-stage radiation primitive inversion; stage metrics still enter the fluxes, geometric sources, stress projections, and CMC spatial operator."
                    if self.use_accepted_metric_for_predictor_radiation_recovery
                    else "Extended-B mode is disabled; the original predictor-stage radiation recovery metric is used."
                ),
                (
                    "Radiation admissibility retries change numerical timestep only; they do not clip conservative variables or alter field equations."
                    if self.use_radiation_admissibility_retry
                    else "No clipping, floor, damping, physical source insertion, or timestep adjustment is introduced."
                ),
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
            "completed_state_budget_count": len(self.completed_state_budgets),
            "completed_state_budget_error_count": sum(
                1 for row in self.stage_snapshots if row.get("completed_state_budget_error")
            ),
            "geometry_predictor_budget_count": len(self.geometry_predictor_budgets),
            "geometry_predictor_budget_error_count": sum(
                1 for row in self.stage_snapshots if row.get("geometry_predictor_budget_error")
            ),
            "progress_sample_count": len(self.progress_samples),
            "rejected_step_attempt_count": len(self.rejected_step_attempts),
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
    use_admissibility_retry = (
        os.environ.get("RADIATION_ADMISSIBILITY_RETRY", "0") == "1"
    )
    trace = RadiationStageTrace(
        output,
        use_accepted_metric_for_predictor_radiation_recovery=use_accepted_recovery,
        use_radiation_admissibility_retry=use_admissibility_retry,
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
