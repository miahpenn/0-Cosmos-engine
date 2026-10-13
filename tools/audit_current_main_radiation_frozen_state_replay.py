#!/usr/bin/env python3
"""Replay one checksum-pinned archived radiation boundary state through today's main campaign runner.

This is a bounded, diagnostics-only provenance/retry/admission check. It does
not replay the lead-in trajectory, alter physical variables, or claim scientific
admission from a short probe.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import platform
import subprocess
import shutil
import time
from types import SimpleNamespace
from pathlib import Path

import numpy as np

import engine.matter_rhs as matter_rhs_module
import engine.production_kernel as production_kernel_module
import engine.v55_pirk_adapter as adapter_module
from engine import matter_system as matter_system_module
from engine.campaign import CampaignConfig
from engine.production_kernel import V55ProductionKernel
from engine.run_production import run_campaign
from engine.scalar_system import ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import metric_slice_from_q

BASE_MAIN_COMMIT = "29d828b4d72b017c12a20b6dc90549f2511f025c"
ARCHIVED_SOURCE_COMMIT = "cdcec8282bff0443f2fe4674b15929822ff38736"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
EXPECTED_NPZ_SHA256 = "ccd1fa0fbabe65c51c3dd7ac2c9227d6d25aa680c169eef1afd3fca4326294cf"
SOURCE_RUN = 38322337250
SOURCE_ARTIFACT = 11778409394
RESOLUTION = 80
R_MAX = 80.0
REQUESTED_DT = 0.03
ACCEPTED_FIRST_RETRY_DT = 0.015
GEOMETRY_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")
SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES_FIELDS = ("rest", "energy_t", "momentum_r")
SPECIES_NAMES = ("dark_matter", "baryons", "radiation")
EXPECTED_DIAGNOSTIC_PATHS = {
    "tools/audit_current_main_radiation_frozen_state_replay.py",
    ".github/workflows/current_main_radiation_frozen_state_replay.yml",
    "tests/test_radiation_replay_rhs_capture_identity.py",
    "tests/test_radiation_l3_stress_projection.py",
    "tests/test_fluid_projection_kappa_contract.py",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_value(*args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], text=True).strip()
    except Exception:
        return "unavailable"


def strict_json(path: Path):
    def reject_constant(value: str):
        raise ValueError(f"non-standard JSON constant {value} in {path}")
    return json.loads(path.read_text(), parse_constant=reject_constant)


def numeric_state_digest(state) -> str:
    """Fingerprint every archived numerical profile and physical time label."""
    digest = hashlib.sha256()
    arrays = [("r", state.grid.centers)]
    arrays.extend((f"geometry.{name}", getattr(state.geometry, name))
                  for name in GEOMETRY_FIELDS)
    arrays.extend((f"scalars.{name}", getattr(state.scalars, name))
                  for name in SCALAR_FIELDS)
    for species_name in SPECIES_NAMES:
        species = getattr(state.matter, species_name)
        arrays.extend((f"matter.{species_name}.{name}", getattr(species, name))
                      for name in SPECIES_FIELDS)
    for name, values in arrays:
        array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
        digest.update(name.encode("utf-8") + b"\0")
        digest.update(str(array.shape).encode("ascii") + b"\0")
        digest.update(array.tobytes())
    for name in ("t", "tau", "e_folds"):
        digest.update(name.encode("ascii") + b"\0")
        digest.update(float(getattr(state, name)).hex().encode("ascii") + b"\0")
    return digest.hexdigest()


def l3_matter_terms_from_projection(alpha, rho, pr, pt) -> dict[str, np.ndarray]:
    """Compute explicit L3 matter terms from independently assembled projections."""
    alpha = np.asarray(alpha, dtype=float)
    rho = np.asarray(rho, dtype=float)
    pr = np.asarray(pr, dtype=float)
    pt = np.asarray(pt, dtype=float)
    return {
        "Aa": -(16.0 * np.pi / 3.0) * alpha * (pr - pt),
        "K": 4.0 * np.pi * alpha * (rho + pr + 2.0 * pt),
    }


def source_increment_roundoff_witness(total, vacuum, expected):
    """Compare a source increment with an independent formula using operand-scaled roundoff.

    No absolute floor or fitted tolerance is used. Returns the signed residual,
    a per-cell floating-point roundoff bound, all-cell pass/fail, and maximum
    residual-to-bound ratio.
    """
    total = np.asarray(total, dtype=float)
    vacuum = np.asarray(vacuum, dtype=float)
    expected = np.asarray(expected, dtype=float)
    residual = (total - vacuum) - expected
    operand_scale = np.abs(total) + np.abs(vacuum) + np.abs(expected)
    bound = 32.0 * np.finfo(np.float64).eps * operand_scale
    residual_abs = np.abs(residual)
    within = np.isfinite(residual_abs) & np.isfinite(bound) & (residual_abs <= bound)
    ratios = np.zeros_like(residual_abs)
    nonzero_bound = bound > 0.0
    ratios[nonzero_bound] = residual_abs[nonzero_bound] / bound[nonzero_bound]
    impossible = (~nonzero_bound) & (residual_abs > 0.0)
    ratios[impossible] = np.finfo(np.float64).max
    ratios[~np.isfinite(residual_abs) | ~np.isfinite(bound)] = np.finfo(np.float64).max
    return residual, bound, bool(np.all(within)), float(np.max(ratios))


def load_frozen_state(npz_path: Path):
    """Build current-main state containers, then replace all numeric profiles."""
    with np.load(npz_path, allow_pickle=False) as data:
        kernel = V55ProductionKernel()
        state = kernel.initialize(
            resolution=RESOLUTION,
            r_max=R_MAX,
            amplitude=0.01,
            width=7.0,
            D_amplitude=1.0e-10,
            include_radiation=True,
        )
        if not np.array_equal(np.asarray(state.grid.centers), data["accepted_r"]):
            raise RuntimeError("frozen-state radial grid differs from N=80/r_max=80")
        for name in GEOMETRY_FIELDS:
            key = f"accepted_geometry_{name}"
            if key not in data:
                raise RuntimeError(f"accepted-state snapshot lacks {key}")
            setattr(state.geometry, name, np.asarray(data[key], dtype=float).copy())
        state.scalars = ScalarFields(*(
            np.asarray(data[f"accepted_scalar_{name}"], dtype=float).copy()
            for name in SCALAR_FIELDS
        ))
        for species_name in SPECIES_NAMES:
            species = getattr(state.matter, species_name)
            for field in SPECIES_FIELDS:
                key = f"accepted_matter_{species_name}_{field}"
                if key not in data:
                    raise RuntimeError(f"accepted-state snapshot lacks {key}")
                setattr(species, field, np.asarray(data[key], dtype=float).copy())
        state.t = float(data["accepted_t"])
        state.tau = float(data["accepted_tau"])
        state.e_folds = float(data["accepted_e_folds"])
    if abs(float(state.t) - 47.130000000001026) > 1.0e-10:
        raise RuntimeError(f"unexpected frozen-state time {state.t:.17g}")
    return kernel, state


def cone_witness(state) -> dict:
    metric = matter_system_module._metric_arrays(
        metric_slice_from_q(state.grid, state.geometry)
    )
    gamma_inv = np.asarray([item.gamma_rr_inv for item in metric], dtype=float)
    energy = np.asarray(state.matter.radiation.energy_t, dtype=float)
    momentum = np.asarray(state.matter.radiation.momentum_r, dtype=float)
    c = energy - np.sqrt(gamma_inv) * np.abs(momentum)
    ratio = np.full_like(energy, np.nan)
    good = energy > 0.0
    ratio[good] = np.sqrt(gamma_inv[good]) * np.abs(momentum[good]) / energy[good]
    i_ratio = int(np.nanargmax(ratio))
    i_c = int(np.nanargmin(c))
    return {
        "max_abs_S_over_E": float(ratio[i_ratio]),
        "max_ratio_cell": i_ratio,
        "max_ratio_radius": float(state.grid.centers[i_ratio]),
        "min_C": float(c[i_c]),
        "min_C_cell": i_c,
        "min_C_radius": float(state.grid.centers[i_c]),
        "all_finite": bool(
            np.all(np.isfinite(energy))
            and np.all(np.isfinite(momentum))
            and np.all(np.isfinite(gamma_inv))
            and np.all(np.isfinite(c))
            and np.all(np.isfinite(ratio))
        ),
    }




class RadiationStageCapture:
    """Read-only tap on real RHS/validator calls during the bounded replay."""

    def __init__(self, output_dir: Path, start_time: float):
        self.output_dir = output_dir
        self.start_time = float(start_time)
        self.arrays: dict[str, np.ndarray] = {}
        self.state_records: dict[str, dict] = {}
        self.rhs_records: dict[tuple[str, int, int], dict] = {}
        self.rhs_geometry: dict[tuple[str, int, int], object] = {}
        self.failed_trials: list[dict] = []
        self.active_state = None
        self.active_time = None
        self.active_dt = None
        self.rhs_counter = 0
        self.current_rhs_record = None
        self.originals = {}
        self._next_state_id = 0
        self._next_rhs_id = 0
        self._next_failed_id = 0
        self.accepted_step_outputs: list[dict] = []
        self.kernel = None
        self.radiation_rhs_component_capture_calls = 0
        self.geometry_records: dict[tuple[str, int, int], dict] = {}
        self.geometry_source_records: list[dict] = []
        self.geometry_counter = 0
        self.step_attempt_counter = 0
        self.active_attempt_id = None
        self._next_geometry_id = 0
        self.boundary_records: list[dict] = []
        self._next_boundary_id = 0

    @staticmethod
    def _copy_array(value):
        return np.ascontiguousarray(np.asarray(value, dtype=np.float64)).copy()

    def _save_state(self, prefix: str, state):
        self.arrays[f"{prefix}__r"] = self._copy_array(state.grid.centers)
        for name in GEOMETRY_FIELDS:
            self.arrays[f"{prefix}__geometry__{name}"] = self._copy_array(
                getattr(state.geometry, name)
            )
        for name in SCALAR_FIELDS:
            self.arrays[f"{prefix}__scalars__{name}"] = self._copy_array(
                getattr(state.scalars, name)
            )
        for species_name in SPECIES_NAMES:
            species = getattr(state.matter, species_name)
            for field in SPECIES_FIELDS:
                self.arrays[f"{prefix}__matter__{species_name}__{field}"] = (
                    self._copy_array(getattr(species, field))
                )

    def _save_metric(self, prefix: str, metric):
        for name in ("r", "a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B"):
            self.arrays[f"{prefix}__metric__{name}"] = self._copy_array(
                getattr(metric, name)
            )

    def _save_conserved(self, prefix: str, conserved):
        for field in SPECIES_FIELDS:
            self.arrays[f"{prefix}__radiation__{field}"] = self._copy_array(
                getattr(conserved, field)
            )

    @staticmethod
    def _geometry_metric(grid, geometry):
        return matter_system_module.BSSNMetricSlice(
            r=np.asarray(grid.centers, dtype=float).copy(),
            a=np.asarray(geometry.a, dtype=float).copy(),
            b=np.asarray(geometry.b, dtype=float).copy(),
            X=np.asarray(geometry.X, dtype=float).copy(),
            alpha=np.asarray(geometry.alpha, dtype=float).copy(),
            beta=np.asarray(geometry.beta, dtype=float).copy(),
            Aa=np.asarray(geometry.Aa, dtype=float).copy(),
            K=np.asarray(geometry.K, dtype=float).copy(),
            Lambda=np.asarray(geometry.Lambda, dtype=float).copy(),
            B=np.asarray(geometry.B, dtype=float).copy(),
        )

    @staticmethod
    def _ratio_profile(metric, conserved):
        # This is exactly the cone ratio used by validate_species_admissibility:
        # sqrt(gamma^rr) * |U_r| / U_E; sqrt(gamma) cancels.
        a = np.asarray(metric.a, dtype=float)
        X = np.asarray(metric.X, dtype=float)
        energy_t = np.asarray(conserved.energy_t, dtype=float)
        momentum_r = np.asarray(conserved.momentum_r, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.sqrt(X * X / a) * np.abs(momentum_r) / energy_t
        margin = energy_t - np.sqrt(X * X / a) * np.abs(momentum_r)
        physical_energy = energy_t / np.asarray([
            m.sqrt_gamma for m in matter_system_module._metric_arrays(metric)
        ], dtype=float)
        physical_abs_momentum = (
            np.sqrt(X * X / a) * np.abs(momentum_r) /
            np.asarray([m.sqrt_gamma for m in matter_system_module._metric_arrays(metric)], dtype=float)
        )
        return ratio, margin, physical_energy, physical_abs_momentum

    def _profile_summary(self, metric, conserved):
        ratio, margin, energy, abs_momentum = self._ratio_profile(metric, conserved)
        lo = max(0, ratio.size - 5)
        outer = []
        for i in range(lo, ratio.size):
            outer.append({
                "i": int(i),
                "r": float(metric.r[i]),
                "U_E": float(conserved.energy_t[i]),
                "U_r": float(conserved.momentum_r[i]),
                "a": float(metric.a[i]),
                "X": float(metric.X[i]),
                "gamma_rr_inv": float(metric.X[i] ** 2 / metric.a[i]),
                "E": float(energy[i]),
                "abs_S": float(abs_momentum[i]),
                "ratio_absS_over_E": float(ratio[i]),
                "conserved_cone_margin": float(margin[i]),
            })
        finite_ratio = np.isfinite(ratio)
        peak = int(np.nanargmax(ratio)) if np.any(finite_ratio) else None
        return {
            "outer_cells": outer,
            "max_ratio_cell": peak,
            "max_ratio": None if peak is None else float(ratio[peak]),
            "min_conserved_cone_margin_cell": int(np.nanargmin(margin)),
            "min_conserved_cone_margin": float(np.nanmin(margin)),
        }

    def begin_step(self, state, dt):
        self.active_state = state
        self.active_time = float(state.t)
        self.active_dt = float(dt)
        self.step_attempt_counter += 1
        self.active_attempt_id = self.step_attempt_counter
        self.rhs_counter = 0
        self.geometry_counter = 0
        self.current_rhs_record = None
        time_key = self.active_time.hex()
        if self.active_time >= self.start_time - 1.0e-10 and time_key not in self.state_records:
            prefix = f"accepted_{self._next_state_id:03d}"
            self._next_state_id += 1
            self._save_state(prefix, state)
            metric = metric_slice_from_q(state.grid, state.geometry)
            self.state_records[time_key] = {
                "prefix": prefix,
                "t": float(state.t),
                "tau": float(state.tau),
                "input_cone": self._profile_summary(metric, state.matter.radiation),
            }

    def finish_step(self, state, dt, advanced):
        if float(advanced.t) >= self.start_time - 1.0e-10:
            metric = metric_slice_from_q(advanced.grid, advanced.geometry)
            summary = {
                "t_before": float(state.t),
                "t": float(advanced.t),
                "tau": float(advanced.tau),
                "dt": float(dt),
                "cone": self._profile_summary(metric, advanced.matter.radiation),
            }
            self.accepted_step_outputs.append(summary)

    def begin_rhs(self, stage_state, ordinal):
        self.current_rhs_record = None
        if self.active_time is None or self.active_attempt_id is None:
            return None
        key = (
            float(self.active_time).hex(),
            int(self.active_attempt_id),
            int(ordinal),
        )
        if self.active_time < self.start_time - 1.0e-10 or key in self.rhs_records:
            return key
        prefix = f"rhs_{self._next_rhs_id:03d}"
        self._next_rhs_id += 1
        self.current_rhs_record = {
            "key": key,
            "prefix": prefix,
            "input_t": float(self.active_time),
            "attempt_id": int(self.active_attempt_id),
            "trial_dt": float(self.active_dt),
            "rhs_ordinal": int(ordinal),
            "stage_geometry": self._geometry_metric(stage_state.grid, stage_state.geometry),
            "stage_state": stage_state,
        }
        return key

    def capture_rhs_components(self, metric, conserved, species, components):
        if self.current_rhs_record is None:
            return
        if species != matter_system_module.Species.RADIATION:
            return
        face_flux, source_e, source_s, q_e, q_s = components
        rec = self.current_rhs_record
        prefix = rec["prefix"]
        self._save_metric(prefix, metric)
        self._save_conserved(prefix, conserved)
        component_values = {
            "face_flux": face_flux,
            "source_energy": source_e,
            "source_momentum": source_s,
            "coupling_energy": q_e,
            "coupling_momentum": q_s,
        }
        rec["components"] = {}
        for name, values in component_values.items():
            copied = self._copy_array(values)
            rec["components"][name] = copied
            self.arrays[f"{prefix}__component__{name}"] = copied

    def finish_rhs(self, stage_state, result, key):
        rec = self.current_rhs_record
        self.current_rhs_record = None
        if rec is None:
            return
        srhs, mrhs, _ = result
        radiation_rhs = mrhs["radiation"]
        prefix = rec["prefix"]
        self._save_conserved(f"{prefix}__rhs", radiation_rhs)
        self._save_state_rhs_geometry(prefix, stage_state)
        rec["radiation_rhs"] = {
            field: self._copy_array(getattr(radiation_rhs, field))
            for field in SPECIES_FIELDS
        }
        rec["stage_matter_radiation"] = {
            field: self._copy_array(getattr(stage_state.matter.radiation, field))
            for field in SPECIES_FIELDS
        }
        self.rhs_records[key] = rec
        self.rhs_geometry[key] = rec["stage_geometry"]

    def _save_state_rhs_geometry(self, prefix, stage_state):
        for name in GEOMETRY_FIELDS:
            self.arrays[f"{prefix}__stage_geometry__{name}"] = self._copy_array(
                getattr(stage_state.geometry, name)
            )
        for name in SCALAR_FIELDS:
            self.arrays[f"{prefix}__stage_scalar__{name}"] = self._copy_array(
                getattr(stage_state.scalars, name)
            )

    @staticmethod
    def _cone_ratio_values(a_values, x_values, conserved):
        a_values = np.asarray(a_values, dtype=float)
        x_values = np.asarray(x_values, dtype=float)
        energy_t = np.asarray(conserved.energy_t, dtype=float)
        momentum_r = np.asarray(conserved.momentum_r, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            root_gamma_inv = np.sqrt(x_values * x_values / a_values)
            ratio = root_gamma_inv * np.abs(momentum_r) / energy_t
            margin = energy_t - root_gamma_inv * np.abs(momentum_r)
        return ratio, margin

    def _explicit_metric_rate_budget(self, predictor_metric, predictor_radiation):
        """Exact read-only split of the explicit a/X update at the failing guard."""
        if self.active_attempt_id is None:
            return {"available": False, "reason": "no active step-attempt identity"}
        key = (float(self.active_time).hex(), int(self.active_attempt_id), 0)
        rec = self.geometry_records.get(key)
        if rec is None:
            return {"available": False, "reason": "no captured geometry-stage RHS at ordinal 0"}
        g0 = rec["stage_metric"]
        explicit = rec["explicit_rates"]
        a0 = np.asarray(g0.a, dtype=float)
        x0 = np.asarray(g0.X, dtype=float)
        b0 = np.asarray(g0.b, dtype=float)
        a_rate = explicit["a"]
        x_rate = explicit["X"]
        b_rate = explicit["b"]
        dt = float(self.active_dt)
        ap = a0 + dt * a_rate
        xp = x0 + dt * x_rate
        bp = b0 + dt * b_rate
        actual_a = np.asarray(predictor_metric.a, dtype=float)
        actual_x = np.asarray(predictor_metric.X, dtype=float)
        actual_b = np.asarray(predictor_metric.b, dtype=float)
        outer = int(a0.size - 1)
        U0 = self.active_state.matter.radiation
        U1 = predictor_radiation

        # Algebraic cross-pairs are not separately evolved counterfactuals.
        metric_pairs = {
            "base": (a0, x0),
            "a_only": (ap, x0),
            "X_only": (a0, xp),
            "a_and_X": (ap, xp),
            "actual_predictor_metric": (actual_a, actual_x),
        }
        ratio_budget = {}
        for u_name, u_state in (("accepted_U", U0), ("predictor_U", U1)):
            ratios = {}
            for metric_name, (aa, xx) in metric_pairs.items():
                ratio, margin = self._cone_ratio_values(aa, xx, u_state)
                ratios[metric_name] = {
                    "outer_cell": outer,
                    "outer_radius": float(self.active_state.grid.centers[outer]),
                    "outer_ratio": float(ratio[outer]),
                    "outer_margin_conserved": float(margin[outer]),
                    "outer_band_max_ratio": float(np.nanmax(ratio[-5:])),
                    "outer_band_max_ratio_cell": int(
                        outer - 4 + np.nanargmax(ratio[-5:])
                    ),
                }
            ratio_budget[u_name] = ratios

        q0 = x0 * x0 / a0
        qa = x0 * x0 / ap
        qx = xp * xp / a0
        qax = xp * xp / ap
        # Symmetric Shapley split shares the nonlinear a-X interaction equally.
        a_contribution = 0.5 * ((qa - q0) + (qax - qx))
        x_contribution = 0.5 * ((qx - q0) + (qax - qa))
        closure = (qax - q0) - (a_contribution + x_contribution)

        def at_outer(values):
            return float(np.asarray(values, dtype=float)[outer])

        return {
            "available": True,
            "method": "read-only explicit RHS; exact finite update; symmetric split of q=X^2/a interaction",
            "input_t": float(self.active_time),
            "trial_dt": dt,
            "geometry_stage_ordinal": 0,
            "geometry_stage_prefix": rec["prefix"],
            "outer_cell": outer,
            "outer_radius": float(self.active_state.grid.centers[outer]),
            "outer_cell_fields_before": {
                "a": at_outer(a0), "X": at_outer(x0), "b": at_outer(b0),
            },
            "outer_cell_explicit_rates": {
                "a_dot": at_outer(a_rate), "X_dot": at_outer(x_rate),
                "b_dot": at_outer(b_rate),
            },
            "outer_cell_fields_expected_predictor": {
                "a": at_outer(ap), "X": at_outer(xp), "b": at_outer(bp),
            },
            "outer_cell_fields_actual_predictor": {
                "a": at_outer(actual_a), "X": at_outer(actual_x),
                "b": at_outer(actual_b),
            },
            "predictor_update_closure_max_abs_all_cells": {
                "a": float(np.max(np.abs(actual_a - ap))),
                "X": float(np.max(np.abs(actual_x - xp))),
                "b": float(np.max(np.abs(actual_b - bp))),
            },
            "predictor_update_closure_outer_cell": {
                "a": at_outer(actual_a - ap),
                "X": at_outer(actual_x - xp),
                "b": at_outer(actual_b - bp),
            },
            "outer_cell_gamma_rr_inv": {
                "base": at_outer(q0), "a_only": at_outer(qa),
                "X_only": at_outer(qx), "a_and_X": at_outer(qax),
                "actual_predictor": at_outer(actual_x * actual_x / actual_a),
                "delta_total": at_outer(qax - q0),
                "delta_attributed_to_a": at_outer(a_contribution),
                "delta_attributed_to_X": at_outer(x_contribution),
                "shapley_closure": at_outer(closure),
            },
            "radiation_cone_cross_pairs": ratio_budget,
            "interpretation_limit": "algebraic decomposition of the immediate predictor update only; not a claim that either field equation or boundary law is wrong",
        }

    def capture_failed_predictor(self, metric, conserved, stage, exc):
        if self.active_state is None or stage != "explicit-predictor":
            return
        prefix = f"failed_predictor_{self._next_failed_id:03d}"
        self._next_failed_id += 1
        self._save_metric(prefix, metric)
        self._save_conserved(prefix, conserved)
        ratio, margin, energy, abs_momentum = self._ratio_profile(metric, conserved)
        lo = max(0, ratio.size - 5)
        trial = {
            "prefix": prefix,
            "input_t": float(self.active_time),
            "attempt_id": int(self.active_attempt_id),
            "requested_trial_dt": float(self.active_dt),
            "stage": str(stage),
            "exception_class": type(exc).__name__,
            "exception_message": str(exc),
            "outer_cell_79_ratio": float(ratio[-1]),
            "outer_cell_79_conserved_margin": float(margin[-1]),
            "outer_cell_79_E": float(energy[-1]),
            "outer_cell_79_abs_S": float(abs_momentum[-1]),
            "outer_band_max_ratio": float(np.nanmax(ratio[lo:])),
            "outer_band_max_ratio_cell": int(lo + np.nanargmax(ratio[lo:])),
            "outer_cells": self._profile_summary(metric, conserved)["outer_cells"],
            "ratio_counterfactuals": {},
        }

        active = self.active_state
        rhs_key = (
            float(self.active_time).hex(),
            int(self.active_attempt_id),
            0,
        )
        rhs_record = self.rhs_records.get(rhs_key)
        input_metric = metric_slice_from_q(active.grid, active.geometry)
        input_radiation = active.matter.radiation
        cmc_metric = self.rhs_geometry.get(rhs_key)
        if cmc_metric is not None:
            metric_options = {
                "accepted_metric": input_metric,
                "cmc_resolved_metric": cmc_metric,
                "explicit_predictor_metric": metric,
            }
            conserved_options = {
                "accepted_conserved_state": input_radiation,
                "explicit_predictor_conserved_state": conserved,
            }
            for u_name, u_state in conserved_options.items():
                for g_name, g_state in metric_options.items():
                    cross_ratio, cross_margin, _, _ = self._ratio_profile(g_state, u_state)
                    trial["ratio_counterfactuals"][f"{u_name}_on_{g_name}"] = {
                        "outer_cell_79_ratio": float(cross_ratio[-1]),
                        "outer_cell_79_conserved_margin": float(cross_margin[-1]),
                        "outer_band_max_ratio": float(np.nanmax(cross_ratio[lo:])),
                        "outer_band_max_ratio_cell": int(lo + np.nanargmax(cross_ratio[lo:])),
                    }

        # Exact algebraic check of the explicit matter predictor against the
        # observed first RHS. This is diagnostic verification, not a projection.
        if rhs_record is not None:
            trial["rhs_capture_key"] = list(rhs_key)
            trial["rhs_capture_trial_dt"] = float(rhs_record["trial_dt"])
            trial["rhs_capture_prefix"] = rhs_record["prefix"]
            rhs = rhs_record["radiation_rhs"]
            dt = float(self.active_dt)
            cells = [int(i) for i in range(lo, ratio.size)]
            update = {}
            for field in ("energy_t", "momentum_r"):
                before = np.asarray(getattr(input_radiation, field), dtype=float)
                after = np.asarray(getattr(conserved, field), dtype=float)
                expected = before + dt * rhs[field]
                update[field] = {
                    "max_abs_predictor_update_closure_all_cells": float(np.max(np.abs(after - expected))),
                    "max_abs_predictor_update_closure_outer_5": float(np.max(np.abs(after[lo:] - expected[lo:]))),
                    "cell_79_actual_delta": float(after[-1] - before[-1]),
                    "cell_79_dt_times_rhs": float(dt * rhs[field][-1]),
                    "cell_79_closure_error": float(after[-1] - expected[-1]),
                }
            trial["explicit_predictor_U_update_closure"] = update
        trial["explicit_metric_rate_budget"] = self._explicit_metric_rate_budget(metric, conserved)
        self.failed_trials.append(trial)

    def _boundary_directional_summary(self, grid, geometry, scalars, matter, original_rhs, original_geometry_stage_terms):
        """Evaluate the exact local radiation-cone directional derivative without evolving state."""
        stage = SimpleNamespace(
            grid=grid,
            geometry=geometry,
            scalars=scalars,
            matter=matter,
            t=float(self.active_time),
            tau=float(self.active_state.tau),
            e_folds=float(self.active_state.e_folds),
        )
        _, matter_rhs, _ = original_rhs(stage)
        radiation_rhs = matter_rhs["radiation"]
        terms = original_geometry_stage_terms(
            grid, geometry, scalars, matter,
            lambda_m=production_kernel_module.LAMBDA_M,
        )
        explicit = terms["explicit"]
        a = np.asarray(geometry.a, dtype=float)
        x = np.asarray(geometry.X, dtype=float)
        energy = np.asarray(matter.radiation.energy_t, dtype=float)
        momentum = np.asarray(matter.radiation.momentum_r, dtype=float)
        denergy = np.asarray(radiation_rhs.energy_t, dtype=float)
        dmomentum = np.asarray(radiation_rhs.momentum_r, dtype=float)
        adot = np.asarray(explicit["a"], dtype=float)
        xdot = np.asarray(explicit["X"], dtype=float)
        q = x * x / a
        sqrt_q = np.sqrt(q)
        ratio = sqrt_q * np.abs(momentum) / energy
        margin = energy - sqrt_q * np.abs(momentum)
        q_log_rate = 2.0 * xdot / x - adot / a
        matter_rate = denergy - sqrt_q * np.sign(momentum) * dmomentum
        metric_rate = -0.5 * sqrt_q * np.abs(momentum) * q_log_rate
        total_rate = matter_rate + metric_rate
        finite = (
            np.isfinite(margin) & np.isfinite(total_rate)
            & np.isfinite(ratio) & (margin > 0.0) & (total_rate < 0.0)
        )
        times = np.full_like(margin, np.inf)
        times[finite] = margin[finite] / (-total_rate[finite])
        fastest = int(np.argmin(times))
        outer = int(len(q) - 1)
        low = max(0, outer - 4)
        def cell(i):
            return {
                "i": int(i),
                "r": float(grid.centers[i]),
                "a": float(a[i]),
                "X": float(x[i]),
                "Aa": float(geometry.Aa[i]),
                "K": float(geometry.K[i]),
                "alpha": float(geometry.alpha[i]),
                "gamma_rr_inv": float(q[i]),
                "ratio_absS_over_E": float(ratio[i]),
                "cone_margin_C": float(margin[i]),
                "Cdot_matter": float(matter_rate[i]),
                "Cdot_metric": float(metric_rate[i]),
                "Cdot_total": float(total_rate[i]),
                "linearized_time_to_boundary": (
                    float(times[i]) if np.isfinite(times[i]) else None
                ),
            }
        return {
            "formula": "Cdot = dU_E - sqrt(q)*sign(U_r)*dU_r - 0.5*sqrt(q)*abs(U_r)*(2*Xdot/X-adot/a), q=X^2/a",
            "outer_cell": cell(outer),
            "outer_band": [cell(i) for i in range(low, outer + 1)],
            "fastest_positive_margin_outward_cell": cell(fastest) if np.isfinite(times[fastest]) else None,
            "radiation_rhs_finite": bool(
                np.all(np.isfinite(denergy)) and np.all(np.isfinite(dmomentum))
            ),
            "geometry_rates_finite": bool(
                np.all(np.isfinite(adot)) and np.all(np.isfinite(xdot))
            ),
        }

    def install(self, kernel):
        self.originals["kernel_rhs"] = kernel._rhs
        original_rhs = kernel._rhs

        def observed_rhs(stage_state):
            ordinal = int(self.rhs_counter)
            self.rhs_counter += 1
            key = self.begin_rhs(stage_state, ordinal)
            try:
                result = original_rhs(stage_state)
            except Exception:
                self.current_rhs_record = None
                raise
            self.finish_rhs(stage_state, result, key)
            return result

        kernel._rhs = observed_rhs
        self.kernel = kernel

        self.originals["validate"] = production_kernel_module.validate_species_admissibility
        original_validate = production_kernel_module.validate_species_admissibility

        def observed_validate(metric, state, species, *, stage):
            try:
                return original_validate(metric, state, species, stage=stage)
            except Exception as exc:
                if species == matter_system_module.Species.RADIATION and stage == "explicit-predictor":
                    self.capture_failed_predictor(metric, state, stage, exc)
                raise

        production_kernel_module.validate_species_admissibility = observed_validate

        # Tap the exact production_kernel species_rhs call. For radiation
        # only, reevaluate the same deterministic component helper with exactly
        # the same stage inputs, then check its flux+source reconstruction
        # against the actual returned RHS. Neither observer changes a value
        # passed to the evolution.
        self.originals["species_rhs"] = production_kernel_module.species_rhs
        original_species_rhs = production_kernel_module.species_rhs
        original_components = matter_rhs_module._species_rhs_components

        def observed_species_rhs(
            metric, metric_derivatives, state, species,
            dphi_t=None, dphi_r=None, beta_dm=-0.04,
        ):
            result = original_species_rhs(
                metric, metric_derivatives, state, species,
                dphi_t=dphi_t, dphi_r=dphi_r, beta_dm=beta_dm,
            )
            if self.current_rhs_record is not None and species == matter_system_module.Species.RADIATION:
                components = original_components(
                    metric, metric_derivatives, state, species,
                    dphi_t=dphi_t, dphi_r=dphi_r, beta_dm=beta_dm,
                )
                self.capture_rhs_components(metric, state, species, components)
                self.radiation_rhs_component_capture_calls += 1
            return result

        production_kernel_module.species_rhs = observed_species_rhs

        self.originals["geometry_stage_terms"] = adapter_module.geometry_stage_terms
        original_geometry_stage_terms = adapter_module.geometry_stage_terms

        def observed_geometry_stage_terms(grid, geometry, fields, matter, lambda_m=2.0):
            ordinal = int(self.geometry_counter)
            self.geometry_counter += 1
            result = original_geometry_stage_terms(
                grid, geometry, fields, matter, lambda_m=lambda_m
            )
            if self.active_time is not None and self.active_time >= self.start_time - 1.0e-10:
                if self.active_attempt_id is None:
                    raise RuntimeError("geometry source capture has no active step-attempt identity")
                key = (float(self.active_time).hex(), int(self.active_attempt_id), ordinal)
                if key not in self.geometry_records:
                    prefix = f"geometry_{self._next_geometry_id:03d}"
                    self._next_geometry_id += 1
                    stage_metric = self._geometry_metric(grid, geometry)
                    self._save_metric(prefix, stage_metric)
                    explicit_rates = {}
                    for field in ("a", "b", "X"):
                        rate = self._copy_array(result["explicit"][field])
                        explicit_rates[field] = rate
                        self.arrays[f"{prefix}__explicit_rate__{field}"] = rate

                    # Exact production L2 and L3 source split at this real PIRK
                    # stage. The extra evaluations below are read-only observers:
                    # their outputs are never passed back to the evolution.
                    _, source_vacuum, _ = adapter_module.vendor_modules()
                    raw_geometry_terms = source_vacuum.geometry_terms(grid, geometry)
                    vacuum_l3 = source_vacuum.primary_l3_rhs(grid, geometry)
                    # Independent source assembly: obtain rho/pr/pt through the
                    # dedicated TotalStressEnergy path on this exact observed stage,
                    # rather than infer the matter term as total-minus-vacuum alone.
                    assembled_stress = assemble_total_stress_energy(
                        grid, geometry, fields, matter
                    )
                    l3_expected = l3_matter_terms_from_projection(
                        geometry.alpha,
                        assembled_stress.rho,
                        assembled_stress.pr,
                        assembled_stress.pt,
                    )
                    l3_Aa_residual, l3_Aa_bound, l3_Aa_pass, l3_Aa_ratio = (
                        source_increment_roundoff_witness(
                            result["primary_l3"]["Aa"],
                            vacuum_l3["Aa"],
                            l3_expected["Aa"],
                        )
                    )
                    l3_K_residual, l3_K_bound, l3_K_pass, l3_K_ratio = (
                        source_increment_roundoff_witness(
                            result["primary_l3"]["K"],
                            vacuum_l3["K"],
                            l3_expected["K"],
                        )
                    )
                    lapse_hessian_Aa = -(
                        np.asarray(raw_geometry_terms["DrDr_alpha"], dtype=float)
                        - np.asarray(raw_geometry_terms["lap_alpha"], dtype=float) / 3.0
                    )
                    ricci_tracefree_Aa = np.asarray(geometry.alpha, dtype=float) * (
                        np.asarray(raw_geometry_terms["Rr_r"], dtype=float)
                        - np.asarray(raw_geometry_terms["R"], dtype=float) / 3.0
                    )
                    source_arrays = {
                        "l2_Aa_lapse_hessian": lapse_hessian_Aa,
                        "l2_Aa_ricci_tracefree": ricci_tracefree_Aa,
                        "l2_Aa_actual": np.asarray(result["primary_l2"]["Aa"], dtype=float),
                        "l2_K_minus_lap_alpha": -np.asarray(raw_geometry_terms["lap_alpha"], dtype=float),
                        "l2_K_actual": np.asarray(result["primary_l2"]["K"], dtype=float),
                        "l3_Aa_vacuum": np.asarray(vacuum_l3["Aa"], dtype=float),
                        "l3_Aa_matter_difference": (
                            np.asarray(result["primary_l3"]["Aa"], dtype=float)
                            - np.asarray(vacuum_l3["Aa"], dtype=float)
                        ),
                        "l3_Aa_actual": np.asarray(result["primary_l3"]["Aa"], dtype=float),
                        "l3_Aa_matter_from_TotalStressEnergy": l3_expected["Aa"],
                        "l3_TotalStressEnergy_rho": np.asarray(assembled_stress.rho, dtype=float),
                        "l3_TotalStressEnergy_pr": np.asarray(assembled_stress.pr, dtype=float),
                        "l3_TotalStressEnergy_pt": np.asarray(assembled_stress.pt, dtype=float),
                        "l3_Aa_matter_identity_residual": l3_Aa_residual,
                        "l3_Aa_matter_operand_scaled_roundoff_bound": l3_Aa_bound,
                        "l3_K_vacuum": np.asarray(vacuum_l3["K"], dtype=float),
                        "l3_K_matter_difference": (
                            np.asarray(result["primary_l3"]["K"], dtype=float)
                            - np.asarray(vacuum_l3["K"], dtype=float)
                        ),
                        "l3_K_actual": np.asarray(result["primary_l3"]["K"], dtype=float),
                        "l3_K_matter_from_TotalStressEnergy": l3_expected["K"],
                        "l3_K_matter_identity_residual": l3_K_residual,
                        "l3_K_matter_operand_scaled_roundoff_bound": l3_K_bound,
                    }
                    for name, values in source_arrays.items():
                        copied = self._copy_array(values)
                        source_arrays[name] = copied
                        self.arrays[f"{prefix}__source__{name}"] = copied

                    n = int(stage_metric.a.size)
                    outer = n - 1
                    low = max(0, n - 5)
                    def at_cell(values, i):
                        return float(np.asarray(values, dtype=float)[i])
                    stage_state = SimpleNamespace(
                        grid=grid,
                        geometry=geometry,
                        scalars=fields,
                        matter=matter,
                        t=float(self.active_time),
                        tau=float(self.active_state.tau),
                        e_folds=float(self.active_state.e_folds),
                    )
                    stage_state_sha256 = numeric_state_digest(stage_state)
                    source_record = {
                        "prefix": prefix,
                        "input_t": float(self.active_time),
                        "trial_dt": float(self.active_dt),
                        "attempt_id": int(self.active_attempt_id),
                        "stage_state_sha256": stage_state_sha256,
                        "ordinal": int(ordinal),
                        "stage_label": {
                            0: "accepted_input_after_current_CMC_solve",
                            1: "explicit_predictor_before_primary_Aa_K_update",
                            2: "primary_predictor_after_outer_boundary_and_CMC",
                        }.get(int(ordinal), "other_geometry_stage"),
                        "outer_cell": outer,
                        "outer_radius": float(grid.centers[outer]),
                        "outer_band": [],
                        "closure_max_abs_all_cells": {
                            "L2_Aa_lapse_hessian_plus_Ricci_minus_actual": float(np.max(np.abs(
                                source_arrays["l2_Aa_lapse_hessian"]
                                + source_arrays["l2_Aa_ricci_tracefree"]
                                - source_arrays["l2_Aa_actual"]
                            ))),
                            "L2_K_minus_lap_alpha_minus_actual": float(np.max(np.abs(
                                source_arrays["l2_K_minus_lap_alpha"]
                                - source_arrays["l2_K_actual"]
                            ))),
                            "L3_Aa_vacuum_plus_matter_difference_minus_actual": float(np.max(np.abs(
                                source_arrays["l3_Aa_vacuum"]
                                + source_arrays["l3_Aa_matter_difference"]
                                - source_arrays["l3_Aa_actual"]
                            ))),
                            "L3_K_vacuum_plus_matter_difference_minus_actual": float(np.max(np.abs(
                                source_arrays["l3_K_vacuum"]
                                + source_arrays["l3_K_matter_difference"]
                                - source_arrays["l3_K_actual"]
                            ))),
                        },
                        "outer_cell_terms": {},
                        "l3_total_stress_projection_oracle": {
                            "assembly": "engine.stress_energy.assemble_total_stress_energy",
                            "projection_input_profile_keys": {
                                "alpha": f"{prefix}__metric__alpha",
                                "rho": f"{prefix}__source__l3_TotalStressEnergy_rho",
                                "pr": f"{prefix}__source__l3_TotalStressEnergy_pr",
                                "pt": f"{prefix}__source__l3_TotalStressEnergy_pt",
                            },
                            "Aa": {
                                "residual_max_abs_all_cells": float(np.max(np.abs(l3_Aa_residual))),
                                "operand_scaled_roundoff_bound_max_abs": float(np.max(l3_Aa_bound)),
                                "max_residual_over_bound": l3_Aa_ratio,
                                "within_bound_all_cells": l3_Aa_pass,
                            },
                            "K": {
                                "residual_max_abs_all_cells": float(np.max(np.abs(l3_K_residual))),
                                "operand_scaled_roundoff_bound_max_abs": float(np.max(l3_K_bound)),
                                "max_residual_over_bound": l3_K_ratio,
                                "within_bound_all_cells": l3_K_pass,
                            },
                        },
                    }
                    for i in range(low, n):
                        row = {
                            "i": int(i),
                            "r": float(grid.centers[i]),
                            "Aa": float(stage_metric.Aa[i]),
                            "K": float(stage_metric.K[i]),
                            "alpha": float(stage_metric.alpha[i]),
                        }
                        row.update({name: at_cell(values, i) for name, values in source_arrays.items()})
                        source_record["outer_band"].append(row)
                    source_record["outer_cell_terms"] = source_record["outer_band"][-1]
                    self.geometry_source_records.append(source_record)

                    self.geometry_records[key] = {
                        "prefix": prefix,
                        "input_t": float(self.active_time),
                        "trial_dt": float(self.active_dt),
                        "attempt_id": int(self.active_attempt_id),
                        "stage_state_sha256": stage_state_sha256,
                        "ordinal": ordinal,
                        "stage_metric": stage_metric,
                        "explicit_rates": explicit_rates,
                        "source_terms": source_arrays,
                    }
            return result

        adapter_module.geometry_stage_terms = observed_geometry_stage_terms


        self.originals["outer_boundary"] = kernel._apply_outer_light_boundary
        original_boundary = kernel._apply_outer_light_boundary

        def observed_outer_boundary(grid, geometry, scalars, matter):
            record = {
                "index": self._next_boundary_id,
                "input_t": float(self.active_time),
                "trial_dt": float(self.active_dt),
                "diagnostic_only": True,
            }
            self._next_boundary_id += 1
            fields = ("a", "b", "X", "Aa", "K", "Lambda", "alpha", "beta", "B")
            saved_before = {}
            try:
                for name in fields:
                    saved_before[name] = self._copy_array(getattr(geometry, name))
                radiation_before = {
                    name: self._copy_array(getattr(matter.radiation, name))
                    for name in SPECIES_FIELDS
                }
                record["before"] = self._boundary_directional_summary(
                    grid, geometry, scalars, matter,
                    original_rhs, original_geometry_stage_terms,
                )
                before_summary_ok = True
            except Exception as exc:
                record["before_capture_error"] = f"{type(exc).__name__}: {exc}"
                before_summary_ok = False

            # Call the actual production boundary exactly once, unchanged.
            result = original_boundary(grid, geometry, scalars, matter)

            try:
                for name in fields:
                    after = self._copy_array(getattr(geometry, name))
                    self.arrays[
                        f"boundary_{record['index']:03d}__geometry__{name}__before"
                    ] = saved_before[name]
                    self.arrays[
                        f"boundary_{record['index']:03d}__geometry__{name}__after"
                    ] = after
                radiation_after = {
                    name: self._copy_array(getattr(matter.radiation, name))
                    for name in SPECIES_FIELDS
                }
                for name in SPECIES_FIELDS:
                    self.arrays[
                        f"boundary_{record['index']:03d}__radiation__{name}__before"
                    ] = radiation_before[name]
                    self.arrays[
                        f"boundary_{record['index']:03d}__radiation__{name}__after"
                    ] = radiation_after[name]
                record["after"] = self._boundary_directional_summary(
                    grid, geometry, scalars, matter,
                    original_rhs, original_geometry_stage_terms,
                )
                record["field_deltas"] = {}
                for name in fields:
                    pre = saved_before[name]
                    post = self._copy_array(getattr(geometry, name))
                    record["field_deltas"][name] = {
                        "max_abs_all_cells": float(np.max(np.abs(post - pre))),
                        "outer_cell_before": float(pre[-1]),
                        "outer_cell_after": float(post[-1]),
                        "outer_cell_delta": float(post[-1] - pre[-1]),
                    }
                record["radiation_deltas"] = {
                    name: {
                        "max_abs_all_cells": float(np.max(np.abs(radiation_after[name] - radiation_before[name]))),
                        "outer_cell_before": float(radiation_before[name][-1]),
                        "outer_cell_after": float(radiation_after[name][-1]),
                    }
                    for name in SPECIES_FIELDS
                }
                record["cone_rate_delta"] = {
                    "outer_cell_Cdot_matter": (
                        record["after"]["outer_cell"]["Cdot_matter"]
                        - record["before"]["outer_cell"]["Cdot_matter"]
                    ),
                    "outer_cell_Cdot_metric": (
                        record["after"]["outer_cell"]["Cdot_metric"]
                        - record["before"]["outer_cell"]["Cdot_metric"]
                    ),
                    "outer_cell_Cdot_total": (
                        record["after"]["outer_cell"]["Cdot_total"]
                        - record["before"]["outer_cell"]["Cdot_total"]
                    ),
                    "outer_cell_C_before": record["before"]["outer_cell"]["cone_margin_C"],
                    "outer_cell_C_after": record["after"]["outer_cell"]["cone_margin_C"],
                } if before_summary_ok and "after" in record else None
            except Exception as exc:
                record["after_capture_error"] = f"{type(exc).__name__}: {exc}"
            self.boundary_records.append(record)
            return result

        kernel._apply_outer_light_boundary = observed_outer_boundary
        return self.restore

    def restore(self):
        if "kernel_rhs" in self.originals:
            original_rhs = self.originals.pop("kernel_rhs")
            if self.kernel is not None:
                self.kernel._rhs = original_rhs
        if "validate" in self.originals:
            production_kernel_module.validate_species_admissibility = self.originals.pop("validate")
        if "species_rhs" in self.originals:
            production_kernel_module.species_rhs = self.originals.pop("species_rhs")
        if "geometry_stage_terms" in self.originals:
            adapter_module.geometry_stage_terms = self.originals.pop("geometry_stage_terms")
        if "outer_boundary" in self.originals and self.kernel is not None:
            self.kernel._apply_outer_light_boundary = self.originals.pop("outer_boundary")

    def _rhs_summary(self):
        rows = []
        for (time_key, attempt_id, ordinal), rec in sorted(
            self.rhs_records.items(),
            key=lambda item: (float.fromhex(item[0][0]), item[0][1], item[0][2]),
        ):
            row = {
                "input_t": rec["input_t"],
                "attempt_id": rec["attempt_id"],
                "trial_dt": rec["trial_dt"],
                "rhs_ordinal": rec["rhs_ordinal"],
                "prefix": rec["prefix"],
                "has_components": "components" in rec,
            }
            components = rec.get("components")
            rhs = rec.get("radiation_rhs")
            if components is not None and rhs is not None:
                face = components["face_flux"]
                inv_dr = 1.0 / float(rec["stage_state"].grid.dr)
                flux_e = -inv_dr * (face[1:, 1] - face[:-1, 1])
                flux_s = -inv_dr * (face[1:, 2] - face[:-1, 2])
                source_e = components["source_energy"]
                source_s = components["source_momentum"]
                q_e = components["coupling_energy"]
                q_s = components["coupling_momentum"]
                net_e = flux_e + source_e + q_e
                net_s = flux_s + source_s + q_s
                row["all_cell_component_closure_max_abs"] = {
                    "energy_t": float(np.max(np.abs(net_e - rhs["energy_t"]))),
                    "momentum_r": float(np.max(np.abs(net_s - rhs["momentum_r"]))),
                }
                lo = max(0, rhs["energy_t"].size - 5)
                row["outer_5_cells"] = []
                for i in range(lo, rhs["energy_t"].size):
                    row["outer_5_cells"].append({
                        "i": int(i),
                        "r": float(rec["stage_state"].grid.centers[i]),
                        "flux_div_energy_t": float(flux_e[i]),
                        "geometric_source_energy_t": float(source_e[i]),
                        "coupling_source_energy_t": float(q_e[i]),
                        "net_rhs_energy_t": float(rhs["energy_t"][i]),
                        "flux_div_momentum_r": float(flux_s[i]),
                        "geometric_source_momentum_r": float(source_s[i]),
                        "coupling_source_momentum_r": float(q_s[i]),
                        "net_rhs_momentum_r": float(rhs["momentum_r"][i]),
                    })
            rows.append(row)
        return rows

    def write(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)
        ordered_states = sorted(
            self.state_records.values(), key=lambda row: float(row["t"])
        )
        state_history = []
        for rec in ordered_states:
            item = {k: v for k, v in rec.items() if k != "prefix"}
            state_history.append(item)
        source_record_keys = [
            (float(rec["input_t"]).hex(), int(rec["attempt_id"]), int(rec["ordinal"]))
            for rec in self.geometry_source_records
        ]
        source_keys_unique = len(source_record_keys) == len(set(source_record_keys))
        source_keys_match_index = (
            len(source_record_keys) == len(self.geometry_records)
            and source_keys_unique
            and set(source_record_keys) == set(self.geometry_records.keys())
        )
        rhs_record_keys = [tuple(rec["key"]) for rec in self.rhs_records.values()]
        rhs_keys_unique = len(rhs_record_keys) == len(set(rhs_record_keys))
        rhs_keys_match_index = (
            len(rhs_record_keys) == len(self.rhs_records)
            and rhs_keys_unique
            and set(rhs_record_keys) == set(self.rhs_records.keys())
            and set(self.rhs_geometry.keys()) == set(self.rhs_records.keys())
            and all(
                rec.get("attempt_id") == key[1]
                and rec.get("rhs_ordinal") == key[2]
                and float(rec["input_t"]).hex() == key[0]
                and float(rec["trial_dt"]) > 0.0
                for key, rec in self.rhs_records.items()
            )
        )
        failed_trial_keys = []
        failed_trial_rhs_match_count = 0
        for trial in self.failed_trials:
            try:
                key = (
                    float(trial["input_t"]).hex(),
                    int(trial["attempt_id"]),
                    0,
                )
                requested_dt = float(trial["requested_trial_dt"])
                rhs_record = self.rhs_records.get(key)
                matched = (
                    rhs_record is not None
                    and int(rhs_record.get("attempt_id", -1)) == key[1]
                    and float(rhs_record.get("input_t", float("nan"))).hex() == key[0]
                    and float(rhs_record.get("trial_dt", float("nan"))).hex() == requested_dt.hex()
                    and trial.get("rhs_capture_key") == list(key)
                    and float(trial.get("rhs_capture_trial_dt", float("nan"))).hex() == requested_dt.hex()
                    and trial.get("rhs_capture_prefix") == rhs_record.get("prefix")
                    and "explicit_predictor_U_update_closure" in trial
                )
                failed_trial_keys.append(key)
                failed_trial_rhs_match_count += int(matched)
            except (KeyError, TypeError, ValueError, OverflowError):
                failed_trial_keys.append(None)
        failed_trial_rhs_identity_pass = (
            len(failed_trial_keys) == len(self.failed_trials)
            and None not in failed_trial_keys
            and len(set(failed_trial_keys)) == len(failed_trial_keys)
            and failed_trial_rhs_match_count == len(self.failed_trials)
        )
        l3_projection_profile_error_count = 0
        for rec in self.geometry_source_records:
            prefix = rec["prefix"]
            expected_profile_keys = {
                "alpha": f"{prefix}__metric__alpha",
                "rho": f"{prefix}__source__l3_TotalStressEnergy_rho",
                "pr": f"{prefix}__source__l3_TotalStressEnergy_pr",
                "pt": f"{prefix}__source__l3_TotalStressEnergy_pt",
            }
            actual_profile_keys = rec.get(
                "l3_total_stress_projection_oracle", {}
            ).get("projection_input_profile_keys")
            profile_ok = actual_profile_keys == expected_profile_keys
            if profile_ok:
                try:
                    reference_shape = self.arrays[f"{prefix}__metric__r"].shape
                    for key in expected_profile_keys.values():
                        arr = np.asarray(self.arrays[key], dtype=float)
                        if arr.shape != reference_shape or not np.all(np.isfinite(arr)):
                            profile_ok = False
                            break
                except (KeyError, TypeError, ValueError):
                    profile_ok = False
            if not profile_ok:
                l3_projection_profile_error_count += 1
        l3_projection_profiles_identity_pass = (
            bool(self.geometry_source_records)
            and l3_projection_profile_error_count == 0
        )
        l3_oracle_stages = [
            {
                "input_t": rec["input_t"],
                "attempt_id": rec["attempt_id"],
                "trial_dt": rec["trial_dt"],
                "ordinal": rec["ordinal"],
                "stage_state_sha256": rec["stage_state_sha256"],
                **rec["l3_total_stress_projection_oracle"],
            }
            for rec in self.geometry_source_records
        ]
        l3_oracle_results = {
            "assembly": "engine.stress_energy.assemble_total_stress_energy",
            "expected_Aa_term": "-16*pi*alpha*(TotalStressEnergy.pr-TotalStressEnergy.pt)/3",
            "expected_K_term": "4*pi*alpha*(TotalStressEnergy.rho+pr+2*pt)",
            "captured_stage_count": len(l3_oracle_stages),
            "all_captured_Aa_stages_within_operand_scaled_roundoff_bound": (
                bool(l3_oracle_stages)
                and all(stage["Aa"]["within_bound_all_cells"] for stage in l3_oracle_stages)
            ),
            "all_captured_K_stages_within_operand_scaled_roundoff_bound": (
                bool(l3_oracle_stages)
                and all(stage["K"]["within_bound_all_cells"] for stage in l3_oracle_stages)
            ),
            "all_captured_stages_have_recomputable_TotalStressEnergy_inputs": (
                l3_projection_profiles_identity_pass
            ),
            "projection_input_profile_error_count": l3_projection_profile_error_count,
            "stages": l3_oracle_stages,
            "interpretation_limit": (
                "This checks the adapter L3 increment against separate TotalStressEnergy "
                "assembly at each captured stage; it is not a full trajectory conservation "
                "or equation-validity proof."
            ),
        }
        summary = {
            "schema": "radiation_explicit_predictor_stage_capture_v5_retry_identity_and_recomputable_independent_l3",
            "diagnostic_only": True,
            "rhs_capture_identity_pass": bool(rhs_keys_match_index),
            "failed_predictor_rhs_identity_pass": bool(failed_trial_rhs_identity_pass),
            "failed_predictor_rhs_identity_match_count": failed_trial_rhs_match_count,
            "failed_predictor_rhs_identity_mismatch_count": len(self.failed_trials) - failed_trial_rhs_match_count,
            "rhs_capture_unique_key_count": len(set(rhs_record_keys)),
            "rhs_capture_record_count": len(self.rhs_records),
            "accepted_state_count": len(self.state_records),
            "captured_rhs_evaluations": len(self.rhs_records),
            "failed_explicit_predictor_candidates": len(self.failed_trials),
            "radiation_rhs_component_capture_calls": self.radiation_rhs_component_capture_calls,
            "captured_geometry_stage_terms": len(self.geometry_records),
            "step_attempt_count": self.step_attempt_counter,
            "geometry_source_term_record_count": len(self.geometry_source_records),
            "geometry_source_term_unique_key_count": len(set(source_record_keys)),
            "geometry_source_capture_identity_pass": bool(source_keys_match_index),
            "l3_total_stress_projection_oracle": l3_oracle_results,
            "geometry_source_term_capture_records": self.geometry_source_records,
            "outer_boundary_capture_calls": len(self.boundary_records),
            "outer_boundary_records": self.boundary_records,
            "geometry_stage_explicit_rates": [
                {
                    "input_t": rec["input_t"],
                    "trial_dt": rec["trial_dt"],
                    "attempt_id": rec["attempt_id"],
                    "stage_state_sha256": rec["stage_state_sha256"],
                    "ordinal": rec["ordinal"],
                    "prefix": rec["prefix"],
                    "outer_cell": int(rec["stage_metric"].a.size - 1),
                    "outer_a_dot": float(rec["explicit_rates"]["a"][-1]),
                    "outer_b_dot": float(rec["explicit_rates"]["b"][-1]),
                    "outer_X_dot": float(rec["explicit_rates"]["X"][-1]),
                }
                for _, rec in sorted(
                    self.geometry_records.items(),
                    key=lambda item: (float.fromhex(item[0][0]), item[0][1], item[0][2]),
                )
            ],
            "accepted_state_history": state_history,
            "accepted_step_outputs": self.accepted_step_outputs,
            "rhs_component_decomposition": self._rhs_summary(),
            "failed_predictor_trials": self.failed_trials,
            "interpretation_limits": [
                "All taps observe calls to current production functions and return their original values unchanged.",
                "Failed candidates are diagnostics only and are never written back to the evolving state.",
                "The first RHS explicit predictor is U0 + dt*RHS0; the captured closure is reported to verify exact reconstruction.",
                "Cone ratios compare conserved states and stage metrics; the counterfactual cross-pairs isolate metric-versus-matter contributions algebraically but are not alternative evolutions.",
                "The bounded run restores the archived accepted state and does not replay the long lead-in trajectory.",
                "This evidence does not admit a production horizon or prove the physical correctness of the outer boundary condition.",
            ],
        }
        summary_path = self.output_dir / "stage_capture.json"
        summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False))
        npz_path = self.output_dir / "stage_capture_profiles.npz"
        np.savez_compressed(npz_path, **self.arrays)
        self.summary_path = summary_path
        self.profiles_path = npz_path
        return summary


class FrozenSeedCampaignKernel:
    """Give the real campaign runner a copy of the exact accepted snapshot."""

    def __init__(self, production_kernel, frozen_seed, capture=None):
        self.production_kernel = production_kernel
        self.capabilities = production_kernel.capabilities
        self.frozen_seed = frozen_seed
        self.capture = capture
        self.initializations = 0

    def initialize(self, *, resolution, r_max, include_radiation=True):
        if int(resolution) != RESOLUTION or float(r_max) != R_MAX:
            raise RuntimeError(
                f"unexpected initializer request N={resolution}, r_max={r_max}"
            )
        if not include_radiation:
            raise RuntimeError("radiation must remain enabled")
        self.initializations += 1
        if self.initializations != 1:
            raise RuntimeError("expected exactly one frozen-state initialization")
        return copy.deepcopy(self.frozen_seed)

    def step(self, state, dt):
        if self.capture is not None:
            self.capture.begin_step(state, dt)
        try:
            advanced = self.production_kernel.step(state, dt)
        except Exception as exc:
            raise
        else:
            if self.capture is not None:
                self.capture.finish_step(state, dt, advanced)
            return advanced

    def diagnostics(self, state, *, profiles=False):
        return self.production_kernel.diagnostics(state, profiles=profiles)

    def checkpoint(self, state, path):
        return self.production_kernel.checkpoint(state, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("runs/source-artifact"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/current-main-radiation-frozen-state-replay"),
    )
    parser.add_argument(
        "--campaign-dir",
        type=Path,
        default=Path("runs/current-main-radiation-frozen-state-replay/campaign"),
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    report = {
        "schema": "current_main_radiation_frozen_state_replay_v1",
        "diagnostic_only": True,
        "status": "audit_not_completed",
        "provenance": {
            "repository": "miahpenn/0-Cosmos-engine",
            "candidate_head": git_value("rev-parse", "HEAD"),
            "expected_base_main_commit": BASE_MAIN_COMMIT,
            "vendor_pin": git_value("rev-parse", "HEAD:vendor/bb-palatini-unified-r0"),
            "expected_vendor_pin": VENDOR_PIN,
            "archived_numerical_source_commit": ARCHIVED_SOURCE_COMMIT,
            "input_artifact_run": SOURCE_RUN,
            "input_artifact_id": SOURCE_ARTIFACT,
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "input_state": {},
        "campaign_runner": {},
        "interpretation_limits": [
            "The exact archived numerical profiles are restored; the approximately 1571-step lead-in is not replayed.",
            "The stored profiles were accepted under archived source cdcec828...; one bounded next-step/short-horizon probe uses current main source.",
            "This checks current runner/retry/reporting behavior on a frozen state, not current-main trajectory equivalence, long-horizon stability, convergence, or conservation.",
            "The first accepted step can be followed by retries and a fail-closed noncompletion; this is not an invitation to relax radiation admissibility or alter state variables.",
            "No production physics equations, coefficients, source terms, boundary conditions, floors, clamps, projections, or damping are changed.",
        ],
    }

    try:
        current_head = git_value("rev-parse", "HEAD")
        if current_head == "unavailable":
            raise RuntimeError("could not identify the checked-out commit")
        vendor = git_value("rev-parse", "HEAD:vendor/bb-palatini-unified-r0")
        if vendor != VENDOR_PIN:
            raise RuntimeError(f"vendor pin mismatch: {vendor}")
        changed_text = subprocess.check_output(
            ["git", "diff", "--name-only", BASE_MAIN_COMMIT, "HEAD"], text=True
        )
        changed = {line for line in changed_text.splitlines() if line}
        if changed != EXPECTED_DIAGNOSTIC_PATHS:
            raise RuntimeError(
                f"unexpected tree delta from pinned main: {sorted(changed)}"
            )

        npz_path = args.input_dir / "boundary_failure_accepted_state.npz"
        metadata_path = args.input_dir / "boundary_failure_accepted_state.json"
        if not npz_path.is_file() or not metadata_path.is_file():
            raise FileNotFoundError("required archived accepted-state files are missing")
        npz_hash = sha256(npz_path)
        if npz_hash != EXPECTED_NPZ_SHA256:
            raise RuntimeError(f"accepted-state NPZ hash mismatch: {npz_hash}")
        metadata = strict_json(metadata_path)
        if metadata.get("rate_fix_source_commit") != ARCHIVED_SOURCE_COMMIT:
            raise RuntimeError("archived accepted-state source metadata mismatch")
        if metadata.get("pinned_vendor_commit") != VENDOR_PIN:
            raise RuntimeError("archived accepted-state vendor metadata mismatch")

        kernel, seed = load_frozen_state(npz_path)
        digest_before = numeric_state_digest(seed)
        cone = cone_witness(seed)
        if not cone["all_finite"]:
            raise RuntimeError("frozen-state radiation cone witness is non-finite")
        report["provenance"].update({
            "frozen_state_npz_sha256": npz_hash,
            "frozen_state_metadata_sha256": sha256(metadata_path),
            "metadata_archived_source_commit": metadata.get("rate_fix_source_commit"),
            "metadata_vendor_pin": metadata.get("pinned_vendor_commit"),
        })
        report["input_state"] = {
            "resolution": RESOLUTION,
            "r_max": R_MAX,
            "t": float(seed.t),
            "tau": float(seed.tau),
            "numeric_state_sha256_before": digest_before,
            "radiation_cone_witness": cone,
            "accepted_state_metadata": {
                "failure_stage": metadata.get("failure", {}).get("stage"),
                "failure_cell": metadata.get("failure", {}).get("grid_index"),
                "archived_t_attempted": metadata.get("attempted_time"),
                "archived_abs_S_over_E": metadata.get("failure", {}).get(
                    "failure_capture", {}
                ).get("ratio_absS_over_E"),
            },
        }

        target_time = math.nextafter(float(seed.t) + REQUESTED_DT, math.inf)
        if not target_time - float(seed.t) > REQUESTED_DT:
            raise RuntimeError("could not construct target beyond nominal dt")
        args.campaign_dir.mkdir(parents=True, exist_ok=True)
        config = CampaignConfig(
            resolutions=(RESOLUTION,),
            r_max=R_MAX,
            cfl=REQUESTED_DT,  # dr=1 for the archived N=80, r_max=80 grid
            final_time=target_time,
            checkpoint_interval=100.0,
            output_dir=str(args.campaign_dir),
        )
        capture = RadiationStageCapture(args.output_dir, float(seed.t))
        wrapped = FrozenSeedCampaignKernel(kernel, seed, capture=capture)
        capture.install(kernel)
        runner_exception = None
        try:
            try:
                returned = run_campaign(kernel=wrapped, config=config)
            except RuntimeError as exc:
                # The strict current runner persists ledgers and summary before
                # raising on an incomplete campaign. Read them below and verify.
                runner_exception = f"{type(exc).__name__}: {exc}"
                returned = None
        finally:
            capture.restore()
            capture_summary = capture.write()

        if capture_summary.get("rhs_capture_identity_pass") is not True:
            raise RuntimeError("RHS records do not map one-to-one to unique attempt/time/ordinal keys")
        if capture_summary.get("failed_predictor_rhs_identity_pass") is not True:
            raise RuntimeError("failed predictor records do not map one-to-one to the RHS for the same attempt and dt")
        if capture_summary.get("geometry_source_capture_identity_pass") is not True:
            raise RuntimeError("geometry source records do not map one-to-one to unique attempt/time/stage keys")
        l3_oracle = capture_summary.get("l3_total_stress_projection_oracle")
        if not isinstance(l3_oracle, dict) or int(l3_oracle.get("captured_stage_count", 0)) <= 0:
            raise RuntimeError("independent L3 TotalStressEnergy witness has no captured stages")
        if l3_oracle.get("all_captured_Aa_stages_within_operand_scaled_roundoff_bound") is not True:
            raise RuntimeError("independent L3 Aa source increment exceeded its operand-scaled roundoff bound")
        if l3_oracle.get("all_captured_K_stages_within_operand_scaled_roundoff_bound") is not True:
            raise RuntimeError("independent L3 K source increment exceeded its operand-scaled roundoff bound")
        if l3_oracle.get("all_captured_stages_have_recomputable_TotalStressEnergy_inputs") is not True:
            raise RuntimeError("independent L3 source projection inputs are missing, non-finite, or mis-keyed in the profile artifact")

        run_dir = args.campaign_dir / f"N{RESOLUTION}"
        retries_path = run_dir / "radiation_step_retries.json"
        ledgers_path = run_dir / "ledgers.json"
        summary_path = args.campaign_dir / "campaign_summary.json"
        for path in (retries_path, ledgers_path, summary_path):
            if not path.is_file():
                raise RuntimeError(f"current runner did not persist {path}")

        retries = strict_json(retries_path)
        ledger = strict_json(ledgers_path)
        summary = strict_json(summary_path)
        if not isinstance(retries, list):
            raise RuntimeError("radiation retry ledger is not a JSON array")
        if ledger.get("radiation_step_retries") != retries:
            raise RuntimeError("combined and standalone retry ledgers differ")
        results = summary.get("results", [])
        if len(results) != 1:
            raise RuntimeError("campaign summary does not contain exactly one result")
        result = results[0]
        admission = summary.get("admission")
        if not isinstance(admission, dict) or not isinstance(admission.get("passed"), bool):
            raise RuntimeError("current strict admission object is missing or malformed")
        if wrapped.initializations != 1:
            raise RuntimeError("runner did not consume exactly one frozen seed")
        digest_after = numeric_state_digest(seed)
        if digest_after != digest_before:
            raise RuntimeError("campaign runner mutated the preserved archived seed")

        retry_check = {
            "entry_count": len(retries),
            "first_entry": retries[0] if retries else None,
            "first_nominal_retry_matches": False,
            "combined_ledger_matches": ledger.get("radiation_step_retries") == retries,
        }
        if retries:
            first = retries[0]
            statuses = [item.get("status") for item in first.get("attempts", [])]
            retry_check["first_nominal_retry_matches"] = (
                math.isclose(float(first.get("requested_dt", math.nan)), REQUESTED_DT,
                             rel_tol=0.0, abs_tol=2.0e-12)
                and math.isclose(float(first.get("accepted_dt", math.nan)),
                                 ACCEPTED_FIRST_RETRY_DT, rel_tol=0.0, abs_tol=2.0e-12)
                and int(first.get("retry_count", -1)) == 1
                and statuses == ["rejected_radiation_admissibility", "accepted"]
            )
        final_time = float(result.get("t_final", math.nan))
        horizon_shortfall = target_time - final_time
        completed = (
            result.get("status") == "completed"
            and admission.get("passed") is True
            and math.isfinite(final_time)
            and abs(final_time - target_time) <= 1.0e-10
        )
        expected_noncompletion = (
            admission.get("passed") is False
            and isinstance(result.get("status"), str)
            and result["status"].startswith("numerical_failure:")
            and "radiation" in result["status"].lower()
            and math.isfinite(final_time)
            and horizon_shortfall > 0.0
            and retry_check["entry_count"] >= 1
            and retry_check["first_nominal_retry_matches"]
        )
        if completed:
            status = "bounded_probe_completed_with_strict_admission"
        elif expected_noncompletion:
            status = "bounded_probe_records_radiation_noncompletion_fail_closed"
        else:
            status = "unexpected_outcome"

        report["stage_capture"] = {
            "summary_schema": capture_summary.get("schema"),
            "accepted_state_count": capture_summary.get("accepted_state_count"),
            "captured_rhs_evaluations": capture_summary.get("captured_rhs_evaluations"),
            "rhs_capture_identity_pass": capture_summary.get("rhs_capture_identity_pass"),
            "rhs_capture_unique_key_count": capture_summary.get("rhs_capture_unique_key_count"),
            "failed_predictor_rhs_identity_pass": capture_summary.get("failed_predictor_rhs_identity_pass"),
            "failed_predictor_rhs_identity_match_count": capture_summary.get("failed_predictor_rhs_identity_match_count"),
            "failed_predictor_rhs_identity_mismatch_count": capture_summary.get("failed_predictor_rhs_identity_mismatch_count"),
            "rhs_capture_record_count": capture_summary.get("rhs_capture_record_count"),
            "failed_explicit_predictor_candidates": capture_summary.get("failed_explicit_predictor_candidates"),
            "radiation_rhs_component_capture_calls": capture_summary.get("radiation_rhs_component_capture_calls"),
            "captured_geometry_stage_terms": capture_summary.get("captured_geometry_stage_terms"),
            "geometry_source_term_capture_records": len(capture_summary.get("geometry_source_term_capture_records", [])),
            "geometry_source_capture_identity_pass": capture_summary.get("geometry_source_capture_identity_pass"),
            "geometry_source_term_unique_key_count": capture_summary.get("geometry_source_term_unique_key_count"),
            "l3_total_stress_projection_oracle": capture_summary.get("l3_total_stress_projection_oracle"),
            "step_attempt_count": capture_summary.get("step_attempt_count"),
            "outer_boundary_capture_calls": capture_summary.get("outer_boundary_capture_calls"),
            "summary_sha256": sha256(capture.summary_path),
            "profiles_sha256": sha256(capture.profiles_path),
        }
        report["campaign_runner"] = {
            "outcome_classification": status,
            "returned_results": returned is not None,
            "runner_exception": runner_exception,
            "resolution": result.get("resolution"),
            "status": result.get("status"),
            "t_final": final_time,
            "target_time": target_time,
            "horizon_shortfall": horizon_shortfall,
            "tau_final": result.get("tau_final"),
            "steps_in_history": result.get("steps"),
            "strict_admission_present": True,
            "admission_passed": admission.get("passed"),
            "admission_errors": admission.get("errors"),
            "retry_ledger": retry_check,
            "numeric_state_sha256_after": digest_after,
            "archived_seed_unchanged": digest_after == digest_before,
            "campaign_summary_keys": sorted(summary.keys()),
        }
        report["status"] = status
        report["elapsed_seconds"] = time.time() - started
        # Bundle the exact frozen input alongside the report so every
        # checksum entry in this result artifact can be independently checked.
        shutil.copy2(npz_path, args.output_dir / "source_boundary_failure_accepted_state.npz")
        shutil.copy2(metadata_path, args.output_dir / "source_boundary_failure_accepted_state.json")
        report_path = args.output_dir / "report.json"
        report_path.write_text(json.dumps(report, indent=2, allow_nan=False))
        report_hash = sha256(report_path)
        metadata_hash = sha256(metadata_path)
        capture_summary_hash = sha256(capture.summary_path)
        capture_profiles_hash = sha256(capture.profiles_path)
        (args.output_dir / "SHA256SUMS.txt").write_text(
            f"{report_hash}  report.json\n"
            f"{npz_hash}  source_boundary_failure_accepted_state.npz\n"
            f"{metadata_hash}  source_boundary_failure_accepted_state.json\n"
            f"{capture_summary_hash}  stage_capture.json\n"
            f"{capture_profiles_hash}  stage_capture_profiles.npz\n"
        )
        print(json.dumps({
            "status": report["status"],
            "candidate_head": current_head,
            "runner_status": result.get("status"),
            "t_final": final_time,
            "target_time": target_time,
            "horizon_shortfall": horizon_shortfall,
            "strict_admission_present": True,
            "admission_passed": admission.get("passed"),
            "retry_entries": len(retries),
            "first_nominal_retry_matches": retry_check["first_nominal_retry_matches"],
            "archived_seed_unchanged": digest_after == digest_before,
            "captured_rhs_evaluations": capture_summary["captured_rhs_evaluations"],
            "rhs_capture_identity_pass": capture_summary["rhs_capture_identity_pass"],
            "rhs_capture_unique_key_count": capture_summary["rhs_capture_unique_key_count"],
            "rhs_capture_record_count": capture_summary["rhs_capture_record_count"],
            "failed_predictor_rhs_identity_pass": capture_summary["failed_predictor_rhs_identity_pass"],
            "failed_predictor_rhs_identity_match_count": capture_summary["failed_predictor_rhs_identity_match_count"],
            "failed_predictor_rhs_identity_mismatch_count": capture_summary["failed_predictor_rhs_identity_mismatch_count"],
            "failed_explicit_predictor_candidates": capture_summary["failed_explicit_predictor_candidates"],
            "radiation_rhs_component_capture_calls": capture_summary["radiation_rhs_component_capture_calls"],
            "captured_geometry_stage_terms": capture_summary["captured_geometry_stage_terms"],
            "geometry_source_term_capture_records": len(capture_summary["geometry_source_term_capture_records"]),
            "l3_total_stress_projection_oracle": capture_summary["l3_total_stress_projection_oracle"],
            "outer_boundary_capture_calls": capture_summary["outer_boundary_capture_calls"],
            "stage_capture_sha256": capture_summary_hash,
            "stage_capture_profiles_sha256": capture_profiles_hash,
            "report_sha256": report_hash,
        }, indent=2, allow_nan=False), flush=True)
        return 0 if status in {
            "bounded_probe_completed_with_strict_admission",
            "bounded_probe_records_radiation_noncompletion_fail_closed",
        } else 2

    except Exception as exc:
        report["status"] = "audit_failed"
        report["error"] = {"class": type(exc).__name__, "message": str(exc)}
        report["elapsed_seconds"] = time.time() - started
        report_path = args.output_dir / "report.json"
        report_path.write_text(json.dumps(report, indent=2, allow_nan=False))
        report_hash = sha256(report_path)
        (args.output_dir / "SHA256SUMS.txt").write_text(f"{report_hash}  report.json\n")
        print(json.dumps({
            "status": report["status"],
            "error": report["error"],
            "report_sha256": report_hash,
        }, indent=2, allow_nan=False), flush=True)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
