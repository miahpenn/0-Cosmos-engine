#!/usr/bin/env python3
"""Canonical zero-step full-Bianchi evaluator for archived Issue #4 slices.

The evaluator consumes checksum-verified frozen-state fixtures. It imports only
the pinned source/vendor components needed for geometry and stress projection;
it does not construct a kernel, evaluate an RHS, call a step, or mutate state.
"""
from __future__ import annotations

import importlib.util
import math
import os
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from tools.issue4_accepted_state_fixture import load_accepted_fixture
from tools.issue4_zero_step_full_bianchi_observer import verify_provenance, require_unchanged


def _finite(name: str, value: Any) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return result


def _load_reviewed_total_stress_helper(diagnostic_root: Path, baseline_root: Path):
    """Load the accepted-anchor helper while preserving its separate baseline pin."""
    path = diagnostic_root / "tools" / "audit_total_stress_energy_divergence.py"
    if not path.is_file():
        raise RuntimeError(f"pinned total-stress helper is missing: {path}")
    old = os.environ.get("COSMOS_PINNED_SOURCE_ROOT")
    os.environ["COSMOS_PINNED_SOURCE_ROOT"] = str(baseline_root.resolve())
    try:
        spec = importlib.util.spec_from_file_location(
            "_issue4_accepted_anchor_total_stress_helper", path
        )
        if spec is None or spec.loader is None:
            raise ImportError(f"cannot load pinned total-stress helper: {path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        return module
    finally:
        if old is None:
            os.environ.pop("COSMOS_PINNED_SOURCE_ROOT", None)
        else:
            os.environ["COSMOS_PINNED_SOURCE_ROOT"] = old


def _load_vendor(diagnostic_root: Path):
    vendor_root = diagnostic_root / "vendor" / "bb-palatini-unified-r0"
    if not vendor_root.is_dir():
        raise RuntimeError(f"pinned vendor checkout is missing: {vendor_root}")
    path = str(vendor_root)
    if path not in sys.path:
        sys.path.insert(0, path)
    import r0_operators
    import vacuum_bssn
    import moving_puncture_pirk
    return r0_operators, vacuum_bssn, moving_puncture_pirk


def _state_from_fixture(fixture: dict[str, Any], r0_operators, vacuum_bssn):
    arrays = fixture["state_arrays"]
    n = fixture["resolution"]
    r = _finite("archived grid centres", arrays["grid.centers"])
    if r.ndim != 1 or r.size != n or n < 5:
        raise ValueError("accepted fixture grid shape is invalid")
    dr = float(r[0] * 2.0)
    r_max = dr * n
    grid = r0_operators.SphericalCellGrid(n, r_max)
    if not np.array_equal(np.asarray(grid.centers), r):
        raise RuntimeError("pinned vendor grid does not exactly reproduce archived centres")

    geometry = vacuum_bssn.VacuumState(**{
        name: _finite(f"geometry.{name}", arrays[f"geometry.{name}"]).copy()
        for name in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
    })
    scalars = SimpleNamespace(**{
        name: _finite(f"scalars.{name}", arrays[f"scalars.{name}"]).copy()
        for name in ("S", "PS", "D", "PD", "phi", "Pi")
    })
    matter = SimpleNamespace(**{
        species: SimpleNamespace(**{
            field: _finite(
                f"matter.{species}.{field}",
                arrays[f"matter.{species}.{field}"],
            ).copy()
            for field in ("rest", "energy_t", "momentum_r")
        })
        for species in ("dark_matter", "baryons", "radiation")
    })
    state = SimpleNamespace(
        grid=grid,
        geometry=geometry,
        scalars=scalars,
        matter=matter,
        t=float(fixture["accepted_time"]),
        electromagnetic=None,
    )
    return state


def _metric_connection_lambda(grid, geometry) -> np.ndarray:
    d1 = grid.cell_derivative_fourth
    r = np.asarray(grid.centers, dtype=np.float64)
    a = _finite("a", geometry.a)
    b = _finite("b", geometry.b)
    da = d1(a, parity=1)
    db = d1(b, parity=1)
    # Exact reference-metric connection formula used by the pinned vendor.
    return _finite(
        "Lambda_metric",
        da / (2.0 * a**2) - db / (a * b) + 2.0 * (1.0 / b - 1.0 / a) / r,
    )


def _summary(values: np.ndarray, r: np.ndarray) -> dict[str, Any]:
    x = _finite("summary profile", values)
    inner = r <= 10.0
    annulus = (r >= 5.0) & (r <= 10.0)
    if not np.any(inner) or not np.any(annulus):
        raise ValueError("accepted comparison masks select no cells")
    idx = int(np.flatnonzero(inner)[np.argmax(np.abs(x[inner]))])
    return {
        "max_abs_r_le_10": float(abs(x[idx])),
        "signed_at_max_r_le_10": float(x[idx]),
        "index_at_max_r_le_10": idx,
        "radius_at_max_r_le_10": float(r[idx]),
        "max_abs_annulus_5_to_10": float(np.max(np.abs(x[annulus]))),
        "rms_annulus_5_to_10": float(np.sqrt(np.mean(x[annulus] ** 2))),
        "signed_cells_0_1_2": [float(x[i]) for i in range(3)],
        "fixed_radius_samples": {
            str(target): float(np.interp(target, r, x))
            for target in (2.0, 3.0, 5.0, 8.0)
        },
    }


def _array_for_path(fixture: dict[str, Any], path: str) -> np.ndarray:
    records = fixture["report"]["raw_arrays"]["array_records"]
    record = next((item for item in records if item["path"] == path), None)
    if record is None:
        raise KeyError(f"accepted fixture lacks required archived array: {path}")
    return _finite(path, fixture["arrays"][record["key"]])


def evaluate_accepted_fixture(
    archive_path: str | Path,
    resolution: int,
    diagnostic_root: str | Path,
    helper_baseline_root: str | Path,
) -> dict[str, Any]:
    """Evaluate one accepted frozen slice under the single locked predictor."""
    diagnostic_root = Path(diagnostic_root).resolve()
    helper_baseline_root = Path(helper_baseline_root).resolve()

    # Provenance is checked before fixture loading, imports, or state evaluation.
    provenance = verify_provenance(diagnostic_root, helper_baseline_root)
    fixture = load_accepted_fixture(archive_path, resolution)
    r0_operators, vacuum_bssn, moving_puncture_pirk = _load_vendor(diagnostic_root)
    total_helper = _load_reviewed_total_stress_helper(
        diagnostic_root, helper_baseline_root
    )

    # In-memory tripwires: any accidental production/PIRK advance fails closed.
    step_calls = []
    kernel_class = total_helper.production_kernel_module.V55ProductionKernel
    if not hasattr(kernel_class, "step"):
        raise RuntimeError("pinned production kernel exposes no step method to guard")

    def forbidden_production_step(*args, **kwargs):
        step_calls.append("ProductionKernel.step")
        raise RuntimeError("zero-step observer forbids ProductionKernel.step")

    kernel_class.step = forbidden_production_step
    if not hasattr(moving_puncture_pirk, "damped_moving_puncture_pirk2_step"):
        raise RuntimeError("pinned vendor has no damped_moving_puncture_pirk2_step method to guard")
    def forbidden_vendor_pirk(*args, **kwargs):
        step_calls.append("vendor.damped_moving_puncture_pirk2_step")
        raise RuntimeError("zero-step observer forbids vendor PIRK advance")
    moving_puncture_pirk.damped_moving_puncture_pirk2_step = forbidden_vendor_pirk

    state = _state_from_fixture(fixture, r0_operators, vacuum_bssn)
    helper_digest_before = total_helper.state_digest(state)
    require_unchanged(fixture["helper_state_digest"], helper_digest_before)

    grid = state.grid
    geom = state.geometry
    d1 = grid.cell_derivative_fourth
    r = np.asarray(grid.centers, dtype=np.float64)
    alpha = _finite("alpha", geom.alpha)
    beta = _finite("beta", geom.beta)
    a = _finite("a", geom.a)
    b = _finite("b", geom.b)
    X = _finite("X", geom.X)
    if np.any(beta != 0.0):
        raise RuntimeError("this accepted-state predictor fixture requires exactly zero shift")

    # Metric-only geometric constraints: replace evolved Lambda only on a
    # disposable geometry copy. The accepted input state remains untouched.
    lambda_metric = _metric_connection_lambda(grid, geom)
    metric_geometry = geom.copy()
    metric_geometry.Lambda = lambda_metric.copy()
    metric_constraints = vacuum_bssn.constraints(grid, metric_geometry)
    H_geom = _finite("metric-only Hamiltonian constraint", metric_constraints["hamiltonian"])
    M_geom = _finite("geometric lower-index momentum constraint", metric_constraints["momentum"])

    # Independently reconstruct all six stress sectors using the reviewed helper.
    projection = total_helper.independent_total_projections(state)
    if not projection.get("projection_crosscheck_passed", False):
        raise RuntimeError("independent total-stress projection cross-check failed")
    rho_total = _finite("rho_total", projection["rho"])
    j_total = _finite("j_total", projection["j"])
    H_total = _finite("H_total", H_geom - 16.0 * math.pi * rho_total)
    M_total = _finite("M_total", M_geom - 8.0 * math.pi * j_total)

    # Keep chi=X^2 distinct from chi_r=-D_1 X/(2X).
    chi = X**2
    chi_r = -d1(X, parity=1) / (2.0 * X)
    C_r = _finite("odd connection difference C^r", geom.Lambda - lambda_metric)
    dC_odd = d1(C_r, parity=-1)
    da_even = d1(a, parity=1)
    W_rr = _finite(
        "tracefree evolved-connection response W^r_r",
        chi * ((2.0 / 3.0) * dC_odd + (da_even / (2.0 * a)) * C_r),
    )
    F_C = _finite("F_C=alpha W^r_r", alpha * W_rr)
    connection_divergence = _finite(
        "spherical divergence D_j(alpha W^j_r)",
        d1(F_C, parity=1)
        + (6.0 * chi_r + 3.0 / r
           + 3.0 * d1(b, parity=1) / (2.0 * b)) * F_C,
    )

    # The accepted slice has beta^r identically zero; retain its measured,
    # separately named transport contribution rather than omitting the term.
    shift_transport = np.zeros_like(r)
    normal_alpha_K_M = _finite("alpha K M_total", alpha * geom.K * M_total)
    normal_lapse_gradient_H = _finite(
        "-(D_r alpha) H_total / 3",
        -(1.0 / 3.0) * d1(alpha, parity=1) * H_total,
    )
    normal_H_gradient = _finite(
        "alpha D_r H_total / 6",
        (1.0 / 6.0) * alpha * d1(H_total, parity=1),
    )

    terms_by_channel = {}
    for epsilon_label, epsilon_value in (("eps2e-05", 2.0e-5), ("eps1e-05", 1.0e-5)):
        tangent = _array_for_path(fixture, f"tangent_{epsilon_label}.Mdot")
        div_pointwise = _array_for_path(
            fixture, f"total_T_{epsilon_label}.total.div_T_r"
        )
        div_fv = _array_for_path(
            fixture, f"total_T_{epsilon_label}.total.div_T_r_fv"
        )
        common = (
            shift_transport + normal_alpha_K_M + normal_lapse_gradient_H
            + normal_H_gradient + connection_divergence
        )
        pointwise_source = _finite(
            "pointwise total-stress Bianchi source",
            -8.0 * math.pi * alpha * div_pointwise,
        )
        fv_source = _finite(
            "finite-volume total-stress Bianchi source",
            -8.0 * math.pi * alpha * div_fv,
        )
        pointwise_predictor = _finite("pointwise full Bianchi predictor", common + pointwise_source)
        fv_predictor = _finite("finite-volume full Bianchi predictor", common + fv_source)
        terms_by_channel[epsilon_label] = {
            "epsilon": epsilon_value,
            "measured_tangent": tangent,
            "pointwise": {
                "div_T_r": div_pointwise,
                "stress_source": pointwise_source,
                "predictor": pointwise_predictor,
                "residual_tangent_minus_predictor": _finite(
                    "pointwise tangent minus predictor", tangent - pointwise_predictor
                ),
            },
            "finite_volume": {
                "div_T_r_fv": div_fv,
                "stress_source": fv_source,
                "predictor": fv_predictor,
                "residual_tangent_minus_predictor": _finite(
                    "finite-volume tangent minus predictor", tangent - fv_predictor
                ),
            },
        }

    helper_digest_after = total_helper.state_digest(state)
    require_unchanged(helper_digest_before, helper_digest_after)
    if helper_digest_after != fixture["helper_state_digest"]:
        raise RuntimeError("accepted input state changed during predictor evaluation")
    if step_calls:
        raise RuntimeError(f"zero-step tripwire was activated: {step_calls}")

    summaries = {}
    for epsilon_label, channels in terms_by_channel.items():
        summaries[epsilon_label] = {}
        for channel_name in ("pointwise", "finite_volume"):
            channel = channels[channel_name]
            summaries[epsilon_label][channel_name] = {
                "predictor": _summary(channel["predictor"], r),
                "residual_tangent_minus_predictor": _summary(
                    channel["residual_tangent_minus_predictor"], r
                ),
            }

    return {
        "schema": "issue4-full-bianchi-accepted-fixture-evaluation/v1",
        "status": "diagnostic_evaluation_complete",
        "state_class": "accepted_t_0p60",
        "resolution": resolution,
        "accepted_time": fixture["accepted_time"],
        "provenance": provenance,
        "input_digests": {
            "upstream_archive_state": fixture["archive_state_digest"],
            "tangent_helper_state": fixture["helper_state_digest"],
        },
        "projection_crosscheck": projection["projection_crosscheck"],
        "structural_inventory": fixture["structural_inventory"],
        "signed_terms": {
            "shift_transport": shift_transport,
            "normal_alpha_K_M": normal_alpha_K_M,
            "normal_lapse_gradient_H": normal_lapse_gradient_H,
            "normal_H_gradient": normal_H_gradient,
            "connection_divergence": connection_divergence,
            "chi": chi,
            "chi_r": chi_r,
            "C_r": C_r,
            "W_rr": W_rr,
            "F_C": F_C,
            "H_geom": H_geom,
            "H_total": H_total,
            "M_geom": M_geom,
            "M_total": M_total,
            "rho_total": rho_total,
            "j_total": j_total,
        },
        "channels": terms_by_channel,
        "summaries": summaries,
        "input_state_digest_before": helper_digest_before,
        "input_state_digest_after": helper_digest_after,
        "production_steps_run": len(step_calls),
        "state_advanced": False,
        "zero_step_tripwire": "armed_and_not_triggered",
        "inherited_assumptions": [
            "locked projected-Bianchi continuum identity and sign/index conventions",
            "pinned vendor BSSN constraint mapping",
            "metric-only Ricci with evolved-Lambda response retained separately",
            "six-sector stress projections and accepted total-stress divergence definitions",
            "pointwise and finite-volume total-stress divergences are separate comparators",
        ],
    }
