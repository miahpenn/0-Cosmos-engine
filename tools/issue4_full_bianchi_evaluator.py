#!/usr/bin/env python3
"""Canonical zero-step full-Bianchi evaluator for archived Issue #4 slices.

The evaluator consumes checksum-verified frozen-state fixtures. It imports only
the pinned source/vendor components needed for geometry and stress projection;
it does not construct a kernel, evaluate an RHS, call a step, or mutate state.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import importlib.util
import json
import math
import os
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np

from tools.issue4_accepted_state_fixture import load_accepted_fixture, sector_accounting_summary
from tools.issue4_manufactured_fixture import load_manufactured_fixture
from tools.issue4_manufactured_regular_family import STATE_ARRAY_NAMES
from tools.issue4_bianchi_reference_regression import assert_reference_regression
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


@contextmanager
def _zero_step_guard(kernel_class, vendor_module, step_calls):
    """Arm both evolution tripwires and restore the original methods on every exit."""
    vendor_name = "damped_moving_puncture_pirk2_step"
    if not hasattr(kernel_class, "step"):
        raise RuntimeError("pinned production kernel exposes no step method to guard")
    if not hasattr(vendor_module, vendor_name):
        raise RuntimeError(f"pinned vendor has no {vendor_name} method to guard")
    original_kernel_step = kernel_class.step
    original_vendor_pirk = getattr(vendor_module, vendor_name)

    def forbidden_production_step(*args, **kwargs):
        step_calls.append("ProductionKernel.step")
        raise RuntimeError("zero-step observer forbids ProductionKernel.step")

    def forbidden_vendor_pirk(*args, **kwargs):
        step_calls.append("vendor.damped_moving_puncture_pirk2_step")
        raise RuntimeError("zero-step observer forbids vendor PIRK advance")

    try:
        kernel_class.step = forbidden_production_step
        setattr(vendor_module, vendor_name, forbidden_vendor_pirk)
        yield
    finally:
        kernel_class.step = original_kernel_step
        setattr(vendor_module, vendor_name, original_vendor_pirk)


def _load_vendor(diagnostic_root: Path):
    """Import vendor modules only when their actual origins match the pinned tree."""
    vendor_root = (diagnostic_root / "vendor" / "bb-palatini-unified-r0").resolve()
    if not vendor_root.is_dir():
        raise RuntimeError(f"pinned vendor checkout is missing: {vendor_root}")
    path = str(vendor_root)
    if path not in sys.path:
        sys.path.insert(0, path)

    def require_pinned_origin(module):
        origin = getattr(module, "__file__", None)
        if not origin:
            raise RuntimeError(
                f"pinned vendor module {module.__name__!r} has no verifiable __file__"
            )
        resolved = Path(origin).resolve()
        try:
            resolved.relative_to(vendor_root)
        except ValueError as exc:
            raise RuntimeError(
                f"vendor module origin mismatch for {module.__name__}: "
                f"{resolved} is outside pinned root {vendor_root}"
            ) from exc
        return module

    modules = []
    for name in ("r0_operators", "vacuum_bssn", "moving_puncture_pirk"):
        cached = sys.modules.get(name)
        if cached is not None:
            # Refuse rather than silently trusting or replacing a preloaded module.
            modules.append(require_pinned_origin(cached))
            continue
        module = __import__(name)
        modules.append(require_pinned_origin(module))
    return tuple(modules)


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
        **{
            species: SimpleNamespace(**{
                field: _finite(
                    f"matter.{species}.{field}",
                    arrays[f"matter.{species}.{field}"],
                ).copy()
                for field in ("rest", "energy_t", "momentum_r")
            })
            for species in ("dark_matter", "baryons", "radiation")
        },
        # The pinned V55MatterState contract includes this opt-in carrier.
        # The accepted archive has no charged-dust arrays, so preserve absence
        # explicitly; do not synthesize a zero-density species or EM source.
        "charged_dust": None,
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
    """Evaluate one accepted frozen slice; accepted loader/regression unchanged."""
    diagnostic_root = Path(diagnostic_root).resolve()
    helper_baseline_root = Path(helper_baseline_root).resolve()
    provenance = verify_provenance(diagnostic_root, helper_baseline_root)
    fixture = load_accepted_fixture(archive_path, resolution)
    return _evaluate_loaded_fixture(
        fixture, resolution, diagnostic_root, helper_baseline_root, provenance,
        state_class="accepted_t_0p60", tangent_semantics="production_rhs_directional",
    )


def manufactured_sector_accounting_summary(report: dict[str, Any]) -> list[dict[str, Any]]:
    """Compact the manufactured helper's raw and time-density-accounted sector ledger."""
    rows = []
    for label, item in report["divergence_reports"].items():
        channels = {
            name: {
                "raw_max_abs_gap": value["max_abs_gap"],
                "max_abs_unexplained_gap": value["max_abs_unexplained_gap"],
                "raw_sum_passed": value["passed"],
                "component_corrected_passed": value["component_corrected_passed"],
            }
            for name, value in item["channel_accounting"].items()
        }
        rows.append({
            "epsilon_label": label,
            "epsilon": 2.0e-5 if label == "eps2e-05" else 1.0e-5,
            "sector_names": item["sector_names"],
            "raw_six_sector_sum_reconstructs_total": item["raw_six_sector_reconstruction_passed"],
            "time_density_component_accounting_passed": item["time_density_component_accounting_passed"],
            "six_sector_reconstruction_after_accounting": item["six_sector_reconstruction_after_time_density_accounting"],
            "raw_sum_max_gap": item["raw_sum_max_gap"],
            "max_abs_unexplained_gap_after_accounting": item["unexplained_max_gap_after_time_density_accounting"],
            "channels": channels,
        })
    return rows


def evaluate_manufactured_fixture(
    archive_path: str | Path,
    resolution: int,
    diagnostic_root: str | Path,
    helper_baseline_root: str | Path,
) -> dict[str, Any]:
    """Evaluate the manufactured class through the same locked predictor."""
    if resolution != 60:
        raise ValueError("manufactured regular-family fixture is pinned to N=60")
    diagnostic_root = Path(diagnostic_root).resolve()
    helper_baseline_root = Path(helper_baseline_root).resolve()
    provenance = verify_provenance(diagnostic_root, helper_baseline_root)
    loaded = load_manufactured_fixture(archive_path)
    report = loaded["report"]
    arrays = loaded["arrays"]
    state_arrays = {
        name: np.asarray(arrays[f"state.base.{name}"], dtype=np.float64)
        for name in STATE_ARRAY_NAMES
    }
    records = report["raw_arrays"]["array_records"]
    paths = {record["path"] for record in records}
    required = {
        "frozen_input_state": {f"state.base.{name}" for name in STATE_ARRAY_NAMES},
        "analytic_family_tangent": {
            "tangent_eps2e-05.Mdot", "tangent_eps1e-05.Mdot",
        },
        "pointwise_total_stress_divergence": {
            "total_T_eps2e-05.total.div_T_r", "total_T_eps1e-05.total.div_T_r",
        },
        "finite_volume_total_stress_divergence": {
            "total_T_eps2e-05.total.div_T_r_fv", "total_T_eps1e-05.total.div_T_r_fv",
        },
    }
    structural_inventory = {
        "lens": "manufactured_fixture_structural_inventory",
        "stored_operand_families": {
            name: sorted(expected.intersection(paths)) for name, expected in required.items()
        },
        "missing_required_witnesses": {
            name: sorted(expected.difference(paths)) for name, expected in required.items()
        },
        "interpretation": (
            "The five analytic family states, six-sector projections, manufactured "
            "central-difference tangent, and pointwise/FV divergence channels are "
            "stored in the verified fixture archive."
        ),
    }
    missing = {
        name: values for name, values in structural_inventory["missing_required_witnesses"].items()
        if values
    }
    if missing:
        raise RuntimeError(f"manufactured fixture structural inventory is incomplete: {missing}")
    fixture = {
        "resolution": resolution,
        "accepted_time": 0.0,
        "state_arrays": state_arrays,
        "arrays": arrays,
        "report": report,
        "helper_state_digest": report["state_digests"]["base"],
        "archive_state_digest": report["state_digests"]["base"],
        "structural_inventory": structural_inventory,
        "state_class": "manufactured_regular_family",
        "tangent_semantics": report["tangent_semantics"],
    }
    return _evaluate_loaded_fixture(
        fixture, resolution, diagnostic_root, helper_baseline_root, provenance,
        state_class="manufactured_regular_family",
        tangent_semantics="analytic_family_central_difference",
    )


def _evaluate_loaded_fixture(
    fixture: dict[str, Any],
    resolution: int,
    diagnostic_root: Path,
    helper_baseline_root: Path,
    provenance: dict[str, Any],
    *,
    state_class: str,
    tangent_semantics: str,
) -> dict[str, Any]:
    r0_operators, vacuum_bssn, moving_puncture_pirk = _load_vendor(diagnostic_root)
    total_helper = _load_reviewed_total_stress_helper(
        diagnostic_root, helper_baseline_root
    )

    # In-memory tripwires are restored even if any later evaluation raises.
    step_calls = []
    kernel_class = total_helper.production_kernel_module.V55ProductionKernel
    with _zero_step_guard(kernel_class, moving_puncture_pirk, step_calls):
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
        dC_odd = _finite("odd-parity derivative of C^r", d1(C_r, parity=-1))
        da_even = _finite("even-parity derivative of a", d1(a, parity=1))
        db_even = _finite("even-parity derivative of b", d1(b, parity=1))
        d_alpha_even = _finite("even-parity derivative of alpha", d1(alpha, parity=1))
        d_H_total_even = _finite("even-parity derivative of H_total", d1(H_total, parity=1))
        W_rr = _finite(
            "tracefree evolved-connection response W^r_r",
            chi * ((2.0 / 3.0) * dC_odd + (da_even / (2.0 * a)) * C_r),
        )
        ricci_evolved = vacuum_bssn.geometry_terms(grid, geom)
        ricci_metric = vacuum_bssn.geometry_terms(grid, metric_geometry)
        delta_R_scalar = _finite("evolved-minus-metric scalar Ricci", ricci_evolved["R"] - ricci_metric["R"])
        delta_Rr_r = _finite("evolved-minus-metric radial Ricci", ricci_evolved["Rr_r"] - ricci_metric["Rr_r"])
        expected_delta_R_scalar = _finite("expected scalar Ricci connection offset", chi * dC_odd)
        expected_delta_Rr_r = _finite("expected radial Ricci connection offset", chi * (dC_odd + (da_even / (2.0 * a)) * C_r))
        W_rr_direct = _finite("direct tracefree Ricci connection offset", delta_Rr_r - delta_R_scalar / 3.0)
        ricci_offset_crosscheck = {
            "scalar_max_abs_gap": float(np.max(np.abs(delta_R_scalar - expected_delta_R_scalar))),
            "radial_max_abs_gap": float(np.max(np.abs(delta_Rr_r - expected_delta_Rr_r))),
            "tracefree_max_abs_gap": float(np.max(np.abs(W_rr - W_rr_direct))),
        }
        if max(ricci_offset_crosscheck.values()) > 1.0e-17:
            raise RuntimeError(f"same-metric Ricci connection-offset cross-check failed: {ricci_offset_crosscheck}")
        F_C = _finite("F_C=alpha W^r_r", alpha * W_rr)
        connection_divergence = _finite(
            "spherical divergence D_j(alpha W^j_r)",
            d1(F_C, parity=1)
            + (6.0 * chi_r + 3.0 / r
               + 3.0 * db_even / (2.0 * b)) * F_C,
        )

        # The accepted slice has beta^r identically zero; retain its measured,
        # separately named transport contribution rather than omitting the term.
        shift_transport = np.zeros_like(r)
        normal_alpha_K_M = _finite("alpha K M_total", alpha * geom.K * M_total)
        normal_lapse_gradient_H = _finite(
            "-(D_r alpha) H_total / 3",
            -(1.0 / 3.0) * d_alpha_even * H_total,
        )
        normal_H_gradient = _finite(
            "alpha D_r H_total / 6",
            (1.0 / 6.0) * alpha * d_H_total_even,
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

        result = {
            "schema": f"issue4-full-bianchi-{state_class}-evaluation/v1",
            "status": "diagnostic_evaluation_complete",
            "state_class": state_class,
            "tangent_semantics": tangent_semantics,
            "resolution": resolution,
            "grid_centers": r.copy(),
            "state_time": fixture["accepted_time"],
            "accepted_time": fixture["accepted_time"] if state_class == "accepted_t_0p60" else None,
            "provenance": provenance,
            "input_digests": {
                "upstream_archive_state": fixture["archive_state_digest"],
                "tangent_helper_state": fixture["helper_state_digest"],
            },
            "projection_crosscheck": projection["projection_crosscheck"],
            "ricci_offset_crosscheck": ricci_offset_crosscheck,
            "structural_inventory": fixture["structural_inventory"],
            "six_sector_divergence_accounting": (
                sector_accounting_summary(fixture["report"])
                if state_class == "accepted_t_0p60"
                else manufactured_sector_accounting_summary(fixture["report"])
            ),
            "signed_terms": {
                "shift_transport": shift_transport,
                "normal_alpha_K_M": normal_alpha_K_M,
                "normal_lapse_gradient_H": normal_lapse_gradient_H,
                "normal_H_gradient": normal_H_gradient,
                "connection_divergence": connection_divergence,
                "chi": chi,
                "chi_r": chi_r,
                "C_r": C_r,
                "dC_odd": dC_odd,
                "da_even": da_even,
                "db_even": db_even,
                "d_alpha_even": d_alpha_even,
                "d_H_total_even": d_H_total_even,
                "Lambda_metric": lambda_metric,
                "beta": beta,
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
        if state_class == "accepted_t_0p60":
            result["accepted_reference_regression"] = assert_reference_regression(result)
        else:
            result["accepted_reference_regression"] = {
                "executed": False,
                "reason": "manufactured class has a separate report; accepted-state regression remains unchanged",
            }
        return result



def _jsonable(value: Any):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _write_json_payload(path: Path, payload: str) -> None:
    """Write one complete JSON document with a real trailing line terminator."""
    path.write_text(payload + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--resolution", type=int, choices=(60, 120, 240), required=True)
    parser.add_argument("--state-class", choices=("accepted", "manufactured"), default="accepted")
    parser.add_argument("--diagnostic-source", type=Path, required=True)
    parser.add_argument("--helper-baseline", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        if args.state_class == "manufactured":
            result = evaluate_manufactured_fixture(
                args.archive, args.resolution, args.diagnostic_source, args.helper_baseline
            )
        else:
            result = evaluate_accepted_fixture(
                args.archive, args.resolution, args.diagnostic_source, args.helper_baseline
            )
        payload = json.dumps(_jsonable(result), indent=2, sort_keys=True, allow_nan=False)
    except Exception as exc:
        payload = json.dumps({
            "schema": f"issue4-full-bianchi-{args.state_class}-evaluation/v1",
            "state_class": "manufactured_regular_family" if args.state_class == "manufactured" else "accepted_t_0p60",
            "status": "refused",
            "reason": f"{type(exc).__name__}: {exc}",
            "production_steps_run": 0,
            "state_advanced": False,
        }, indent=2, sort_keys=True)
        if args.output is None:
            print(payload)
        else:
            _write_json_payload(args.output, payload)
        return 2
    if args.output is None:
        print(payload)
    else:
            _write_json_payload(args.output, payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
