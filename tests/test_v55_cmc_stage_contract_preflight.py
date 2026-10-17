"""Observation-only production CMC stage-contract preflight.

This exercises copied Stage-0, explicit-predictor, and primary-predictor slices
without calling V55ProductionKernel.step. The production step's ordering is
also checked structurally, so this test can detect contract drift without
reimplementing a full production step or duplicating the existing resolution
sweeps in #339/#340. Native-minus-reference residuals are descriptive only.
"""

import ast
import hashlib
import inspect
import json
import textwrap

import numpy as np

from engine import v55_matter
from engine.matter_system import Species, validate_species_admissibility
from engine.native_geometry import (
    NativeSphericalCellGrid,
    native_bssn_geometry_terms,
    native_explicit_rhs,
    native_lambda_l2_rhs,
    native_lambda_l3_rhs,
    native_primary_l2_rhs,
    native_primary_l3_rhs,
)
from engine.production_kernel import LAMBDA_M, ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.v55_pirk_adapter import vendor_modules


_GEOMETRY_FIELDS = (
    "a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B",
)
_SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
_MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
_MATTER_SPECIES = ("dark_matter", "baryons", "radiation", "charged_dust")
_EXPLICIT_PREDICTOR_FIELDS = ("a", "b", "X")
_PRIMARY_FIELDS = ("Aa", "K")
_LAMBDA_M = LAMBDA_M


def _slice_sha256(geometry, scalars, matter):
    """Hash copied/accepted profile arrays with dtype and shape provenance."""
    digest = hashlib.sha256()
    for label, obj, names in (
        ("geometry", geometry, _GEOMETRY_FIELDS),
        ("scalars", scalars, _SCALAR_FIELDS),
    ):
        for name in names:
            values = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            digest.update(label.encode("ascii") + b"\0" + name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    for species in _MATTER_SPECIES:
        obj = getattr(matter, species, None)
        if obj is None:
            digest.update(b"matter\0" + species.encode("ascii") + b"\0none\0")
            continue
        for name in _MATTER_FIELDS:
            values = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            digest.update(b"matter\0" + species.encode("ascii") + b"\0")
            digest.update(name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _production_stage_call_order():
    """Lock the sequence of solve/RHS/boundary calls without executing step()."""
    source = textwrap.dedent(inspect.getsource(V55ProductionKernel.step))
    tree = ast.parse(source)
    observed = []
    relevant = {
        "_solve_lapse_for_state",
        "_geometry_stage_terms",
        "_apply_outer_boundary_for_state",
    }
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr in relevant:
            observed.append((node.lineno, node.col_offset, node.func.attr))
    observed.sort()
    names = [name for _, _, name in observed]
    expected = [
        "_solve_lapse_for_state",
        "_geometry_stage_terms",
        "_solve_lapse_for_state",
        "_geometry_stage_terms",
        "_apply_outer_boundary_for_state",
        "_solve_lapse_for_state",
        "_geometry_stage_terms",
        "_apply_outer_boundary_for_state",
        "_solve_lapse_for_state",
    ]
    assert names == expected, (
        "production CMC stage order changed; reconcile the stage contract "
        f"before substituting a native component: {names}"
    )

    # Also lock gauge reset and boundary placement around those stage calls.
    fragments = (
        "g0.alpha = self._solve_lapse_for_state(",
        "g0.beta.fill(0.0)",
        "g0.B.fill(0.0)",
        'gterms0 = self._geometry_stage_terms(',
        "g_explicit1.beta.fill(0.0)",
        "g_explicit1.B.fill(0.0)",
        "g_explicit1.alpha = self._solve_lapse_for_state(",
        "gterms_pred = self._geometry_stage_terms(",
        "self._apply_outer_boundary_for_state(grid, g1, s1, m1, em1)",
        "g1.alpha = self._solve_lapse_for_state(",
        "g1.beta.fill(0.0)",
        "g1.B.fill(0.0)",
        "gterms1 = self._geometry_stage_terms(",
        "self._apply_outer_boundary_for_state(grid, gnew, snew, mnew, emnew)",
        "gnew.alpha = self._solve_lapse_for_state(",
        "gnew.beta.fill(0.0)",
        "gnew.B.fill(0.0)",
    )
    cursor = 0
    for fragment in fragments:
        position = source.find(fragment, cursor)
        assert position >= 0, (
            "production stage contract fragment moved or changed: "
            f"{fragment!r}"
        )
        cursor = position + len(fragment)
    return names


def _regional_profile_summary(values, masks):
    values = np.asarray(values, dtype=float)
    assert values.ndim == 1 and np.all(np.isfinite(values))
    summary = {}
    for name, mask in masks.items():
        selected = values[mask]
        # A few cells are needed to make the descriptive profile meaningful;
        # this is not a residual acceptance threshold.
        assert selected.size >= 4, f"{name}: insufficient cells for a profile"
        summary[name] = {
            "cells": int(selected.size),
            "rms": float(np.sqrt(np.mean(selected**2))),
            "max_abs": float(np.max(np.abs(selected))),
        }
    return summary


def _field_comparison(native, reference, field_names, masks):
    summary = {}
    for name in field_names:
        lhs = np.asarray(native[name], dtype=float)
        rhs = np.asarray(reference[name], dtype=float)
        assert lhs.shape == rhs.shape, f"{name}: native/reference shape mismatch"
        residual = lhs - rhs
        assert np.all(np.isfinite(residual)), f"{name}: non-finite signed residual"
        summary[name] = _regional_profile_summary(residual, masks)
    return summary


def _compare_stage_rhs(grid, geometry, scalars, matter, native_grid, terms):
    """Compare only same-input compatible geometry blocks; report source increments apart."""
    _, vacuum, _ = vendor_modules()
    native = {
        "explicit": native_explicit_rhs(native_grid, geometry),
        "primary_l2": native_primary_l2_rhs(native_grid, geometry),
        "primary_l3": native_primary_l3_rhs(native_grid, geometry),
        "lambda_l2": native_lambda_l2_rhs(
            native_grid, geometry, lambda_m=_LAMBDA_M
        ),
        "lambda_l3": native_lambda_l3_rhs(native_grid, geometry),
    }
    reference_primary_l3_geometry = vacuum.primary_l3_rhs(grid, geometry)
    reference_lambda_l3_geometry = np.asarray(
        vacuum.lambda_l3_rhs(grid, geometry), dtype=float
    )
    projection = v55_matter.total_matter_projection(
        grid, geometry, scalars, matter
    )
    primary_source = {
        "Aa": (
            -(16.0 * np.pi / 3.0)
            * geometry.alpha
            * (projection["pr"] - projection["pt"])
        ),
        "K": (
            4.0 * np.pi
            * geometry.alpha
            * (projection["rho"] + projection["pr"] + 2.0 * projection["pt"])
        ),
    }
    lambda_current_source = (
        -8.0 * np.pi * _LAMBDA_M * geometry.alpha * projection["j"] / geometry.a
    )

    # These are implementation/source-ownership checks on the new stage-local
    # copies. The established formula is reused; it is not independently re-derived.
    np.testing.assert_allclose(
        terms["primary_l3"]["Aa"],
        np.asarray(reference_primary_l3_geometry["Aa"]) + primary_source["Aa"],
        rtol=2.0e-14, atol=1.0e-22,
    )
    np.testing.assert_allclose(
        terms["primary_l3"]["K"],
        np.asarray(reference_primary_l3_geometry["K"]) + primary_source["K"],
        rtol=2.0e-14, atol=1.0e-22,
    )
    np.testing.assert_allclose(
        terms["lambda_l3"],
        reference_lambda_l3_geometry + lambda_current_source,
        rtol=2.0e-14, atol=1.0e-24,
    )

    areal = native_bssn_geometry_terms(native_grid, geometry)["areal_radius"]
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    for name, mask in masks.items():
        assert int(np.count_nonzero(mask)) >= 4, name

    comparison = {
        # Alpha/beta RHSs are deliberately excluded: production CMC solves
        # alpha and pins beta/B instead of advancing moving-puncture gauge RHSs.
        "explicit_kinematic_a_b_X": _field_comparison(
            native["explicit"], terms["explicit"], _EXPLICIT_PREDICTOR_FIELDS, masks
        ),
        "primary_l2_Aa_K": _field_comparison(
            native["primary_l2"], terms["primary_l2"], _PRIMARY_FIELDS, masks
        ),
        "primary_l3_geometry_only_Aa_K": _field_comparison(
            native["primary_l3"], reference_primary_l3_geometry,
            _PRIMARY_FIELDS, masks,
        ),
        "lambda_l2_geometry_only": _regional_profile_summary(
            native["lambda_l2"] - np.asarray(terms["lambda_l2"]), masks
        ),
        "lambda_l3_geometry_only": _regional_profile_summary(
            native["lambda_l3"] - reference_lambda_l3_geometry, masks
        ),
    }
    source_profiles = {
        "primary_l3_Aa_matter_stress": _regional_profile_summary(
            primary_source["Aa"], masks
        ),
        "primary_l3_K_matter_stress": _regional_profile_summary(
            primary_source["K"], masks
        ),
        "lambda_l3_current": _regional_profile_summary(
            lambda_current_source, masks
        ),
        "projected_current_j": _regional_profile_summary(projection["j"], masks),
    }
    return comparison, source_profiles


def _observe_stage(
    kernel, grid, geometry, scalars, matter, native_grid, stage_name, *, compare_native
):
    before = _slice_sha256(geometry, scalars, matter)
    stage = ProductionState(
        grid=grid, geometry=geometry, scalars=scalars, matter=matter
    )
    scalar_rhs, matter_rhs, _ = kernel._rhs(stage)
    terms = kernel._geometry_stage_terms(
        grid, geometry, scalars, matter, None
    )
    comparison = None
    source_profiles = None
    if compare_native:
        comparison, source_profiles = _compare_stage_rhs(
            grid, geometry, scalars, matter, native_grid, terms
        )
    after = _slice_sha256(geometry, scalars, matter)
    assert after == before, (
        f"{stage_name}: RHS evaluation mutated its copied stage input: "
        f"before={before}, after={after}"
    )
    return scalar_rhs, matter_rhs, terms, {
        "stage": stage_name,
        "geometry_gauge": "CMC-resolved alpha; beta=0; B=0",
        "stage_input_sha256_before_rhs": before,
        "stage_input_sha256_after_rhs": after,
        "input_unchanged_by_rhs_evaluation": True,
        "rhs_comparison": (
            comparison if compare_native else
            "reuse existing #340 frozen-slice comparison; avoid duplicate matrix"
        ),
        "matter_source_profiles": source_profiles,
        "new_residual_thresholds": False,
    }


def test_v55_cmc_stage_contract_preflight_observes_copied_predictor_slices(
    monkeypatch, capsys
):
    """Observe the CMC predictor contract without executing the production step."""
    contract_events = _production_stage_call_order()
    kernel = V55ProductionKernel()
    accepted = kernel.initialize(
        resolution=48, r_max=12.0, include_radiation=True
    )
    accepted_hash_before = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )

    def forbidden_step(*_args, **_kwargs):
        raise AssertionError("CMC stage preflight must never call production .step")

    monkeypatch.setattr(kernel, "step", forbidden_step)

    grid = accepted.grid
    g0 = accepted.geometry.copy()
    s0 = accepted.scalars.copy()
    m0 = accepted.matter.copy()
    native_grid = NativeSphericalCellGrid(grid.n, grid.r_max)
    np.testing.assert_allclose(
        native_grid.centers, grid.centers, rtol=0.0, atol=1.0e-14
    )
    dt = 1.0e-5

    active_stage = ["stage0_cmc_refresh"]
    lapse_trace = []
    original_solve = kernel._solve_lapse_for_state

    def traced_solve(grid_arg, geometry_arg, scalars_arg, matter_arg, em_arg):
        result = original_solve(
            grid_arg, geometry_arg, scalars_arg, matter_arg, em_arg
        )
        alpha = np.asarray(result[0], dtype=float)
        lapse_trace.append({
            "stage": active_stage[0],
            "alpha_shape": list(alpha.shape),
            "alpha_finite": bool(np.all(np.isfinite(alpha))),
            "alpha_min": float(np.min(alpha)),
            "alpha_max": float(np.max(alpha)),
        })
        return result

    monkeypatch.setattr(kernel, "_solve_lapse_for_state", traced_solve)

    # Stage 0: resolve CMC lapse on the copied accepted slice and pin the
    # production spherical-gauge variables before the first RHS evaluation.
    g0.alpha = kernel._solve_lapse_for_state(grid, g0, s0, m0, None)[0]
    g0.beta.fill(0.0)
    g0.B.fill(0.0)
    g0.assert_finite_positive()
    validate_species_admissibility(
        metric_slice_from_q(grid, g0), m0.radiation, Species.RADIATION,
        stage="stage0-preflight",
    )
    scalar_rhs0, matter_rhs0, terms0, stage0_record = _observe_stage(
        kernel, grid, g0, s0, m0, native_grid, "stage0",
        compare_native=False,
    )
    assert np.all(g0.beta == 0.0) and np.all(g0.B == 0.0)

    # Explicit predictor: follow the production partition on copies only.
    g_explicit1 = g0.copy()
    for name in ("a", "b", "X"):
        setattr(
            g_explicit1, name,
            getattr(g0, name) + dt * terms0["explicit"][name],
        )
    g_explicit1.beta.fill(0.0)
    g_explicit1.B.fill(0.0)
    s1 = ScalarFields(*(
        getattr(s0, name) + dt * getattr(scalar_rhs0, name)
        for name in _SCALAR_FIELDS
    ))
    m1_base = tuple(
        kernel._add_matter(getattr(m0, name), matter_rhs0[name], dt)
        for name in ("dark_matter", "baryons", "radiation")
    )
    m1 = V55MatterState(*m1_base, charged_dust=None)
    validate_species_admissibility(
        metric_slice_from_q(grid, g_explicit1), m1.radiation, Species.RADIATION,
        stage="explicit-predictor-preflight",
    )
    active_stage[0] = "explicit_predictor_cmc_refresh"
    g_explicit1.alpha = kernel._solve_lapse_for_state(
        grid, g_explicit1, s1, m1, None
    )[0]
    validate_species_admissibility(
        metric_slice_from_q(grid, g_explicit1), m1.radiation, Species.RADIATION,
        stage="explicit-predictor-preflight",
    )
    _, _, terms_pred, predictor_record = _observe_stage(
        kernel, grid, g_explicit1, s1, m1, native_grid, "explicit_predictor",
        compare_native=True,
    )
    assert np.all(g_explicit1.beta == 0.0) and np.all(g_explicit1.B == 0.0)

    # Primary predictor: update native-candidate primary/Lambda blocks using
    # their established stage partition, then apply the existing boundary and
    # re-solve CMC before the second-stage RHS observation.
    g1 = g_explicit1.copy()
    g1.Aa = g0.Aa + dt * (
        0.5 * terms0["primary_l2"]["Aa"]
        + 0.5 * terms_pred["primary_l2"]["Aa"]
        + terms0["primary_l3"]["Aa"]
    )
    g1.K = g0.K + dt * (
        0.5 * terms0["primary_l2"]["K"]
        + 0.5 * terms_pred["primary_l2"]["K"]
        + terms0["primary_l3"]["K"]
    )
    g1.Lambda = g0.Lambda + dt * (
        0.5 * terms0["lambda_l2"]
        + 0.5 * terms_pred["lambda_l2"]
        + terms0["lambda_l3"]
    )
    validate_species_admissibility(
        metric_slice_from_q(grid, g1), m1.radiation, Species.RADIATION,
        stage="primary-predictor-preflight",
    )
    kernel._apply_outer_boundary_for_state(grid, g1, s1, m1, None)
    active_stage[0] = "primary_predictor_cmc_refresh"
    g1.alpha = kernel._solve_lapse_for_state(grid, g1, s1, m1, None)[0]
    g1.beta.fill(0.0)
    g1.B.fill(0.0)
    g1.assert_finite_positive()
    validate_species_admissibility(
        metric_slice_from_q(grid, g1), m1.radiation, Species.RADIATION,
        stage="primary-predictor-preflight",
    )
    _, _, _, primary_record = _observe_stage(
        kernel, grid, g1, s1, m1, native_grid, "primary_predictor",
        compare_native=True,
    )

    assert [item["stage"] for item in lapse_trace] == [
        "stage0_cmc_refresh",
        "explicit_predictor_cmc_refresh",
        "primary_predictor_cmc_refresh",
    ]
    assert all(item["alpha_finite"] for item in lapse_trace)
    assert all(item["alpha_min"] > 0.0 for item in lapse_trace)
    accepted_hash_after = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )
    assert accepted_hash_after == accepted_hash_before, (
        "CMC stage preflight modified the caller-owned accepted input"
    )

    record = {
        "schema": "v55-cmc-stage-contract-preflight/v1",
        "scope": (
            "copied-state observation through primary predictor only; production "
            "V55ProductionKernel.step is structurally inspected but never called"
        ),
        "resolution": int(grid.n),
        "r_max": float(grid.r_max),
        "dt": dt,
        "native_and_production_grid_centres_match": True,
        "production_step_called": False,
        "production_stage_call_order": contract_events,
        "cmc_lapse_refresh_trace": lapse_trace,
        "observed_stages": [stage0_record, predictor_record, primary_record],
        "accepted_state_sha256_before": accepted_hash_before,
        "accepted_state_sha256_after": accepted_hash_after,
        "accepted_state_unchanged": True,
        "source_ownership": {
            "primary_l3_matter_stress": "reported separately at predictor/primary-predictor",
            "lambda_l3_current": "reported separately at predictor/primary-predictor",
            "geometry_only_rhs": "compared only at identical stage inputs",
        },
        "new_acceptance_thresholds": False,
        "scientific_admission": False,
        "interpretation_limit": (
            "single-resolution stage-contract/software evidence only; not a full-step, "
            "full-trajectory, conservation, or scientific-admission result"
        ),
    }
    with capsys.disabled():
        print(
            "NATIVE_V55_CMC_STAGE_CONTRACT_PREFLIGHT "
            + json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False),
            flush=True,
        )
