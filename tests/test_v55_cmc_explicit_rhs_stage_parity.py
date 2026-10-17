"""Stage-local parity witness for the native explicit RHS on production CMC slices.

The current production step is run once on a deep-copied initialized state solely
to expose its real CMC-owned stage inputs. The original adapter terms are
returned unchanged. Native output is evaluated read-only on an identical copy;
only a, b, and X are compared because the production CMC path solves alpha and
pins beta/B rather than evolving the native moving-puncture gauge RHS.
"""

from __future__ import annotations

import copy
import hashlib
import json

import numpy as np

from engine.native_geometry import NativeSphericalCellGrid, native_explicit_rhs
from engine.production_kernel import V55ProductionKernel


_GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
_SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
_MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
_MATTER_SPECIES = ("dark_matter", "baryons", "radiation", "charged_dust")
_COMPARED_FIELDS = ("a", "b", "X")
_ANNULUS_RMS_LIMIT = 1.0e-12  # Reuse the existing native explicit-RHS gate from #323.


def _slice_sha256(geometry, scalars, matter):
    """Hash numerical inputs with field names, dtype, and shape provenance."""
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


def _profile_hash(geometry):
    digest = hashlib.sha256()
    for name in _GEOMETRY_FIELDS:
        values = np.ascontiguousarray(np.asarray(getattr(geometry, name)))
        digest.update(name.encode("ascii") + b"\0")
        digest.update(values.dtype.str.encode("ascii") + b"\0")
        digest.update(repr(values.shape).encode("ascii") + b"\0")
        digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _regional_summary(residual, areal_radius):
    residual = np.asarray(residual, dtype=float)
    areal_radius = np.asarray(areal_radius, dtype=float)
    assert residual.ndim == 1 and residual.shape == areal_radius.shape
    assert np.all(np.isfinite(residual))
    masks = {
        "centre_R_areal_lt_2": areal_radius < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (
            (areal_radius >= 2.0) & (areal_radius <= 6.0)
        ),
        "outer_R_areal_gt_6": areal_radius > 6.0,
    }
    regions = {}
    for name, mask in masks.items():
        selected = residual[mask]
        assert selected.size >= 4, f"{name}: insufficient cells for observation"
        regions[name] = {
            "cells": int(selected.size),
            "rms": float(np.sqrt(np.mean(selected**2))),
            "max_abs": float(np.max(np.abs(selected))),
        }
    return regions


def test_native_explicit_rhs_matches_actual_cmc_stage_inputs(capsys):
    """Compare only a/b/X rates on real stage inputs; never substitute outputs."""
    kernel = V55ProductionKernel()
    accepted = kernel.initialize(
        resolution=64, r_max=16.0, include_radiation=True
    )
    accepted_sha_before = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )

    # Isolate the one bounded production step. Captured stage inputs are the
    # actual inputs reaching the adapter, not a separately reconstructed step.
    step_input = copy.deepcopy(accepted)
    step_input_sha_before = _slice_sha256(
        step_input.geometry, step_input.scalars, step_input.matter
    )
    assert step_input_sha_before == accepted_sha_before

    native_grid = NativeSphericalCellGrid(step_input.grid.n, step_input.grid.r_max)
    np.testing.assert_allclose(
        native_grid.centers, step_input.grid.centers, rtol=0.0, atol=1.0e-14
    )

    stage_names = ("stage0_cmc_slice", "explicit_predictor_cmc_slice",
                   "primary_predictor_cmc_slice")
    stage_records = []
    original_stage_terms = kernel._geometry_stage_terms

    def observe_real_stage(grid, geometry, scalars, matter, electromagnetic):
        index = len(stage_records)
        assert index < len(stage_names), (
            "production CMC stage call count changed; update the stage contract first"
        )
        stage_name = stage_names[index]
        live_sha_before = _slice_sha256(geometry, scalars, matter)
        snapshot = geometry.copy()
        snapshot_sha_before = _profile_hash(snapshot)

        # Obtain the authoritative adapter RHS on the real production stage.
        # Return it unchanged after read-only native observation.
        reference_terms = original_stage_terms(
            grid, geometry, scalars, matter, electromagnetic
        )
        live_sha_after_adapter = _slice_sha256(geometry, scalars, matter)
        assert live_sha_after_adapter == live_sha_before, (
            f"{stage_name}: adapter mutated its stage input"
        )

        snapshot_sha_after_adapter = _profile_hash(snapshot)
        assert snapshot_sha_after_adapter == snapshot_sha_before

        native_terms = native_explicit_rhs(native_grid, snapshot)
        snapshot_sha_after_native = _profile_hash(snapshot)
        live_sha_after_native = _slice_sha256(geometry, scalars, matter)
        assert snapshot_sha_after_native == snapshot_sha_before, (
            f"{stage_name}: native explicit RHS mutated its copied geometry"
        )
        assert live_sha_after_native == live_sha_before, (
            f"{stage_name}: native evaluation changed the live production stage"
        )

        areal_radius = (
            np.asarray(native_grid.centers, dtype=float)
            * np.sqrt(np.asarray(snapshot.b, dtype=float))
            / np.asarray(snapshot.X, dtype=float)
        )
        field_records = {}
        for field in _COMPARED_FIELDS:
            native_values = np.asarray(native_terms[field], dtype=float)
            reference_values = np.asarray(reference_terms["explicit"][field], dtype=float)
            assert native_values.shape == reference_values.shape == (grid.n,)
            assert np.all(np.isfinite(native_values)), f"{stage_name}/{field}: native non-finite"
            assert np.all(np.isfinite(reference_values)), f"{stage_name}/{field}: adapter non-finite"
            residual = native_values - reference_values
            regional = _regional_summary(residual, areal_radius)
            annulus_rms = regional["fixed_annulus_2_le_R_areal_le_6"]["rms"]
            assert annulus_rms < _ANNULUS_RMS_LIMIT, (
                f"{stage_name}/{field}: native-vs-adapter annulus RMS "
                f"{annulus_rms:.6e} exceeds the existing #323 gate "
                f"{_ANNULUS_RMS_LIMIT:.1e}"
            )
            field_records[field] = {
                "native_minus_adapter": regional,
                "existing_annulus_rms_limit": _ANNULUS_RMS_LIMIT,
                "annulus_gate_passed": True,
            }

        stage_records.append({
            "stage": stage_name,
            "input_gauge": "production CMC stage; alpha refreshed as required; beta/B pinned to zero",
            "stage_input_sha256_before_adapter": live_sha_before,
            "stage_input_sha256_after_adapter": live_sha_after_adapter,
            "stage_input_sha256_after_native": live_sha_after_native,
            "geometry_snapshot_sha256_before": snapshot_sha_before,
            "geometry_snapshot_sha256_after": snapshot_sha_after_native,
            "native_and_adapter_used_identical_geometry": True,
            "compared_fields": list(_COMPARED_FIELDS),
            "native_moving_puncture_alpha_beta_rhs_compared": False,
            "field_residuals": field_records,
        })
        return reference_terms

    # The observer is instrumentation only. It does not replace a production RHS
    # or mutate the kernel's stage state; the original adapter terms are returned.
    kernel._geometry_stage_terms = observe_real_stage
    result = kernel.step(step_input, 1.0e-5)

    assert len(stage_records) == 3, (
        "expected Stage-0, explicit-predictor, and primary-predictor CMC RHS inputs"
    )
    assert result is not None
    step_input_sha_after = _slice_sha256(
        step_input.geometry, step_input.scalars, step_input.matter
    )
    accepted_sha_after = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )
    assert step_input_sha_after == step_input_sha_before, (
        "production step mutated numerical arrays on its isolated input copy"
    )
    assert accepted_sha_after == accepted_sha_before, (
        "instrumented copied step mutated the caller-owned accepted state"
    )

    report = {
        "schema": "v55-cmc-explicit-rhs-stage-parity/v1",
        "resolution": int(step_input.grid.n),
        "r_max": float(step_input.grid.r_max),
        "dt": 1.0e-5,
        "production_step_used_for_stage_capture": True,
        "step_input_is_deep_copy": True,
        "production_adapter_terms_substituted": False,
        "native_runtime_wired": False,
        "native_grid_centres_match": True,
        "accepted_state_sha256_before": accepted_sha_before,
        "accepted_state_sha256_after": accepted_sha_after,
        "step_copy_sha256_before": step_input_sha_before,
        "step_copy_sha256_after": step_input_sha_after,
        "accepted_and_copy_numerical_arrays_unchanged": True,
        "captured_stage_count": len(stage_records),
        "stage_records": stage_records,
        "existing_annulus_rms_gate": _ANNULUS_RMS_LIMIT,
        "centre_and_outer_regions": "reported for diagnosis only; no added threshold",
        "scientific_admission": False,
        "interpretation_limit": (
            "one-step, stage-local explicit-RHS parity only; not whole-step native "
            "equivalence, trajectory stability/conservation, production wiring, or admission"
        ),
    }
    print(
        "NATIVE_V55_CMC_EXPLICIT_RHS_STAGE_PARITY "
        + json.dumps(report, sort_keys=True, separators=(",", ":"), allow_nan=False),
        flush=True,
    )
