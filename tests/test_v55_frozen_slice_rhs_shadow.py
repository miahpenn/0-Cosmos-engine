"""Observation-only V5.5 frozen-slice RHS/source decomposition witness.

This evaluates source-free native RHS blocks and the canonical adapter on the
same copied geometry slice. It is not a full-step comparison: the initialized
production state uses CMC gauge, while the isolated native step evolves the
moving-puncture gauge. No production state is advanced or modified here.
"""

import hashlib
import json

import numpy as np

from engine import v55_matter
from engine.native_geometry import (
    NativeSphericalCellGrid,
    native_bssn_geometry_terms,
    native_explicit_rhs,
    native_lambda_l2_rhs,
    native_lambda_l3_rhs,
    native_primary_l2_rhs,
    native_primary_l3_rhs,
)
from engine.production_kernel import V55ProductionKernel
from engine.v55_pirk_adapter import geometry_stage_terms, vendor_modules


_GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
_SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
_MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
_MATTER_SPECIES = ("dark_matter", "baryons", "radiation")
_LAMBDA_M = 2.0


def _slice_sha256(geometry, scalars, matter):
    digest = hashlib.sha256()
    groups = (
        ("geometry", geometry, _GEOMETRY_FIELDS),
        ("scalars", scalars, _SCALAR_FIELDS),
    )
    for label, obj, names in groups:
        for name in names:
            values = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            digest.update(label.encode("ascii") + b"\0" + name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    for species in _MATTER_SPECIES:
        obj = getattr(matter, species)
        for name in _MATTER_FIELDS:
            values = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            digest.update(b"matter\0" + species.encode("ascii") + b"\0")
            digest.update(name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _regional_profile_summary(values, masks):
    values = np.asarray(values, dtype=float)
    assert values.ndim == 1 and np.all(np.isfinite(values))
    result = {}
    for name, mask in masks.items():
        selected = values[mask]
        assert selected.size >= 4, f"{name}: insufficient cells for a descriptive profile"
        result[name] = {
            "cells": int(selected.size),
            "rms": float(np.sqrt(np.mean(selected**2))),
            "max_abs": float(np.max(np.abs(selected))),
        }
    return result


def _regional_comparison_summary(native, reference, masks):
    residual = np.asarray(native, dtype=float) - np.asarray(reference, dtype=float)
    assert residual.shape == np.asarray(reference).shape
    assert np.all(np.isfinite(residual))
    return _regional_profile_summary(residual, masks)


def test_v55_frozen_slice_reports_native_geometry_rhs_and_matter_source_decomposition(capsys):
    # Keep the established production initializer/gauge, but only evaluate RHSs
    # on copies: this test never calls V55ProductionKernel.step.
    accepted = V55ProductionKernel().initialize(
        resolution=80, r_max=40.0, include_radiation=True
    )
    accepted_hash_before = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )
    geometry = accepted.geometry.copy()
    scalars = accepted.scalars.copy()
    matter = accepted.matter.copy()

    # A tiny, regular radial current on the copied dust profile exercises the
    # canonical Lambda-L3 source without asserting a new matter law.
    r = np.asarray(accepted.grid.centers, dtype=float)
    regular_velocity = 1.0e-4 * r / (1.0 + r * r)
    matter.dark_matter.momentum_r = (
        np.asarray(matter.dark_matter.momentum_r, dtype=float)
        + np.asarray(matter.dark_matter.rest, dtype=float) * regular_velocity
    )
    copy_hash_before = _slice_sha256(geometry, scalars, matter)

    native_grid = NativeSphericalCellGrid(accepted.grid.n, accepted.grid.r_max)
    np.testing.assert_allclose(
        native_grid.centers, accepted.grid.centers, rtol=0.0, atol=1.0e-14
    )
    areal = native_bssn_geometry_terms(native_grid, geometry)["areal_radius"]
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    for name, mask in masks.items():
        assert int(np.count_nonzero(mask)) >= 4, name

    projection = v55_matter.total_matter_projection(
        accepted.grid, geometry, scalars, matter
    )
    assert np.all(np.isfinite(projection["rho"]))
    assert np.all(np.isfinite(projection["pr"]))
    assert np.all(np.isfinite(projection["pt"]))
    assert np.all(np.isfinite(projection["j"]))
    assert np.any(projection["j"] != 0.0), "copied-state fixture must exercise Lambda-L3 current"

    reference = geometry_stage_terms(
        accepted.grid, geometry, scalars, matter, lambda_m=_LAMBDA_M
    )
    native = {
        "explicit": native_explicit_rhs(native_grid, geometry),
        "primary_l2": native_primary_l2_rhs(native_grid, geometry),
        "primary_l3": native_primary_l3_rhs(native_grid, geometry),
        "lambda_l2": native_lambda_l2_rhs(native_grid, geometry, lambda_m=_LAMBDA_M),
        "lambda_l3": native_lambda_l3_rhs(native_grid, geometry),
    }
    _, vacuum, _ = vendor_modules()
    vacuum_primary_l3 = vacuum.primary_l3_rhs(accepted.grid, geometry)
    vacuum_lambda_l3 = np.asarray(
        vacuum.lambda_l3_rhs(accepted.grid, geometry), dtype=float
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
    lambda_source = (
        -8.0 * np.pi
        * _LAMBDA_M
        * geometry.alpha
        * projection["j"]
        / geometry.a
    )

    # Lock the adapter source split on this same initialized/copied slice.
    np.testing.assert_allclose(
        reference["primary_l3"]["Aa"],
        np.asarray(vacuum_primary_l3["Aa"]) + primary_source["Aa"],
        rtol=2.0e-14,
        atol=1.0e-22,
    )
    np.testing.assert_allclose(
        reference["primary_l3"]["K"],
        np.asarray(vacuum_primary_l3["K"]) + primary_source["K"],
        rtol=2.0e-14,
        atol=1.0e-22,
    )
    np.testing.assert_allclose(
        reference["lambda_l3"],
        vacuum_lambda_l3 + lambda_source,
        rtol=2.0e-14,
        atol=1.0e-24,
    )

    # Compare only same-slice geometry-stage RHS profiles. Remove the explicit
    # source increments before comparing L3; do not compare post-step gauge data.
    adjusted_reference = {
        "explicit": reference["explicit"],
        "primary_l2": reference["primary_l2"],
        "primary_l3": {
            "Aa": np.asarray(reference["primary_l3"]["Aa"]) - primary_source["Aa"],
            "K": np.asarray(reference["primary_l3"]["K"]) - primary_source["K"],
        },
        "lambda_l2": reference["lambda_l2"],
        "lambda_l3": np.asarray(reference["lambda_l3"]) - lambda_source,
    }
    compared = {
        "explicit": {
            name: _regional_comparison_summary(
                native["explicit"][name], adjusted_reference["explicit"][name], masks
            )
            for name in ("a", "b", "X", "alpha", "beta")
        },
        "primary_l2": {
            name: _regional_comparison_summary(
                native["primary_l2"][name], adjusted_reference["primary_l2"][name], masks
            )
            for name in ("Aa", "K")
        },
        "primary_l3_geometry_only": {
            name: _regional_comparison_summary(
                native["primary_l3"][name], adjusted_reference["primary_l3"][name], masks
            )
            for name in ("Aa", "K")
        },
        "lambda_l2": _regional_comparison_summary(
            native["lambda_l2"], adjusted_reference["lambda_l2"], masks
        ),
        "lambda_l3_geometry_only": _regional_comparison_summary(
            native["lambda_l3"], adjusted_reference["lambda_l3"], masks
        ),
    }
    source_profiles = {
        "primary_l3_Aa": _regional_profile_summary(primary_source["Aa"], masks),
        "primary_l3_K": _regional_profile_summary(primary_source["K"], masks),
        "lambda_l3_current": _regional_profile_summary(lambda_source, masks),
        "projected_current_j": _regional_profile_summary(projection["j"], masks),
    }

    copy_hash_after = _slice_sha256(geometry, scalars, matter)
    accepted_hash_after = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )
    assert copy_hash_after == copy_hash_before, (
        f"frozen copied input mutated during RHS evaluation: "
        f"before={copy_hash_before}, after={copy_hash_after}"
    )
    assert accepted_hash_after == accepted_hash_before, (
        f"accepted production state mutated: "
        f"before={accepted_hash_before}, after={accepted_hash_after}"
    )

    record = {
        "schema": "v55-frozen-slice-rhs-shadow/v1",
        "resolution": int(accepted.grid.n),
        "r_max": float(accepted.grid.r_max),
        "input_gauge": "CMC-solved alpha; beta and B zero; frozen slice only",
        "comparison_scope": "RHS profiles only; not a CMC-versus-moving-puncture full-step comparison",
        "native_and_reference_grid_centres_match": True,
        "accepted_state_sha256_before": accepted_hash_before,
        "accepted_state_sha256_after": accepted_hash_after,
        "copied_state_sha256_before": copy_hash_before,
        "copied_state_sha256_after": copy_hash_after,
        "source_profiles": source_profiles,
        "native_minus_adapter_geometry_only_rhs": compared,
        "new_acceptance_thresholds": False,
    }
    with capsys.disabled():
        print(
            "NATIVE_V55_FROZEN_SLICE_SHADOW "
            + json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False),
            flush=True,
        )
