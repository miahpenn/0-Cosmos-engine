"""Refinement witness for native/pinned primary-L2 on initialized CMC slices.

This checks convergence of the independent Aa_L2 discrete residual on the
actual production initializer's smooth geometry profiles. It does not promote
the pinned discretization to formula authority or change production routing.
"""

import json

import numpy as np

from engine.native_geometry import (
    NativeSphericalCellGrid,
    native_bssn_geometry_terms,
    native_primary_l2_rhs,
)
from engine.production_kernel import V55ProductionKernel
from engine.v55_pirk_adapter import vendor_modules


def test_primary_l2_aa_residual_refines_on_initialized_cmc_slices(capsys):
    _, vacuum, _ = vendor_modules()
    records = []
    aa_errors = []
    k_errors = []

    for n in (80, 160, 320):
        accepted = V55ProductionKernel().initialize(
            resolution=n, r_max=40.0, include_radiation=True
        )
        geometry = accepted.geometry.copy()
        native_grid = NativeSphericalCellGrid(n, 40.0)
        np.testing.assert_allclose(
            native_grid.centers, accepted.grid.centers, rtol=0.0, atol=1.0e-14
        )

        fields = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
        before = {name: np.asarray(getattr(geometry, name)).copy() for name in fields}
        native_geometry = native_bssn_geometry_terms(native_grid, geometry)
        native = native_primary_l2_rhs(native_grid, geometry)
        pinned = vacuum.primary_l2_rhs(accepted.grid, geometry)
        areal = np.asarray(native_geometry["areal_radius"], dtype=float)
        mask = (areal >= 2.0) & (areal <= 6.0)
        assert int(np.count_nonzero(mask)) >= 8, f"N={n}: insufficient annulus cells"

        aa_delta = np.asarray(native["Aa"]) - np.asarray(pinned["Aa"])
        k_delta = np.asarray(native["K"]) - np.asarray(pinned["K"])
        assert np.all(np.isfinite(aa_delta)), f"N={n}: non-finite Aa difference"
        assert np.all(np.isfinite(k_delta)), f"N={n}: non-finite K difference"
        aa_rms = float(np.sqrt(np.mean(aa_delta[mask] ** 2)))
        k_rms = float(np.sqrt(np.mean(k_delta[mask] ** 2)))
        aa_errors.append(aa_rms)
        k_errors.append(k_rms)
        records.append({
            "resolution": n,
            "r_max": 40.0,
            "annulus_cells": int(np.count_nonzero(mask)),
            "aa_rms_native_minus_pinned": aa_rms,
            "aa_max_abs_native_minus_pinned": float(np.max(np.abs(aa_delta[mask]))),
            "k_rms_native_minus_pinned": k_rms,
        })
        for name, values in before.items():
            np.testing.assert_array_equal(getattr(geometry, name), values)

    # The K lapse-Laplacian identity agrees to roundoff on this pinned grid.
    assert max(k_errors) < 1.0e-10, f"K_L2 annulus RMS: {k_errors}"

    # Aa uses the independently formed Ricci/Hessian trace-free combination.
    # Require refinement on these initialized profiles, not bitwise identity
    # or any fitted coefficient/absolute vendor-agreement threshold.
    assert aa_errors[1] < aa_errors[0], f"Aa_L2 annulus RMS did not refine: {aa_errors}"
    assert aa_errors[2] < aa_errors[1], f"Aa_L2 annulus RMS did not refine: {aa_errors}"

    with capsys.disabled():
        print(
            "NATIVE_V55_PRIMARY_L2_CMC_REFINEMENT "
            + json.dumps(
                {
                    "schema": "v55-primary-l2-cmc-refinement/v1",
                    "scope": "fixed-physical-areal-annulus; diagnostic RHS only",
                    "records": records,
                    "aa_monotonic_refinement": True,
                    "k_roundoff_contract": True,
                    "new_absolute_error_threshold": False,
                },
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ),
            flush=True,
        )
