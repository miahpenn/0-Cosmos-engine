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
        # Reuse the existing pinned geometry decomposition to attribute the
        # trace-free Aa RHS residual to its discrete Hessian and Ricci pieces.
        pinned_geometry = vacuum.geometry_terms(accepted.grid, geometry)
        alpha = np.asarray(geometry.alpha, dtype=float)
        native_hessian_tf = -(
            native_geometry["DrDr_alpha"] - native_geometry["lap_alpha"] / 3.0
        )
        pinned_hessian_tf = -(
            pinned_geometry["DrDr_alpha"] - pinned_geometry["lap_alpha"] / 3.0
        )
        native_ricci_tf = alpha * (
            native_geometry["Rr_r"] - native_geometry["R"] / 3.0
        )
        pinned_ricci_tf = alpha * (
            pinned_geometry["Rr_r"] - pinned_geometry["R"] / 3.0
        )
        hessian_delta = native_hessian_tf - pinned_hessian_tf
        ricci_delta = native_ricci_tf - pinned_ricci_tf
        areal = np.asarray(native_geometry["areal_radius"], dtype=float)
        region_masks = {
            "centre_R_areal_lt_2": areal < 2.0,
            "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
            "outer_R_areal_gt_6": areal > 6.0,
        }
        mask = region_masks["fixed_annulus_2_le_R_areal_le_6"]
        assert int(np.count_nonzero(mask)) >= 8, f"N={n}: insufficient annulus cells"

        # Form full signed residual profiles first; region reductions below are
        # descriptive for centre and outer cells, not additional pass/fail limits.
        aa_delta = np.asarray(native["Aa"], dtype=float) - np.asarray(pinned["Aa"], dtype=float)
        k_delta = np.asarray(native["K"], dtype=float) - np.asarray(pinned["K"], dtype=float)
        assert np.all(np.isfinite(aa_delta)), f"N={n}: non-finite Aa difference"
        assert np.all(np.isfinite(k_delta)), f"N={n}: non-finite K difference"
        assert np.all(np.isfinite(hessian_delta)), f"N={n}: non-finite Hessian-piece difference"
        assert np.all(np.isfinite(ricci_delta)), f"N={n}: non-finite Ricci-piece difference"
        # This is an algebraic bookkeeping identity for the published RHS form,
        # not a convergence or scientific-admission threshold.
        decomposition_remainder = aa_delta - (hessian_delta + ricci_delta)
        assert np.all(np.isfinite(decomposition_remainder)), f"N={n}: invalid decomposition remainder"
        region_summaries = {}
        for region, region_mask in region_masks.items():
            cells = int(np.count_nonzero(region_mask))
            if cells == 0:
                raise AssertionError(f"N={n}: empty diagnostic region {region}")
            aa_region = aa_delta[region_mask]
            k_region = k_delta[region_mask]
            hessian_region = hessian_delta[region_mask]
            ricci_region = ricci_delta[region_mask]
            remainder_region = decomposition_remainder[region_mask]
            region_summaries[region] = {
                "cells": cells,
                "r_areal_min": float(np.min(areal[region_mask])),
                "r_areal_max": float(np.max(areal[region_mask])),
                "aa_rms_native_minus_pinned": float(np.sqrt(np.mean(aa_region**2))),
                "aa_max_abs_native_minus_pinned": float(np.max(np.abs(aa_region))),
                "aa_hessian_tracefree_rms": float(np.sqrt(np.mean(hessian_region**2))),
                "aa_hessian_tracefree_max_abs": float(np.max(np.abs(hessian_region))),
                "aa_ricci_tracefree_rms": float(np.sqrt(np.mean(ricci_region**2))),
                "aa_ricci_tracefree_max_abs": float(np.max(np.abs(ricci_region))),
                "aa_decomposition_remainder_max_abs": float(np.max(np.abs(remainder_region))),
                "k_rms_native_minus_pinned": float(np.sqrt(np.mean(k_region**2))),
                "k_max_abs_native_minus_pinned": float(np.max(np.abs(k_region))),
            }
        aa_rms = region_summaries["fixed_annulus_2_le_R_areal_le_6"]["aa_rms_native_minus_pinned"]
        k_rms = region_summaries["fixed_annulus_2_le_R_areal_le_6"]["k_rms_native_minus_pinned"]
        aa_errors.append(aa_rms)
        k_errors.append(k_rms)
        records.append({
            "resolution": n,
            "r_max": 40.0,
            "annulus_cells": int(np.count_nonzero(mask)),
            "aa_rms_native_minus_pinned": aa_rms,
            "aa_max_abs_native_minus_pinned": float(np.max(np.abs(aa_delta[mask]))),
            "k_rms_native_minus_pinned": k_rms,
            "regional_residuals": region_summaries,
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
                    "schema": "v55-primary-l2-cmc-refinement/v3",
                    "scope": "regional residuals plus trace-free Aa term decomposition; diagnostic RHS only",
                    "records": records,
                    "aa_monotonic_refinement": True,
                    "k_roundoff_contract": True,
                    "new_absolute_error_threshold": False,
                    "termwise_decomposition_is_diagnostic_only": True,
                },
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ),
            flush=True,
        )
