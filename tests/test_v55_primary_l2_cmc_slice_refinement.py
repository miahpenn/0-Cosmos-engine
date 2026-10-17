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
        X = np.asarray(geometry.X, dtype=float)
        # The pinned BSSN Ricci has explicit evolved-Lambda pieces in addition
        # to its metric terms. Split those pieces using the pinned RHS's own
        # derivative conventions, then compare them with the existing native
        # metric Ricci and reference-connection correction separately.
        pinned_lambda_p = vacuum._d(accepted.grid, geometry.Lambda, -1)
        pinned_a_p = vacuum._d(accepted.grid, geometry.a, 1)
        pinned_b_p = vacuum._d(accepted.grid, geometry.b, 1)
        # Expose the evolved-versus-metric connection split under each
        # implementation's own derivative closure. These are diagnostic
        # primitives only; neither derivative convention is silently swapped.
        native_lambda_metric = np.asarray(native_geometry["Lambda_metric"], dtype=float)
        native_delta_lambda = np.asarray(
            native_geometry["connection_constraint"], dtype=float
        )
        native_lambda_p = native_grid.cell_derivative_fourth(
            np.asarray(geometry.Lambda, dtype=float), parity=-1
        )
        native_lambda_metric_p = native_grid.cell_derivative_fourth(
            native_lambda_metric, parity=-1
        )
        native_delta_lambda_p = native_grid.cell_derivative_fourth(
            native_delta_lambda, parity=-1
        )
        native_a_p = native_grid.cell_derivative_fourth(
            np.asarray(geometry.a, dtype=float), parity=1
        )
        native_b_p = native_grid.cell_derivative_fourth(
            np.asarray(geometry.b, dtype=float), parity=1
        )
        r_centres = np.asarray(native_grid.centers, dtype=float)
        a_profile = np.asarray(geometry.a, dtype=float)
        b_profile = np.asarray(geometry.b, dtype=float)
        X_profile = np.asarray(geometry.X, dtype=float)
        Lambda_profile = np.asarray(geometry.Lambda, dtype=float)
        pinned_lambda_metric = (
            pinned_a_p / (2.0 * a_profile**2)
            - pinned_b_p / (a_profile * b_profile)
            + 2.0 * (1.0 / b_profile - 1.0 / a_profile) / r_centres
        )
        pinned_delta_lambda = Lambda_profile - pinned_lambda_metric
        pinned_lambda_metric_p = vacuum._d(
            accepted.grid, pinned_lambda_metric, -1
        )
        pinned_delta_lambda_p = vacuum._d(
            accepted.grid, pinned_delta_lambda, -1
        )
        pinned_metric_lambda_aux = (
            1.0 - a_profile / b_profile
        ) / r_centres**2
        pinned_connection_Rr_r = X**2 * (
            pinned_lambda_p + 0.5 * (pinned_a_p / geometry.a) * geometry.Lambda
        )
        pinned_connection_R = X**2 * pinned_lambda_p
        pinned_metric_Rr_r = pinned_geometry["Rr_r"] - pinned_connection_Rr_r
        pinned_metric_R = pinned_geometry["R"] - pinned_connection_R
        native_metric_ricci_tf = alpha * (
            native_geometry["Rr_r_metric"] - native_geometry["R_metric"] / 3.0
        )
        pinned_metric_ricci_tf = alpha * (
            pinned_metric_Rr_r - pinned_metric_R / 3.0
        )
        native_connection_ricci_tf = alpha * (
            native_geometry["delta_Rr_r_connection"]
            - native_geometry["delta_R_scalar_connection"] / 3.0
        )
        pinned_connection_ricci_tf = alpha * (
            pinned_connection_Rr_r - pinned_connection_R / 3.0
        )
        metric_ricci_delta = native_metric_ricci_tf - pinned_metric_ricci_tf
        connection_ricci_delta = native_connection_ricci_tf - pinned_connection_ricci_tf
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
        assert np.all(np.isfinite(metric_ricci_delta)), f"N={n}: non-finite metric-Ricci difference"
        assert np.all(np.isfinite(connection_ricci_delta)), f"N={n}: non-finite connection-Ricci difference"
        # These are algebraic bookkeeping identities of existing RHS terms,
        # not convergence or scientific-admission thresholds.
        decomposition_remainder = aa_delta - (hessian_delta + ricci_delta)
        ricci_split_remainder = ricci_delta - (
            metric_ricci_delta + connection_ricci_delta
        )
        assert np.all(np.isfinite(decomposition_remainder)), f"N={n}: invalid Aa decomposition remainder"
        assert np.all(np.isfinite(ricci_split_remainder)), f"N={n}: invalid Ricci split remainder"
        region_summaries = {}
        for region, region_mask in region_masks.items():
            cells = int(np.count_nonzero(region_mask))
            if cells == 0:
                raise AssertionError(f"N={n}: empty diagnostic region {region}")
            aa_region = aa_delta[region_mask]
            k_region = k_delta[region_mask]
            hessian_region = hessian_delta[region_mask]
            ricci_region = ricci_delta[region_mask]
            metric_ricci_region = metric_ricci_delta[region_mask]
            connection_ricci_region = connection_ricci_delta[region_mask]
            remainder_region = decomposition_remainder[region_mask]
            ricci_remainder_region = ricci_split_remainder[region_mask]
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
                "aa_ricci_metric_tracefree_rms": float(np.sqrt(np.mean(metric_ricci_region**2))),
                "aa_ricci_metric_tracefree_max_abs": float(np.max(np.abs(metric_ricci_region))),
                "aa_ricci_connection_tracefree_rms": float(np.sqrt(np.mean(connection_ricci_region**2))),
                "aa_ricci_connection_tracefree_max_abs": float(np.max(np.abs(connection_ricci_region))),
                "aa_ricci_split_remainder_max_abs": float(np.max(np.abs(ricci_remainder_region))),
                "aa_decomposition_remainder_max_abs": float(np.max(np.abs(remainder_region))),
                "k_rms_native_minus_pinned": float(np.sqrt(np.mean(k_region**2))),
                "k_max_abs_native_minus_pinned": float(np.max(np.abs(k_region))),
            }
        aa_rms = region_summaries["fixed_annulus_2_le_R_areal_le_6"]["aa_rms_native_minus_pinned"]
        k_rms = region_summaries["fixed_annulus_2_le_R_areal_le_6"]["k_rms_native_minus_pinned"]
        aa_errors.append(aa_rms)
        k_errors.append(k_rms)
        # Keep signed samples for the first regular-origin cells in addition
        # to regional reductions. These are descriptive only and deliberately
        # have no per-cell acceptance threshold.
        origin_indices = np.flatnonzero(areal < 2.0)[:8]
        origin_cell_samples = []
        for index in origin_indices:
            origin_cell_samples.append({
                "index": int(index),
                "coordinate_r": float(native_grid.centers[index]),
                "areal_radius": float(areal[index]),
                "aa_total_native_minus_pinned": float(aa_delta[index]),
                "aa_hessian_tracefree_native_minus_pinned": float(hessian_delta[index]),
                "aa_ricci_tracefree_native_minus_pinned": float(ricci_delta[index]),
                "aa_ricci_metric_piece_native_minus_pinned": float(metric_ricci_delta[index]),
                "aa_ricci_connection_piece_native_minus_pinned": float(connection_ricci_delta[index]),
                "aa_hessian_plus_ricci_remainder": float(decomposition_remainder[index]),
                "ricci_metric_plus_connection_remainder": float(ricci_split_remainder[index]),
                "native": {
                    "Rr_r_metric": float(native_geometry["Rr_r_metric"][index]),
                    "R_metric": float(native_geometry["R_metric"][index]),
                    "delta_Rr_r_connection": float(native_geometry["delta_Rr_r_connection"][index]),
                    "delta_R_scalar_connection": float(native_geometry["delta_R_scalar_connection"][index]),
                    "tracefree_metric_piece": float(native_metric_ricci_tf[index]),
                    "tracefree_connection_piece": float(native_connection_ricci_tf[index]),
                },
                "pinned": {
                    "Rr_r_total": float(pinned_geometry["Rr_r"][index]),
                    "R_scalar_total": float(pinned_geometry["R"][index]),
                    "Rr_r_metric_remainder": float(pinned_metric_Rr_r[index]),
                    "R_scalar_metric_remainder": float(pinned_metric_R[index]),
                    "Rr_r_explicit_Lambda_piece": float(pinned_connection_Rr_r[index]),
                    "R_scalar_explicit_Lambda_piece": float(pinned_connection_R[index]),
                    "tracefree_metric_piece": float(pinned_metric_ricci_tf[index]),
                    "tracefree_Lambda_piece": float(pinned_connection_ricci_tf[index]),
                },
            })
        # Emit the first four regular-origin cells' connection/Ricci
        # primitives side-by-side. This exposes derivative-convention and
        # metric-vs-evolved-Lambda differences without introducing tolerances.
        origin_primitive_samples = []
        for index in origin_indices[:4]:
            primitive_arrays = {
                "native_Lambda_metric": native_lambda_metric,
                "native_deltaLambda": native_delta_lambda,
                "native_Lambda_prime": native_lambda_p,
                "native_Lambda_metric_prime": native_lambda_metric_p,
                "native_deltaLambda_prime": native_delta_lambda_p,
                "native_a_prime": native_a_p,
                "native_b_prime": native_b_p,
                "pinned_Lambda_metric": pinned_lambda_metric,
                "pinned_deltaLambda": pinned_delta_lambda,
                "pinned_Lambda_prime": pinned_lambda_p,
                "pinned_Lambda_metric_prime": pinned_lambda_metric_p,
                "pinned_deltaLambda_prime": pinned_delta_lambda_p,
                "pinned_a_prime": pinned_a_p,
                "pinned_b_prime": pinned_b_p,
                "pinned_metric_lambda_aux": pinned_metric_lambda_aux,
                "native_Rr_r_metric": native_geometry["Rr_r_metric"],
                "native_R_scalar_metric": native_geometry["R_metric"],
                "native_delta_Rr_r_connection": native_geometry["delta_Rr_r_connection"],
                "native_delta_R_scalar_connection": native_geometry["delta_R_scalar_connection"],
                "native_Rr_r_total": native_geometry["Rr_r"],
                "native_R_scalar_total": native_geometry["R"],
                "pinned_Rr_r_total": pinned_geometry["Rr_r"],
                "pinned_R_scalar_total": pinned_geometry["R"],
                "pinned_Rr_r_metric_remainder": pinned_metric_Rr_r,
                "pinned_R_scalar_metric_remainder": pinned_metric_R,
                "pinned_Rr_r_explicit_Lambda_piece": pinned_connection_Rr_r,
                "pinned_R_scalar_explicit_Lambda_piece": pinned_connection_R,
            }
            for name, values in primitive_arrays.items():
                assert np.all(np.isfinite(np.asarray(values))), (
                    f"N={n}, cell={index}: non-finite primitive {name}"
                )
            origin_primitive_samples.append({
                "index": int(index),
                "coordinate_r": float(r_centres[index]),
                "areal_radius": float(areal[index]),
                "evolved_Lambda": float(Lambda_profile[index]),
                "native": {
                    "Lambda_metric": float(native_lambda_metric[index]),
                    "deltaLambda": float(native_delta_lambda[index]),
                    "Lambda_prime": float(native_lambda_p[index]),
                    "Lambda_metric_prime": float(native_lambda_metric_p[index]),
                    "deltaLambda_prime": float(native_delta_lambda_p[index]),
                    "a_prime": float(native_a_p[index]),
                    "b_prime": float(native_b_p[index]),
                    "Rr_r_metric": float(native_geometry["Rr_r_metric"][index]),
                    "R_scalar_metric": float(native_geometry["R_metric"][index]),
                    "delta_Rr_r_connection": float(native_geometry["delta_Rr_r_connection"][index]),
                    "delta_R_scalar_connection": float(native_geometry["delta_R_scalar_connection"][index]),
                    "Rr_r_total": float(native_geometry["Rr_r"][index]),
                    "R_scalar_total": float(native_geometry["R"][index]),
                },
                "pinned": {
                    "Lambda_metric": float(pinned_lambda_metric[index]),
                    "deltaLambda": float(pinned_delta_lambda[index]),
                    "Lambda_prime": float(pinned_lambda_p[index]),
                    "Lambda_metric_prime": float(pinned_lambda_metric_p[index]),
                    "deltaLambda_prime": float(pinned_delta_lambda_p[index]),
                    "a_prime": float(pinned_a_p[index]),
                    "b_prime": float(pinned_b_p[index]),
                    "metric_lambda_aux": float(pinned_metric_lambda_aux[index]),
                    "Rr_r_total": float(pinned_geometry["Rr_r"][index]),
                    "R_scalar_total": float(pinned_geometry["R"][index]),
                    "Rr_r_metric_remainder": float(pinned_metric_Rr_r[index]),
                    "R_scalar_metric_remainder": float(pinned_metric_R[index]),
                    "Rr_r_explicit_Lambda_piece": float(pinned_connection_Rr_r[index]),
                    "R_scalar_explicit_Lambda_piece": float(pinned_connection_R[index]),
                },
                "native_minus_pinned": {
                    "Lambda_metric": float(native_lambda_metric[index] - pinned_lambda_metric[index]),
                    "deltaLambda": float(native_delta_lambda[index] - pinned_delta_lambda[index]),
                    "deltaLambda_prime": float(native_delta_lambda_p[index] - pinned_delta_lambda_p[index]),
                    "Rr_r_total": float(native_geometry["Rr_r"][index] - pinned_geometry["Rr_r"][index]),
                    "R_scalar_total": float(native_geometry["R"][index] - pinned_geometry["R"][index]),
                },
            })
        records.append({
            "resolution": n,
            "r_max": 40.0,
            "annulus_cells": int(np.count_nonzero(mask)),
            "aa_rms_native_minus_pinned": aa_rms,
            "aa_max_abs_native_minus_pinned": float(np.max(np.abs(aa_delta[mask]))),
            "k_rms_native_minus_pinned": k_rms,
            "regional_residuals": region_summaries,
            "origin_cell_signed_samples": origin_cell_samples,
            "origin_primitive_samples_first_four": origin_primitive_samples,
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
                    "schema": "v55-primary-l2-cmc-refinement/v6",
                    "scope": "regional residuals, first-eight signed term samples, and first-four Lambda/Ricci primitives; diagnostic RHS only",
                    "records": records,
                    "aa_monotonic_refinement": True,
                    "k_roundoff_contract": True,
                    "new_absolute_error_threshold": False,
                    "termwise_decomposition_is_diagnostic_only": True,
                    "ricci_metric_connection_split_is_diagnostic_only": True,
                    "origin_cell_samples_are_diagnostic_only": True,
                },
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ),
            flush=True,
        )
