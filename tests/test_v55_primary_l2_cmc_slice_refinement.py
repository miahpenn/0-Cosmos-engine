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

        # Expose raw metric-Ricci primitives under each implementation's own
        # derivative closure. Native geometry differentiates areal radius in
        # proper-radial form; the pinned BSSN path differentiates a, b and X
        # in coordinate form. This block is diagnostic only.
        native_A_profile = a_profile / X_profile**2
        native_areal_profile = r_centres * np.sqrt(b_profile) / X_profile
        native_A_p = native_grid.cell_derivative_fourth(native_A_profile, parity=1)
        native_A_pp = native_grid.cell_second_derivative_fourth(native_A_profile, parity=1)
        native_R_p = native_grid.cell_derivative_fourth(native_areal_profile, parity=-1)
        native_R_pp = native_grid.cell_second_derivative_fourth(native_areal_profile, parity=-1)
        native_a_pp = native_grid.cell_second_derivative_fourth(a_profile, parity=1)
        native_b_pp = native_grid.cell_second_derivative_fourth(b_profile, parity=1)
        native_X_p = native_grid.cell_derivative_fourth(X_profile, parity=1)
        native_X_pp = native_grid.cell_second_derivative_fourth(X_profile, parity=1)
        native_R_proper_p = native_R_p / np.sqrt(native_A_profile)
        native_R_proper_pp = (
            native_R_pp - 0.5 * (native_A_p / native_A_profile) * native_R_p
        ) / native_A_profile

        pinned_A_profile = a_profile / X_profile**2
        pinned_areal_profile = r_centres * np.sqrt(b_profile) / X_profile
        pinned_A_p = vacuum._d(accepted.grid, pinned_A_profile, 1)
        pinned_A_pp = vacuum._d2(accepted.grid, pinned_A_profile, 1)
        pinned_R_p = vacuum._d(accepted.grid, pinned_areal_profile, -1)
        pinned_R_pp = vacuum._d2(accepted.grid, pinned_areal_profile, -1)
        pinned_a_pp = vacuum._d2(accepted.grid, a_profile, 1)
        pinned_b_pp = vacuum._d2(accepted.grid, b_profile, 1)
        pinned_X_p = vacuum._d(accepted.grid, X_profile, 1)
        pinned_X_pp = vacuum._d2(accepted.grid, X_profile, 1)
        pinned_R_proper_p = pinned_R_p / np.sqrt(pinned_A_profile)
        pinned_R_proper_pp = (
            pinned_R_pp - 0.5 * (pinned_A_p / pinned_A_profile) * pinned_R_p
        ) / pinned_A_profile

        # Literal expansion of the pinned radial-Ricci bracket, to localize
        # the origin discrepancy without replacing either implementation.
        pinned_chi_p = -0.5 * pinned_X_p / X_profile
        pinned_chi_pp = -0.5 * (
            pinned_X_pp / X_profile - (pinned_X_p / X_profile) ** 2
        )
        pinned_radial_bracket_terms = {
            "half_a_second_over_a": 0.5 * pinned_a_pp / a_profile,
            "minus_a_lambda_prime": -a_profile * pinned_lambda_p,
            "minus_three_quarters_a_logprime_sq": -0.75 * (pinned_a_p / a_profile) ** 2,
            "half_b_logprime_sq": 0.5 * (pinned_b_p / b_profile) ** 2,
            "minus_half_lambda_a_prime": -0.5 * Lambda_profile * pinned_a_p,
            "a_prime_over_r_b": pinned_a_p / (r_centres * b_profile),
            "metric_lambda_aux": 2.0 * pinned_metric_lambda_aux * (
                1.0 + r_centres * pinned_b_p / b_profile
            ),
            "four_chi_second": 4.0 * pinned_chi_pp,
            "chi_prime_coupling": -2.0 * pinned_chi_p * (
                pinned_a_p / a_profile - pinned_b_p / b_profile - 2.0 / r_centres
            ),
        }
        pinned_radial_bracket_total = np.sum(
            np.stack(tuple(pinned_radial_bracket_terms.values())), axis=0
        )
        pinned_radial_bracket_metric = np.sum(
            np.stack(tuple(
                values for key, values in pinned_radial_bracket_terms.items()
                if key not in {"minus_a_lambda_prime", "minus_half_lambda_a_prime"}
            )), axis=0
        )
        pinned_Rr_r_from_bracket = -X_profile**2 / a_profile * pinned_radial_bracket_total
        pinned_Rr_r_metric_from_bracket = -X_profile**2 / a_profile * pinned_radial_bracket_metric
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
                "native_A_profile": native_A_profile,
                "native_A_prime": native_A_p,
                "native_A_second": native_A_pp,
                "native_areal_radius": native_areal_profile,
                "native_areal_radius_prime": native_R_p,
                "native_areal_radius_second": native_R_pp,
                "native_a_second": native_a_pp,
                "native_b_second": native_b_pp,
                "native_X_prime": native_X_p,
                "native_X_second": native_X_pp,
                "native_R_proper_prime": native_R_proper_p,
                "native_R_proper_second": native_R_proper_pp,
                "pinned_A_profile": pinned_A_profile,
                "pinned_A_prime": pinned_A_p,
                "pinned_A_second": pinned_A_pp,
                "pinned_areal_radius": pinned_areal_profile,
                "pinned_areal_radius_prime": pinned_R_p,
                "pinned_areal_radius_second": pinned_R_pp,
                "pinned_a_second": pinned_a_pp,
                "pinned_b_second": pinned_b_pp,
                "pinned_X_prime": pinned_X_p,
                "pinned_X_second": pinned_X_pp,
                "pinned_chi_prime": pinned_chi_p,
                "pinned_chi_second": pinned_chi_pp,
                "pinned_radial_bracket_total": pinned_radial_bracket_total,
                "pinned_radial_bracket_metric": pinned_radial_bracket_metric,
                "pinned_Rr_r_from_bracket": pinned_Rr_r_from_bracket,
                "pinned_Rr_r_metric_from_bracket": pinned_Rr_r_metric_from_bracket,
                **{f"pinned_radial_term_{key}": values for key, values in pinned_radial_bracket_terms.items()},
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
                    "a_second": float(native_a_pp[index]),
                    "b_prime": float(native_b_p[index]),
                    "b_second": float(native_b_pp[index]),
                    "X_prime": float(native_X_p[index]),
                    "X_second": float(native_X_pp[index]),
                    "A_profile": float(native_A_profile[index]),
                    "A_prime": float(native_A_p[index]),
                    "A_second": float(native_A_pp[index]),
                    "areal_radius_profile": float(native_areal_profile[index]),
                    "areal_radius_prime": float(native_R_p[index]),
                    "areal_radius_second": float(native_R_pp[index]),
                    "areal_radius_proper_prime": float(native_R_proper_p[index]),
                    "areal_radius_proper_second": float(native_R_proper_pp[index]),
                    "Rr_r_metric": float(native_geometry["Rr_r_metric"][index]),
                    "Rtheta_theta_metric": float(native_geometry["Rtheta_theta_metric"][index]),
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
                    "a_second": float(pinned_a_pp[index]),
                    "b_prime": float(pinned_b_p[index]),
                    "b_second": float(pinned_b_pp[index]),
                    "X_prime": float(pinned_X_p[index]),
                    "X_second": float(pinned_X_pp[index]),
                    "A_profile": float(pinned_A_profile[index]),
                    "A_prime": float(pinned_A_p[index]),
                    "A_second": float(pinned_A_pp[index]),
                    "areal_radius_profile": float(pinned_areal_profile[index]),
                    "areal_radius_prime": float(pinned_R_p[index]),
                    "areal_radius_second": float(pinned_R_pp[index]),
                    "areal_radius_proper_prime": float(pinned_R_proper_p[index]),
                    "areal_radius_proper_second": float(pinned_R_proper_pp[index]),
                    "chi_prime": float(pinned_chi_p[index]),
                    "chi_second": float(pinned_chi_pp[index]),
                    "radial_bracket_total": float(pinned_radial_bracket_total[index]),
                    "radial_bracket_metric": float(pinned_radial_bracket_metric[index]),
                    "radial_bracket_term_values": {
                        key: float(values[index]) for key, values in pinned_radial_bracket_terms.items()
                    },
                    "Rr_r_from_bracket": float(pinned_Rr_r_from_bracket[index]),
                    "Rr_r_metric_from_bracket": float(pinned_Rr_r_metric_from_bracket[index]),
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
                    "schema": "v55-primary-l2-cmc-refinement/v7",
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
