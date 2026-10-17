"""Observation-only native CMC stage RHS shadow on copied stage slices.

The stage-transition ordering is pinned by test_native_cmc_stage_contract_preflight.py.
This witness builds only stage 0, the explicit predictor, and the primary
predictor from copies using the existing production-owned CMC/RHS/boundary
helpers. It never calls V55ProductionKernel.step. The comparison is restricted
to the explicit metric fields a/b/X; alpha/beta RHSs belong to the moving-
puncture gauge and are not compared against production CMC gauge evolution.

PR #340 remains the initialized frozen-slice/source decomposition baseline.
This is one fixed-grid stage observation, not a repeated resolution matrix.
"""
from __future__ import annotations

import hashlib
import inspect
import json

import numpy as np

from engine.matter_system import Species, validate_species_admissibility
from engine.native_geometry import (
    NativeSphericalCellGrid,
    native_bssn_geometry_terms,
    native_cmc_explicit_metric_rhs,
    native_lambda_l2_rhs,
    native_primary_l2_rhs,
)
from engine.production_kernel import LAMBDA_M, ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState, metric_slice_from_q
from engine.v55_pirk_adapter import vendor_modules


_GEOMETRY_FIELDS = ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B")
_SCALAR_FIELDS = ("S", "PS", "D", "PD", "phi", "Pi")
_MATTER_SPECIES = ("dark_matter", "baryons", "radiation", "charged_dust")
_MATTER_FIELDS = ("rest", "energy_t", "momentum_r")
_EXPLICIT_METRIC_FIELDS = ("a", "b", "X")
_DT = 1.0e-3


def _slice_sha256(geometry, scalars, matter):
    digest = hashlib.sha256()
    for label, obj, fields in (
        ("geometry", geometry, _GEOMETRY_FIELDS),
        ("scalars", scalars, _SCALAR_FIELDS),
    ):
        for name in fields:
            values = np.ascontiguousarray(np.asarray(getattr(obj, name)))
            digest.update(label.encode("ascii") + b"\0" + name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))

    for species_name in _MATTER_SPECIES:
        species = getattr(matter, species_name, None)
        if species is None:
            continue
        for name in _MATTER_FIELDS:
            values = np.ascontiguousarray(np.asarray(getattr(species, name)))
            digest.update(b"matter\0" + species_name.encode("ascii") + b"\0")
            digest.update(name.encode("ascii") + b"\0")
            digest.update(values.dtype.str.encode("ascii") + b"\0")
            digest.update(repr(values.shape).encode("ascii") + b"\0")
            digest.update(values.tobytes(order="C"))
    return digest.hexdigest()


def _regional_rhs_residuals(native_rhs, adapter_rhs, geometry, native_grid):
    areal = np.asarray(
        native_bssn_geometry_terms(native_grid, geometry)["areal_radius"],
        dtype=float,
    )
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    summaries = {}
    for region, mask in masks.items():
        selected_count = int(np.count_nonzero(mask))
        assert selected_count >= 4, (
            f"{region}: expected at least four descriptive cells, got {selected_count}"
        )
        fields = {}
        for name in _EXPLICIT_METRIC_FIELDS:
            native = np.asarray(native_rhs[name], dtype=float)
            reference = np.asarray(adapter_rhs[name], dtype=float)
            assert native.shape == reference.shape == areal.shape, (
                f"{region}/{name}: incompatible profile shapes"
            )
            residual = native - reference
            assert np.all(np.isfinite(residual)), (
                f"{region}/{name}: non-finite native-minus-adapter residual"
            )
            # Exact same-input parity was observed on all three CMC-aligned
            # stages and all radial regions. Enforce it here rather than letting
            # bitwise equality remain descriptive-only JSON metadata.
            assert np.array_equal(native, reference), (
                f"{region}/{name}: native explicit RHS differs bitwise from the "
                "pinned adapter on identical CMC stage inputs"
            )
            selected = residual[mask]
            fields[name] = {
                "signed_mean_native_minus_adapter": float(np.mean(selected)),
                "rms_native_minus_adapter": float(np.sqrt(np.mean(selected**2))),
                "max_abs_native_minus_adapter": float(np.max(np.abs(selected))),
                "bitwise_equal_on_region": bool(np.array_equal(native[mask], reference[mask])),
            }
        summaries[region] = {"cells": selected_count, "fields": fields}
    return summaries


def _regional_primary_l2_residuals(native_rhs, adapter_rhs, geometry, native_grid):
    """Describe native-minus-pinned Aa/K L2 residuals on identical CMC slices.

    This is a diagnostic report only: unlike the established explicit a/b/X
    parity contract, the primary-L2 residuals have no pass/fail tolerance here.
    """
    areal = np.asarray(
        native_bssn_geometry_terms(native_grid, geometry)["areal_radius"],
        dtype=float,
    )
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    summaries = {}
    for region, mask in masks.items():
        selected_count = int(np.count_nonzero(mask))
        assert selected_count >= 4, (
            f"{region}: expected at least four descriptive cells, got {selected_count}"
        )
        fields = {}
        for name in ("Aa", "K"):
            native = np.asarray(native_rhs[name], dtype=float)
            reference = np.asarray(adapter_rhs[name], dtype=float)
            assert native.shape == reference.shape == areal.shape, (
                f"{region}/{name}: incompatible primary-L2 profile shapes"
            )
            residual = native - reference
            assert np.all(np.isfinite(residual)), (
                f"{region}/{name}: non-finite native-minus-adapter primary-L2 residual"
            )
            selected = residual[mask]
            fields[name] = {
                "signed_mean_native_minus_adapter": float(np.mean(selected)),
                "rms_native_minus_adapter": float(np.sqrt(np.mean(selected**2))),
                "max_abs_native_minus_adapter": float(np.max(np.abs(selected))),
            }
        summaries[region] = {"cells": selected_count, "fields": fields}
    return {
        "regions": summaries,
        "first_four_signed_native_minus_adapter": {
            name: [
                float(value)
                for value in (
                    np.asarray(native_rhs[name], dtype=float)
                    - np.asarray(adapter_rhs[name], dtype=float)
                )[:4]
            ]
            for name in ("Aa", "K")
        },
        "interpretation": (
            "same-input CMC-stage software diagnostic; descriptive signed residuals "
            "only, no added centre/annulus/outer acceptance threshold"
        ),
    }


def _regional_lambda_l2_residuals(native_rhs, adapter_rhs, geometry, native_grid):
    """Report same-input geometry-only Lambda-L2 residuals by physical region.

    Lambda-L3 advection and the local matter-current source are separate blocks
    and deliberately excluded. No numerical acceptance threshold is added.
    """
    areal = np.asarray(
        native_bssn_geometry_terms(native_grid, geometry)["areal_radius"],
        dtype=float,
    )
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    native = np.asarray(native_rhs, dtype=float)
    reference = np.asarray(adapter_rhs, dtype=float)
    assert native.shape == reference.shape == areal.shape, (
        "Lambda-L2: incompatible native/adapter profile shapes"
    )
    residual = native - reference
    assert np.all(np.isfinite(residual)), (
        "Lambda-L2: non-finite native-minus-adapter residual"
    )
    summaries = {}
    for region, mask in masks.items():
        selected_count = int(np.count_nonzero(mask))
        assert selected_count >= 4, (
            f"{region}: expected at least four descriptive cells, got {selected_count}"
        )
        selected = residual[mask]
        summaries[region] = {
            "cells": selected_count,
            "signed_mean_native_minus_adapter": float(np.mean(selected)),
            "rms_native_minus_adapter": float(np.sqrt(np.mean(selected**2))),
            "max_abs_native_minus_adapter": float(np.max(np.abs(selected))),
        }
    return {
        "regions": summaries,
        "first_four_signed_native_minus_adapter": [
            float(value) for value in residual[:4]
        ],
        "interpretation": (
            "same-input geometry-only Lambda-L2 software diagnostic; Lambda-L3 "
            "advection and matter-current source are excluded; descriptive only"
        ),
    }


def _primary_l2_term_split_shadow(
    native_geometry_terms,
    pinned_geometry_terms,
    native_primary_l2,
    pinned_primary_l2,
    geometry,
    native_grid,
):
    """Describe lapse-Hessian and trace-free Ricci contributions on one CMC slice.

    The term definitions are reused from the established #346/#348 primary-L2
    diagnostic contract. This is report-only; no centre or region threshold is
    added and the pinned routine is a comparator rather than formula authority.
    """
    alpha = np.asarray(geometry.alpha, dtype=float)

    def split(terms):
        lapse = -(
            np.asarray(terms["DrDr_alpha"], dtype=float)
            - np.asarray(terms["lap_alpha"], dtype=float) / 3.0
        )
        ricci = alpha * (
            np.asarray(terms["Rr_r"], dtype=float)
            - np.asarray(terms["R"], dtype=float) / 3.0
        )
        return lapse, ricci, lapse + ricci

    native_lapse, native_ricci, native_total = split(native_geometry_terms)
    pinned_lapse, pinned_ricci, pinned_total = split(pinned_geometry_terms)
    terms = {
        "aa_lapse_hessian_tracefree": native_lapse - pinned_lapse,
        "aa_tracefree_ricci": native_ricci - pinned_ricci,
        "aa_reconstructed_from_terms": native_total - pinned_total,
        "aa_native_rhs_closure": np.asarray(native_primary_l2["Aa"], dtype=float) - native_total,
        "aa_pinned_rhs_closure": np.asarray(pinned_primary_l2["Aa"], dtype=float) - pinned_total,
    }
    areal = np.asarray(
        native_bssn_geometry_terms(native_grid, geometry)["areal_radius"],
        dtype=float,
    )
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    regional = {}
    for term_name, residual in terms.items():
        residual = np.asarray(residual, dtype=float)
        assert residual.shape == areal.shape, f"{term_name}: incompatible profile shape"
        assert np.all(np.isfinite(residual)), f"{term_name}: non-finite term residual"
        regions = {}
        for region_name, mask in masks.items():
            count = int(np.count_nonzero(mask))
            assert count >= 4, f"{region_name}: expected at least four cells"
            selected = residual[mask]
            regions[region_name] = {
                "cells": count,
                "signed_mean": float(np.mean(selected)),
                "rms": float(np.sqrt(np.mean(selected**2))),
                "max_abs": float(np.max(np.abs(selected))),
            }
        regional[term_name] = {
            "regions": regions,
            "first_four_signed_values": [float(value) for value in residual[:4]],
        }
    return {
        "component_residuals": regional,
        "term_definition": (
            "Aa_L2 = -(D_r D_r alpha - lap(alpha)/3) "
            "+ alpha*(R^r_r - R/3); use the same CMC slice and local alpha"
        ),
        "interpretation": (
            "test-only term attribution on an already reconstructed CMC stage; "
            "signed residuals and closure are descriptive, no new acceptance threshold"
        ),
    }


def _metric_lambda_chain_shadow(
    native_geometry,
    pinned_geometry,
    geometry,
    native_grid,
    pinned_grid,
    vacuum,
):
    """Measure metric-Lambda derivative chain residual on one copied CMC slice.

    This compares D1(Lambda_metric) with the derivative of the explicit
    metric-to-Lambda relation built from the same D1/D2 operands. Its predicted
    TF-Ricci response is offline-only and keeps evolved deltaLambda correction
    separate. No chain-rule substitution is made in the native or pinned RHS.
    """
    r = np.asarray(native_grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    alpha = np.asarray(geometry.alpha, dtype=float)
    Lambda = np.asarray(geometry.Lambda, dtype=float)

    def channels(d1, d2):
        ap, bp = d1(a, 1), d1(b, 1)
        app, bpp = d2(a, 1), d2(b, 1)
        term_a = ap / (2.0 * a**2)
        term_b = -bp / (a * b)
        delta = 1.0 / b - 1.0 / a
        term_radial = 2.0 * delta / r
        lambda_metric = term_a + term_b + term_radial

        lambda_prime_discrete = d1(lambda_metric, -1)
        chain_terms = {
            "a_prime_over_2a_squared": app / (2.0 * a**2) - ap**2 / a**3,
            "minus_b_prime_over_ab": (
                -bpp / (a * b)
                + bp * (ap * b + a * bp) / (a**2 * b**2)
            ),
            "radial_inverse_metric_term": (
                -2.0 * delta / r**2
                + 2.0 * (-bp / b**2 + ap / a**2) / r
            ),
        }
        discrete_terms = {
            "a_prime_over_2a_squared": d1(term_a, -1),
            "minus_b_prime_over_ab": d1(term_b, -1),
            "radial_inverse_metric_term": d1(term_radial, -1),
        }
        residual_terms = {
            name: chain_terms[name] - discrete_terms[name] for name in chain_terms
        }
        chain_total = sum(chain_terms.values())
        delta_total = chain_total - lambda_prime_discrete
        term_recombined = sum(residual_terms.values())

        delta_lambda = Lambda - lambda_metric
        delta_lambda_prime = d1(delta_lambda, -1)
        evolved_lambda_prime = d1(Lambda, -1)

        return {
            "lambda_metric": lambda_metric,
            "lambda_prime_discrete": lambda_prime_discrete,
            "lambda_prime_chain": chain_total,
            "lambda_prime_chain_residual": delta_total,
            "chain_term_residuals": residual_terms,
            "chain_term_recombined": term_recombined,
            "lambda_prime_chain_reconstruction_closure": chain_total - sum(chain_terms.values()),
            "lambda_prime_term_residual_closure": term_recombined - delta_total,
            "lambda_metric_reconstruction_difference": (
                lambda_metric - np.asarray(
                    native_geometry["Lambda_metric"], dtype=float
                )
            ),
            "evolved_lambda_prime": evolved_lambda_prime,
            "delta_lambda_prime": delta_lambda_prime,
            "evolved_lambda_derivative_identity_closure": (
                evolved_lambda_prime - (lambda_prime_discrete + delta_lambda_prime)
            ),
        }

    native_d1 = lambda values, parity=1: np.asarray(
        native_grid.cell_derivative_fourth(values, parity=parity), dtype=float
    )
    native_d2 = lambda values, parity=1: np.asarray(
        native_grid.cell_second_derivative_fourth(values, parity=parity), dtype=float
    )
    pinned_d1 = lambda values, parity=1: np.asarray(
        vacuum._d(pinned_grid, values, parity), dtype=float
    )
    pinned_d2 = lambda values, parity=1: np.asarray(
        vacuum._d2(pinned_grid, values, parity), dtype=float
    )

    native = channels(native_d1, native_d2)
    pinned = channels(pinned_d1, pinned_d2)

    # The TF-Ricci response to changing the metric-Lambda prime operand is
    # +(2/3)*alpha*X^2*delta(Lambda_metric'). Build this only on the pinned
    # metric-channel remainder after the separately evolved deltaLambda
    # connection correction has been removed.
    response_factor = (2.0 / 3.0) * alpha * X**2
    pinned_response_terms = {
        name: response_factor * values
        for name, values in pinned["chain_term_residuals"].items()
    }
    pinned_response = sum(pinned_response_terms.values())

    pinned_total_tf = alpha * (
        np.asarray(pinned_geometry["Rr_r"], dtype=float)
        - np.asarray(pinned_geometry["R"], dtype=float) / 3.0
    )
    pinned_a_prime = pinned_d1(a, 1)
    pinned_delta_lambda = Lambda - pinned["lambda_metric"]
    pinned_delta_lambda_prime = pinned_d1(pinned_delta_lambda, -1)
    pinned_connection_tf = alpha * X**2 * (
        (2.0 / 3.0) * pinned_delta_lambda_prime
        + 0.5 * (pinned_a_prime / a) * pinned_delta_lambda
    )
    pinned_metric_tf = pinned_total_tf - pinned_connection_tf
    native_metric_tf = alpha * (
        np.asarray(native_geometry["Rr_r_metric"], dtype=float)
        - np.asarray(native_geometry["R_metric"], dtype=float) / 3.0
    )
    baseline_metric_residual = native_metric_tf - pinned_metric_tf
    counterfactual_metric_residual = baseline_metric_residual - pinned_response

    profiles = {
        "native_lambda_prime_chain_residual": native["lambda_prime_chain_residual"],
        "pinned_lambda_prime_chain_residual": pinned["lambda_prime_chain_residual"],
        "native_chain_residual_a_term": native["chain_term_residuals"]["a_prime_over_2a_squared"],
        "native_chain_residual_b_term": native["chain_term_residuals"]["minus_b_prime_over_ab"],
        "native_chain_residual_radial_term": native["chain_term_residuals"]["radial_inverse_metric_term"],
        "pinned_chain_residual_a_term": pinned["chain_term_residuals"]["a_prime_over_2a_squared"],
        "pinned_chain_residual_b_term": pinned["chain_term_residuals"]["minus_b_prime_over_ab"],
        "pinned_chain_residual_radial_term": pinned["chain_term_residuals"]["radial_inverse_metric_term"],
        "pinned_chain_response_a_term": pinned_response_terms["a_prime_over_2a_squared"],
        "pinned_chain_response_b_term": pinned_response_terms["minus_b_prime_over_ab"],
        "pinned_chain_response_radial_term": pinned_response_terms["radial_inverse_metric_term"],
        "pinned_chain_response_total": pinned_response,
        "metric_tf_residual_before_counterfactual": baseline_metric_residual,
        "metric_tf_residual_after_chain_counterfactual": counterfactual_metric_residual,
        "metric_tf_counterfactual_delta": (
            counterfactual_metric_residual - baseline_metric_residual
        ),
        "pinned_lambda_prime_term_reconstruction_closure": (
            pinned["lambda_prime_chain_residual"] - sum(
                pinned["chain_term_residuals"].values()
            )
        ),
        "pinned_evolved_lambda_derivative_identity_closure": (
            pinned["evolved_lambda_derivative_identity_closure"]
        ),
        "pinned_lambda_metric_reconstruction_difference": (
            pinned["lambda_metric_reconstruction_difference"]
        ),
        "native_lambda_metric_reconstruction_difference": (
            native["lambda_metric_reconstruction_difference"]
        ),
    }

    areal = np.asarray(native_geometry["areal_radius"], dtype=float)
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    report = {}
    for name, raw in profiles.items():
        values = np.asarray(raw, dtype=float)
        assert values.shape == areal.shape, f"{name}: incompatible profile shape"
        assert np.all(np.isfinite(values)), f"{name}: non-finite Lambda-chain profile"
        regions = {}
        for region_name, mask in masks.items():
            cells = int(np.count_nonzero(mask))
            assert cells >= 4, f"{region_name}: expected at least four cells"
            selected = values[mask]
            regions[region_name] = {
                "cells": cells,
                "signed_mean": float(np.mean(selected)),
                "rms": float(np.sqrt(np.mean(selected**2))),
                "max_abs": float(np.max(np.abs(selected))),
            }
        report[name] = {
            "regions": regions,
            "first_four_signed_values": [float(value) for value in values[:4]],
        }

    return {
        "profiles": report,
        "definitions": {
            "metric_Lambda": "D1(a)/(2*a^2) - D1(b)/(a*b) + 2/r*(1/b-1/a)",
            "chain_derivative": (
                "D2(a)/(2*a^2)-D1(a)^2/a^3 - D2(b)/(a*b) "
                "+ D1(b)*(D1(a)*b+a*D1(b))/(a^2*b^2) "
                "- 2/r^2*(1/b-1/a) + 2/r*(-D1(b)/b^2+D1(a)/a^2)"
            ),
            "TF_Ricci_response": "(2/3)*alpha*X^2*(chain_prime - discrete_D1(Lambda_metric))",
            "counterfactual": (
                "subtract this response from the baseline native-minus-pinned "
                "deltaLambda-aligned metric-channel residual; diagnostic only"
            ),
        },
        "interpretation": (
            "same-input discrete derivative attribution. An offline chain-rule "
            "counterfactual is not an evolved BSSN operator and is not a production "
            "substitution recommendation."
        ),
    }


def _primary_l2_ricci_term_contributions(
    native_geometry,
    pinned_geometry,
    geometry,
    native_grid,
    pinned_grid,
    vacuum,
):
    """Compare signed native warped-product and pinned BSSN TF-Ricci contributions.

    The decompositions are deliberately representation-specific: native uses
    proper-radial / areal-radius curvature terms; pinned uses its expanded
    radial-minus-one-third-scalar bracket. This is not a one-to-one mapping
    between physical terms, and component magnitudes are not error diagnoses.
    """
    r = np.asarray(native_grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    alpha = np.asarray(geometry.alpha, dtype=float)
    Lambda = np.asarray(geometry.Lambda, dtype=float)
    areal = np.asarray(native_geometry["areal_radius"], dtype=float)

    # Native warped-product representation, on its own direct composite
    # profile derivatives.
    A = a / X**2
    R_p = np.asarray(
        native_grid.cell_derivative_fourth(areal, parity=-1), dtype=float
    )
    R_pp = np.asarray(
        native_grid.cell_second_derivative_fourth(areal, parity=-1), dtype=float
    )
    A_p = np.asarray(
        native_grid.cell_derivative_fourth(A, parity=1), dtype=float
    )
    R_l = R_p / np.sqrt(A)
    R_ll = (R_pp - 0.5 * (A_p / A) * R_p) / A
    native_components = {
        "proper_radial_curvature": -(2.0 / 3.0) * alpha * R_ll / areal,
        "angular_curvature": (2.0 / 3.0) * alpha * (R_l**2 - 1.0) / areal**2,
    }
    native_reconstructed = sum(native_components.values())
    native_direct = alpha * (
        np.asarray(native_geometry["Rr_r_metric"], dtype=float)
        - np.asarray(native_geometry["R_metric"], dtype=float) / 3.0
    )

    # Pinned BSSN representation. Build terms of (radial bracket -
    # scalar bracket/3), retaining each pinned derivative operation.
    d1 = lambda values, parity=1: np.asarray(
        vacuum._d(pinned_grid, values, parity), dtype=float
    )
    d2 = lambda values, parity=1: np.asarray(
        vacuum._d2(pinned_grid, values, parity), dtype=float
    )
    ap, bp = d1(a, 1), d1(b, 1)
    app, bpp = d2(a, 1), d2(b, 1)
    Xp = d1(X, 1)
    Xpp = d2(X, 1)
    chip = -0.5 * Xp / X
    chipp = -0.5 * (Xpp / X - (Xp / X) ** 2)
    Lambda_p = d1(Lambda, -1)
    metric_lambda = (1.0 - a / b) / r**2

    bracket_terms = {
        "second_derivatives_a_b": app / (3.0 * a) - bpp / (3.0 * b),
        "evolved_Lambda_prime": -(2.0 / 3.0) * a * Lambda_p,
        "quadratic_a_b_first_derivatives": (
            -(5.0 / 12.0) * (ap / a) ** 2
            + (1.0 / 3.0) * (bp / b) ** 2
        ),
        "undifferentiated_Lambda_times_a_prime": -0.5 * Lambda * ap,
        "radial_a_prime": ap / (r * b),
        "reference_metric_lambda": (
            (2.0 / 3.0) * metric_lambda
            + 2.0 * metric_lambda * r * bp / b
        ),
        "scalar_angular_b_prime": (
            -(2.0 / 3.0) * (3.0 - a / b) * bp / (r * b)
        ),
        "chi_second_derivative": (4.0 / 3.0) * chipp,
        "chi_prime_squared": -(8.0 / 3.0) * chip**2,
        "chi_prime_first_derivative_coupling": (
            -(2.0 / 3.0) * chip * (ap / a + bp / b + 2.0 / r)
        ),
    }
    inverse_physical_radial = X**2 / a
    pinned_components = {
        name: -alpha * inverse_physical_radial * contribution
        for name, contribution in bracket_terms.items()
    }
    pinned_reconstructed = sum(pinned_components.values())
    pinned_direct = alpha * (
        np.asarray(pinned_geometry["Rr_r"], dtype=float)
        - np.asarray(pinned_geometry["R"], dtype=float) / 3.0
    )

    profiles = {
        "native_component_reconstruction_closure": native_reconstructed - native_direct,
        "pinned_component_reconstruction_closure": pinned_reconstructed - pinned_direct,
        "total_residual_from_components": native_reconstructed - pinned_reconstructed,
        "total_residual_direct": native_direct - pinned_direct,
        "total_residual_component_closure": (
            (native_reconstructed - pinned_reconstructed)
            - (native_direct - pinned_direct)
        ),
    }
    profiles.update({
        "native_" + name: values for name, values in native_components.items()
    })
    profiles.update({
        "pinned_" + name: values for name, values in pinned_components.items()
    })

    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    report = {}
    for name, raw in profiles.items():
        values = np.asarray(raw, dtype=float)
        assert values.shape == areal.shape, f"{name}: incompatible profile shape"
        assert np.all(np.isfinite(values)), f"{name}: non-finite term contribution"
        regions = {}
        for region_name, mask in masks.items():
            cells = int(np.count_nonzero(mask))
            assert cells >= 4, f"{region_name}: expected at least four cells"
            selected = values[mask]
            regions[region_name] = {
                "cells": cells,
                "signed_mean": float(np.mean(selected)),
                "rms": float(np.sqrt(np.mean(selected**2))),
                "max_abs": float(np.max(np.abs(selected))),
            }
        report[name] = {
            "regions": regions,
            "first_four_signed_values": [float(value) for value in values[:4]],
        }

    return {
        "profiles": report,
        "native_component_definition": (
            "TF_Ricci = -(2/3)*alpha*R_ll/R + "
            "(2/3)*alpha*(R_l^2-1)/R^2 from warped-product curvature"
        ),
        "pinned_component_definition": (
            "TF_Ricci = -(alpha*X^2/a)*(radial_bracket - scalar_bracket/3), "
            "expanded into signed BSSN bracket contributions"
        ),
        "pinned_bracket_terms": {
            name: "alpha*X^2/a with overall minus sign applied to: " + name
            for name in bracket_terms
        },
        "interpretation": (
            "same-input copied CMC-stage arithmetic attribution only; native and "
            "pinned components are different representations, not termwise physical "
            "counterparts. Closures test reconstruction only; no term is designated "
            "wrong and no scientific-admission or production-correction claim follows."
        ),
    }


def _primary_l2_ricci_primitive_shadow(
    native_geometry,
    pinned_geometry,
    geometry,
    native_grid,
    pinned_grid,
    vacuum,
):
    """Attribute trace-free Ricci residuals to metric and connection pieces.

    This reuses the initialized-slice #346/#348 primitive convention, now on
    the identical already-copied CMC stage. Native and pinned derivative
    channels use their own existing stencils; no formula or threshold is tuned.
    """
    r = np.asarray(native_grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    alpha = np.asarray(geometry.alpha, dtype=float)
    Lambda = np.asarray(geometry.Lambda, dtype=float)

    native_lambda_metric = np.asarray(native_geometry["Lambda_metric"], dtype=float)
    native_delta_lambda = np.asarray(
        native_geometry["connection_constraint"], dtype=float
    )
    native_lambda_p = np.asarray(
        native_grid.cell_derivative_fourth(Lambda, parity=-1), dtype=float
    )
    native_lambda_metric_p = np.asarray(
        native_grid.cell_derivative_fourth(native_lambda_metric, parity=-1),
        dtype=float,
    )
    native_delta_lambda_p = np.asarray(
        native_grid.cell_derivative_fourth(native_delta_lambda, parity=-1),
        dtype=float,
    )

    # Match the pinned convention already used by the initialized-slice
    # primitive witness: vendor derivative helpers, not native stencils.
    pinned_a_p = np.asarray(vacuum._d(pinned_grid, a, 1), dtype=float)
    pinned_b_p = np.asarray(vacuum._d(pinned_grid, b, 1), dtype=float)
    pinned_lambda_metric = (
        pinned_a_p / (2.0 * a**2)
        - pinned_b_p / (a * b)
        + 2.0 * (1.0 / b - 1.0 / a) / r
    )
    pinned_delta_lambda = Lambda - pinned_lambda_metric
    pinned_lambda_p = np.asarray(vacuum._d(pinned_grid, Lambda, -1), dtype=float)
    pinned_lambda_metric_p = np.asarray(
        vacuum._d(pinned_grid, pinned_lambda_metric, -1), dtype=float
    )
    pinned_delta_lambda_p = np.asarray(
        vacuum._d(pinned_grid, pinned_delta_lambda, -1), dtype=float
    )

    native_metric_tf = alpha * (
        np.asarray(native_geometry["Rr_r_metric"], dtype=float)
        - np.asarray(native_geometry["R_metric"], dtype=float) / 3.0
    )
    native_connection_tf = alpha * (
        np.asarray(native_geometry["delta_Rr_r_connection"], dtype=float)
        - np.asarray(native_geometry["delta_R_scalar_connection"], dtype=float) / 3.0
    )

    # The pinned geometry terms expose total Ricci components. Isolate the
    # existing BSSN connection-form correction with the same primitive formula
    # documented in #348, then recover the metric-derived remainder by subtraction.
    pinned_connection_Rr_r = X**2 * (
        pinned_lambda_p + 0.5 * (pinned_a_p / a) * Lambda
    )
    pinned_connection_R = X**2 * pinned_lambda_p
    pinned_metric_Rr_r = (
        np.asarray(pinned_geometry["Rr_r"], dtype=float) - pinned_connection_Rr_r
    )
    pinned_metric_R = (
        np.asarray(pinned_geometry["R"], dtype=float) - pinned_connection_R
    )
    pinned_metric_tf = alpha * (pinned_metric_Rr_r - pinned_metric_R / 3.0)
    pinned_connection_tf = alpha * (
        pinned_connection_Rr_r - pinned_connection_R / 3.0
    )

    native_total_tf = alpha * (
        np.asarray(native_geometry["Rr_r"], dtype=float)
        - np.asarray(native_geometry["R"], dtype=float) / 3.0
    )
    pinned_total_tf = alpha * (
        np.asarray(pinned_geometry["Rr_r"], dtype=float)
        - np.asarray(pinned_geometry["R"], dtype=float) / 3.0
    )

    # Convention-aligned split: native geometry adds the delta-Lambda
    # connection correction to its metric Ricci, while the pinned BSSN
    # expression uses full Lambda inside its principal brackets. Subtracting
    # the full-Lambda term (the legacy split below) is a different allocation.
    # Here subtract only the pinned delta-Lambda correction, matching native.
    aligned_pinned_connection_Rr_r = X**2 * (
        pinned_delta_lambda_p
        + 0.5 * (pinned_a_p / a) * pinned_delta_lambda
    )
    aligned_pinned_connection_R = X**2 * pinned_delta_lambda_p
    aligned_pinned_metric_Rr_r = (
        np.asarray(pinned_geometry["Rr_r"], dtype=float)
        - aligned_pinned_connection_Rr_r
    )
    aligned_pinned_metric_R = (
        np.asarray(pinned_geometry["R"], dtype=float)
        - aligned_pinned_connection_R
    )
    aligned_pinned_metric_tf = alpha * (
        aligned_pinned_metric_Rr_r - aligned_pinned_metric_R / 3.0
    )
    aligned_pinned_connection_tf = alpha * (
        aligned_pinned_connection_Rr_r - aligned_pinned_connection_R / 3.0
    )
    aligned_residual_profiles = {
        "metric_tracefree_ricci": native_metric_tf - aligned_pinned_metric_tf,
        "connection_tracefree_ricci": native_connection_tf - aligned_pinned_connection_tf,
        "total_tracefree_ricci": native_total_tf - pinned_total_tf,
        "ricci_component_closure": (
            (native_metric_tf - aligned_pinned_metric_tf)
            + (native_connection_tf - aligned_pinned_connection_tf)
            - (native_total_tf - pinned_total_tf)
        ),
        "pinned_delta_lambda_split_closure": pinned_total_tf - (
            aligned_pinned_metric_tf + aligned_pinned_connection_tf
        ),
        "connection_correction_native_minus_pinned_delta_lambda": (
            native_connection_tf - aligned_pinned_connection_tf
        ),
    }

    residual_profiles = {
        "metric_tracefree_ricci": native_metric_tf - pinned_metric_tf,
        "connection_tracefree_ricci": native_connection_tf - pinned_connection_tf,
        "total_tracefree_ricci": native_total_tf - pinned_total_tf,
        "ricci_component_closure": (
            (native_metric_tf - pinned_metric_tf)
            + (native_connection_tf - pinned_connection_tf)
            - (native_total_tf - pinned_total_tf)
        ),
        "lambda_metric": native_lambda_metric - pinned_lambda_metric,
        "delta_lambda": native_delta_lambda - pinned_delta_lambda,
        "lambda_prime": native_lambda_p - pinned_lambda_p,
        "lambda_metric_prime": native_lambda_metric_p - pinned_lambda_metric_p,
        "delta_lambda_prime": native_delta_lambda_p - pinned_delta_lambda_p,
    }

    areal = np.asarray(native_geometry["areal_radius"], dtype=float)
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    report = {}
    for name, residual_raw in residual_profiles.items():
        residual = np.asarray(residual_raw, dtype=float)
        assert residual.shape == areal.shape, f"{name}: incompatible profile shape"
        assert np.all(np.isfinite(residual)), f"{name}: non-finite primitive residual"
        regions = {}
        for region_name, mask in masks.items():
            cells = int(np.count_nonzero(mask))
            assert cells >= 4, f"{region_name}: expected at least four cells"
            selected = residual[mask]
            regions[region_name] = {
                "cells": cells,
                "signed_mean": float(np.mean(selected)),
                "rms": float(np.sqrt(np.mean(selected**2))),
                "max_abs": float(np.max(np.abs(selected))),
            }
        report[name] = {
            "regions": regions,
            "first_four_signed_native_minus_pinned": [
                float(value) for value in residual[:4]
            ],
        }

    convention_aligned_report = {}
    for name, residual_raw in aligned_residual_profiles.items():
        residual = np.asarray(residual_raw, dtype=float)
        assert residual.shape == areal.shape, f"aligned {name}: incompatible profile shape"
        assert np.all(np.isfinite(residual)), f"aligned {name}: non-finite primitive residual"
        regions = {}
        for region_name, mask in masks.items():
            cells = int(np.count_nonzero(mask))
            assert cells >= 4, f"{region_name}: expected at least four cells"
            selected = residual[mask]
            regions[region_name] = {
                "cells": cells,
                "signed_mean": float(np.mean(selected)),
                "rms": float(np.sqrt(np.mean(selected**2))),
                "max_abs": float(np.max(np.abs(selected))),
            }
        convention_aligned_report[name] = {
            "regions": regions,
            "first_four_signed_native_minus_pinned": [
                float(value) for value in residual[:4]
            ],
        }
    return {
        "profiles": report,
        "convention_aligned_profiles": convention_aligned_report,
        "metric_component_definition": (
            "alpha*(Rr_r_metric - R_metric/3), using each implementation's own "
            "metric Ricci primitives"
        ),
        "connection_component_definition": (
            "legacy PR #367 allocation: alpha*(delta_Rr_r_connection - "
            "delta_R_scalar_connection/3); pinned side extracts a full-Lambda term"
        ),
        "convention_aligned_split_definition": (
            "native and pinned sides both extract only the deltaLambda connection "
            "correction: X^2*(D(deltaLambda) + 0.5*(D(a)/a)*deltaLambda) "
            "for Rr_r and X^2*D(deltaLambda) for R; the pinned metric remainder "
            "is then obtained by subtraction from pinned total Ricci"
        ),
        "interpretation": (
            "same-input CMC-stage primitive attribution; descriptive residuals only, "
            "no new centre/annulus/outer acceptance threshold"
        ),
    }


def _metric_product_chain_rule_shadow(geometry, native_grid):
    """Report discrete product/chain-rule residuals on one already-copied CMC slice.

    Composite A=a/X^2 and R=r*sqrt(b)/X are differentiated directly with
    native D1/D2 operators, then compared against product/chain expressions
    built from the same operators applied to a, b, and X. These residuals
    measure finite-discrete noncommutation, not a continuum-law defect.
    """
    r = np.asarray(native_grid.centers, dtype=float)
    a = np.asarray(geometry.a, dtype=float)
    b = np.asarray(geometry.b, dtype=float)
    X = np.asarray(geometry.X, dtype=float)
    A = a / X**2
    R = r * np.sqrt(b) / X

    d1 = native_grid.cell_derivative_fourth
    d2 = native_grid.cell_second_derivative_fourth

    # Direct derivatives of the composite profiles used by native geometry.
    A_p_direct = np.asarray(d1(A, parity=1), dtype=float)
    A_pp_direct = np.asarray(d2(A, parity=1), dtype=float)
    R_p_direct = np.asarray(d1(R, parity=-1), dtype=float)
    R_pp_direct = np.asarray(d2(R, parity=-1), dtype=float)

    # Derivatives of operand profiles on the same grid and with their native
    # even parity. The expressions below are the continuum product rules,
    # evaluated discretely using these same operand derivatives.
    a_p = np.asarray(d1(a, parity=1), dtype=float)
    a_pp = np.asarray(d2(a, parity=1), dtype=float)
    b_p = np.asarray(d1(b, parity=1), dtype=float)
    b_pp = np.asarray(d2(b, parity=1), dtype=float)
    X_p = np.asarray(d1(X, parity=1), dtype=float)
    X_pp = np.asarray(d2(X, parity=1), dtype=float)

    log_A_p = a_p / a - 2.0 * X_p / X
    log_A_pp = (
        a_pp / a
        - (a_p / a) ** 2
        - 2.0 * (X_pp / X - (X_p / X) ** 2)
    )
    A_p_chain = A * log_A_p
    A_pp_chain = A * (log_A_p**2 + log_A_pp)

    q = 1.0 / r + 0.5 * b_p / b - X_p / X
    q_p = (
        -1.0 / r**2
        + 0.5 * (b_pp / b - (b_p / b) ** 2)
        - (X_pp / X - (X_p / X) ** 2)
    )
    R_p_chain = R * q
    R_pp_chain = R * (q**2 + q_p)

    # Propagate the primitive derivative residuals through the warped-product
    # conversion to proper radial distance: d/dl=A^(-1/2) d/dr.
    R_l_direct = R_p_direct / np.sqrt(A)
    R_l_chain = R_p_chain / np.sqrt(A)
    R_ll_direct = (
        R_pp_direct - 0.5 * (A_p_direct / A) * R_p_direct
    ) / A
    R_ll_chain = (
        R_pp_chain - 0.5 * (A_p_chain / A) * R_p_chain
    ) / A

    residual_profiles = {
        "A_first_derivative_product_rule": A_p_direct - A_p_chain,
        "A_second_derivative_product_rule": A_pp_direct - A_pp_chain,
        "areal_R_first_derivative_chain_rule": R_p_direct - R_p_chain,
        "areal_R_second_derivative_chain_rule": R_pp_direct - R_pp_chain,
        "proper_R_l_chain_propagation": R_l_direct - R_l_chain,
        "proper_R_ll_chain_propagation": R_ll_direct - R_ll_chain,
    }
    # The following identity decomposes the R_ll residual into the direct
    # R'' residual and the A',R' product residuals. It is an algebraic closure
    # check, not an independent numerical-accuracy threshold.
    R_ll_from_primitives = (
        (R_pp_direct - R_pp_chain)
        - 0.5 * (
            (A_p_direct / A) * R_p_direct
            - (A_p_chain / A) * R_p_chain
        )
    ) / A
    residual_profiles["R_ll_primitive_decomposition_closure"] = (
        residual_profiles["proper_R_ll_chain_propagation"] - R_ll_from_primitives
    )

    areal = R
    masks = {
        "centre_R_areal_lt_2": areal < 2.0,
        "fixed_annulus_2_le_R_areal_le_6": (areal >= 2.0) & (areal <= 6.0),
        "outer_R_areal_gt_6": areal > 6.0,
    }
    profiles = {}
    for name, residual_raw in residual_profiles.items():
        residual = np.asarray(residual_raw, dtype=float)
        assert residual.shape == areal.shape, f"{name}: incompatible profile shape"
        assert np.all(np.isfinite(residual)), f"{name}: non-finite chain-rule residual"
        regions = {}
        for region_name, mask in masks.items():
            cells = int(np.count_nonzero(mask))
            assert cells >= 4, f"{region_name}: expected at least four cells"
            selected = residual[mask]
            regions[region_name] = {
                "cells": cells,
                "signed_mean": float(np.mean(selected)),
                "rms": float(np.sqrt(np.mean(selected**2))),
                "max_abs": float(np.max(np.abs(selected))),
            }
        profiles[name] = {
            "regions": regions,
            "first_four_signed_values": [
                float(value) for value in residual[:4]
            ],
        }
    return {
        "profiles": profiles,
        "definitions": {
            "A": "a/X^2; direct D1/D2 minus product-rule expressions from operand D1/D2",
            "areal_R": "r*sqrt(b)/X; direct D1/D2 minus logarithmic chain-rule expressions",
            "proper_R_l": "D(R)/sqrt(A)",
            "proper_R_ll": "(D2(R) - 0.5*(D(A)/A)*D(R))/A",
            "operator": "native grid fourth-order parity-aware cell D1/D2; even a,b,X,A and odd R",
        },
        "interpretation": (
            "same-input copied CMC-stage discrete derivative attribution only; "
            "these residuals are not a continuum-law defect, new acceptance "
            "threshold, physical correction, or scientific-admission witness"
        ),
    }


def _observe_stage(kernel, native_grid, stage, *, evaluate_coupled_rhs):
    geometry = stage.geometry
    scalars = stage.scalars
    matter = stage.matter
    digest_before = _slice_sha256(geometry, scalars, matter)

    scalar_rhs = None
    matter_rhs = None
    if evaluate_coupled_rhs:
        scalar_rhs, matter_rhs, _ = kernel._rhs(stage)

    blocks = kernel._geometry_stage_terms(
        stage.grid, geometry, scalars, matter, stage.electromagnetic
    )
    # The CMC-only explicit metric API returns a/b/X; keep primary-L2
    # comparison separate because its known centre residual is not a parity gate.
    native_rhs = native_cmc_explicit_metric_rhs(native_grid, geometry)
    native_primary_l2 = native_primary_l2_rhs(native_grid, geometry)
    primary_l2_shadow = _regional_primary_l2_residuals(
        native_primary_l2, blocks["primary_l2"], geometry, native_grid
    )
    native_lambda_l2 = native_lambda_l2_rhs(
        native_grid, geometry, lambda_m=LAMBDA_M
    )
    lambda_l2_shadow = _regional_lambda_l2_residuals(
        native_lambda_l2, blocks["lambda_l2"], geometry, native_grid
    )
    _, vacuum, _ = vendor_modules()
    native_geometry_terms = native_bssn_geometry_terms(native_grid, geometry)
    pinned_geometry_terms = vacuum.geometry_terms(stage.grid, geometry)
    metric_product_chain_rule_shadow = _metric_product_chain_rule_shadow(
        geometry, native_grid
    )
    primary_l2_terms_shadow = _primary_l2_term_split_shadow(
        native_geometry_terms,
        pinned_geometry_terms,
        native_primary_l2,
        blocks["primary_l2"],
        geometry,
        native_grid,
    )
    primary_l2_ricci_primitives = _primary_l2_ricci_primitive_shadow(
        native_geometry_terms,
        pinned_geometry_terms,
        geometry,
        native_grid,
        stage.grid,
        vacuum,
    )
    primary_l2_ricci_term_contributions = _primary_l2_ricci_term_contributions(
        native_geometry_terms,
        pinned_geometry_terms,
        geometry,
        native_grid,
        stage.grid,
        vacuum,
    )
    metric_lambda_chain_shadow = _metric_lambda_chain_shadow(
        native_geometry_terms,
        pinned_geometry_terms,
        geometry,
        native_grid,
        stage.grid,
        vacuum,
    )
    digest_after = _slice_sha256(geometry, scalars, matter)
    assert digest_after == digest_before, (
        f"{stage.geometry!r}: production helper or RHS evaluation mutated its "
        "copied input slice"
    )

    record = {
        "stage": stage.stage_name if hasattr(stage, "stage_name") else "unnamed",
        "geometry_input_sha256_before": digest_before,
        "geometry_input_sha256_after": digest_after,
        "native_minus_adapter_explicit_metric_rhs": _regional_rhs_residuals(
            native_rhs, blocks["explicit"], geometry, native_grid
        ),
        "native_minus_adapter_primary_l2": primary_l2_shadow,
        "native_minus_adapter_lambda_l2": lambda_l2_shadow,
        "primary_l2_term_split_native_minus_pinned": primary_l2_terms_shadow,
        "primary_l2_ricci_primitive_shadow": primary_l2_ricci_primitives,
        "primary_l2_ricci_term_contributions": primary_l2_ricci_term_contributions,
        "metric_lambda_chain_shadow": metric_lambda_chain_shadow,
        "metric_product_chain_rule_shadow": metric_product_chain_rule_shadow,
        "compared_fields": list(_EXPLICIT_METRIC_FIELDS),
        "primary_l2_fields_compared": ["Aa", "K"],
        "primary_l2_stage_observation": (
            "same-input geometry-only native-vs-adapter primary-L2 diagnostics; "
            "signed regional residuals only, no equality tolerance or admission gate"
        ),
        "lambda_l2_stage_observation": (
            "same-input geometry-only native-vs-adapter Lambda-L2 diagnostics; "
            "Lambda-L3 matter-current/advection remains separately owned"
        ),
        "native_interface": "native_cmc_explicit_metric_rhs",
        "native_interface_contract": "returns exactly a/b/X; requires positive CMC lapse and beta=0; has no B-driver dependency",
        "excluded_gauge_rhs_fields": {
            "alpha": "production CMC lapse is re-solved, not advanced by moving-puncture 1+log RHS",
            "beta": "production CMC branch pins beta/B to zero instead of advancing moving-puncture shift",
        },
    }
    return blocks, scalar_rhs, matter_rhs, record


def test_native_explicit_rhs_on_cmc_stage_slices_without_calling_full_step(capsys):
    kernel = V55ProductionKernel()
    accepted = kernel.initialize(
        resolution=80, r_max=40.0, include_radiation=True
    )
    accepted_hash_before = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )

    grid = accepted.grid
    native_grid = NativeSphericalCellGrid(grid.n, grid.r_max)
    np.testing.assert_allclose(
        native_grid.centers, grid.centers, rtol=0.0, atol=1.0e-14
    )

    # Stage 0: only copies are refreshed onto the current CMC gauge.
    g0 = accepted.geometry.copy()
    s0 = accepted.scalars.copy()
    m0 = accepted.matter.copy()
    em0 = accepted.electromagnetic
    validate_species_admissibility(
        metric_slice_from_q(grid, g0), m0.radiation, Species.RADIATION,
        stage="initial",
    )
    g0.alpha = kernel._solve_lapse_for_state(
        grid, g0, s0, m0, em0
    )[0]
    g0.beta.fill(0.0)
    g0.B.fill(0.0)
    validate_species_admissibility(
        metric_slice_from_q(grid, g0), m0.radiation, Species.RADIATION,
        stage="initial",
    )

    stage0 = ProductionState(
        grid=grid,
        geometry=g0,
        scalars=s0,
        matter=m0,
        electromagnetic=em0,
        t=accepted.t,
        tau=accepted.tau,
        e_folds=accepted.e_folds,
    )
    blocks0, srhs0, mrhs0, record0 = _observe_stage(
        kernel, native_grid, stage0, evaluate_coupled_rhs=True
    )
    record0["stage"] = "stage0_cmc"

    # Explicit predictor: copy production's source-locked a/b/X update and
    # scalar/matter predictor, then refresh the lapse on that CMC slice.
    g_explicit1 = g0.copy()
    for name in ("a", "b", "X"):
        setattr(
            g_explicit1,
            name,
            getattr(g0, name) + _DT * blocks0["explicit"][name],
        )
    g_explicit1.beta.fill(0.0)
    g_explicit1.B.fill(0.0)

    s1 = ScalarFields(*(
        getattr(s0, name) + _DT * getattr(srhs0, name)
        for name in _SCALAR_FIELDS
    ))
    m1_base = tuple(
        kernel._add_matter(getattr(m0, name), mrhs0[name], _DT)
        for name in ("dark_matter", "baryons", "radiation")
    )
    m1_charged = (
        None if m0.charged_dust is None
        else kernel._add_matter(m0.charged_dust, mrhs0["charged_dust"], _DT)
    )
    m1 = V55MatterState(*m1_base, charged_dust=m1_charged)

    validate_species_admissibility(
        metric_slice_from_q(grid, g_explicit1), m1.radiation,
        Species.RADIATION, stage="explicit-predictor",
    )
    g_explicit1.alpha = kernel._solve_lapse_for_state(
        grid, g_explicit1, s1, m1, em0
    )[0]
    validate_species_admissibility(
        metric_slice_from_q(grid, g_explicit1), m1.radiation,
        Species.RADIATION, stage="explicit-predictor",
    )
    predictor = ProductionState(
        grid=grid,
        geometry=g_explicit1,
        scalars=s1,
        matter=m1,
        electromagnetic=em0,
        t=accepted.t,
        tau=accepted.tau,
        e_folds=accepted.e_folds,
    )
    blocks_pred, _, _, record_pred = _observe_stage(
        kernel, native_grid, predictor, evaluate_coupled_rhs=False
    )
    record_pred["stage"] = "explicit_predictor_cmc"

    # Primary predictor uses the existing CMC L2/L3 split and source channels.
    # Apply the established outer boundary only after Lambda is current, then
    # refresh CMC lapse and pin beta/B before observing the next RHS slice.
    g1 = g_explicit1.copy()
    g1.Aa = g0.Aa + _DT * (
        0.5 * blocks0["primary_l2"]["Aa"]
        + 0.5 * blocks_pred["primary_l2"]["Aa"]
        + blocks0["primary_l3"]["Aa"]
    )
    g1.K = g0.K + _DT * (
        0.5 * blocks0["primary_l2"]["K"]
        + 0.5 * blocks_pred["primary_l2"]["K"]
        + blocks0["primary_l3"]["K"]
    )
    g1.Lambda = g0.Lambda + _DT * (
        0.5 * blocks0["lambda_l2"]
        + 0.5 * blocks_pred["lambda_l2"]
        + blocks0["lambda_l3"]
    )
    validate_species_admissibility(
        metric_slice_from_q(grid, g1), m1.radiation,
        Species.RADIATION, stage="primary-predictor",
    )
    kernel._apply_outer_boundary_for_state(grid, g1, s1, m1, em0)
    g1.alpha = kernel._solve_lapse_for_state(
        grid, g1, s1, m1, em0
    )[0]
    g1.beta.fill(0.0)
    g1.B.fill(0.0)
    g1.assert_finite_positive()
    validate_species_admissibility(
        metric_slice_from_q(grid, g1), m1.radiation,
        Species.RADIATION, stage="primary-predictor",
    )

    primary_predictor = ProductionState(
        grid=grid,
        geometry=g1,
        scalars=s1,
        matter=m1,
        electromagnetic=em0,
        t=accepted.t,
        tau=accepted.tau,
        e_folds=accepted.e_folds,
    )
    blocks1, _, _, record1 = _observe_stage(
        kernel, native_grid, primary_predictor, evaluate_coupled_rhs=True
    )
    record1["stage"] = "primary_predictor_cmc"

    accepted_hash_after = _slice_sha256(
        accepted.geometry, accepted.scalars, accepted.matter
    )
    assert accepted_hash_after == accepted_hash_before, (
        "accepted production input changed while observing copied CMC stages"
    )

    source = inspect.getsource(V55ProductionKernel.step).encode("utf-8")
    report = {
        "schema": "native-cmc-stage-rhs-shadow/v9",
        "production_step_called": False,
        "stage_build": (
            "copied-state reconstruction using current production-owned lapse/RHS/"
            "boundary helpers and the merged source-level stage-order contract"
        ),
        "production_step_source_sha256": hashlib.sha256(source).hexdigest(),
        "resolution": int(grid.n),
        "r_max": float(grid.r_max),
        "dt": _DT,
        "native_and_production_grid_centres_match": True,
        "accepted_state_sha256_before": accepted_hash_before,
        "accepted_state_sha256_after": accepted_hash_after,
        "stage_observations": [record0, record_pred, record1],
        "reuse": [
            "PR #340 initialized frozen CMC RHS/source decomposition",
            "PR #351 CMC stage-order/boundary/reset source contract",
        ],
        "new_numeric_thresholds": False,
        "production_runtime_selector_added": False,
        "scientific_admission": False,
        "interpretation_limit": (
            "This compares a/b/X explicit RHS on three source-reconstructed CMC "
            "stage slices; it is not a full production-step capture, full-step "
            "equivalence result, or scientific-admission witness."
        ),
    }
    with capsys.disabled():
        print(
            "NATIVE_CMC_STAGE_RHS_SHADOW "
            + json.dumps(report, sort_keys=True, separators=(",", ":"), allow_nan=False),
            flush=True,
        )
