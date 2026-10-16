"""Analytic controls for the pinned production HLL/source pairing observer."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))
SPEC = importlib.util.spec_from_file_location(
    "_issue4_production_hll_operator_pairing",
    ROOT / "tools" / "audit_issue4_production_hll_operator_pairing.py",
)
if SPEC is None or SPEC.loader is None:
    raise ImportError("cannot load production HLL operator-pairing diagnostic")
DIAG = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = DIAG
SPEC.loader.exec_module(DIAG)

from engine.matter_system import (  # noqa: E402
    BSSNMetricSlice,
    Species,
    _metric_arrays,
    _species_rhs_components,
    geometry_metric_derivatives,
    initialize_radiation,
)


def test_flat_uniform_radiation_hll_momentum_flux_matches_exact_face_flux():
    n = 48
    dr = 0.25
    centers = (np.arange(n, dtype=np.float64) + 0.5) * dr
    ones = np.ones(n, dtype=np.float64)
    zeros = np.zeros(n, dtype=np.float64)
    metric = BSSNMetricSlice(
        r=centers,
        a=ones.copy(),
        b=ones.copy(),
        X=ones.copy(),
        alpha=ones.copy(),
        beta=zeros.copy(),
        Aa=zeros.copy(),
        K=zeros.copy(),
        Lambda=zeros.copy(),
        B=zeros.copy(),
    )
    metrics = _metric_arrays(metric)
    energy = 0.6
    radiation = initialize_radiation(metrics, np.full(n, energy, dtype=np.float64))
    # a,b,X,alpha,beta are constant; the production helper supplies the exact
    # spherical derivative of gamma_theta-theta through its explicit 2/r term.
    metric_derivatives = geometry_metric_derivatives(
        metric, lambda values, parity: np.zeros_like(values, dtype=np.float64)
    )
    face_flux, _, _, _, _ = _species_rhs_components(
        metric,
        metric_derivatives,
        radiation,
        Species.RADIATION,
        dphi_t=zeros,
        dphi_r=zeros,
        beta_dm=0.0,
    )
    faces = np.arange(n + 1, dtype=np.float64) * dr
    exact = (energy / 3.0) * faces**2
    assert np.allclose(face_flux[:, 2], exact, rtol=3.0e-13, atol=2.0e-14)
    # A reversed-flux control must be rejected by the exact analytic fixture.
    correct_gap = float(np.sqrt(np.mean((face_flux[:, 2] - exact) ** 2)))
    wrong_sign_gap = float(np.sqrt(np.mean((-face_flux[:, 2] - exact) ** 2)))
    assert correct_gap < 2.0e-14
    assert wrong_sign_gap > 1.0e-3


def test_independent_geometric_source_matches_pinned_covariant_contraction():
    kernel = DIAG.AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=48, r_max=24.0, include_radiation=True)
    n = state.grid.n
    # A smooth manufactured non-unit lapse and regular even scalar profile
    # ensure the conditional Q_r source test is not vacuous. This is a copied
    # analytic/control state, not a production trajectory or step.
    r = np.asarray(state.grid.centers, dtype=np.float64)
    state.geometry.alpha[:] = 0.75 + 0.20 * np.exp(-(r / 6.0) ** 2)
    state.scalars.phi[:] = 0.01 + 0.001 * r**2
    rates = {"scalars": {"phi": np.zeros(n, dtype=np.float64)}}
    snapshot = DIAG.operator_snapshot(state, -0.04, rates)
    assert snapshot["max_abs_geometric_source_gap_over_kappa"] < 2.0e-14


def test_conditional_qr_source_discrepancy_is_recorded_not_silently_relabelled():
    kernel = DIAG.AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=48, r_max=24.0, include_radiation=True)
    n = state.grid.n
    r = np.asarray(state.grid.centers, dtype=np.float64)
    alpha = 0.75 + 0.20 * np.exp(-(r / 6.0) ** 2)
    state.geometry.alpha[:] = alpha
    state.scalars.phi[:] = 0.01 + 0.001 * r**2
    rates = {"scalars": {"phi": np.zeros(n, dtype=np.float64)}}
    snapshot = DIAG.operator_snapshot(state, -0.04, rates)
    raw = snapshot["raw_arrays"]["dark_matter"]
    actual = np.asarray(raw["interaction_momentum_source_actual_over_kappa"])
    expected = np.asarray(
        raw["interaction_momentum_source_conditional_covariant_over_kappa"]
    )
    gap = np.asarray(raw["interaction_momentum_source_gap_over_kappa"])
    # Under the explicitly conditional covector premise, the required source
    # is sqrt(-g) Q_r. The pinned implementation's measured q_s is compared
    # against it and must retain the known alpha-weighting discrepancy.
    assert np.allclose(gap, actual * (1.0 - alpha), rtol=2.0e-12, atol=2.0e-20)
    assert np.max(np.abs(gap)) > 1.0e-12
    assert np.max(np.abs(expected)) > 1.0e-12
