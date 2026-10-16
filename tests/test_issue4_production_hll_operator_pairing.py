"""Analytic controls for the pinned production HLL/source pairing observer."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pytest
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


def test_optional_charged_dust_and_em_lanes_fail_closed():
    kernel = DIAG.AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=24, r_max=24.0, include_radiation=True)
    n = len(state.grid.centers)
    rates = {"scalars": {"phi": np.zeros(n, dtype=np.float64)}}

    # These sentinel objects deliberately exercise the scope guard before
    # projection or helper capture; they are not an EM physics test.
    state.matter.charged_dust = object()
    with pytest.raises(ValueError, match="optional charged-dust/Maxwell"):
        DIAG.operator_snapshot(state, -0.04, rates)

    state.matter.charged_dust = None
    state.electromagnetic = object()
    with pytest.raises(ValueError, match="optional charged-dust/Maxwell"):
        DIAG.operator_snapshot(state, -0.04, rates)


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
    assert face_flux.shape == (n + 1, 3)
    assert np.array_equal(face_flux[0], np.zeros(3, dtype=np.float64))
    assert np.allclose(face_flux[:, 2], exact, rtol=3.0e-13, atol=2.0e-14)
    # A reversed-flux control must be rejected by the exact analytic fixture.
    correct_gap = float(np.sqrt(np.mean((face_flux[:, 2] - exact) ** 2)))
    wrong_sign_gap = float(np.sqrt(np.mean((-face_flux[:, 2] - exact) ** 2)))
    assert correct_gap < 2.0e-14
    assert wrong_sign_gap > 1.0e-3


def test_independent_geometric_source_matches_pinned_covariant_contraction():
    kernel = DIAG.AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=48, r_max=24.0, include_radiation=True)
    n = len(state.grid.centers)
    # A smooth manufactured non-unit lapse and regular even scalar profile
    # ensure the conditional Q_r source test is not vacuous. This is a copied
    # analytic/control state, not a production trajectory or step.
    r = np.asarray(state.grid.centers, dtype=np.float64)
    state.geometry.alpha[:] = 0.75 + 0.20 * np.exp(-(r / 6.0) ** 2)
    state.scalars.phi[:] = 0.01 + 0.001 * r**2
    rates = {"scalars": {"phi": np.zeros(n, dtype=np.float64)}}
    snapshot = DIAG.operator_snapshot(state, -0.04, rates)
    assert snapshot["max_abs_geometric_source_gap_over_kappa"] < 2.0e-14
    assert snapshot["origin_face_flux_exact_zero"]
    assert snapshot["native_component_arrays_preserved"]
    raw_radiation = snapshot["raw_arrays"]["radiation"]
    native_flux = np.asarray(raw_radiation["native_hll_face_flux"], dtype=np.float64)
    native_source_s = np.asarray(
        raw_radiation["native_geometric_momentum_source"], dtype=np.float64
    )
    assert native_flux.shape == (n + 1, 3)
    assert np.array_equal(
        native_flux[:, 2] / DIAG.KAPPA,
        np.asarray(raw_radiation["hll_face_flux_momentum_over_kappa"]),
    )
    assert np.array_equal(
        native_source_s / DIAG.KAPPA,
        np.asarray(raw_radiation["geometric_momentum_source_actual_over_kappa"]),
    )
    assert all(
        np.all(np.isfinite(np.asarray(raw_radiation[key], dtype=np.float64)))
        for key in (
            "native_hll_face_flux",
            "native_geometric_energy_source",
            "native_geometric_momentum_source",
            "native_interaction_energy_source",
            "native_interaction_momentum_source",
        )
    )


def test_conditional_qr_source_discrepancy_is_recorded_not_silently_relabelled():
    kernel = DIAG.AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=48, r_max=24.0, include_radiation=True)
    n = len(state.grid.centers)
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


def test_manufactured_probe_rejects_sign_flipped_flux_and_omitted_sources():
    kernel = DIAG.AUDIT.V55ProductionKernel()
    state = DIAG.make_manufactured_operator_probe(kernel, n=48, r_max=24.0)
    n = len(state.grid.centers)
    snapshot = DIAG.operator_snapshot(
        state, -0.04, {"scalars": {"phi": np.zeros(n, dtype=np.float64)}}
    )
    assert snapshot["origin_face_flux_exact_zero"]
    assert snapshot["native_component_arrays_preserved"]
    dm = snapshot["species"]["dark_matter"]
    controls = dm["negative_controls"]

    assert controls["sign_flipped_hll_minus_projected_momentum_flux"]["max_abs"] > (
        dm["hll_momentum_flux_minus_projected_stress_face_flux"]["max_abs"]
    )
    assert controls["sign_flipped_hll_minus_projected_momentum_flux"]["max_abs"] > 1.0e-12
    assert controls["omitted_geometric_momentum_source_residual"]["max_abs"] > 1.0e-12
    assert controls["omitted_interaction_momentum_source_residual"]["max_abs"] > 1.0e-12

    assert set(snapshot["scalar_sector_projections"]) == {"S", "D", "phi"}
    assert set(snapshot["raw_scalar_arrays"]["phi"]) == {"rho", "pr", "pt", "j"}
    assert snapshot["raw_scalar_arrays"]["phi"]["rho"]


def test_combined_discrete_budget_matches_native_flux_source_signs_and_controls():
    kernel = DIAG.AUDIT.V55ProductionKernel()
    state = DIAG.make_manufactured_operator_probe(kernel, n=64, r_max=24.0)
    n = len(state.grid.centers)
    snapshot = DIAG.operator_snapshot(
        state, -0.04, {"scalars": {"phi": np.zeros(n, dtype=np.float64)}}
    )
    dr = float(state.grid.dr)

    for species_name, raw in snapshot["raw_arrays"].items():
        rest = -np.diff(np.asarray(raw["hll_face_flux_rest_over_kappa"], dtype=np.float64)) / dr
        energy = (
            -np.diff(np.asarray(raw["hll_face_flux_energy_over_kappa"], dtype=np.float64)) / dr
            + np.asarray(raw["geometric_energy_source_actual_over_kappa"], dtype=np.float64)
            + np.asarray(raw["interaction_energy_source_actual_over_kappa"], dtype=np.float64)
        )
        momentum = (
            -np.diff(np.asarray(raw["hll_face_flux_momentum_over_kappa"], dtype=np.float64)) / dr
            + np.asarray(raw["geometric_momentum_source_actual_over_kappa"], dtype=np.float64)
            + np.asarray(raw["interaction_momentum_source_actual_over_kappa"], dtype=np.float64)
        )
        np.testing.assert_allclose(raw["production_rhs_rest_over_kappa"], rest, rtol=2e-13, atol=2e-14)
        np.testing.assert_allclose(raw["production_rhs_energy_over_kappa"], energy, rtol=2e-13, atol=2e-14)
        np.testing.assert_allclose(raw["production_rhs_momentum_over_kappa"], momentum, rtol=2e-13, atol=2e-14)

        flipped = (
            np.diff(np.asarray(raw["hll_face_flux_momentum_over_kappa"], dtype=np.float64)) / dr
            + np.asarray(raw["geometric_momentum_source_actual_over_kappa"], dtype=np.float64)
            + np.asarray(raw["interaction_momentum_source_actual_over_kappa"], dtype=np.float64)
        )
        omit_geometric = (
            -np.diff(np.asarray(raw["hll_face_flux_momentum_over_kappa"], dtype=np.float64)) / dr
            + np.asarray(raw["interaction_momentum_source_actual_over_kappa"], dtype=np.float64)
        )
        omit_interaction = (
            -np.diff(np.asarray(raw["hll_face_flux_momentum_over_kappa"], dtype=np.float64)) / dr
            + np.asarray(raw["geometric_momentum_source_actual_over_kappa"], dtype=np.float64)
        )
        np.testing.assert_allclose(raw["production_rhs_momentum_sign_flipped_flux_control_over_kappa"], flipped, rtol=2e-13, atol=2e-14)
        np.testing.assert_allclose(raw["production_rhs_momentum_omitted_geometric_source_control_over_kappa"], omit_geometric, rtol=2e-13, atol=2e-14)
        np.testing.assert_allclose(raw["production_rhs_momentum_omitted_interaction_source_control_over_kappa"], omit_interaction, rtol=2e-13, atol=2e-14)

        report = snapshot["species"][species_name]["production_rhs_budget"]
        for channel in ("rest", "energy", "radial_momentum"):
            assert report[channel]["cell_0"]["index"] == 0
            assert report[channel]["cell_1"]["index"] == 1
            assert report[channel]["annulus_2_to_10"]["cell_count"] > 0
            assert "signed_rms" in report[channel]["annulus_2_to_10"]
        assert report["negative_controls"]["radial_momentum_sign_flipped_flux"]["cell_0"]["index"] == 0
        assert report["negative_controls"]["radial_momentum_omitted_geometric_source"]["cell_1"]["index"] == 1
        assert report["negative_controls"]["radial_momentum_omitted_interaction_source"]["annulus_2_to_10"]["cell_count"] > 0


def test_base_state_total_t_audit_keeps_admission_closed():
    kernel = DIAG.AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=24, r_max=24.0, include_radiation=True)
    initial_hash = DIAG.AUDIT.state_digest(state)
    epsilons = (2.0e-5, 1.0e-5)

    witness = DIAG.matched_total_t_witness(
        kernel, state, -0.04, epsilons, initial_hash
    )

    contract = witness["matched_state_contract"]
    assert contract["same_in_memory_base_state_as_hll_operator_snapshot"] is True
    assert contract["initial_state_sha256"] == initial_hash
    assert contract["state_digest_unchanged_after_total_t_audit"] is True
    assert contract["matched_epsilon_values"] == list(epsilons)
    assert contract["all_total_t_directional_arrays_finite"] is True
    assert contract["production_steps_advanced"] == 0
    assert contract["conditional_exchange_law_accepted"] is False
    assert contract["independent_total_stress_energy_divergence_closed"] is False
    assert contract["scientific_admission"] is False
    assert len(witness["epsilon_rows"]) == len(epsilons)
    assert len(witness["raw_arrays"]["epsilon_0_div_T_t"]) == len(state.grid.centers)
    assert len(witness["raw_arrays"]["epsilon_1_div_T_r_fv"]) == len(state.grid.centers)
    assert witness["conditional_rest_current_source"]["model_law_accepted"] is False
    assert DIAG.AUDIT.state_digest(state) == initial_hash


def test_hll_and_total_t_witness_share_exact_epsilon_state_pair():
    row = DIAG.run_resolution(
        n=24, r_max=24.0, epsilons=(2.0e-5, 1.0e-5)
    )
    assert row["production_steps_advanced"] == 0
    assert row["base_state_unchanged"] is True
    all_reconstruction_contracts_passed = True
    all_reconstruction_diagnostics = []
    for coupling_name, coupling_witnesses in row["directional_epsilon_witnesses"].items():
        assert set(coupling_witnesses) == {"epsilon_2e-05", "epsilon_1e-05"}
        for witness in coupling_witnesses.values():
            paired = witness["same_state_total_T_directional_witness"]
            assert paired["state_pairing_verified"] is True
            assert paired["base_state_sha256"] == row["initial_state_sha256"]
            assert paired["plus_state_sha256"] == witness["plus_state_sha256"]
            assert paired["minus_state_sha256"] == witness["minus_state_sha256"]
            assert paired["epsilon"] == witness["epsilon"]
            assert paired["projection_crosscheck_passed"] is True
            decomposition = paired["sector_decomposition"]
            expected_sectors = {"dark_matter", "baryons", "radiation", "S", "D", "phi"}
            assert set(decomposition["sectors"]) == expected_sectors
            assert decomposition["state_pairing"]["all_sector_arrays_finite"] is True
            reconstruction_debug = {}
            for channel_name, channel in decomposition["reconstruction"]["channels"].items():
                sum_values = np.asarray(channel["sum_of_sectors"], dtype=np.float64)
                total_values = np.asarray(paired["raw_signed_arrays"][channel_name], dtype=np.float64)
                reconstruction_debug[channel_name] = {
                    "passed": channel["passed"],
                    "max_abs_gap": channel["max_abs_gap"],
                    "max_abs_sum": float(np.max(np.abs(sum_values))),
                    "max_abs_total": float(np.max(np.abs(total_values))),
                    "atol": channel["atol"],
                    "rtol": channel["rtol"],
                    "max_gap_allowed_by_global_scale": channel["atol"] + channel["rtol"] * max(
                        float(np.max(np.abs(sum_values))), float(np.max(np.abs(total_values)))
                    ),
                    "max_abs_sector_l1_scale": float(np.max(np.sum(np.stack([
                        np.abs(np.asarray(sector["raw_signed_arrays"][channel_name], dtype=np.float64))
                        for sector in decomposition["sectors"].values()
                    ], axis=0), axis=0))),
                    "max_gap_over_sector_l1": float(np.max(
                        np.abs(np.asarray(channel["sum_minus_total"], dtype=np.float64)) /
                        np.maximum(
                            np.sum(np.stack([
                                np.abs(np.asarray(sector["raw_signed_arrays"][channel_name], dtype=np.float64))
                                for sector in decomposition["sectors"].values()
                            ], axis=0), axis=0),
                            np.finfo(np.float64).tiny,
                        )
                    )),
                    "max_gap_index": int(np.argmax(np.abs(np.asarray(channel["sum_minus_total"], dtype=np.float64)))),
                    "value_at_max_gap": {
                        "sum_of_sectors": float(sum_values[int(np.argmax(np.abs(np.asarray(channel["sum_minus_total"], dtype=np.float64))))]),
                        "total": float(total_values[int(np.argmax(np.abs(np.asarray(channel["sum_minus_total"], dtype=np.float64))))]),
                        "gap": float(np.asarray(channel["sum_minus_total"], dtype=np.float64)[int(np.argmax(np.abs(np.asarray(channel["sum_minus_total"], dtype=np.float64))))]),
                        "sector_terms": {
                            name: float(np.asarray(sector["raw_signed_arrays"][channel_name], dtype=np.float64)[int(np.argmax(np.abs(np.asarray(channel["sum_minus_total"], dtype=np.float64))))])
                            for name, sector in decomposition["sectors"].items()
                        },
                    },
                }
            print("SECTOR_RECON_DEBUG=", reconstruction_debug, flush=True)
            reconstruction_contract = bool(
                decomposition["reconstruction"]["all_sectors_reconstruct_total"]
            )
            all_reconstruction_contracts_passed &= reconstruction_contract
            all_reconstruction_diagnostics.append({
                "coupling_arm": coupling_name,
                "epsilon": witness["epsilon"],
                "reconstruction_contract": reconstruction_contract,
                "channels": reconstruction_debug,
            })
            for key in ("div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv"):
                values = np.asarray(paired["raw_signed_arrays"][key], dtype=np.float64)
                assert values.shape == (24,)
                assert np.all(np.isfinite(values))
                summary = paired[f"{key}_summary"]
                assert "cell_0" in summary and "cell_1" in summary
                assert "annulus_2_to_10_rms" in summary
                rec = decomposition["reconstruction"]["channels"][key]
                all_reconstruction_contracts_passed &= bool(rec["passed"])
                all_reconstruction_contracts_passed &= bool(rec["max_abs_gap"] <= rec["atol"] + rec["rtol"] * max(
                    np.max(np.abs(np.asarray(rec["sum_of_sectors"], dtype=np.float64))),
                    np.max(np.abs(values)),
                ))
                np.testing.assert_allclose(
                    np.asarray(rec["sum_of_sectors"], dtype=np.float64),
                    values, rtol=rec["rtol"], atol=rec["atol"],
                )
                np.testing.assert_allclose(
                    np.asarray(rec["sum_minus_total"], dtype=np.float64),
                    np.asarray(rec["sum_of_sectors"], dtype=np.float64) - values,
                    rtol=0.0, atol=5.0e-25,
                )
                for sector in decomposition["sectors"].values():
                    sector_values = np.asarray(sector["raw_signed_arrays"][key], dtype=np.float64)
                    assert sector_values.shape == (24,)
                    assert np.all(np.isfinite(sector_values))
                    sector_summary = sector["summaries"][key]
                    assert "cell_0" in sector_summary and "cell_1" in sector_summary
                    assert "annulus_2_to_10_rms" in sector_summary
    print("SECTOR_RECON_ALL_DEBUG=", all_reconstruction_diagnostics, flush=True)
    assert all_reconstruction_contracts_passed, all_reconstruction_diagnostics
    assert row["matched_total_t_state_hashes_match"] is True
