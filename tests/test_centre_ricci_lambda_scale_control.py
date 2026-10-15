from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from analyze_centre_ricci_lambda_scale_control import analyze_profile, analyze_scales


def test_lambda_radial_scale_control_reports_signed_centre_and_fixed_annulus_witnesses():
    report = analyze_scales()
    assert report["schema"] == "centre_ricci_lambda_radial_scale_control_v1"
    assert report["resolutions"] == [80, 160, 320]
    assert report["length_scales"] == [4.0, 8.0, 16.0]
    assert report["fixed_annulus"] == [1.0, 4.0]

    for scale in report["scales"]:
        rows = scale["resolutions"]
        assert [row["resolution"] for row in rows] == [80, 160, 320]
        for row in rows:
            d4 = row["D4_sampled_odd_Lambda_minus_exact"]
            chain = row["D1_chain_rule_minus_exact"]
            gap = row["D4_minus_chain_rule"]
            for channel in (d4, chain, gap):
                assert np.isfinite(channel["cell0"]["signed_error"])
                assert np.isfinite(channel["cell0"]["error_over_r2"])
                assert np.isfinite(channel["cell1"]["signed_error"])
                assert np.isfinite(channel["cell1"]["error_over_r2"])
                assert channel["fixed_annulus_max_abs"] > 0.0
            assert row["channel_response"]["full_three_term_recombination_max_abs"] < 2e-12
            assert row["channel_response"]["radial_term_identity_max_abs"] < 2e-12


def test_both_derivative_representations_converge_at_fixed_physical_coordinates():
    report = analyze_scales()
    for scale in report["scales"]:
        for channel in (
            "D4_sampled_odd_Lambda_minus_exact",
            "D1_chain_rule_minus_exact",
            "D4_minus_chain_rule",
        ):
            orders = scale["fixed_annulus_observed_orders"][channel]
            # The first scale is the most demanding on this grid; retain margin
            # below the nominal fourth-order values without hiding a failed slope.
            assert len(orders) == 2
            assert min(row["observed_order"] for row in orders) > 3.7
            assert all(row["error_ratio"] > 1.0 for row in orders)


def test_longer_length_scales_reduce_absolute_error_without_changing_domain():
    report = analyze_scales()
    for n_index in range(3):
        for channel in (
            "D4_sampled_odd_Lambda_minus_exact",
            "D1_chain_rule_minus_exact",
            "D4_minus_chain_rule",
        ):
            errors = [
                scale["resolutions"][n_index][channel]["fixed_annulus_max_abs"]
                for scale in report["scales"]
            ]
            assert errors[0] > errors[1] > errors[2]


def test_determinant_one_manufactured_metric_identities_are_exact():
    row = analyze_profile(160, 8.0)
    assert row["all_profiles_finite"]
    # These identities follow from a=q^(-4/3), b=q^(2/3), X=q^(-2/3).
    # Unit test exercises their numerically evaluated form on the fixed grid.
    r = (np.arange(160, dtype=float) + 0.5) * (32.0 / 160)
    q = 1.0 + 0.1 * (1.0 - np.exp(-((r / 8.0) ** 2)))
    a, b, X = q ** (-4.0 / 3.0), q ** (2.0 / 3.0), q ** (-2.0 / 3.0)
    np.testing.assert_allclose(a * b * b, 1.0, rtol=2e-15, atol=2e-15)
    np.testing.assert_allclose(np.sqrt(a) / X, 1.0, rtol=2e-15, atol=2e-15)
    np.testing.assert_allclose(r * np.sqrt(b) / X, r * q, rtol=2e-15, atol=2e-15)
