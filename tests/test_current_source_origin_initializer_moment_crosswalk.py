"""Static integrity controls for the current-source origin audit wrapper."""
from tools import audit_current_main_centre_proper_radial_ricci as archived
from tools import audit_current_source_origin_initializer_moment_crosswalk as live


def test_current_source_pins_are_explicit_and_distinct_from_archive():
    assert live.SOURCE_COMMIT == "ea2dabc17426b714eb445871c3c6b880638eee38"
    assert live.VENDOR_COMMIT == "d6052d605673ce9d82cdc99b0df79855fa2ea215"
    assert live.PINNED_OBSERVER_BLOB == "da485d220664269cd52563728acfc9441678cda0"
    assert set(live.CURRENT_SOURCE_BLOBS) == set(archived.SOURCE_BLOBS)
    assert live.CURRENT_SOURCE_BLOBS["engine/matter_system.py"] == (
        "e52a117690253ddad7b24d116145a65474a351d2"
    )
    assert archived.SOURCE_COMMIT == "285a2ee0c1a3b97cddeb4519d8e749109b115aab"
    assert live.SOURCE_COMMIT != archived.SOURCE_COMMIT


def test_crosswalk_is_one_no_evolution_grid_and_no_exact_moment_substitution():
    assert live.RESOLUTIONS == (80,)
    assert live.R_MAX == 40.0
    assert live.OUTPUT_DIR.name == (
        "current-source-origin-initializer-moment-crosswalk-2026-10-17"
    )


def test_delta_over_r2_uses_the_declared_metric_delta():
    import numpy as np

    r = np.array([0.25, 0.75, 1.25], dtype=float)
    a = 1.0 + 0.02 * r**2
    b = 1.0 - 0.01 * r**2
    X = np.ones_like(r)
    Lambda = 0.003 * r
    source = 1.0 + 0.2 * r**2
    B_reconstructed = np.ones_like(r)
    B_from_X6 = np.ones_like(r)
    rows = live.first_six_metric_profiles(
        r, a, b, X, Lambda, source, B_reconstructed, B_from_X6
    )

    expected = (1.0 / b - 1.0 / a) / (r * r)
    np.testing.assert_allclose(
        [row["delta_over_r2"] for row in rows],
        expected,
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        [row["Lambda_over_r"] for row in rows],
        np.full_like(r, 0.003),
        rtol=0.0,
        atol=1.0e-18,
    )
