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


def test_authorized_runtime_source_exception_is_exact_and_narrow():
    assert live.AUTHORIZED_RUNTIME_HEAD_BLOBS == {
        "engine/production_kernel.py": "cc4ae0f58dd5df4e05f596f6e677526586956f26",
        "engine/v55_pirk_adapter.py": "be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3",
    }
    assert live.AUTHORIZED_LANE_B_HEAD_BLOBS == {
        "engine/matter_system.py": "0969f4d996e4f8e0fc5b70a7dff147f024907302",
        "engine/matter_rhs.py": "0153c232ee15b5a9d3dbac37ea858e4d2d730011",
        "engine/scalar_system.py": "e2c2dae9a016b4d6454cf74a9eae8b0167d5733d",
        "engine/production_kernel.py": "1720615628a98a8d7cfc9fcadcc984fe28c88be1",
    }
    assert live.AUTHORIZED_NATIVE_GEOMETRY_TRANSITIONS == {
        "engine/native_geometry/geometry.py": {
            "historical_blob": "9e916b562a137ada0a8ba0268a5d2140162654d3",
            "approved_head_blob": "1bf12c1139572ff4597eaec0570aa43f2563fd91",
        },
    }
    assert live._native_geometry_transition_matches(
        "engine/native_geometry/geometry.py",
        "9e916b562a137ada0a8ba0268a5d2140162654d3",
        "1bf12c1139572ff4597eaec0570aa43f2563fd91",
    )
    assert not live._native_geometry_transition_matches(
        "engine/native_geometry/geometry.py",
        "f" * 40,
        "1bf12c1139572ff4597eaec0570aa43f2563fd91",
    )
    assert not live._native_geometry_transition_matches(
        "engine/native_geometry/geometry.py",
        "9e916b562a137ada0a8ba0268a5d2140162654d3",
        "f" * 40,
    )
    for path, approved_head in live.AUTHORIZED_RUNTIME_HEAD_BLOBS.items():
        expected = live.CURRENT_SOURCE_BLOBS[path]
        assert live._source_blob_matches(path, expected, expected, expected)
        assert live._source_blob_matches(path, expected, expected, approved_head)
        assert not live._source_blob_matches(path, expected, expected, "f" * 40)
        assert not live._source_blob_matches(path, expected, "f" * 40, approved_head)
    assert live._unexpected_runtime_source_paths(
        ["engine/production_kernel.py", "engine/v55_pirk_adapter.py",
         "engine/native_geometry/geometry.py"]
    ) == []
    assert live._unexpected_runtime_source_paths(
        ["engine/production_kernel.py", "engine/matter_system.py", "engine/matter_rhs.py", "engine/scalar_system.py"]
    ) == []
    assert live._unexpected_runtime_source_paths(["engine/production_kernel.py", "engine/unapproved.py"]) == ["engine/unapproved.py"]


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
