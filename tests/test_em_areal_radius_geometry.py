"""Geometry edge-case contracts for the research-only radial-Maxwell state."""

from types import SimpleNamespace

import numpy as np
import pytest

from engine.radial_maxwell import RadialMaxwellState


def _nonmonotone_geometry(*, centers=(0.5, 1.5, 2.5)):
    # Positive metric samples whose face-interpolated areal radius first
    # decreases and then increases. This is intentional: R(r) need not be
    # monotone on every valid strong-field spherical slice.
    grid = SimpleNamespace(
        n=3,
        faces=np.array([0.0, 1.0, 2.0, 3.0]),
        centers=np.asarray(centers, dtype=float),
    )
    geometry = SimpleNamespace(
        b=np.array([9.0, 0.25, 1.0]),
        X=np.ones(3),
    )
    state = RadialMaxwellState(
        charge_to_rest_mass=1.0,
        enclosed_charge_faces=np.array([0.0, 1.0, 2.0, 3.0]),
    )
    return state, grid, geometry


def test_nonmonotone_positive_areal_radius_faces_and_gauss_field_are_valid():
    state, grid, geometry = _nonmonotone_geometry()

    radius_faces = state.areal_radius_faces(grid, geometry)
    expected_faces = np.array(
        [0.0, np.sqrt(4.625), 2.0 * np.sqrt(0.625), 3.0]
    )

    # The centre face is regular, all other face radii are positive, and the
    # non-monotonicity is preserved rather than clamped or rejected.
    np.testing.assert_allclose(radius_faces, expected_faces, rtol=0.0, atol=1e-15)
    assert radius_faces[0] == 0.0
    assert np.all(radius_faces[1:] > 0.0)
    assert radius_faces[1] > radius_faces[2]

    radius_centres = grid.centers * np.sqrt(geometry.b) / geometry.X
    assert radius_centres[0] > radius_centres[1] > 0.0
    q_centres = 0.5 * (
        state.enclosed_charge_faces[:-1] + state.enclosed_charge_faces[1:]
    )
    field = state.electric_field_hat(grid, geometry)
    np.testing.assert_allclose(
        field, q_centres / np.square(radius_centres), rtol=0.0, atol=1e-15
    )
    assert np.all(np.isfinite(field))


def test_true_zero_cell_centred_areal_radius_fails_before_field_division():
    # Unlike the regular origin *face*, a zero-radius cell centre would make
    # Q/R**2 singular. The reconstruction must fail closed rather than emit
    # inf/NaN or silently regularize the input.
    state, grid, geometry = _nonmonotone_geometry(centers=(0.0, 1.5, 2.5))

    with pytest.raises(
        ValueError, match="cell-centred areal radius must be positive"
    ):
        state.electric_field_hat(grid, geometry)


@pytest.mark.parametrize(
    "metric",
    [
        {"b": np.array([9.0, 0.0, 1.0]), "X": np.ones(3)},
        {"b": np.array([9.0, 0.25, 1.0]), "X": np.array([1.0, -1.0, 1.0])},
    ],
)
def test_areal_radius_face_reconstruction_rejects_nonpositive_metric(metric):
    state, grid, _ = _nonmonotone_geometry()
    geometry = SimpleNamespace(**metric)

    with pytest.raises(ValueError, match="strictly positive"):
        state.areal_radius_faces(grid, geometry)
