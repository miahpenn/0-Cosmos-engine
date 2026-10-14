"""Manufactured discrete contracts for EM current/boundary exchange.

These tests validate the algebra of the selected face-charge energy functional.
They do not select a physical carrier, current law, or boundary convention.
Physical campaign residuals must remain unthresholded.
"""
import numpy as np
import pytest


def _roundoff_bound(*operands):
    """64-epsilon bound scaled only by the actual algebraic operands."""
    scale = sum(abs(float(value)) for value in operands)
    return 64.0 * np.finfo(np.float64).eps * max(
        scale, np.finfo(np.float64).tiny
    )


def _manufactured_terms(charge_faces, current_faces, weight, dr):
    """Independently differentiate U=dr/2 sum_i w_i q_i**2.

    q_i is the arithmetic mean of neighboring face-enclosed charges and
    Qdot_face=-I_face. The interior-face and outer-face contributions are
    separated only after each discrete term has been formed.
    """
    Q = np.asarray(charge_faces, dtype=np.float64)
    I = np.asarray(current_faces, dtype=np.float64)
    w = np.asarray(weight, dtype=np.float64)
    assert Q.ndim == I.ndim == w.ndim == 1
    assert Q.size == I.size == w.size + 1
    assert np.all(np.isfinite(Q))
    assert np.all(np.isfinite(I))
    assert np.all(np.isfinite(w))
    assert float(dr) > 0.0
    assert I[0] == 0.0  # regular-centre current

    q_cell = 0.5 * (Q[:-1] + Q[1:])
    qdot_cell = -0.5 * (I[:-1] + I[1:])
    full_charge_tangent = float(dr) * float(
        np.sum(w * q_cell * qdot_cell)
    )

    # Face j between cells j-1 and j contributes to both neighboring
    # cell energies. The right boundary contributes to the last cell only.
    interior_power = -0.5 * float(dr) * float(np.sum(
        (w[:-1] * q_cell[:-1] + w[1:] * q_cell[1:]) * I[1:-1]
    ))
    outer_power = (
        -0.5 * float(dr) * float(w[-1] * q_cell[-1] * I[-1])
    )

    # The full cell-adjoint recipient projection uses the same averaged
    # face current as qdot. The interior-only recipient projection excludes
    # the outer face and therefore pairs with an explicit boundary term.
    full_cell_recipient = float(dr) * float(np.sum(
        w * q_cell * 0.5 * (I[:-1] + I[1:])
    ))
    interior_recipient = -interior_power

    return {
        "q_cell": q_cell,
        "full_charge_tangent": full_charge_tangent,
        "interior_power": interior_power,
        "outer_power": outer_power,
        "full_cell_recipient": full_cell_recipient,
        "interior_recipient": interior_recipient,
    }


@pytest.mark.parametrize(
    "current_faces",
    [
        np.array([0.0, 0.2, -0.15, 0.07, 0.0]),   # closed outer face
        np.array([0.0, 0.0, 0.0, 0.0, 0.31]),     # outer face only
        np.array([0.0, 0.2, -0.15, 0.07, 0.31]),  # mixed interior + boundary
    ],
)
def test_manufactured_field_tangent_splits_into_interior_and_outer_terms(
    current_faces,
):
    charge_faces = np.array([0.0, 0.8, 1.7, 2.0, 2.4])
    weight = np.array([1.4, 0.9, 0.7, 0.45])
    dr = 0.25
    terms = _manufactured_terms(charge_faces, current_faces, weight, dr)

    tangent = terms["full_charge_tangent"]
    interior = terms["interior_power"]
    outer = terms["outer_power"]
    # This is an algebraic identity, so only an operand-scaled roundoff
    # allowance is permitted.
    assert abs(tangent - (interior + outer)) <= _roundoff_bound(
        tangent, interior, outer
    )


@pytest.mark.parametrize(
    "current_faces",
    [
        np.array([0.0, 0.2, -0.15, 0.07, 0.0]),
        np.array([0.0, 0.0, 0.0, 0.0, 0.31]),
        np.array([0.0, 0.2, -0.15, 0.07, 0.31]),
    ],
)
def test_manufactured_donor_recipient_closes_under_each_consistent_convention(
    current_faces,
):
    charge_faces = np.array([0.0, 0.8, 1.7, 2.0, 2.4])
    weight = np.array([1.4, 0.9, 0.7, 0.45])
    dr = 0.25
    terms = _manufactured_terms(charge_faces, current_faces, weight, dr)

    tangent = terms["full_charge_tangent"]
    outer = terms["outer_power"]
    full_cell_recipient = terms["full_cell_recipient"]
    interior_recipient = terms["interior_recipient"]

    # Convention A: the cell recipient uses the full face-adjoint current,
    # so the full functional tangent is paired directly with it.
    full_balance = tangent + full_cell_recipient
    assert abs(full_balance) <= _roundoff_bound(
        tangent, full_cell_recipient
    )

    # Convention B: recipient work represents interior-face exchange only;
    # the outer-face term is then present explicitly, exactly once.
    interior_balance = tangent + interior_recipient - outer
    assert abs(interior_balance) <= _roundoff_bound(
        tangent, interior_recipient, outer
    )


def test_manufactured_case_exposes_mixing_full_cell_work_with_separate_boundary():
    """A boundary term must not be added to a full-cell balance by default."""
    charge_faces = np.array([0.0, 0.8, 1.7, 2.0, 2.4])
    current_faces = np.array([0.0, 0.0, 0.0, 0.0, 0.31])
    weight = np.array([1.4, 0.9, 0.7, 0.45])
    dr = 0.25
    terms = _manufactured_terms(charge_faces, current_faces, weight, dr)

    # With boundary current alone, the full charge tangent is the outer term.
    assert abs(
        terms["full_charge_tangent"] - terms["outer_power"]
    ) <= _roundoff_bound(
        terms["full_charge_tangent"], terms["outer_power"]
    )

    # The full cell-adjoint recipient already cancels that full tangent.
    # Subtracting the boundary term again creates exactly -P_out, not zero.
    mixed_convention_residual = (
        terms["full_charge_tangent"]
        + terms["full_cell_recipient"]
        - terms["outer_power"]
    )
    expected = -terms["outer_power"]
    assert abs(mixed_convention_residual - expected) <= _roundoff_bound(
        mixed_convention_residual, expected
    )
    assert abs(expected) > 0.0
