"""Regression for the archived GEAR-10 REV2 rung-6 drag-sign discriminator."""
from __future__ import annotations

import sys
from pathlib import Path

ARCHIVE_DIR = (
    Path(__file__).resolve().parents[1]
    / "archive"
    / "early_solvers"
    / "GEAR10"
)
if str(ARCHIVE_DIR) not in sys.path:
    sys.path.insert(0, str(ARCHIVE_DIR))

import validate_darksector_REV2 as validation  # noqa: E402


def test_rung6_identity_rejects_test_only_drag_sign_flip() -> None:
    """A wrong drag sign plateaus while the corrected integrated identity refines."""
    model = validation.solver.Model(beta=-0.2, lam=1.0, D=0.0, V0=2.07)
    sample_counts = (501, 2001, 8001)
    correct = validation._compute_rung6_rows(
        model, k=1.0, sample_counts=sample_counts, flip_drag_sign=False
    )
    faulty = validation._compute_rung6_rows(
        model, k=1.0, sample_counts=sample_counts, flip_drag_sign=True
    )

    correct_max = [row["normalized_residual_max_side"] for row in correct]
    faulty_max = [row["normalized_residual_max_side"] for row in faulty]
    assert correct_max[0] > correct_max[1] > correct_max[2], correct_max
    assert correct_max[-1] < correct_max[0] / 100.0, correct_max
    assert min(faulty_max) > 1.0e-3, faulty_max
    assert 0.5 < faulty_max[-1] / faulty_max[0] < 2.0, faulty_max

    for row in [*correct, *faulty]:
        # The sum-of-sides denominator is never larger than the max-side denominator.
        assert row["normalized_residual_max_side"] >= row["normalized_residual_sum_sides"]
        assert row["normalized_residual"] == row["normalized_residual_max_side"]
