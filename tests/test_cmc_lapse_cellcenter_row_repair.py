"""Regression for the cell-centred CMC lapse origin row.

The spatial grid stores cell centres at (i+1/2)dr. The CMC solve must use
its parity-aware derivative stencil at cell 0; alpha[0] == alpha[1] is not
the centre Neumann condition for this grid.
"""
import numpy as np

from tools.audit_cmc_lapse_cellcenter_row_repair import audit_one


def test_parity_aware_cmc_row_closes_same_stencil_as_production_rhs():
    row = audit_one(80)
    assert row["finite"]
    assert row["alpha_min"] > 0.0

    # PDE equations are enforced on cells 0..N-2. Cell N-1 is the alpha=1
    # normalization row and must not be counted as a PDE row.
    assert row["row_equation_residual"]["max_abs_all_pde_rows_i0_to_nminus2"] < 2.0e-11
    assert row["production_Kdot_minus_target"]["max_abs_all_pde_rows_i0_to_nminus2"] < 2.0e-10

    # A nonflat regular profile should not be artificially flattened across
    # the first two cell centres by the old alpha[0] == alpha[1] row.
    assert abs(row["alpha0_minus_alpha1"]) > 1.0e-10
    assert row["evolution_advanced"] is False
