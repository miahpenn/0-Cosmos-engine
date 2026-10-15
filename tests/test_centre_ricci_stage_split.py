import numpy as np

from engine.production_kernel import V55ProductionKernel
from engine.v55_pirk_adapter import vendor_modules
from tools.audit_centre_ricci_first_stage_split import decompose_primary_l2


def test_primary_l2_signed_split_closes_on_pinned_initializer():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=32,
        r_max=80.0,
        amplitude=0.01,
        width=7.0,
        D_amplitude=1.0e-10,
        include_radiation=True,
    )
    _, vacuum, _ = vendor_modules()
    actual = vacuum.primary_l2_rhs(state.grid, state.geometry)
    split = decompose_primary_l2(state.grid, state.geometry, actual["Aa"], vacuum)

    assert all(np.all(np.isfinite(values)) for values in split.values())
    assert np.all(
        np.abs(split["primary_l2_Aa_closure"])
        <= split["primary_l2_Aa_closure_bound"]
    )
    np.testing.assert_allclose(
        split["vendor_tf_ricci"],
        split["lambda_D4_ricci_piece_unweighted"]
        + split["lambda_algebraic_ricci_piece_unweighted"]
        + split["metric_tf_ricci_remainder_unweighted"],
        rtol=0.0,
        atol=8.0 * np.finfo(float).eps
        * max(1.0, float(np.max(np.abs(split["vendor_tf_ricci"])))),
    )
