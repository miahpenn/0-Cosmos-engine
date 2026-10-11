"""Fast schema/identity checks for the diagnostic-only outer constraint audit."""
import numpy as np

from engine.production_kernel import V55ProductionKernel
from engine.outer_constraint_audit import record_sample, DECOMPOSITION_TOLERANCE


def test_full_grid_constraint_profile_is_finite_and_reconstructs_hamiltonian():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=16,
        r_max=80.0,
        D_amplitude=1.0e-4,
        include_radiation=True,
    )
    sample = record_sample(state, sample_index=0, expected_target=0.0)

    profile = sample["profile"]
    assert len(profile) == 16
    assert [row["cell"] for row in profile] == list(range(16))
    assert set(sample["regional_summary"]) == {
        "full_grid", "inner_r_le_20", "middle_20_lt_r_lt_64", "outer_r_ge_64"
    }
    assert sum(
        sample["regional_summary"][name]["H"]["cell_count"]
        for name in ("inner_r_le_20", "middle_20_lt_r_lt_64", "outer_r_ge_64")
    ) == 16
    assert sample["decomposition_error_max_abs_full_grid"] <= DECOMPOSITION_TOLERANCE
    numeric = [
        value
        for row in profile
        for key, value in row.items()
        if key != "cell"
    ]
    assert np.all(np.isfinite(np.asarray(numeric, dtype=float)))
    assert sample["regional_summary"]["full_grid"]["H"]["max_cell"] in range(16)
    assert sample["regional_summary"]["full_grid"]["M"]["max_cell"] in range(16)
