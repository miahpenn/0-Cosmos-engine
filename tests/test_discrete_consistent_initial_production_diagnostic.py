import math
import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy


def test_opt_in_discrete_consistent_initial_state_closes_discrete_hamiltonian():
    """The opt-in solve must close the exact production H definition at t=0."""
    kernel = V55ProductionKernel()
    state, B, history, info = discrete_consistent_state(
        kernel,
        resolution=40,
        r_max=80.0,
        max_iter=12,
        tol=1.0e-13,
        max_cond=1.0e12,
        return_info=True,
    )
    _, vacuum, _ = adapter.vendor_modules()
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    H = np.asarray(vacuum.constraints(state.grid, state.geometry)["hamiltonian"], dtype=float) \
        - 16.0 * math.pi * np.asarray(total.rho, dtype=float)

    assert info["final_max_residual"] <= info["tolerance"]
    assert np.all(np.isfinite(B)) and np.all(B > 0.0)
    assert np.all(np.isfinite(H))
    assert float(np.max(np.abs(H))) < 1.0e-12
    assert float(np.max(np.abs(state.geometry.a * state.geometry.b**2 - 1.0))) < 1.0e-12
    assert history[-1] == info["final_max_residual"]
