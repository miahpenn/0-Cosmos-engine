import numpy as np

from engine.reference_v55_geodesic import H_constraint, initial, NG


def test_v55_reference_initial_constraints_are_finite():
    U, grid, bg = initial(nphys=40, rmax=20.0, Aamp=0.003, width=7.0)
    H, M, R, em = H_constraint(U, grid, bg)
    sl = slice(NG + 4, -NG - 4)
    for arr in (H, M, R, em.rho):
        assert np.all(np.isfinite(arr[sl]))
