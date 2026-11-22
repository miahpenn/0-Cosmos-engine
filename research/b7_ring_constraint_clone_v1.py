"""RESEARCH ONLY. B7 W2: research clone of the engine's constrained initial-data reconstruction (engine/v55_initial.py:258-296).

Same arithmetic and operation order as the engine loop, with the momentum density j = -(PD D_r) extended by -(PS S_r) (flag),
and K_t = 0 as a labelled research-only isolated option (flag). Engine functions are imported read-only; nothing in
engine/ or vendor/ is changed. Design and expected verdicts: Project doc claude/b7-w2-constraint-clone-design.md.
"""
import math

import numpy as np

from engine import scalar_system as _ss
from engine.scalar_system import cosmos_potential, sd_potential
from engine.v55_initial import CorrectedInitialData, _cumulative_origin_integral, _d1
from engine.v55_matter import RHO_B_PRESENT, RHO_DM_PRESENT, RHO_R_PRESENT, initialize_from_archive

PHI0, PI0 = 0.179055, 0.00341913          # the engine's literals (v55_initial.py ~212)


def background_H0(include_radiation=True, sea_floor=True):
    """H0 exactly as v55_initial.py:~215-225 (sea_floor False = the hand_floor_reference branch)."""
    rho_r = RHO_R_PRESENT if include_radiation else 0.0
    rho0 = 0.5 * PI0 ** 2 + float(cosmos_potential(np.array([PHI0]))[0]) + RHO_DM_PRESENT + RHO_B_PRESENT + rho_r
    if sea_floor:
        rho0 += float(_ss.sea_floor_density(PHI0))
    return math.sqrt(max(rho0 / 3.0, 0.0))


def reconstruct_slice(grid, vacuum_state_factory, S, PS, D, PD, phi, Pi, *, include_ps_momentum=True,
                      include_ps_kinetic=True, Kt_mode="hubble", include_radiation=True, sea_floor=True,
                      n_iter=6):
    r = np.asarray(grid.centers)
    n = len(r)
    dr = float(grid.dr)
    rho_r = RHO_R_PRESENT if include_radiation else 0.0
    H0 = background_H0(include_radiation, sea_floor)
    rho_charged = np.zeros(n, dtype=float)
    dm0, b0 = RHO_DM_PRESENT, RHO_B_PRESENT
    if Kt_mode == "hubble":
        Kt = np.full(n, -H0)
    elif Kt_mode == "zero":
        Kt = np.zeros(n)                      # LABELLED research-only isolated start; the engine has none
    else:
        raise ValueError("Kt_mode must be 'hubble' or 'zero'")
    Sp = _d1(grid, S, 1)
    Dp = _d1(grid, D, 1)
    Ktp = _d1(grid, Kt, 1)
    j = -(PD * Dp)
    if include_ps_momentum:
        j = j - PS * Sp
    kin = 1.0 if include_ps_kinetic else 0.0
    B = np.ones_like(r)
    for _ in range(n_iter):
        if _ss.coupling_active():
            rhoS = 0.5 * (kin * PS ** 2 + B * Sp ** 2)
            rhoD = 0.5 * (PD ** 2 + B * Dp ** 2) + sd_potential(S, D, phi=phi)
        else:
            rhoS = 0.5 * (kin * PS ** 2 + B * Sp ** 2) + 0.5 * S ** 2
            rhoD = 0.5 * (PD ** 2 + B * Dp ** 2) - 0.5 * D ** 2
            if _ss.SD_QUARTIC_LAMBDA != 0.0:
                rhoD = rhoD + sd_potential(S, D) - (0.5 * S ** 2 - 0.5 * D ** 2)
        rhoPhi = (0.5 * Pi ** 2 + cosmos_potential(phi)) / (8.0 * math.pi)
        rhoMatter = (dm0 + b0 + rho_r + rho_charged) / (8.0 * math.pi)
        rho = rhoS + rhoD + rhoPhi + rhoMatter
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
        src = 1.0 - 8.0 * math.pi * r * r * rho + 2.0 * r * r * Kr * Kt + r * r * Kt * Kt
        integ = _cumulative_origin_integral(grid, r, dr, src)
        B = integ / r
    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError("invalid constrained initial metric")
    A = 1.0 / np.sqrt(B)
    Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
    K = Kr + 2.0 * Kt
    Aa = Kr - K / 3.0
    a = A ** (4.0 / 3.0)
    b = A ** (-2.0 / 3.0)
    X = A ** (-1.0 / 3.0)
    Lambda = (_d1(grid, a, 1) / (2.0 * a * a) - _d1(grid, b, 1) / (a * b) + 2.0 / r * (1.0 / b - 1.0 / a))
    state = vacuum_state_factory(grid)
    state.a[:] = a
    state.b[:] = b
    state.X[:] = X
    state.alpha[:] = 1.0
    state.beta[:] = 0.0
    state.Aa[:] = Aa
    state.K[:] = K
    state.Lambda[:] = Lambda
    state.B[:] = 0.75 * Lambda
    scalars = __import__("engine.scalar_system", fromlist=["ScalarFields"]).ScalarFields(S, PS, D, PD, phi, Pi)
    matter = initialize_from_archive(grid, state, include_radiation=include_radiation, charged_dust_density=None)
    return CorrectedInitialData(geometry=state, scalars=scalars, matter=matter, H0=H0, include_radiation=include_radiation), {"j": j, "B": B}


def research_state(kernel, grid, init):
    """Production state on the open edge as V55ProductionKernel.initialize builds it (production_kernel.py:~233-262)."""
    from engine.production_kernel import ProductionState
    init.geometry.alpha = kernel._solve_lapse_for_state(grid, init.geometry, init.scalars, init.matter, None)[0]
    init.geometry.beta.fill(0.0)
    init.geometry.B.fill(0.0)
    return ProductionState(grid=grid, geometry=init.geometry, scalars=init.scalars, matter=init.matter,
                           electromagnetic=None, e_folds=0.0)
