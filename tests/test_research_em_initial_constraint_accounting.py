"""Initial-slice constraint accounting for the opt-in radial-EM candidate.

This is a research-only cross-component audit. It compares otherwise identical
separate-carrier initial data with q/m=0 and q/m!=0, then checks the actual
initial geometry/constraint diagnostic against the canonical stress projection.
It tests inclusion and observability, not constraint convergence, field units,
a GEAR-derived carrier, or scientific admission.
"""
import math

import numpy as np
import pytest

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy


ETA = 0.25
R_MAX = 12.0
DENSITY = 2.0e-5


def _initial_state(resolution, charge_to_mass):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=R_MAX,
        include_radiation=True,
        radial_maxwell_charge_to_mass=charge_to_mass,
        charged_dust_density=np.full(resolution, DENSITY, dtype=float),
    )
    assert state.electromagnetic is not None
    assert state.matter.charged_dust is not None
    return kernel, state


@pytest.mark.parametrize("resolution", [24, 48])
def test_nonzero_em_initial_data_reaches_metric_and_constraint_accounting(resolution):
    # A zero-q/m control keeps the *additional dust carrier* and its profile
    # fixed, isolating the initial electric-field contribution.
    kernel, zero = _initial_state(resolution, 0.0)
    _, charged = _initial_state(resolution, ETA)

    assert zero.electromagnetic.charge_to_rest_mass == 0.0
    assert np.all(zero.electromagnetic.enclosed_charge_faces == 0.0)
    assert charged.electromagnetic.charge_to_rest_mass == ETA
    assert charged.electromagnetic.enclosed_charge_faces[0] == 0.0

    em_zero = zero.electromagnetic.stress_projection(zero.grid, zero.geometry)
    em = charged.electromagnetic.stress_projection(charged.grid, charged.geometry)
    assert np.all(em_zero["rho"] == 0.0)
    assert np.all(np.isfinite(em["rho"]))
    assert np.max(np.abs(em["rho"])) > 0.0
    assert all(np.all(np.isfinite(v)) for v in em.values())

    # The initial radial Hamiltonian reconstruction claims to include u_EM.
    # On identical density/profile data, changing only q/m must therefore
    # change the initialized spatial metric (X is independent of the later
    # CMC lapse solve). This is a measured source-to-initial-data interaction,
    # not a proof that the residual converges to zero.
    metric_delta = float(np.max(np.abs(charged.geometry.X - zero.geometry.X)))
    assert math.isfinite(metric_delta)
    assert metric_delta > 0.0

    # Recompute the actual constraint residual on this initialized slice with
    # and without the EM projection, using the same raw geometric operators.
    # Their difference must be exactly the single canonical 16*pi*rho_EM
    # term; EM has j=0, so its addition cannot change the momentum constraint.
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(charged.grid, charged.geometry)
    matter_only = assemble_total_stress_energy(
        charged.grid, charged.geometry, charged.scalars, charged.matter
    )
    total = assemble_total_stress_energy(
        charged.grid,
        charged.geometry,
        charged.scalars,
        charged.matter,
        electromagnetic_projection=em,
    )
    H_without_em = np.asarray(raw["hamiltonian"], dtype=float) - (
        16.0 * math.pi * matter_only.rho
    )
    H_with_em = np.asarray(raw["hamiltonian"], dtype=float) - (
        16.0 * math.pi * total.rho
    )
    M_without_em = np.asarray(raw["momentum"], dtype=float) - (
        8.0 * math.pi * matter_only.j
    )
    M_with_em = np.asarray(raw["momentum"], dtype=float) - (
        8.0 * math.pi * total.j
    )

    np.testing.assert_allclose(
        H_without_em - H_with_em,
        16.0 * math.pi * em["rho"],
        rtol=2.0e-11,
        atol=2.0e-13,
    )
    np.testing.assert_array_equal(M_without_em, M_with_em)

    # Verify the diagnostics surfaced by the initialized machine are computed
    # from the same stress-inclusive residual, while recording their actual
    # (not assumed-zero) magnitude for follow-up convergence work.
    observed = kernel.diagnostics(charged)
    for name, values in (
        ("hamiltonian_full_absmax", H_with_em),
        ("momentum_full_absmax", M_with_em),
    ):
        expected = float(np.max(np.abs(values)))
        assert math.isfinite(expected)
        assert math.isfinite(observed[name])
        assert observed[name] == pytest.approx(expected, rel=2.0e-11, abs=2.0e-13)

    # Explicitly ensure this fixture really is a nonzero field/stress case:
    # a zero projection or an accidentally omitted EM term must not pass.
    assert np.max(np.abs(H_without_em - H_with_em)) > 0.0
