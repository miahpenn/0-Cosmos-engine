"""Frozen constraint-consumer source contract for Issue #4.

These tests validate only the factors/signs implemented by the pinned consumer:
    H_total = H_geom - 16*pi*rho_total
    M_total = M_geom -  8*pi*j_total
They do NOT decide the physical normalization of archived S/D fields, prove
the full Bianchi propagation law, or grant scientific admission.  The archive
action-to-T_mu_nu mapping remains a separate reviewed research obligation.
"""
import numpy as np
import pytest

from engine import v55_pirk_adapter as adapter
from engine.matter_system import (
    ConservedSpecies,
    Species,
    _metric_arrays,
    initialize_dust,
    initialize_radiation,
)
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import KAPPA, ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_matter import V55MatterState, metric_slice_from_q
from tools.audit_current_source_momentum_constraint_tangent import (
    assert_supported_baseline_source_scope,
    observe_constraint,
    state_digest,
)


SECTORS = ("S", "D", "COSMOS_phi", "dark_matter", "baryons", "radiation")


def _zero_conserved(n: int) -> ConservedSpecies:
    z = np.zeros(n, dtype=np.float64)
    return ConservedSpecies(z.copy(), z.copy(), z.copy())


def _zero_state(template: ProductionState) -> ProductionState:
    n = template.grid.n
    z = np.zeros(n, dtype=np.float64)
    scalars = ScalarFields(z.copy(), z.copy(), z.copy(), z.copy(), z.copy(), z.copy())
    metrics = _metric_arrays(metric_slice_from_q(template.grid, template.geometry))
    # Valencia dust recovery rejects exactly zero rest-current. Represent an
    # inactive dust sector by a tiny positive, zero-velocity state instead.
    rho_floor = np.full(n, 1.0e-20, dtype=np.float64)
    matter = V55MatterState(
        initialize_dust(metrics, rho_floor),
        initialize_dust(metrics, rho_floor),
        initialize_radiation(metrics, np.zeros(n, dtype=np.float64)),
    )
    return ProductionState(
        grid=template.grid,
        geometry=template.geometry.copy(),
        scalars=scalars,
        matter=matter,
        t=template.t,
        tau=template.tau,
        e_folds=template.e_folds,
    )


def _manufactured_state(
    template: ProductionState,
    sector: str,
    sign: float = 1.0,
    *,
    homogeneous: bool = False,
) -> ProductionState:
    state = _zero_state(template)
    grid = state.grid
    r = np.asarray(grid.centers, dtype=np.float64)
    metrics = _metric_arrays(metric_slice_from_q(grid, state.geometry))
    n = grid.n
    z = np.zeros(n, dtype=np.float64)

    if homogeneous:
        # Nonzero homogeneous sources, but zero canonical radial current.
        fields = ScalarFields(
            np.full(n, 0.01), z.copy(),
            np.full(n, 0.02), z.copy(),
            np.full(n, 0.1), z.copy(),
        )
        rho_dm = np.full(n, 2.0e-5)
        rho_b = np.full(n, 4.0e-6)
        rho_rad = np.full(n, 8.0e-9)
        matter = V55MatterState(
            initialize_dust(metrics, rho_dm),
            initialize_dust(metrics, rho_b),
            initialize_radiation(metrics, rho_rad),
        )
        state.scalars = fields
        state.matter = matter
        return state

    fields = {
        "S": z.copy(),
        "PS": z.copy(),
        "D": z.copy(),
        "PD": z.copy(),
        "phi": z.copy(),
        "Pi": z.copy(),
    }
    if sector == "S":
        fields["S"] = 0.03 * np.exp(-(r / 4.0) ** 2)
        fields["PS"] = sign * 0.004 * np.exp(-(r / 5.0) ** 2)
    elif sector == "D":
        fields["D"] = 0.02 * np.exp(-(r / 5.0) ** 2)
        fields["PD"] = sign * 0.003 * np.exp(-(r / 6.0) ** 2)
    elif sector == "COSMOS_phi":
        fields["phi"] = 0.025 * np.exp(-(r / 6.0) ** 2)
        fields["Pi"] = sign * 0.004 * np.exp(-(r / 7.0) ** 2)
    state.scalars = ScalarFields(*(np.asarray(fields[name]) for name in
                                  ("S", "PS", "D", "PD", "phi", "Pi")))

    rho = 1.0e-4 * (1.0 + 0.2 * np.exp(-(r / 5.0) ** 2))
    velocity = sign * 0.025 * r / (1.0 + r) * np.exp(-(r / 10.0) ** 2)
    rho_floor = np.full(n, 1.0e-20, dtype=np.float64)
    dm = initialize_dust(metrics, rho, velocity) if sector == "dark_matter" else initialize_dust(metrics, rho_floor)
    baryons = initialize_dust(metrics, rho, velocity) if sector == "baryons" else initialize_dust(metrics, rho_floor)
    radiation = initialize_radiation(metrics, rho, velocity) if sector == "radiation" else initialize_radiation(metrics, np.zeros(n, dtype=np.float64))
    state.matter = V55MatterState(dm, baryons, radiation)
    assert sector in SECTORS
    return state


def _consumer_arrays(state: ProductionState):
    """Reconstruct the arrays whose absolute maxima are reported by diagnostics."""
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(state.grid, state.geometry)
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    h_total = np.asarray(raw["hamiltonian"]) - 16.0 * np.pi * total.rho
    m_total = np.asarray(raw["momentum"]) - 8.0 * np.pi * total.j
    return h_total, m_total, total


@pytest.fixture(scope="module")
def pinned_template():
    # Geometry is held fixed across the sign-reversal pair so changes in H/M
    # isolate the explicit source projection, not an evolved metric response.
    return V55ProductionKernel().initialize(
        resolution=32, r_max=16.0, include_radiation=True
    )


@pytest.mark.parametrize("sector", SECTORS)
def test_sector_current_sign_reversal_and_constraint_consumer_coefficients(
    pinned_template, sector
):
    plus = _manufactured_state(pinned_template, sector, +1.0)
    minus = _manufactured_state(pinned_template, sector, -1.0)
    h_plus, m_plus, total_plus = _consumer_arrays(plus)
    h_minus, m_minus, total_minus = _consumer_arrays(minus)

    # A momentum sign reversal leaves rho, p_r, p_t invariant and flips j.
    np.testing.assert_allclose(total_plus.rho, total_minus.rho, rtol=3e-13, atol=1e-24)
    np.testing.assert_allclose(total_plus.pr, total_minus.pr, rtol=3e-13, atol=1e-24)
    np.testing.assert_allclose(total_plus.pt, total_minus.pt, rtol=3e-13, atol=1e-24)
    np.testing.assert_allclose(total_plus.j, -total_minus.j, rtol=3e-13, atol=1e-24)
    assert float(np.max(np.abs(total_plus.j))) > 0.0, f"{sector} control has zero current"

    # These are the exact factors currently wired into the consumer.  They
    # are intentionally not interpreted as proof of the archived S/D physics.
    _, vacuum, _ = adapter.vendor_modules()
    raw = vacuum.constraints(plus.grid, plus.geometry)
    np.testing.assert_allclose(
        h_plus - np.asarray(raw["hamiltonian"]),
        -16.0 * np.pi * total_plus.rho,
        rtol=3e-13, atol=1e-23,
    )
    np.testing.assert_allclose(
        m_plus - np.asarray(raw["momentum"]),
        -8.0 * np.pi * total_plus.j,
        rtol=3e-13, atol=1e-23,
    )
    np.testing.assert_allclose(
        h_plus, h_minus, rtol=3e-13, atol=1e-22,
    )
    np.testing.assert_allclose(
        m_plus - m_minus, -16.0 * np.pi * total_plus.j,
        rtol=3e-13, atol=1e-22,
    )

    # Exercise the real diagnostic consumer too: its reported peak constraints
    # must match the independently reconstructed arrays above.
    kernel = V55ProductionKernel()
    for state, h_expected, m_expected, total in (
        (plus, h_plus, m_plus, total_plus),
        (minus, h_minus, m_minus, total_minus),
    ):
        observed = kernel.diagnostics(state, profiles=True)
        np.testing.assert_allclose(observed["total_rho"], total.rho, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(observed["total_j"], total.j, rtol=0.0, atol=0.0)
        assert observed["hamiltonian_full_absmax"] == pytest.approx(
            float(np.max(np.abs(h_expected))), rel=3e-13, abs=1e-23
        )
        assert observed["momentum_full_absmax"] == pytest.approx(
            float(np.max(np.abs(m_expected))), rel=3e-13, abs=1e-23
        )


def test_homogeneous_nonzero_sources_are_a_zero_current_null_control(pinned_template):
    state = _manufactured_state(pinned_template, "all", homogeneous=True)
    _h, _m, total = _consumer_arrays(state)
    assert np.any(np.abs(total.rho) > 0.0)
    assert np.all(total.j == 0.0)
    pieces = observe_constraint(state)
    assert np.all(pieces["J_total"] == 0.0)


def test_baseline_six_sector_audit_rejects_opt_in_source_states(pinned_template):
    from engine.matter_system import initialize_dust

    charged = _zero_state(pinned_template)
    metrics = _metric_arrays(metric_slice_from_q(charged.grid, charged.geometry))
    density = np.full(charged.grid.n, 1.0e-8)
    charged.matter.charged_dust = initialize_dust(metrics, density)
    with pytest.raises(RuntimeError, match="charged_dust"):
        assert_supported_baseline_source_scope(charged)
    with pytest.raises(RuntimeError, match="charged_dust"):
        state_digest(charged)
    with pytest.raises(RuntimeError, match="charged_dust"):
        observe_constraint(charged)

    em = _zero_state(pinned_template)
    # The scope guard must reject the feature's presence before trying to
    # inspect or interpret its concrete class/source-normalization convention.
    em.electromagnetic = object()
    with pytest.raises(RuntimeError, match="radial EM"):
        assert_supported_baseline_source_scope(em)
    with pytest.raises(RuntimeError, match="radial EM"):
        state_digest(em)
    with pytest.raises(RuntimeError, match="radial EM"):
        observe_constraint(em)
