"""Curved-slice, live-stage charge/Lorentz work balance audit.

Research-only discriminator for the opt-in radial-Maxwell candidate. It tests
whether the *actual* face-current charge update and cell-centred Lorentz work
approach the same integrated work after subtracting the explicit metric work
of the field-energy functional and exposing the finite outer-face contribution.

This is not a full donor+field+gravity conservation proof: the candidate uses
the zero-shift CMC branch, and a complete global budget must also include the
matter transport flux and the gravitational/geometry source budget. No model,
normalization, or production capability is authorized by this test.
"""
import json
import math

import numpy as np

from engine.matter_system import _metric_arrays, initialize_dust
from engine.production_kernel import V55ProductionKernel
from engine.radial_maxwell import RadialMaxwellState
from engine.v55_matter import V55MatterState, metric_slice_from_q


FOUR_PI = 4.0 * math.pi
ETA = 0.25
R_MAX = 12.0
DENSITY = 2.0e-5


def _geometry_snapshot(geometry):
    return {
        name: np.asarray(getattr(geometry, name), dtype=float).copy()
        for name in ("a", "b", "X", "alpha", "beta")
    }


def _weight(grid, geometry):
    # sqrt(gamma)/R_areal**4 for
    # R_areal=r*sqrt(b)/X and sqrt(gamma)=r**2*sqrt(a)*b/X**3.
    r = np.asarray(grid.centers, dtype=float)
    return (
        np.asarray(geometry["X"], dtype=float)
        * np.sqrt(np.asarray(geometry["a"], dtype=float))
        / (np.square(r) * np.asarray(geometry["b"], dtype=float))
    )


def _field_energy(grid, geometry, electromagnetic):
    q_cell = 0.5 * (
        electromagnetic.enclosed_charge_faces[:-1]
        + electromagnetic.enclosed_charge_faces[1:]
    )
    return 0.5 * grid.dr * float(np.sum(_weight(grid, geometry) * q_cell**2))


def _stage_powers(grid, snapshot, rest_mass_flux, explicit_terms):
    em = snapshot["electromagnetic"]
    geometry = snapshot["geometry"]
    q_cell = 0.5 * (
        em.enclosed_charge_faces[:-1] + em.enclosed_charge_faces[1:]
    )
    weight = _weight(grid, geometry)
    current_faces = em.charge_flux_faces(rest_mass_flux)
    q_cell_dot = -0.5 * (current_faces[:-1] + current_faces[1:])

    charge_power = grid.dr * float(np.sum(weight * q_cell * q_cell_dot))
    log_weight_dot = (
        explicit_terms["X"] / geometry["X"]
        + 0.5 * explicit_terms["a"] / geometry["a"]
        - explicit_terms["b"] / geometry["b"]
    )
    metric_power = 0.5 * grid.dr * float(
        np.sum(weight * q_cell**2 * log_weight_dot)
    )

    # source_energy is the same stage-local alpha*sqrt(gamma)*E_r*J^r
    # returned by RadialMaxwellState.lorentz_source. Its radial array is per
    # unit solid angle, so integrate over 4*pi and the coordinate cell width.
    matter_lorentz_power = (
        FOUR_PI * grid.dr * float(np.sum(snapshot["source_energy"]))
    )

    # In the discrete U_EM functional, the outer face Q_N enters only through
    # Q_cell,N-1=(Q_{N-1}+Q_N)/2. Its charge-rate contribution is therefore
    # -dr/2 * (sqrt(gamma_rr)*Ehat)_last * I_outer.
    outer_current_power = (
        -0.5 * grid.dr * weight[-1] * q_cell[-1] * current_faces[-1]
    )
    return {
        "charge_power": charge_power,
        "metric_power": metric_power,
        "matter_lorentz_power": matter_lorentz_power,
        "outer_current_power": outer_current_power,
    }


def _run_probe(resolution, dt):
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=resolution,
        r_max=R_MAX,
        include_radiation=True,
        radial_maxwell_charge_to_mass=ETA,
        charged_dust_density=np.full(resolution, DENSITY, dtype=float),
    )

    # The analytic velocity profile below vanishes at both radial endpoints.
    # The finite-volume outer-face flux is nevertheless measured from the actual
    # last-cell/stage reconstruction; do not force it to zero by assumption.
    # Rebuild Q from the exact rest state so Gauss starts consistent; no
    # existing GEAR species is charged or relabelled.
    r = np.asarray(state.grid.centers, dtype=float)
    x = r / R_MAX
    velocity = 0.08 * x * (1.0 - x)
    carrier = initialize_dust(
        _metric_arrays(metric_slice_from_q(state.grid, state.geometry)),
        np.full(resolution, DENSITY, dtype=float),
        velocity,
    )
    state.matter = V55MatterState(
        state.matter.dark_matter,
        state.matter.baryons,
        state.matter.radiation,
        charged_dust=carrier,
    )
    state.electromagnetic = RadialMaxwellState.from_charged_dust(
        state.grid, carrier, ETA
    )

    stage_snapshots = []
    face_fluxes = []
    explicit_terms = []
    rhs = kernel._rhs
    flux = kernel._charged_dust_mass_flux_faces
    geometry_terms = kernel._geometry_stage_terms

    def record_rhs(stage):
        result = rhs(stage)
        if stage.electromagnetic is not None:
            metric = metric_slice_from_q(stage.grid, stage.geometry)
            _, _, values = stage.electromagnetic.lorentz_source(
                stage.grid, stage.geometry, metric, stage.matter.charged_dust
            )
            stage_snapshots.append({
                "geometry": _geometry_snapshot(stage.geometry),
                "electromagnetic": stage.electromagnetic.copy(),
                "source_energy": np.asarray(
                    values["source_energy"], dtype=float
                ).copy(),
            })
        return result

    def record_flux(stage):
        values = flux(stage)
        face_fluxes.append(np.asarray(values, dtype=float).copy())
        return values

    def record_geometry(grid, geometry, scalars, matter, electromagnetic):
        result = geometry_terms(
            grid, geometry, scalars, matter, electromagnetic
        )
        explicit_terms.append({
            name: np.asarray(result["explicit"][name], dtype=float).copy()
            for name in ("a", "b", "X")
        })
        return result

    kernel._rhs = record_rhs
    kernel._charged_dust_mass_flux_faces = record_flux
    kernel._geometry_stage_terms = record_geometry
    accepted = kernel.step(state, dt)

    # Actual order inside step(): start geometry block, explicit predictor
    # block, then primary-predictor geometry block. The two matter/current
    # RHS stages are start and primary predictor.
    assert len(stage_snapshots) == 2
    assert len(face_fluxes) == 2
    assert len(explicit_terms) == 3

    start_power = _stage_powers(
        state.grid, stage_snapshots[0], face_fluxes[0], explicit_terms[0]
    )
    predictor_power = _stage_powers(
        state.grid, stage_snapshots[1], face_fluxes[1], explicit_terms[2]
    )

    # The candidate explicitly selects zero shift; don't silently generalize
    # this witness to a nonzero-shift evolution.
    assert np.all(stage_snapshots[0]["geometry"]["beta"] == 0.0)
    assert np.all(stage_snapshots[1]["geometry"]["beta"] == 0.0)
    assert np.max(np.abs(face_fluxes[0])) > 0.0
    assert np.max(np.abs(stage_snapshots[0]["source_energy"])) > 0.0

    u_start = _field_energy(
        state.grid,
        stage_snapshots[0]["geometry"],
        stage_snapshots[0]["electromagnetic"],
    )
    u_end = _field_energy(
        accepted.grid,
        _geometry_snapshot(accepted.geometry),
        accepted.electromagnetic,
    )
    avg = lambda key: 0.5 * (
        start_power[key] + predictor_power[key]
    )

    # ΔU_EM contains (i) charge work from the face update and (ii) geometry
    # work through w=X*sqrt(a)/(r^2*b). Subtract the RK2 geometry work, add
    # the reciprocal matter Lorentz work, and remove the separately exposed
    # outer-face contribution. What remains is the spatial discrete
    # face-current/cell-current mismatch plus the bounded time-integration
    # truncation error, suitable for a mesh-convergence discriminator.
    residual_rate = (
        (u_end - u_start) / dt
        - avg("metric_power")
        + avg("matter_lorentz_power")
        - avg("outer_current_power")
    )
    scale = 0.5 * (
        abs(start_power["charge_power"] - start_power["outer_current_power"])
        + abs(predictor_power["charge_power"] - predictor_power["outer_current_power"])
        + abs(start_power["matter_lorentz_power"])
        + abs(predictor_power["matter_lorentz_power"])
    )
    relative_residual = abs(residual_rate) / max(scale, 1.0e-30)

    # Negative control: a reversed reciprocal Lorentz sign should destroy the
    # matter/field exchange cancellation, while keeping geometry and boundary
    # terms unchanged. Keep this value separate from the accepted residual.
    field_energy_rate = (u_end - u_start) / dt
    wrong_lorentz_sign_rate = (
        field_energy_rate
        - avg("metric_power")
        - avg("matter_lorentz_power")
        - avg("outer_current_power")
    )
    wrong_lorentz_sign_relative = (
        abs(wrong_lorentz_sign_rate) / max(scale, 1.0e-30)
    )

    assert math.isfinite(residual_rate)
    assert math.isfinite(relative_residual)
    assert math.isfinite(wrong_lorentz_sign_rate)
    assert math.isfinite(wrong_lorentz_sign_relative)
    return {
        "resolution": resolution,
        "dt": dt,
        "field_energy_rate": field_energy_rate,
        "residual_rate": residual_rate,
        "relative_residual": relative_residual,
        "wrong_lorentz_sign_residual_rate": wrong_lorentz_sign_rate,
        "wrong_lorentz_sign_relative_residual": wrong_lorentz_sign_relative,
        "metric_power": avg("metric_power"),
        "matter_lorentz_power": avg("matter_lorentz_power"),
        "outer_current_power": avg("outer_current_power"),
        "charge_power": avg("charge_power"),
    }


def test_curved_slice_em_exchange_residual_converges_on_live_pirk_stages(record_property, capsys):
    # Keep the same physical domain and reduce dt with dr so both spatial and
    # stage-quadrature errors are bounded. This is one deterministic step per
    # resolution, not a production campaign or trajectory-admission test.
    results = [
        _run_probe(n, dt=0.002 * 24.0 / n)
        for n in (24, 48, 96)
    ]
    errors = [row["relative_residual"] for row in results]

    for row in results:
        # Scale each named contribution against the largest measured power;
        # this rejects a vacuous zero term without imposing an arbitrary
        # absolute unit-dependent cutoff.
        power_scale = max(
            abs(row[name])
            for name in (
                "field_energy_rate", "charge_power", "metric_power",
                "matter_lorentz_power", "outer_current_power",
            )
        )
        roundoff_floor = 128.0 * np.finfo(float).eps * max(power_scale, 1.0e-30)
        assert abs(row["metric_power"]) > roundoff_floor, row
        assert abs(row["outer_current_power"]) > roundoff_floor, row
        assert abs(row["matter_lorentz_power"]) > roundoff_floor, row
        # With the wrong sign, the normalized residual must remain visibly
        # larger than roundoff and a material fraction of the real Lorentz
        # work; otherwise this fixture would not discriminate the sign.
        assert abs(row["wrong_lorentz_sign_residual_rate"]) > max(
            0.5 * abs(row["matter_lorentz_power"]), roundoff_floor
        ), row

    assert all(math.isfinite(value) for value in errors), results
    assert all(value >= 0.0 for value in errors), results
    assert errors[1] < errors[0], results
    assert errors[2] < errors[1], results

    refinement_ratios = {
        "N24_to_N48": errors[0] / errors[1],
        "N48_to_N96": errors[1] / errors[2],
    }
    report = {
        "scope": "research-only; one step; zero shift; not full conservation/admission",
        "results": results,
        "refinement_ratios": refinement_ratios,
    }
    # Structured output is preserved in the Actions log by the targeted -s
    # workflow step; record_property also keeps it machine-readable in JUnit
    # output for runners that enable that report format.
    record_property(
        "em_curved_slice_exchange_evidence_json",
        json.dumps(report, sort_keys=True, separators=(",", ":")),
    )
    with capsys.disabled():
        print("EM_CURVED_SLICE_EXCHANGE_EVIDENCE=" + json.dumps(report, sort_keys=True), flush=True)
