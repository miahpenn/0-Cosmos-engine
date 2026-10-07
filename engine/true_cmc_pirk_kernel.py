"""Stage-aware archive-anchored CMC V5.5 strong-field kernel.

This kernel changes only the foliation layer relative to the existing V5.5
reference-metric PIRK system. The CMC lapse is solved from the ADM K equation
at the beginning of a step, on the explicit predictor, and after the completed
primary state. The spatial shift is held at beta=0 so the lapse/foliation gate
is isolated from shift-driver behavior.

The V5.5 scalar, conservative matter, stress-energy, Misner-Sharp, and cycle
ledger definitions are inherited unchanged. No lapse floor, bounce, reset,
fitted coefficient, interface source, or matter clip is introduced.
"""

from __future__ import annotations

import math
import numpy as np
from .production_kernel import V55ProductionKernel, ProductionState, LAMBDA_M
from .scalar_system import ScalarFields
from .matter_system import ConservedSpecies
from .v55_matter import metric_slice_from_q
from .cosmology_observables import append_efolds
from .handoff import handoff_from_ledger
from . import v55_pirk_adapter as adapter


from .cmc_gauge import solve_cmc_lapse

OUTER_CMC_FRACTION = 0.20


def solve_archive_cmc_lapse(grid, geometry, scalars, matter):
    """Compatibility name for the shared production CMC solver.

    The production and isolated true-CMC kernels must carry the same
    stage-aware gauge operator and the same proper-volume CMC target.
    """
    return solve_cmc_lapse(grid, geometry, scalars, matter)
 
 
class V55TrueCMCPIRKKernel(V55ProductionKernel):
    """V5.5 with stage-aware CMC lapse and zero spatial shift."""

    validation_state = "implementation_smoke_pending"

    @staticmethod
    def _regularize(grid, geometry):
        out = V55ProductionKernel._enforce_center_regularity(
            grid, geometry
        )
        out.beta.fill(0.0)
        out.B.fill(0.0)
        return out

    def initialize(
        self,
        resolution: int = 80,
        r_max: float = 40.0,
        *,
        amplitude: float = 0.01,
        width: float = 7.0,
        D_amplitude: float = 1.0e-10,
        include_radiation: bool = True,
    ) -> ProductionState:
        grid_ops, vacuum, _ = adapter.vendor_modules()
        grid = grid_ops.SphericalCellGrid(resolution, r_max)
        init = __import__(
            "engine.v55_initial",
            fromlist=["build_initial_data"],
        ).build_initial_data(
            grid,
            vacuum.flat_state,
            amplitude=amplitude,
            width=width,
            D_amplitude=D_amplitude,
            include_radiation=include_radiation,
        )
        geometry = self._regularize(grid, init.geometry)
        geometry.alpha = solve_archive_cmc_lapse(
            grid, geometry, init.scalars, init.matter
        )[0]
        geometry.beta.fill(0.0)
        geometry.B.fill(0.0)
        geometry.assert_finite_positive()
        return ProductionState(
            grid=grid,
            geometry=geometry,
            scalars=init.scalars,
            matter=init.matter,
            e_folds=0.0,
        )

    @staticmethod
    def _validate_matter_state(grid, geometry, matter):
        metric = metric_slice_from_q(grid, geometry)
        metrics = __import__(
            "engine.matter_system",
            fromlist=["_metric_arrays", "Species", "primitives"],
        )
        arrays = metrics._metric_arrays(metric)
        for name, species in (
            ("dark_matter", metrics.Species.DARK_MATTER),
            ("baryons", metrics.Species.BARYON),
            ("radiation", metrics.Species.RADIATION),
        ):
            metrics.primitives(arrays, getattr(matter, name), species)

    @staticmethod
    def _matter_euler(base, rhs, dt):
        return type(base)(
            *(
                ConservedSpecies(
                    getattr(base, name).rest + dt * rhs[name].rest,
                    getattr(base, name).energy_t + dt * rhs[name].energy_t,
                    getattr(base, name).momentum_r + dt * rhs[name].momentum_r,
                )
                for name in ("dark_matter", "baryons", "radiation")
            )
        )

    @staticmethod
    def _matter_trapezoid(base, rhs0, rhs1, dt):
        return type(base)(
            *(
                ConservedSpecies(
                    getattr(base, name).rest + 0.5 * dt * (
                        rhs0[name].rest + rhs1[name].rest
                    ),
                    getattr(base, name).energy_t + 0.5 * dt * (
                        rhs0[name].energy_t + rhs1[name].energy_t
                    ),
                    getattr(base, name).momentum_r + 0.5 * dt * (
                        rhs0[name].momentum_r + rhs1[name].momentum_r
                    ),
                )
                for name in ("dark_matter", "baryons", "radiation")
            )
        )

    def step(self, state: ProductionState, dt: float) -> ProductionState:
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")

        grid = state.grid

        g0 = self._regularize(grid, state.geometry.copy())
        s0 = state.scalars
        m0 = state.matter
        g0.alpha = solve_archive_cmc_lapse(
            grid, g0, s0, m0
        )[0]

        srhs0, mrhs0, _ = self._rhs(
            ProductionState(
                grid=grid,
                geometry=g0,
                scalars=s0,
                matter=m0,
                t=state.t,
                tau=state.tau,
                e_folds=state.e_folds,
            )
        )
        gterms0 = adapter.geometry_stage_terms(
            grid, g0, s0, m0, lambda_m=LAMBDA_M
        )

        # PIRK explicit predictor: a,b,X evolve; alpha is a CMC-constrained
        # variable and beta/B are fixed by the isolated zero-shift gauge.
        gpred = g0.copy()
        for name in ("a", "b", "X"):
            setattr(
                gpred,
                name,
                getattr(g0, name) + dt * gterms0["explicit"][name],
            )
        gpred = self._regularize(grid, gpred)

        spred = ScalarFields(
            *(
                getattr(s0, name) + dt * getattr(srhs0, name)
                for name in ("S", "PS", "D", "PD", "phi", "Pi")
            )
        )
        mpred = self._matter_euler(m0, mrhs0, dt)
        self._validate_matter_state(grid, gpred, mpred)
        gpred.alpha = solve_archive_cmc_lapse(
            grid, gpred, spred, mpred
        )[0]

        gterms_pred = adapter.geometry_stage_terms(
            grid, gpred, spred, mpred, lambda_m=LAMBDA_M
        )

        g1 = gpred.copy()
        g1.Aa = g0.Aa + dt * (
            0.5 * gterms0["primary_l2"]["Aa"]
            + 0.5 * gterms_pred["primary_l2"]["Aa"]
            + gterms0["primary_l3"]["Aa"]
        )
        g1.K = g0.K + dt * (
            0.5 * gterms0["primary_l2"]["K"]
            + 0.5 * gterms_pred["primary_l2"]["K"]
            + gterms0["primary_l3"]["K"]
        )
        g1.Lambda = g0.Lambda + dt * (
            0.5 * gterms0["lambda_l2"]
            + 0.5 * gterms_pred["lambda_l2"]
            + gterms0["lambda_l3"]
        )
        g1 = self._regularize(grid, g1)
        g1.alpha = gpred.alpha.copy()
        g1.assert_finite_positive()

        stage1 = ProductionState(
            grid=grid,
            geometry=g1,
            scalars=spred,
            matter=mpred,
            t=state.t,
            tau=state.tau,
            e_folds=state.e_folds,
        )
        srhs1, mrhs1, _ = self._rhs(stage1)
        gterms1 = adapter.geometry_stage_terms(
            grid, g1, spred, mpred, lambda_m=LAMBDA_M
        )

        gnew = g0.copy()
        for name in ("a", "b", "X"):
            setattr(
                gnew,
                name,
                getattr(g0, name) + 0.5 * dt * (
                    gterms0["explicit"][name] + gterms1["explicit"][name]
                ),
            )
        gnew.beta.fill(0.0)
        gnew.B.fill(0.0)
        gnew.alpha = g1.alpha.copy()

        snew = ScalarFields(
            *(
                getattr(s0, name) + 0.5 * dt * (
                    getattr(srhs0, name) + getattr(srhs1, name)
                )
                for name in ("S", "PS", "D", "PD", "phi", "Pi")
            )
        )
        mnew = self._matter_trapezoid(m0, mrhs0, mrhs1, dt)
        self._validate_matter_state(grid, gnew, mnew)

        gnew = self._regularize(grid, gnew)

        # Final primary PIRK block uses the stage CMC lapse as the temporal
        # gauge value, then the completed slice gets its own CMC solve.
        final_u_old_v = g1.copy()
        for name in ("a", "b", "X"):
            setattr(final_u_old_v, name, getattr(gnew, name))

        _, vacuum, _ = adapter.vendor_modules()
        l2_final = vacuum.primary_l2_rhs(grid, final_u_old_v)
        gnew.Aa = g0.Aa + 0.5 * dt * (
            gterms0["primary_l2"]["Aa"]
            + l2_final["Aa"]
            + gterms0["primary_l3"]["Aa"]
            + gterms1["primary_l3"]["Aa"]
        )
        gnew.K = g0.K + 0.5 * dt * (
            gterms0["primary_l2"]["K"]
            + l2_final["K"]
            + gterms0["primary_l3"]["K"]
            + gterms1["primary_l3"]["K"]
        )

        final_primary = final_u_old_v.copy()
        final_primary.Aa = gnew.Aa.copy()
        final_primary.K = gnew.K.copy()
        gnew.Lambda = g0.Lambda + 0.5 * dt * (
            gterms0["lambda_l2"]
            + vacuum.lambda_l2_rhs(
                grid, final_primary, lambda_m=LAMBDA_M
            )
            + gterms0["lambda_l3"]
            + gterms1["lambda_l3"]
        )

        gnew = self._regularize(grid, gnew)
        gnew.alpha = solve_archive_cmc_lapse(
            grid, gnew, snew, mnew
        )[0]
        gnew.beta.fill(0.0)
        gnew.B.fill(0.0)
        gnew.assert_finite_positive()

        candidate = ProductionState(
            grid=grid,
            geometry=gnew,
            scalars=snew,
            matter=mnew,
            t=state.t + dt,
            tau=state.tau,
            e_folds=state.e_folds,
            cycle=state.cycle,
            history=list(state.history),
            handoffs=list(state.handoffs),
        )

        obs = self.diagnostics(candidate, profiles=False)
        candidate.e_folds = append_efolds(
            state.e_folds, obs["H_eff"], dt
        )
        previous_H = (
            state.history[-1]["H_eff"]
            if state.history else obs["H_eff"]
        )
        candidate.tau += dt * obs["tau_rate"]
        event = candidate.cycle.observe(
            len(candidate.history),
            candidate.t,
            candidate.tau,
            previous_H,
            obs["H_eff"],
        )

        hp = None
        if event is not None:
            # Cycle handoffs are event records, not per-step samples. Localize
            # the H_eff crossing inside this step using the solved endpoints.
            h0 = float(previous_H)
            h1 = float(obs["H_eff"])
            denom = h0 - h1
            fraction = 0.5 if abs(denom) <= 1.0e-300 else h0 / denom
            fraction = min(1.0, max(0.0, fraction))

            event_obs = dict(obs)
            if state.history:
                previous_obs = state.history[-1]
                for key in (
                    "phi_outer", "R_sigma", "M_MS",
                    "chi_sigma", "flux_T", "work_pR",
                ):
                    if key in previous_obs and key in obs:
                        event_obs[key] = float(
                            previous_obs[key]
                            + fraction * (obs[key] - previous_obs[key])
                        )
            event_obs["H_eff"] = 0.0
            event_t = state.t + fraction * dt
            event_tau = state.tau + fraction * (candidate.tau - state.tau)

            s_event = float(
                state.scalars.S[0]
                + fraction * (candidate.scalars.S[0] - state.scalars.S[0])
            )
            d_event = float(
                state.scalars.D[0]
                + fraction * (candidate.scalars.D[0] - state.scalars.D[0])
            )
            sdot_event = float(
                srhs0.PS[0]
                + fraction * (srhs1.PS[0] - srhs0.PS[0])
            )
            ddot_event = float(
                srhs0.PD[0]
                + fraction * (srhs1.PD[0] - srhs0.PD[0])
            )
            hp = handoff_from_ledger(
                event_t,
                event_tau,
                event_obs,
                s_event,
                sdot_event,
                d_event,
                ddot_event,
            )

        obs = {
            key: value
            for key, value in obs.items()
            if key not in (
                "total_rho",
                "total_pr",
                "total_pt",
                "total_j",
                "areal_radius",
                "chi",
                "misner_sharp",
            )
        }
        obs["t"] = candidate.t
        obs["tau"] = candidate.tau
        obs["cycle_event"] = event.kind if event else None
        candidate.history.append(obs)
        if hp is not None:
            candidate.handoffs.append(hp)
        return candidate
