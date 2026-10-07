"""Moving-gauge V5.5 strong-field production kernel.

This branch changes only the strong-field evolution/gauge representation:
the pinned reference spherical BSSN moving-puncture/1+log PIRK2 ordering is
used for (a,b,X,alpha,beta,Aa,K,Lambda,B). The V5.5 scalar, conservative
matter, stress-energy, and source equations remain exactly those used by the
existing production kernel.

The moving gauge is the archive-backed horizon-penetrating candidate. No
lapse floor, reset, bounce rule, fitted damping coefficient, or matter clip is
introduced.

The centre regularity projection is retained, but the vendor helper's
algebraic B=3/4 Lambda assignment is not allowed to overwrite the independently
evolved moving-gauge B field. The projection is therefore applied to the
algebraic/centre identities while preserving the PIRK B update.
"""

from __future__ import annotations

import numpy as np

from .production_kernel import (
    V55ProductionKernel,
    ProductionState,
    LAMBDA_M,
)
from .scalar_system import ScalarFields
from .matter_system import ConservedSpecies
from .v55_initial import build_initial_data
from .cosmology_observables import append_efolds
from .handoff import handoff_from_ledger
from . import v55_pirk_adapter as adapter


class V55MovingPIRKKernel(V55ProductionKernel):
    """V5.5 with the pinned horizon-penetrating moving gauge."""

    validation_state = "implementation_smoke_pending"

    @staticmethod
    def _enforce_moving_center_regularity(grid, geometry):
        # The pinned algebraic projection supplies exact conformal determinant
        # and spherical-centre regularity, but also assigns B=3/4 Lambda.
        # In the moving gauge B is an independently evolved PIRK variable, so
        # preserve its solved value while using the projection for the exact
        # regularity identities.
        B_saved = geometry.B.copy()
        projected = V55ProductionKernel._enforce_center_regularity(
            grid, geometry
        )
        projected.B = B_saved
        return projected

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
        _, vacuum, _ = adapter.vendor_modules()
        grid_ops, _, _ = adapter.vendor_modules()
        grid = grid_ops.SphericalCellGrid(resolution, r_max)
        init = build_initial_data(
            grid,
            vacuum.flat_state,
            amplitude=amplitude,
            width=width,
            D_amplitude=D_amplitude,
            include_radiation=include_radiation,
        )
        geometry = self._enforce_moving_center_regularity(
            grid, init.geometry
        )
        # build_initial_data already supplies alpha=1, beta=0 and the moving
        # gauge seed B=3/4 Lambda. Do not re-solve CMC here.
        geometry.assert_finite_positive()
        return ProductionState(
            grid=grid,
            geometry=geometry,
            scalars=init.scalars,
            matter=init.matter,
            e_folds=0.0,
        )

    def step(self, state: ProductionState, dt: float) -> ProductionState:
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")

        grid = state.grid
        g0 = self._enforce_moving_center_regularity(
            grid, state.geometry.copy()
        )
        s0 = state.scalars
        m0 = state.matter

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

        # Vendor moving-puncture PIRK2 explicit block:
        # a,b,X,alpha,beta are explicit; Aa,K are partially implicit.
        u_names = ("a", "b", "X", "alpha", "beta")
        g_explicit1 = g0.copy()
        for name in u_names:
            setattr(
                g_explicit1,
                name,
                getattr(g0, name) + dt * gterms0["explicit"][name],
            )
        g_explicit1 = self._enforce_moving_center_regularity(
            grid, g_explicit1
        )

        s1 = ScalarFields(
            *(
                getattr(s0, name) + dt * getattr(srhs0, name)
                for name in ("S", "PS", "D", "PD", "phi", "Pi")
            )
        )
        m1 = V55MovingPIRKKernel._advance_matter(m0, mrhs0, dt)

        # First PIRK primary stage: use the updated explicit block but the
        # previous Aa/K values, exactly as in the pinned two-stage ordering.
        gterms_pred = adapter.geometry_stage_terms(
            grid, g_explicit1, s1, m1, lambda_m=LAMBDA_M
        )
        g1 = g_explicit1.copy()
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
        g1.B = g0.B + 0.75 * (g1.Lambda - g0.Lambda)
        g1 = self._enforce_moving_center_regularity(grid, g1)
        self._apply_outer_light_boundary(grid, g1, s1, m1)
        g1.assert_finite_positive()

        stage1 = ProductionState(
            grid=grid,
            geometry=g1,
            scalars=s1,
            matter=m1,
            t=state.t,
            tau=state.tau,
            e_folds=state.e_folds,
        )
        srhs1, mrhs1, _ = self._rhs(stage1)
        gterms1 = adapter.geometry_stage_terms(
            grid, g1, s1, m1, lambda_m=LAMBDA_M
        )

        # Second explicit PIRK block.
        gnew = g0.copy()
        for name in u_names:
            setattr(
                gnew,
                name,
                getattr(g0, name) + 0.5 * dt * (
                    gterms0["explicit"][name] + gterms1["explicit"][name]
                ),
            )

        snew = ScalarFields(
            *(
                getattr(s0, name) + 0.5 * dt * (
                    getattr(srhs0, name) + getattr(srhs1, name)
                )
                for name in ("S", "PS", "D", "PD", "phi", "Pi")
            )
        )
        mnew = V55MovingPIRKKernel._advance_matter_trapezoid(
            m0, mrhs0, mrhs1, dt
        )

        # Final implicit primary block, following the pinned PIRK2 sequence.
        final_u_old_v = g1.copy()
        for name in u_names:
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
        ll2_final = vacuum.lambda_l2_rhs(
            grid, final_primary, lambda_m=LAMBDA_M
        )
        gnew.Lambda = g0.Lambda + 0.5 * dt * (
            gterms0["lambda_l2"]
            + ll2_final
            + gterms0["lambda_l3"]
            + gterms1["lambda_l3"]
        )
        gnew.B = g0.B + 0.75 * (gnew.Lambda - g0.Lambda)

        gnew = self._enforce_moving_center_regularity(grid, gnew)
        self._apply_outer_light_boundary(grid, gnew, snew, mnew)
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
            state.history[-1]["H_eff"] if state.history else obs["H_eff"]
        )
        candidate.tau += dt * obs["tau_rate"]
        event = candidate.cycle.observe(
            len(candidate.history),
            candidate.t,
            candidate.tau,
            previous_H,
            obs["H_eff"],
        )

        S_t = float(srhs1.S[0])
        D_t = float(srhs1.D[0])
        hp = handoff_from_ledger(
            candidate.t,
            candidate.tau,
            obs,
            float(candidate.scalars.S[0]),
            S_t,
            float(candidate.scalars.D[0]),
            D_t,
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
        candidate.handoffs.append(hp)
        return candidate

    @staticmethod
    def _advance_matter(base, rhs, dt):
        return type(base)(
            *(
                V55MovingPIRKKernel._add_matter(
                    getattr(base, name), rhs[name], dt
                )
                for name in ("dark_matter", "baryons", "radiation")
            )
        )

    @staticmethod
    def _advance_matter_trapezoid(base, rhs0, rhs1, dt):
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

    @staticmethod
    def _add_matter(base, rhs, factor):
        return ConservedSpecies(
            base.rest + factor * rhs.rest,
            base.energy_t + factor * rhs.energy_t,
            base.momentum_r + factor * rhs.momentum_r,
        )
