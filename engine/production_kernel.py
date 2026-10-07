# Campaign diagnostic build: CMC turnaround response instrumentation remains diagnostic-only.
# Deep worldtube causal audit: radiation transport/source and COSMOS characteristics.
"""Integrated V5.5 production state and synchronized PIRK2 evolution.

The machine carries one global spherical spacetime with:
- reference-metric BSSN geometry;
- local S/D scalar fields;
- corrected COSMOS scalar;
- conservative beta-coupled dark matter;
- conservative baryonic dust;
- conservative radiation p=rho/3;
- derived proper time, trapping, Misner-Sharp, and cycle/handoff ledgers.

The external vendor supplies numerical geometry operators/PIRK ordering only.
All V5.5 matter equations and source terms remain local to this repository.

No bounce, reset, branch flip, shell source, fitted feedback coefficient,
lapse floor, ejection threshold, or physical stop is present.
"""
from dataclasses import dataclass, field
from pathlib import Path
import math
import numpy as np

from .handoff import CycleLedger, handoff_from_ledger
from .cmc_gauge import solve_cmc_lapse
from .production_contract import KernelCapabilities
from .invariant_diagnostics import (
    misner_sharp_mass_from_chi,
    trapping_indicator_from_areal_radius,
)
from .cosmology_observables import effective_hubble, append_efolds
from .stress_energy import assemble_total_stress_energy
from .matter_rhs import species_rhs
from .matter_system import (
    Species,
    ConservedSpecies,
    SphericalMetric,
    geometry_metric_derivatives,
    _metric_arrays,
    primitives as matter_primitives,
    _reconstructed_primitive,
    _hll_flux,
    _valencia_flux,
    _valencia_source,
)
from .scalar_system import ScalarFields, scalar_rhs_arrays
from .v55_initial import build_initial_data
from .v55_matter import (
    V55MatterState,
    dm_density,
    metric_slice_from_q,
    total_matter_projection,
)
from . import v55_pirk_adapter as adapter


BETA_DM = -0.04
LAMBDA_M = 2.0


@dataclass
class ProductionState:
    grid: object
    geometry: object
    scalars: ScalarFields
    matter: V55MatterState
    t: float = 0.0
    tau: float = 0.0
    e_folds: float = 0.0
    cycle: CycleLedger = field(default_factory=CycleLedger)
    history: list[dict] = field(default_factory=list)
    handoffs: list = field(default_factory=list)


class V55ProductionKernel:
    """Unvalidated production-facing implementation of the full state graph."""

    capabilities = KernelCapabilities(
        reference_metric_bssn=True,
        pirk2=True,
        moving_gauge=True,
        dm_momentum=True,
        baryon_momentum=True,
        radiation=True,
        invariant_trapping=True,
        misner_sharp_current=True,
    )
    validation_state = "implemented_not_campaign_validated"

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
        init = build_initial_data(
            grid,
            vacuum.flat_state,
            amplitude=amplitude,
            width=width,
            D_amplitude=D_amplitude,
            include_radiation=include_radiation,
        )
        # Initialize the production slice on the same stage-aware CMC gauge
        # that is carried throughout evolution.
        init.geometry = self._enforce_center_regularity(
            grid, init.geometry
        )
        init.geometry.alpha = self._solve_lapse(
            grid, init.geometry, init.scalars, init.matter
        )[0]
        init.geometry.beta.fill(0.0)
        init.geometry.B.fill(0.0)
        return ProductionState(
            grid=grid,
            geometry=init.geometry,
            scalars=init.scalars,
            matter=init.matter,
            e_folds=0.0,
        )

    @staticmethod
    def _d1(grid, values, parity):
        return grid.cell_derivative_fourth(values, parity=parity)

    @staticmethod
    def _solve_lapse(grid, geometry, scalars, matter):
        return solve_cmc_lapse(grid, geometry, scalars, matter)

    @staticmethod
    def _enforce_center_regularity(grid, geometry):
        """Project only exact algebraic/spherical-centre BSSN identities.

        The pinned reference kernel defines this as a numerical regularity
        projection: a*b**2=1 everywhere, a/b=1+O(r**2), Aa=O(r**2), and odd
        connection/shift variables are O(r) at the first cell. No field
        equation, source, gauge target, or physical outcome is introduced.
        """
        _, vacuum, _ = adapter.vendor_modules()
        return vacuum.enforce_algebraic_regularity(grid, geometry)

    def _rhs(self, state: ProductionState):
        metric = metric_slice_from_q(state.grid, state.geometry)

        rho_dm = dm_density(metric, state.matter)
        srhs = scalar_rhs_arrays(
            state.grid,
            state.geometry,
            state.scalars,
            beta_dm=BETA_DM,
            rho_dm=rho_dm,
        )

        md = geometry_metric_derivatives(
            metric,
            lambda values, parity: self._d1(state.grid, values, parity),
        )
        dphi_t = srhs.phi
        dphi_r = self._d1(state.grid, state.scalars.phi, 1)

        mrhs = {}
        for name, species in (
            ("dark_matter", Species.DARK_MATTER),
            ("baryons", Species.BARYON),
            ("radiation", Species.RADIATION),
        ):
            mrhs[name] = species_rhs(
                metric,
                md,
                getattr(state.matter, name),
                species,
                dphi_t=dphi_t if species is Species.DARK_MATTER else None,
                dphi_r=dphi_r if species is Species.DARK_MATTER else None,
                beta_dm=BETA_DM,
            )
        return srhs, mrhs, md

    @staticmethod
    def _apply_outer_light_boundary(grid, geometry, scalars, matter) -> None:
        """Apply the pinned R0 incoming-light constraint at the finite radius.

        This is a boundary-condition operation only. It uses the already
        evolved total stress-energy projections; no fitted boundary parameter
        or new physical source is introduced.
        """
        adapter.vendor_modules()
        from bssn_characteristic_boundary import apply_light_constraint_boundary

        total = assemble_total_stress_energy(
            grid, geometry, scalars, matter
        )
        apply_light_constraint_boundary(
            grid, geometry, total.rho, total.j
        )

    @staticmethod
    def _add_matter(base: ConservedSpecies, rhs, factor: float) -> ConservedSpecies:
        return ConservedSpecies(
            base.rest + factor * rhs.rest,
            base.energy_t + factor * rhs.energy_t,
            base.momentum_r + factor * rhs.momentum_r,
        )

    def step(self, state: ProductionState, dt: float) -> ProductionState:
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")

        grid = state.grid
        g0 = self._enforce_center_regularity(
            grid, state.geometry.copy()
        )
        s0 = state.scalars
        m0 = state.matter

        # The lapse is elliptically determined by the current CMC slice.
        # Spatial shift is zero in the spherical CMC branch.
        g0.alpha = self._solve_lapse(
            grid, g0, s0, m0
        )[0]
        g0.beta.fill(0.0)
        g0.B.fill(0.0)

        srhs0, mrhs0, _ = self._rhs(
            ProductionState(
                grid=grid, geometry=g0, scalars=s0, matter=m0,
                t=state.t, tau=state.tau, e_folds=state.e_folds,
            )
        )
        gterms0 = adapter.geometry_stage_terms(
            grid, g0, s0, m0, lambda_m=LAMBDA_M
        )

        # TRUE-PIRK predictor with CMC: advance only the metric's
        # explicit variables. Alpha and shift are gauge variables here and are
        # resolved from the predictor slice rather than integrated by 1+log.
        g_explicit1 = g0.copy()
        for name in ("a", "b", "X"):
            setattr(
                g_explicit1,
                name,
                getattr(g0, name) + dt * gterms0["explicit"][name],
            )
        g_explicit1 = self._enforce_center_regularity(
            grid, g_explicit1
        )
        g_explicit1.beta.fill(0.0)
        g_explicit1.B.fill(0.0)

        s1 = ScalarFields(
            *(
                getattr(s0, name) + dt * getattr(srhs0, name)
                for name in ("S", "PS", "D", "PD", "phi", "Pi")
            )
        )
        m1 = V55MatterState(
            *(
                self._add_matter(getattr(m0, name), mrhs0[name], dt)
                for name in ("dark_matter", "baryons", "radiation")
            )
        )

        g_explicit1.alpha = self._solve_lapse(
            grid, g_explicit1, s1, m1
        )[0]
        gterms_pred = adapter.geometry_stage_terms(
            grid, g_explicit1, s1, m1, lambda_m=LAMBDA_M
        )
        _, vacuum, _ = adapter.vendor_modules()

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
        g1 = self._enforce_center_regularity(grid, g1)
        # The pinned R0 radial system supplies a parameter-free incoming-light
        # characteristic reconstruction for Aa at the finite outer worldtube.
        # Apply it only after Lambda is current, as in the reference PIRK/CPBC
        # sequence. The CMC lapse is then re-solved on the conditioned slice.
        self._apply_outer_light_boundary(grid, g1, s1, m1)
        # Resolve the CMC lapse on the full primary predictor, then
        # use that gauge state for the second split evaluation.
        g1.alpha = self._solve_lapse(
            grid, g1, s1, m1
        )[0]
        g1.beta.fill(0.0)
        g1.B.fill(0.0)
        g1.assert_finite_positive()

        # Re-evaluate the complete split on the primary-updated predictor.
        stage1 = ProductionState(
            grid=grid, geometry=g1, scalars=s1, matter=m1,
            t=state.t, tau=state.tau, e_folds=state.e_folds,
        )
        srhs1, mrhs1, _ = self._rhs(stage1)
        gterms1 = adapter.geometry_stage_terms(
            grid, g1, s1, m1, lambda_m=LAMBDA_M
        )

        # Final explicit block uses the TRUE-PIRK trapezoidal pairing with the
        # second split evaluation above. CMC supplies alpha after the primary
        # state is formed; it is not an evolved 1+log variable.
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

        # Regularize the explicit block before the final L2 evaluation: the
        # spherical-center identities are part of the numerical state domain.
        gnew = self._enforce_center_regularity(grid, gnew)

        snew = ScalarFields(
            *(
                getattr(s0, name) + 0.5 * dt * (
                    getattr(srhs0, name) + getattr(srhs1, name)
                )
                for name in ("S", "PS", "D", "PD", "phi", "Pi")
            )
        )
        mnew = V55MatterState(
            *(
                ConservedSpecies(
                    getattr(m0, name).rest + 0.5 * dt * (
                        mrhs0[name].rest + mrhs1[name].rest
                    ),
                    getattr(m0, name).energy_t + 0.5 * dt * (
                        mrhs0[name].energy_t + mrhs1[name].energy_t
                    ),
                    getattr(m0, name).momentum_r + 0.5 * dt * (
                        mrhs0[name].momentum_r + mrhs1[name].momentum_r
                    ),
                )
                for name in ("dark_matter", "baryons", "radiation")
            )
        )

        # Final implicit evaluation is taken on the final explicit state,
        # preserving the archived PIRK partition.
        final_u_old_v = g1.copy()
        for name in gterms0["explicit"]:
            setattr(final_u_old_v, name, getattr(gnew, name))

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
        # Apply the same incoming-light reconstruction after Lambda is current
        # on the completed final slice. This replaces the unconstrained finite-
        # radius boundary mode without altering the interior equations.
        gnew = self._enforce_center_regularity(grid, gnew)
        self._apply_outer_light_boundary(grid, gnew, snew, mnew)
        # Final stage-aware CMC solve on the completed final A/K/matter slice.
        gnew.alpha = self._solve_lapse(
            grid, gnew, snew, mnew
        )[0]
        gnew.beta.fill(0.0)
        gnew.B.fill(0.0)

        candidate = ProductionState(
            grid=grid, geometry=gnew, scalars=snew, matter=mnew,
            t=state.t + dt, tau=state.tau, e_folds=state.e_folds,
            cycle=state.cycle, history=list(state.history),
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
            len(candidate.history), candidate.t, candidate.tau,
            previous_H, obs["H_eff"],
        )

        S_t = float(srhs1.S[0])
        D_t = float(srhs1.D[0])
        hp = handoff_from_ledger(
            candidate.t, candidate.tau, obs,
            float(candidate.scalars.S[0]), S_t,
            float(candidate.scalars.D[0]), D_t,
        )

        obs = {
            key: value for key, value in obs.items()
            if key not in (
                "total_rho", "total_pr", "total_pt", "total_j",
                "areal_radius", "chi", "misner_sharp"
            )
        }
        obs["t"] = candidate.t
        obs["tau"] = candidate.tau
        obs["cycle_event"] = event.kind if event else None
        candidate.history.append(obs)
        candidate.handoffs.append(hp)
        return candidate

    def diagnostics(self, state: ProductionState, *, profiles: bool = False) -> dict:
        grid = state.grid
        geom = state.geometry

        _, vacuum, moving = adapter.vendor_modules()
        metric = metric_slice_from_q(grid, geom)
        total = assemble_total_stress_energy(
            grid, geom, state.scalars, state.matter
        )
        raw = vacuum.constraints(grid, geom)

        H = raw["hamiltonian"] - 16.0 * math.pi * total.rho
        M = raw["momentum"] - 8.0 * math.pi * total.j

        r = np.asarray(grid.centers)
        R = r * np.sqrt(geom.b) / geom.X
        Rr = grid.cell_derivative_fourth(R, parity=1)
        Ktheta = geom.K / 3.0 - geom.Aa / 2.0
        normal_dR = -R * Ktheta
        chi = trapping_indicator_from_areal_radius(
            Rr,
            normal_dR,
            geom.X**2 / geom.a,
        )
        mass = misner_sharp_mass_from_chi(R, chi)

        explicit = moving.moving_puncture_explicit_rhs(grid, geom)
        Rdot = R * (
            0.5 * explicit["b"] / geom.b
            - explicit["X"] / geom.X
        )

        # Turnaround diagnostic: expose the CMC projection and radial lapse
        # response on the same slice. Diagnostic only; no evolution change.
        from .cmc_gauge import target_kdot
        cmc_kdot = target_kdot(grid, geom, state.scalars, state.matter)
        alpha_r = grid.cell_derivative_fourth(geom.alpha, parity=1)
        i_alpha = int(np.argmax(np.abs(alpha_r)))

        # Radiation worldtube audit: expose the actual outer characteristic
        # state and compare the physical outer-face flux to the HLL flux on the
        # last interior face. Diagnostic only; the evolution closure is
        # unchanged.
        rad_metrics = _metric_arrays(metric)
        rad_prim = matter_primitives(
            rad_metrics, state.matter.radiation, Species.RADIATION
        )
        q_rad = rad_prim[-1]
        vhat_rad = math.sqrt(rad_metrics[-1].gamma_rr) * q_rad.v_r
        cs_rad = 1.0 / math.sqrt(3.0)
        s_minus_rad = (
            geom.alpha[-1] * (vhat_rad - cs_rad)
            / (1.0 - vhat_rad * cs_rad)
            - geom.beta[-1]
        )
        s_plus_rad = (
            geom.alpha[-1] * (vhat_rad + cs_rad)
            / (1.0 + vhat_rad * cs_rad)
            - geom.beta[-1]
        )
        rad_outer_flux = _valencia_flux(rad_metrics[-1], q_rad)
        m0_rad, m1_rad = rad_metrics[-2], rad_metrics[-1]
        mf_rad = SphericalMetric(
            alpha=0.5 * (m0_rad.alpha + m1_rad.alpha),
            beta=0.5 * (m0_rad.beta + m1_rad.beta),
            gamma_rr=0.5 * (m0_rad.gamma_rr + m1_rad.gamma_rr),
            gamma_rr_inv=0.5 * (
                m0_rad.gamma_rr_inv + m1_rad.gamma_rr_inv
            ),
            gamma_thth=0.5 * (
                m0_rad.gamma_thth + m1_rad.gamma_thth
            ),
            gamma_thth_inv=0.5 * (
                m0_rad.gamma_thth_inv + m1_rad.gamma_thth_inv
            ),
            sqrt_gamma=0.5 * (
                m0_rad.sqrt_gamma + m1_rad.sqrt_gamma
            ),
        )
        ql_rad = _reconstructed_primitive(
            rad_prim, grid.n - 2, "right", mf_rad.gamma_rr
        )
        qr_rad = _reconstructed_primitive(
            rad_prim, grid.n - 1, "left", mf_rad.gamma_rr
        )
        rad_inner_flux = _hll_flux(
            mf_rad, ql_rad, qr_rad, Species.RADIATION
        )
        # One more interior face lets us distinguish a boundary closure
        # contribution from a genuine outgoing pulse entering the last cell.
        m2_rad, m3_rad = rad_metrics[-3], rad_metrics[-2]
        mf2_rad = SphericalMetric(
            alpha=0.5 * (m2_rad.alpha + m3_rad.alpha),
            beta=0.5 * (m2_rad.beta + m3_rad.beta),
            gamma_rr=0.5 * (m2_rad.gamma_rr + m3_rad.gamma_rr),
            gamma_rr_inv=0.5 * (m2_rad.gamma_rr_inv + m3_rad.gamma_rr_inv),
            gamma_thth=0.5 * (m2_rad.gamma_thth + m3_rad.gamma_thth),
            gamma_thth_inv=0.5 * (m2_rad.gamma_thth_inv + m3_rad.gamma_thth_inv),
            sqrt_gamma=0.5 * (m2_rad.sqrt_gamma + m3_rad.sqrt_gamma),
        )
        ql2_rad = _reconstructed_primitive(
            rad_prim, grid.n - 3, "right", mf2_rad.gamma_rr
        )
        qr2_rad = _reconstructed_primitive(
            rad_prim, grid.n - 2, "left", mf2_rad.gamma_rr
        )
        rad_penultimate_flux = _hll_flux(
            mf2_rad, ql2_rad, qr2_rad, Species.RADIATION
        )
        dr_rad = float(grid.dr)
        rad_outer_div_E = -(
            float(rad_outer_flux[1]) - float(rad_inner_flux[1])
        ) / dr_rad
        rad_outer_div_S = -(
            float(rad_outer_flux[2]) - float(rad_inner_flux[2])
        ) / dr_rad
        rad_penultimate_div_E = -(
            float(rad_inner_flux[1]) - float(rad_penultimate_flux[1])
        ) / dr_rad
        rad_penultimate_div_S = -(
            float(rad_inner_flux[2]) - float(rad_penultimate_flux[2])
        ) / dr_rad
        rad_md = geometry_metric_derivatives(
            metric,
            lambda values, parity: grid.cell_derivative_fourth(values, parity=parity),
        )
        rad_outer_source_E, rad_outer_source_S = _valencia_source(
            rad_metrics[-1], q_rad, float(geom.K[-1]), float(geom.Aa[-1]),
            float(rad_md.radial["alpha"][-1]),
            float(rad_md.radial["beta"][-1]),
            float(rad_md.radial["rr"][-1]),
            float(rad_md.radial["thth"][-1]),
        )
        rad_penultimate_source_E, rad_penultimate_source_S = _valencia_source(
            rad_metrics[-2], rad_prim[-2], float(geom.K[-2]), float(geom.Aa[-2]),
            float(rad_md.radial["alpha"][-2]),
            float(rad_md.radial["beta"][-2]),
            float(rad_md.radial["rr"][-2]),
            float(rad_md.radial["thth"][-2]),
        )
        # Dust worldtube audit: compare outer-cell transport divergence with
        # the geometric metric source. Dark matter additionally has the archived
        # scalar exchange term, which is intentionally kept separate here.
        dust_audits = {}
        for label, species, field in (
            ("dark_matter", Species.DARK_MATTER, state.matter.dark_matter),
            ("baryons", Species.BARYON, state.matter.baryons),
        ):
            prim = matter_primitives(rad_metrics, field, species)
            q_outer = prim[-1]
            outer_flux = _valencia_flux(rad_metrics[-1], q_outer)
            ql_d = _reconstructed_primitive(
                prim, grid.n - 2, "right", mf_rad.gamma_rr
            )
            qr_d = _reconstructed_primitive(
                prim, grid.n - 1, "left", mf_rad.gamma_rr
            )
            inner_flux = _hll_flux(
                mf_rad, ql_d, qr_d, species
            )
            div_d = -(
                float(outer_flux[1]) - float(inner_flux[1])
            ) / dr_rad, -(
                float(outer_flux[2]) - float(inner_flux[2])
            ) / dr_rad
            geom_src = _valencia_source(
                rad_metrics[-1], q_outer, float(geom.K[-1]), float(geom.Aa[-1]),
                float(rad_md.radial["alpha"][-1]),
                float(rad_md.radial["beta"][-1]),
                float(rad_md.radial["rr"][-1]),
                float(rad_md.radial["thth"][-1]),
            )
            dust_audits[label] = {
                "div_E": div_d[0],
                "div_S": div_d[1],
                "geom_source_E": float(geom_src[0]),
                "geom_source_S": float(geom_src[1]),
                "transport_to_geom_source_E": abs(div_d[0]) / max(abs(float(geom_src[0])), 1.0e-300),
                "transport_to_geom_source_S": abs(div_d[1]) / max(abs(float(geom_src[1])), 1.0e-300),
                "vhat": math.sqrt(rad_metrics[-1].gamma_rr) * q_outer.v_r,
            }

        # COSMOS scalar characteristics at the same worldtube.
        # With the evolved equations, W_in = Pi + c*phi_r has speed
        # -beta-alpha*c, while W_out = Pi - c*phi_r has speed
        # -beta+alpha*c. Diagnostic only; no boundary condition is changed.
        phi_r = grid.cell_derivative_fourth(state.scalars.phi, parity=1)
        phi_char_in = state.scalars.Pi[-1] + math.sqrt(rad_metrics[-1].gamma_rr_inv) * phi_r[-1]
        phi_char_out = state.scalars.Pi[-1] - math.sqrt(rad_metrics[-1].gamma_rr_inv) * phi_r[-1]
        phi_speed_in = -rad_metrics[-1].beta - rad_metrics[-1].alpha * math.sqrt(rad_metrics[-1].gamma_rr_inv)
        phi_speed_out = -rad_metrics[-1].beta + rad_metrics[-1].alpha * math.sqrt(rad_metrics[-1].gamma_rr_inv)
        flux_e_scale = max(
            abs(float(rad_inner_flux[1])),
            abs(float(rad_outer_flux[1])),
            1.0e-300,
        )
        flux_s_scale = max(
            abs(float(rad_inner_flux[2])),
            abs(float(rad_outer_flux[2])),
            1.0e-300,
        )
        rad_E_outer = (
            state.matter.radiation.energy_t[-1]
            / rad_metrics[-1].sqrt_gamma
        )
        rad_S_outer = (
            state.matter.radiation.momentum_r[-1]
            / rad_metrics[-1].sqrt_gamma
        )
        rad_Sabs_outer = math.sqrt(
            max(
                rad_metrics[-1].gamma_rr_inv * rad_S_outer * rad_S_outer,
                0.0,
            )
        )
        rad_ratio_outer = rad_Sabs_outer / max(abs(rad_E_outer), 1.0e-300)

        # CMC boundary audit: evaluate the elliptic equation on the solved
        # lapse through the last interior cell. Diagnostic only; no evolution
        # quantity or boundary condition is changed here.
        inv_cmc = geom.X**2 / geom.a
        a_r_cmc = grid.cell_derivative_fourth(geom.a, parity=1)
        b_r_cmc = grid.cell_derivative_fourth(geom.b, parity=1)
        X_r_cmc = grid.cell_derivative_fourth(geom.X, parity=1)
        alpha_rr_cmc = grid.cell_second_derivative_fourth(geom.alpha, parity=1)
        c_cmc = (
            -0.5 * a_r_cmc / geom.a
            + b_r_cmc / geom.b
            - X_r_cmc / geom.X
            + 2.0 / r
        )
        Q_cmc = (
            1.5 * geom.Aa * geom.Aa
            + geom.K * geom.K / 3.0
            + 4.0 * math.pi * (
                total.rho + total.pr + 2.0 * total.pt
            )
        )
        cmc_residual = (
            alpha_rr_cmc
            + c_cmc * alpha_r
            - (Q_cmc / inv_cmc) * geom.alpha
            + cmc_kdot / inv_cmc
        )
        outer_cmc = r >= 0.8 * grid.r_max
        cmc_last_interior = grid.n - 2

        # Geometry-collapse witness.  This is diagnostic only: it does not
        # modify any evolved variable.  In a unit-determinant conformal
        # metric, a*b^2=1 and therefore da/a + 2 db/b must vanish.
        det = geom.a * geom.b**2
        conformal_trace_rhs = (
            explicit["a"] / geom.a
            + 2.0 * explicit["b"] / geom.b
        )

        # Split the primary curvature source into the pinned vacuum operator
        # and the local total stress-energy contribution.  This exposes whether
        # the first departure is already present in vacuum geometry or enters
        # through the matter projection, without introducing a new source.
        vacuum_l3 = vacuum.primary_l3_rhs(grid, geom)
        vacuum_l2 = vacuum.primary_l2_rhs(grid, geom)
        matter_l3 = {
            key: np.asarray(adapter.primary_l3_with_matter(
                grid, geom, state.scalars, state.matter
            )[key]) - np.asarray(vacuum_l3[key])
            for key in ("Aa", "K")
        }

        surface = int(np.argmin(np.abs(r - 10.0)))
        i_sigma = surface
        outer = r >= 0.8 * grid.r_max
        coordinate_energy_flux = 4.0 * math.pi * R[surface]**2 * (
            geom.alpha[surface] * (geom.X[surface]**2 / geom.a[surface])
            * total.j[surface]
            - geom.beta[surface] * total.rho[surface]
        )
        flux_T = 4.0 * math.pi * R[surface]**2 * (
            -geom.alpha[surface] * (geom.X[surface]**2 / geom.a[surface])
            * total.j[surface] * Rr[surface]
        )
        work_pR = -4.0 * math.pi * R[surface]**2 * (
            total.pr[surface] * Rdot[surface]
        )

        # Spatial witnesses for the first runaway source. Diagnostic only:
        # record where the vacuum and explicit geometry RHS are largest and
        # expose the local slice state there so the next run can distinguish
        # a boundary, center, or interior operator failure.
        def _witness(arr):
            values = np.asarray(arr, dtype=float)
            idx = int(np.argmax(np.abs(values)))
            return idx, float(r[idx]), float(values[idx])

        i_vAa, r_vAa, v_vAa = _witness(vacuum_l3["Aa"])
        i_vK, r_vK, v_vK = _witness(vacuum_l3["K"])
        i_l2Aa, r_l2Aa, v_l2Aa = _witness(vacuum_l2["Aa"])
        i_l2K, r_l2K, v_l2K = _witness(vacuum_l2["K"])
        i_ea, r_ea, v_ea = _witness(explicit["a"] / geom.a)
        i_eb, r_eb, v_eb = _witness(explicit["b"] / geom.b)
        i_conn, r_conn, v_conn = _witness(raw["connection"])
        i_H, r_H, v_H = _witness(H)

        roots = []
        sign_change = chi[:-1] * chi[1:] <= 0.0
        for i in np.where(sign_change)[0]:
            if chi[i] != chi[i + 1]:
                roots.append(
                    float(
                        r[i] - chi[i] * (r[i + 1] - r[i])
                        / (chi[i + 1] - chi[i])
                    )
                )

        H_eff = effective_hubble(geom, grid.volumes)
        H_center = -float(geom.K[0]) / 3.0
        H_coord_center = float(geom.alpha[0]) * H_center
        out = {
            "H_eff": H_eff,
            "H_center": H_center,
            "H_coord_center": H_coord_center,
            "H_eff_over_H_center": (
                H_eff / H_center if abs(H_center) > 1.0e-300 else float("nan")
            ),
            "e_folds": float(state.e_folds),
            "tau_rate": float(geom.alpha[0]),
            "R_sigma": float(R[surface]),
            "M_MS": float(mass[surface]),
            "chi_sigma": float(chi[surface]),
            "trapping_min": float(np.min(chi)),
            "trapped_roots": roots,
            "coordinate_energy_flux": float(coordinate_energy_flux),
            "flux_T": float(flux_T),
            "work_pR": float(work_pR),
            "rho_outer": float(np.mean(total.rho[outer])),
            "p_outer": float(np.mean(total.pr[outer])),
            "j_outer": float(np.mean(total.j[outer])),
            "rho_total_max": float(np.max(total.rho)),
            "hamiltonian_max": float(np.max(np.abs(H[2:]))),
            "momentum_max": float(np.max(np.abs(M[2:]))),
            "connection_max": float(
                np.max(np.abs(raw["connection"][2:]))
            ),
            "determinant_min": float(np.min(det)),
            "determinant_constraint_max": float(
                np.max(np.abs(det - 1.0))
            ),
            "conformal_trace_rhs_max": float(np.max(np.abs(conformal_trace_rhs))),
            "explicit_da_over_a_max": float(np.max(np.abs(explicit["a"] / geom.a))),
            "explicit_db_over_b_max": float(np.max(np.abs(explicit["b"] / geom.b))),
            "vacuum_l3_Aa_max": float(np.max(np.abs(vacuum_l3["Aa"]))),
            "matter_l3_Aa_max": float(np.max(np.abs(matter_l3["Aa"]))),
            "vacuum_l3_K_max": float(np.max(np.abs(vacuum_l3["K"]))),
            "matter_l3_K_max": float(np.max(np.abs(matter_l3["K"]))),
            "lapse_min": float(np.min(geom.alpha)),
            "lapse_max": float(np.max(geom.alpha)),
            "Rdot_sigma": float(Rdot[surface]),
            "cmc_kdot": float(cmc_kdot),
            "cmc_residual_outer_max": float(
                np.max(np.abs(cmc_residual[outer_cmc & (np.arange(grid.n) < grid.n - 1)]))
            ),
            "cmc_residual_last_interior": float(cmc_residual[cmc_last_interior]),
            "cmc_residual_penultimate_interior": float(cmc_residual[max(1, grid.n - 3)]),
            "cmc_boundary_alpha": float(geom.alpha[-1]),
            "cmc_boundary_alpha_prev": float(geom.alpha[-2]),
            "cmc_boundary_alpha_drop_per_dr": float(
                (geom.alpha[-1] - geom.alpha[-2]) / grid.dr
            ),
            "cmc_boundary_alpha_r": float(alpha_r[-1]),
            "alpha_r_max": float(np.max(np.abs(alpha_r))),
            "alpha_r_max_r": float(r[i_alpha]),
            "alpha_r_sigma": float(alpha_r[i_sigma]),
            "alpha_sigma": float(geom.alpha[i_sigma]),
            "radiation_outer_ratio": float(rad_ratio_outer),
            "radiation_outer_vhat": float(vhat_rad),
            "radiation_outer_s_minus": float(s_minus_rad),
            "radiation_outer_s_plus": float(s_plus_rad),
            "radiation_outer_all_outgoing": float(s_minus_rad >= 0.0),
            "radiation_boundary_energy_flux": float(rad_outer_flux[1]),
            "radiation_inner_face_energy_flux": float(rad_inner_flux[1]),
            "radiation_boundary_momentum_flux": float(rad_outer_flux[2]),
            "radiation_inner_face_momentum_flux": float(rad_inner_flux[2]),
            "radiation_boundary_energy_flux_rel_mismatch": float(
                abs(float(rad_outer_flux[1] - rad_inner_flux[1]))
                / flux_e_scale
            ),
            "radiation_boundary_momentum_flux_rel_mismatch": float(
                abs(float(rad_outer_flux[2] - rad_inner_flux[2]))
                / flux_s_scale
            ),
            "radiation_outer_div_E": float(rad_outer_div_E),
            "radiation_outer_div_S": float(rad_outer_div_S),
            "radiation_penultimate_div_E": float(rad_penultimate_div_E),
            "radiation_penultimate_div_S": float(rad_penultimate_div_S),
            "radiation_outer_source_E": float(rad_outer_source_E),
            "radiation_outer_source_S": float(rad_outer_source_S),
            "radiation_penultimate_source_E": float(rad_penultimate_source_E),
            "radiation_penultimate_source_S": float(rad_penultimate_source_S),
            "radiation_outer_transport_to_source_E": float(
                abs(rad_outer_div_E) / max(abs(rad_outer_source_E), 1.0e-300)
            ),
            "radiation_outer_transport_to_source_S": float(
                abs(rad_outer_div_S) / max(abs(rad_outer_source_S), 1.0e-300)
            ),
            "cosmos_phi_outer_char_in": float(phi_char_in),
            "cosmos_phi_outer_char_out": float(phi_char_out),
            "cosmos_phi_outer_speed_in": float(phi_speed_in),
            "cosmos_phi_outer_speed_out": float(phi_speed_out),
            "cosmos_phi_outer_incoming": float(phi_speed_in < 0.0),
            "cosmos_phi_outer_outgoing": float(phi_speed_out > 0.0),
            "cosmos_phi_outer_abs_out_in_ratio": float(
                abs(phi_char_out) / max(abs(phi_char_in), 1.0e-300)
            ),
            "dark_matter_outer_div_E": float(dust_audits["dark_matter"]["div_E"]),
            "dark_matter_outer_div_S": float(dust_audits["dark_matter"]["div_S"]),
            "dark_matter_outer_geom_source_E": float(dust_audits["dark_matter"]["geom_source_E"]),
            "dark_matter_outer_geom_source_S": float(dust_audits["dark_matter"]["geom_source_S"]),
            "dark_matter_outer_transport_to_geom_source_E": float(dust_audits["dark_matter"]["transport_to_geom_source_E"]),
            "dark_matter_outer_transport_to_geom_source_S": float(dust_audits["dark_matter"]["transport_to_geom_source_S"]),
            "dark_matter_outer_vhat": float(dust_audits["dark_matter"]["vhat"]),
            "baryons_outer_div_E": float(dust_audits["baryons"]["div_E"]),
            "baryons_outer_div_S": float(dust_audits["baryons"]["div_S"]),
            "baryons_outer_geom_source_E": float(dust_audits["baryons"]["geom_source_E"]),
            "baryons_outer_geom_source_S": float(dust_audits["baryons"]["geom_source_S"]),
            "baryons_outer_transport_to_geom_source_E": float(dust_audits["baryons"]["transport_to_geom_source_E"]),
            "baryons_outer_transport_to_geom_source_S": float(dust_audits["baryons"]["transport_to_geom_source_S"]),
            "baryons_outer_vhat": float(dust_audits["baryons"]["vhat"]),
            "K_sigma": float(geom.K[i_sigma]),
            "Aa_sigma": float(geom.Aa[i_sigma]),
            "phi_outer": float(np.mean(state.scalars.phi[outer])),
            "Pi_outer": float(np.mean(state.scalars.Pi[outer])),
            "vacuum_l3_Aa_max_idx": i_vAa,
            "vacuum_l3_Aa_max_r": r_vAa,
            "vacuum_l3_Aa_at_max": v_vAa,
            "vacuum_l3_Aa_state_a": float(geom.a[i_vAa]),
            "vacuum_l3_Aa_state_b": float(geom.b[i_vAa]),
            "vacuum_l3_Aa_state_X": float(geom.X[i_vAa]),
            "vacuum_l3_Aa_state_Aa": float(geom.Aa[i_vAa]),
            "vacuum_l3_Aa_state_K": float(geom.K[i_vAa]),
            "vacuum_l3_Aa_state_alpha": float(geom.alpha[i_vAa]),
            "vacuum_l3_K_max_idx": i_vK,
            "vacuum_l3_K_max_r": r_vK,
            "vacuum_l3_K_at_max": v_vK,
            "vacuum_l3_K_state_a": float(geom.a[i_vK]),
            "vacuum_l3_K_state_b": float(geom.b[i_vK]),
            "vacuum_l3_K_state_X": float(geom.X[i_vK]),
            "vacuum_l3_K_state_Aa": float(geom.Aa[i_vK]),
            "vacuum_l3_K_state_K": float(geom.K[i_vK]),
            "vacuum_l3_K_state_alpha": float(geom.alpha[i_vK]),
            "vacuum_l2_Aa_max_r": r_l2Aa,
            "vacuum_l2_Aa_at_max": v_l2Aa,
            "vacuum_l2_K_max_r": r_l2K,
            "vacuum_l2_K_at_max": v_l2K,
            "explicit_da_over_a_max_r": r_ea,
            "explicit_da_over_a_at_max": v_ea,
            "explicit_db_over_b_max_r": r_eb,
            "explicit_db_over_b_at_max": v_eb,
            "connection_max_r": r_conn,
            "connection_at_max": v_conn,
            "hamiltonian_max_r": r_H,
            "hamiltonian_at_max": v_H,
        }

        if profiles:
            out.update({
                "total_rho": total.rho,
                "total_pr": total.pr,
                "total_pt": total.pt,
                "total_j": total.j,
                "areal_radius": R,
                "chi": chi,
                "misner_sharp": mass,
            })
        return out

    def checkpoint(self, state: ProductionState, path: str | Path) -> None:
        write_checkpoint_npz(state, path)


def write_checkpoint_npz(
    state: ProductionState,
    path: str | Path,
) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    arrays = {
        "t": np.asarray(state.t),
        "tau": np.asarray(state.tau),
        "e_folds": np.asarray(state.e_folds),
        "r": np.asarray(state.grid.centers),
        "a": state.geometry.a,
        "b": state.geometry.b,
        "X": state.geometry.X,
        "alpha": state.geometry.alpha,
        "beta": state.geometry.beta,
        "Aa": state.geometry.Aa,
        "K": state.geometry.K,
        "Lambda": state.geometry.Lambda,
        "B": state.geometry.B,
        "S": state.scalars.S,
        "PS": state.scalars.PS,
        "D": state.scalars.D,
        "PD": state.scalars.PD,
        "phi": state.scalars.phi,
        "Pi": state.scalars.Pi,
        "dm_rest": state.matter.dark_matter.rest,
        "dm_energy_t": state.matter.dark_matter.energy_t,
        "dm_momentum_r": state.matter.dark_matter.momentum_r,
        "baryon_rest": state.matter.baryons.rest,
        "baryon_energy_t": state.matter.baryons.energy_t,
        "baryon_momentum_r": state.matter.baryons.momentum_r,
        "radiation_rest": state.matter.radiation.rest,
        "radiation_energy_t": state.matter.radiation.energy_t,
        "radiation_momentum_r": state.matter.radiation.momentum_r,
    }
    np.savez_compressed(out, **arrays)
