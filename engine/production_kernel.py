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
    geometry_metric_derivatives,
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
GAUGE_ETA = 2.0  # pinned vendor moving-puncture damping; not a fitted physics parameter


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
        g0 = state.geometry.copy()
        s0 = state.scalars
        m0 = state.matter

        srhs0, mrhs0, _ = self._rhs(state)
        gterms0 = adapter.geometry_stage_terms(
            grid, g0, s0, m0, lambda_m=LAMBDA_M
        )

        # TRUE-PIRK predictor: advance every explicit variable first, then
        # evaluate the implicit terms on that predictor together with the
        # explicitly advanced matter/scalars.
        g_explicit1 = g0.copy()
        for name, value in gterms0["explicit"].items():
            setattr(
                g_explicit1,
                name,
                getattr(g0, name) + dt * value,
            )

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
        g1.B = (
            g0.B + 0.75 * (g1.Lambda - g0.Lambda)
            - dt * GAUGE_ETA * g0.B
        )
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
        # second split evaluation above.
        gnew = g0.copy()
        for name in gterms0["explicit"]:
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
        gnew.B = (
            g0.B
            + 0.75 * (gnew.Lambda - g0.Lambda)
            - 0.5 * dt * GAUGE_ETA * (g0.B + g1.B)
        )

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

        surface = int(np.argmin(np.abs(r - 10.0)))
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

        out = {
            "H_eff": effective_hubble(geom, grid.volumes),
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
            "determinant_min": float(
                np.min(geom.a * geom.b**2)
            ),
            "determinant_constraint_max": float(
                np.max(np.abs(geom.a * geom.b**2 - 1.0))
            ),
            "lapse_min": float(np.min(geom.alpha)),
            "lapse_max": float(np.max(geom.alpha)),
            "Rdot_sigma": float(Rdot[surface]),
            "phi_outer": float(np.mean(state.scalars.phi[outer])),
            "Pi_outer": float(np.mean(state.scalars.Pi[outer])),
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
