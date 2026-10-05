"""Integrated V5.5 production state and synchronized PIRK2 step.

This is the first production-facing machine layer that evolves:
- reference-metric BSSN geometry from the pinned numerical kernel;
- V5.5 local S/D and COSMOS scalars;
- conservative DM/baryon/radiation states;
- the derived proper clock and cycle ledger.

No branch reset, bounce, ejection threshold, interface source, or fitted
feedback coefficient is present. The only hard termination in a host runner
should be non-finite state or an explicitly stated computational resource
limit.
"""
from dataclasses import dataclass, field
import math
from pathlib import Path
import numpy as np

from .handoff import CycleLedger, handoff_from_ledger
from .matter_rhs import species_rhs
from .matter_system import geometry_metric_derivatives
from .production_contract import KernelCapabilities
from .scalar_system import ScalarFields, scalar_rhs_arrays
from .v55_initial import build_initial_data
from .v55_matter import V55MatterState, dm_density, total_matter_projection
from . import v55_pirk_adapter as adapter
from .invariant_diagnostics import (
    misner_sharp_mass_from_chi,
    trapping_indicator_from_areal_radius,
)

BETA_DM = -0.04


@dataclass
class ProductionState:
    grid: object
    geometry: object
    scalars: ScalarFields
    matter: V55MatterState
    t: float = 0.0
    tau: float = 0.0
    cycle: CycleLedger = field(default_factory=CycleLedger)
    history: list[dict] = field(default_factory=list)
    handoffs: list = field(default_factory=list)


@dataclass(frozen=True)
class ProductionKernel:
    """Capability-bearing kernel object used by the final campaign contract."""
    capabilities: KernelCapabilities = KernelCapabilities(
        reference_metric_bssn=True,
        pirk2=True,
        moving_gauge=True,
        dm_momentum=True,
        baryon_momentum=True,
        radiation=True,
        invariant_trapping=True,
        misner_sharp_current=True,
    )

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
        ops, vacuum, _ = adapter.vendor_modules()
        grid = ops.SphericalCellGrid(resolution, r_max)
        init = build_initial_data(
            grid,
            vacuum.VacuumState,
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
            t=0.0,
            tau=0.0,
        )

    @staticmethod
    def _grid_d1(grid, values, parity):
        return grid.cell_derivative_fourth(values, parity=parity)

    def _scalar_and_matter_rhs(self, state: ProductionState):
        rho_dm = dm_density(
            adapter_vm_slice(state),
            state.matter,
        )
        scalar_rhs = scalar_rhs_arrays(
            state.grid,
            state.geometry,
            state.scalars,
            beta_dm=BETA_DM,
            rho_dm=rho_dm,
        )

        md = geometry_metric_derivatives(
            adapter_vm_slice(state),
            self._grid_d1,
        )
        dphi_t = scalar_rhs.phi
        dphi_r = self._grid_d1(
            state.grid, state.scalars.phi, 1
        )

        matter_rhs = {
            name: species_rhs(
                adapter_vm_slice(state),
                md,
                getattr(state.matter, name),
                species,
                dphi_t=dphi_t if species.name == "DARK_MATTER" else None,
                dphi_r=dphi_r if species.name == "DARK_MATTER" else None,
                beta_dm=BETA_DM,
            )
            for name, species in (
                ("dark_matter", __import__("engine.matter_system", fromlist=["Species"]).Species.DARK_MATTER),
                ("baryons", __import__("engine.matter_system", fromlist=["Species"]).Species.BARYON),
                ("radiation", __import__("engine.matter_system", fromlist=["Species"]).Species.RADIATION),
            )
        }
        return scalar_rhs, matter_rhs, md

    def step(self, state: ProductionState, dt: float) -> ProductionState:
        if not dt > 0.0 or not np.isfinite(dt):
            raise ValueError("dt must be finite and positive")

        grid = state.grid
        g0 = state.geometry
        f0 = state.scalars
        m0 = state.matter

        # Stage 0: explicit geometry/scalars/matter plus PIRK primary terms.
        gterms0 = adapter.geometry_stage_terms(
            grid, g0, f0, m0
        )
        srhs0, mrhs0, _ = self._scalar_and_matter_rhs(state)

        g1 = g0.copy()
        for name, rhs in gterms0["explicit"].items():
            setattr(g1, name, getattr(g0, name) + dt * rhs)

        # PIRK primary predictor uses explicit-stage geometry and stage-0
        # primary source.  This is the same variable-role ordering as the
        # pinned kernel, but with V5.5 sources.
        l2_0 = gterms0["primary_l2"]
        l3_0 = gterms0["primary_l3"]
        ll2_0 = gterms0["lambda_l2"]
        ll3_0 = gterms0["lambda_l3"]

        old_primary = g0.copy()
        old_primary.a = g1.a
        old_primary.b = g1.b
        old_primary.X = g1.X
        old_primary.alpha = g1.alpha
        old_primary.beta = g1.beta
        l2_u1 = adapter.vendor_modules()[1].primary_l2_rhs(grid, old_primary)

        g1.Aa = g0.Aa + dt * (
            0.5 * l2_0["Aa"] + 0.5 * l2_u1["Aa"] + l3_0["Aa"]
        )
        g1.K = g0.K + dt * (
            0.5 * l2_0["K"] + 0.5 * l2_u1["K"] + l3_0["K"]
        )
        g1.Lambda = g0.Lambda + dt * (
            0.5 * ll2_0 + 0.5 * adapter.vendor_modules()[1].lambda_l2_rhs(
                grid, g1, lambda_m=2.0
            ) + ll3_0
        )
        g1.B = g0.B + 0.75 * (g1.Lambda - g0.Lambda)

        s1 = ScalarFields(
            *(getattr(f0, name) + dt * getattr(srhs0, name)
              for name in ("S", "PS", "D", "PD", "phi", "Pi"))
        )
        m1 = V55MatterState(
            *(ConservedSpecies(
                getattr(m0, name).rest + dt * getattr(mrhs0, name).rest,
                getattr(m0, name).energy_t + dt * getattr(mrhs0, name).energy_t,
                getattr(m0, name).momentum_r + dt * getattr(mrhs0, name).momentum_r,
            ) for name in ("dark_matter", "baryons", "radiation"))
        )

        stage1 = ProductionState(
            grid=grid,
            geometry=g1,
            scalars=s1,
            matter=m1,
            t=state.t,
            tau=state.tau,
        )

        # Stage 1 explicit RHS uses the complete predicted U state.
        gterms1 = adapter.geometry_stage_terms(
            grid, g1, s1, m1
        )
        srhs1, mrhs1, _ = self._scalar_and_matter_rhs(stage1)

        # Final explicit block.
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
                getattr(f0, name) + 0.5 * dt * (
                    getattr(srhs0, name) + getattr(srhs1, name)
                )
                for name in ("S", "PS", "D", "PD", "phi", "Pi")
            )
        )
        mnew = V55MatterState(
            *(
                ConservedSpecies(
                    getattr(m0, name).rest + 0.5 * dt * (
                        getattr(mrhs0, name).rest + getattr(mrhs1, name).rest
                    ),
                    getattr(m0, name).energy_t + 0.5 * dt * (
                        getattr(mrhs0, name).energy_t + getattr(mrhs1, name).energy_t
                    ),
                    getattr(m0, name).momentum_r + 0.5 * dt * (
                        getattr(mrhs0, name).momentum_r + getattr(mrhs1, name).momentum_r
                    ),
                )
                for name in ("dark_matter", "baryons", "radiation")
            )
        )

        # Final PIRK primary correction using the final explicit block.
        final_u_old_v = g1.copy()
        for name in gterms0["explicit"]:
            setattr(final_u_old_v, name, getattr(gnew, name))
        _, vacuum, _ = adapter.vendor_modules()
        l2_final = vacuum.primary_l2_rhs(grid, final_u_old_v)
        ll2_final = vacuum.lambda_l2_rhs(
            grid, final_u_old_v, lambda_m=2.0
        )

        gnew.Aa = g0.Aa + 0.5 * dt * (
            l2_0["Aa"] + l2_final["Aa"]
            + l3_0["Aa"] + gterms1["primary_l3"]["Aa"]
        )
        gnew.K = g0.K + 0.5 * dt * (
            l2_0["K"] + l2_final["K"]
            + l3_0["K"] + gterms1["primary_l3"]["K"]
        )
        gnew.Lambda = g0.Lambda + 0.5 * dt * (
            ll2_0 + ll2_final
            + ll3_0 + gterms1["lambda_l3"]
        )
        gnew.B = g0.B + 0.75 * (gnew.Lambda - g0.Lambda)

        candidate = ProductionState(
            grid=grid,
            geometry=gnew,
            scalars=snew,
            matter=mnew,
            t=state.t + dt,
            tau=state.tau,
            cycle=state.cycle,
            history=list(state.history),
            handoffs=list(state.handoffs),
        )

        obs = self.diagnostics(candidate)
        previous_H = (
            state.history[-1]["H_eff"] if state.history else obs["H_eff"]
        )
        event = candidate.cycle.observe(
            len(candidate.history),
            candidate.t,
            candidate.tau + dt * obs["tau_rate"],
            previous_H,
            obs["H_eff"],
        )
        candidate.tau += dt * obs["tau_rate"]

        hp = handoff_from_ledger(
            candidate.t,
            candidate.tau,
            obs,
            float(candidate.scalars.S[0]),
            float(srhs1.S[0]),
            float(candidate.scalars.D[0]),
            float(srhs1.D[0]),
        )
        candidate.handoffs.append(hp)
        obs["cycle_event"] = event.kind if event else None
        obs["t"] = candidate.t
        obs["tau"] = candidate.tau
        obs["R_sigma"] = hp.radius
        obs["M_MS_sigma"] = hp.misner_sharp_mass
        candidate.history.append(obs)

        return candidate

    def diagnostics(self, state: ProductionState) -> dict:
        grid = state.grid
        geom = state.geometry
        scalars = state.scalars
        matter = state.matter

        _, vacuum, moving = adapter.vendor_modules()
        metric = adapter_vm_slice(state)
        total = total_matter_projection(
            grid, geom, scalars, matter
        )
        raw_constraints = vacuum.constraints(grid, geom)
        H = raw_constraints["hamiltonian"] - 16.0 * math.pi * total["rho"]
        M = raw_constraints["momentum"] - 8.0 * math.pi * total["j"]

        r = np.asarray(grid.centers)
        R = r * np.sqrt(geom.b) / geom.X
        Rr = grid.cell_derivative_fourth(R, parity=1)
        Ktheta = geom.K / 3.0 - geom.Aa / 2.0
        normal_Rt = -R * Ktheta
        chi = trapping_indicator_from_areal_radius(
            Rr, normal_Rt, geom.X**2 / geom.a
        )
        mass = misner_sharp_mass_from_chi(R, chi)

        eu = moving.moving_puncture_explicit_rhs(grid, geom)
        Rdot = R * (0.5 * eu["b"] / geom.b - eu["X"] / geom.X)

        surface = int(np.argmin(np.abs(r - 10.0)))
        outer = r >= 0.8 * grid.r_max

        E = total["rho"]
        j = total["j"]
        coordinate_energy_flux = 4.0 * math.pi * R[surface]**2 * (
            geom.alpha[surface] * (geom.X[surface]**2 / geom.a[surface])
            * j[surface] - geom.beta[surface] * E[surface]
        )
        work_pR = -4.0 * math.pi * R[surface]**2 * (
            total["pr"][surface] * Rdot[surface]
        )

        H_eff = -float(
            np.sum(grid.volumes * geom.K) / np.sum(grid.volumes)
        )

        roots = []
        signs = chi[:-1] * chi[1:]
        for i in np.where(signs <= 0.0)[0]:
            if chi[i] != chi[i + 1]:
                roots.append(
                    float(
                        r[i] - chi[i] * (r[i + 1] - r[i]) /
                        (chi[i + 1] - chi[i])
                    )
                )

        return {
            "H_eff": H_eff,
            "tau_rate": float(geom.alpha[0]),
            "R_sigma": float(R[surface]),
            "M_MS": float(mass[surface]),
            "chi_sigma": float(chi[surface]),
            "trapping_min": float(np.min(chi)),
            "trapped_roots": roots,
            "coordinate_energy_flux": float(coordinate_energy_flux),
            "flux_T": float(
                4.0 * math.pi * R[surface]**2
                * (-geom.alpha[surface] * (geom.X[surface]**2 / geom.a[surface])
                   * j[surface] * Rr[surface])
            ),
            "work_pR": float(work_pR),
            "rho_outer": float(np.mean(total["rho"][outer])),
            "p_outer": float(np.mean(total["pr"][outer])),
            "j_outer": float(np.mean(total["j"][outer])),
            "rho_total_max": float(np.max(total["rho"])),
            "hamiltonian_max": float(np.max(np.abs(H[2:]))),
            "momentum_max": float(np.max(np.abs(M[2:]))),
            "connection_max": float(np.max(np.abs(raw_constraints["connection"][2:]))),
            "determinant_min": float(np.min(geom.a * geom.b**2)),
            "lapse_min": float(np.min(geom.alpha)),
            "lapse_max": float(np.max(geom.alpha)),
            "Rdot_sigma": float(Rdot[surface]),
            "phi_outer": float(np.mean(scalars.phi[outer])),
            "Pi_outer": float(np.mean(scalars.Pi[outer])),
            "total_rho": total["rho"],
            "total_pr": total["pr"],
            "total_pt": total["pt"],
            "total_j": total["j"],
            "areal_radius": R,
            "chi": chi,
            "misner_sharp": mass,
        }


def adapter_vm_slice(state: ProductionState):
    from .v55_matter import metric_slice_from_q
    return metric_slice_from_q(state.grid, state.geometry)


def write_checkpoint_npz(
    state: ProductionState,
    path: str | Path,
) -> None:
    """Write the solved state without imposing any restart/reset operation."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    z = {
        "t": np.asarray(state.t),
        "tau": np.asarray(state.tau),
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
        "radiation_energy_t": state.matter.radiation.energy_t,
        "radiation_momentum_r": state.matter.radiation.momentum_r,
    }
    np.savez_compressed(out, **z)
