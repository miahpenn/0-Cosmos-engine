"""Diagnostic-only D-mode source audit used by the preregistered campaign."""
import math
import numpy as np
from .production_kernel import ProductionState, V55ProductionKernel
from .v55_initial import build_initial_data
from .v55_matter import metric_slice_from_q, project_species
from .scalar_system import KAPPA
from . import v55_pirk_adapter as adapter

class SourceAuditKernel(V55ProductionKernel):
    def __init__(self, d_amplitude):
        super().__init__()
        self.d_amplitude = float(d_amplitude)

    def initialize(self, resolution=160, r_max=160.0, **kwargs):
        grid_ops, vacuum, _ = adapter.vendor_modules()
        grid = grid_ops.SphericalCellGrid(resolution, r_max)
        init = build_initial_data(
            grid,
            vacuum.flat_state,
            D_amplitude=self.d_amplitude,
            include_radiation=True,
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

    def diagnostics(self, state, *, profiles=False):
        out = super().diagnostics(state, profiles=profiles)
        g = state.grid
        geom = state.geometry
        fields = state.scalars
        r = np.asarray(g.centers)
        invr = geom.X**2 / geom.a

        def d1(x, parity=1):
            return g.cell_derivative_fourth(x, parity=parity)

        # Existing scalar stress projection, decomposed by sector.
        sectors = {}
        for label, f, p, V, scale in (
            ("S", fields.S, fields.PS, 0.5 * fields.S**2, 1.0),
            ("D", fields.D, fields.PD, -0.5 * fields.D**2, 1.0),
            ("COSMOS", fields.phi, fields.Pi,
             8.242522415500654e-5 * (np.exp(-fields.phi) - 0.10), KAPPA),
        ):
            fp = d1(f)
            gr2 = invr * fp * fp
            rho = (0.5 * p*p + 0.5*gr2 + V) / scale
            pr = (0.5 * p*p + 0.5*gr2 - V) / scale
            pt = (0.5 * p*p - 0.5*gr2 - V) / scale
            sectors[label] = {
                "rho": rho,
                "pr": pr,
                "pt": pt,
                "K_source": 4.0 * math.pi * geom.alpha * (rho + pr + 2.0*pt),
            }

        metric = metric_slice_from_q(g, geom)
        fluids = {}
        for label, species, field in (
            ("dark_matter", 0, state.matter.dark_matter),
            ("baryons", 1, state.matter.baryons),
            ("radiation", 2, state.matter.radiation),
        ):
            from engine.matter_system import Species
            sp = (Species.DARK_MATTER, Species.BARYON, Species.RADIATION)[species]
            piece = project_species(metric, field, sp)
            fluids[label] = {
                "rho": piece["rho"],
                "pr": piece["pr"],
                "pt": piece["pt"],
                "K_source": 4.0 * math.pi * geom.alpha * (
                    piece["rho"] + piece["pr"] + 2.0*piece["pt"]
                ),
            }

        active_D = 2.0*fields.PD**2 + fields.D**2
        i_grad = int(np.argmax(np.abs(d1(geom.alpha))))
        i_center = 0
        i_active = int(np.argmax(active_D))

        # Exact existing geometry K source split: vacuum plus total matter.
        _, vacuum, _ = adapter.vendor_modules()
        vacuum_l3 = vacuum.primary_l3_rhs(g, geom)
        total_l3 = adapter.primary_l3_with_matter(
            g, geom, fields, state.matter
        )
        matter_K = np.asarray(total_l3["K"]) - np.asarray(vacuum_l3["K"])

        def at(arr, idx):
            return float(np.asarray(arr)[idx])

        out.update({
            "D_center": at(fields.D, i_center),
            "PD_center": at(fields.PD, i_center),
            "D_growth_rate_center": (
                at(geom.alpha, i_center) * at(fields.PD, i_center) / at(fields.D, i_center)
                if abs(at(fields.D, i_center)) > 1.0e-300 else 0.0
            ),
            "D_active_center": at(active_D, i_center),
            "D_active_max": at(active_D, i_active),
            "D_active_max_r": at(r, i_active),
            "source_probe_r": at(r, i_grad),
            "D_at_source_probe": at(fields.D, i_grad),
            "D_active_at_source_probe": at(active_D, i_grad),
            "PD_at_source_probe": at(fields.PD, i_grad),
            "matter_K_source_center": at(matter_K, i_center),
            "matter_K_source_probe": at(matter_K, i_grad),
            "vacuum_K_source_center": at(vacuum_l3["K"], i_center),
            "vacuum_K_source_probe": at(vacuum_l3["K"], i_grad),
        })
        for label, data in sectors.items():
            for key, arr in data.items():
                out[f"{label}_{key}_center"] = at(arr, i_center)
                out[f"{label}_{key}_probe"] = at(arr, i_grad)
        for label, data in fluids.items():
            for key, arr in data.items():
                out[f"{label}_{key}_center"] = at(arr, i_center)
                out[f"{label}_{key}_probe"] = at(arr, i_grad)
        return out
