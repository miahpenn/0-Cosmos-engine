"""Research-only closed-loop recirculation probe.

The boundary quantity is a spatial integral of the illustrative shell-density
tail, not a conserved physical flux. The synthetic lapse is explicitly labelled
and must not be substituted for the engine metric lapse. This module is not
imported by the production runtime and performs no time evolution.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class StationaryRecirculationProbe:
    r_min: float = 0.0
    r_max: float = 24.0
    resolution: int = 24
    core_scale: float = 5.0
    stream_speed: float = 0.8
    feedback_gain: float = 0.15
    shell_width: float = 3.0

    def __post_init__(self) -> None:
        for name in ("r_min", "r_max", "core_scale", "stream_speed",
                     "feedback_gain", "shell_width"):
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if self.r_max <= self.r_min:
            raise ValueError("r_max must be greater than r_min")
        if isinstance(self.resolution, (bool, np.bool_)) or not isinstance(
            self.resolution, (int, np.integer)
        ):
            raise ValueError("resolution must be an integer")
        if self.resolution < 2:
            raise ValueError("resolution must be at least 2")
        if self.core_scale <= 0.0 or self.shell_width <= 0.0:
            raise ValueError("core_scale and shell_width must be positive")

    @staticmethod
    def _finite_query(name: str, value: float) -> float:
        result = float(value)
        if not np.isfinite(result):
            raise ValueError(f"{name} must be finite")
        return result

    def evaluate_slice(self, *, tau_local: float, t_cosmos: float) -> dict[str, Any]:
        tau = self._finite_query("tau_local", tau_local)
        t = self._finite_query("t_cosmos", t_cosmos)

        scale = 1.0 + 0.05 * tau
        width_factor = 1.0 + 0.02 * t
        shell_amp = 0.5 + 0.1 * t
        if scale <= 0.0:
            raise ValueError("tau_local produces a non-positive profile scale")
        if width_factor <= 0.0:
            raise ValueError("t_cosmos produces a non-positive shell width")
        if shell_amp < 0.0:
            raise ValueError("t_cosmos produces a negative shell amplitude")

        r = np.linspace(self.r_min, self.r_max, self.resolution)
        effective_r = r / scale
        core = 10.0 * np.exp(-((effective_r / self.core_scale) ** 2))
        shell_center = 15.0 + self.stream_speed * t
        effective_width = self.shell_width / width_factor
        shell = shell_amp * np.exp(
            -0.5 * ((effective_r - shell_center) / effective_width) ** 2
        )

        # Include the exact physical cutoff by interpolation so the integral
        # does not jump merely because a grid node crosses the boundary.
        cutoff = self.r_min + 0.75 * (self.r_max - self.r_min)
        if cutoff < self.r_min or cutoff > self.r_max:
            raise ValueError("boundary cutoff lies outside the spatial grid")
        tail_mask = r > cutoff
        tail_r = np.concatenate(([cutoff], r[tail_mask]))
        tail_shell = np.concatenate(([np.interp(cutoff, r, shell)], shell[tail_mask]))
        boundary_tail_integral = float(
            np.sum(0.5 * (tail_shell[1:] + tail_shell[:-1]) * np.diff(tail_r))
        )
        recirculation = (
            self.feedback_gain
            * boundary_tail_integral
            * np.exp(-((effective_r / (self.core_scale * 1.5)) ** 2))
        )
        lapse_synthetic = (
            1.0
            - 0.1 * np.exp(-((r / self.core_scale) ** 2)) * np.cos(0.1 * t)
        )

        return {
            "r": r.copy(),
            "tau_local": tau,
            "t_cosmos": t,
            "effective_r": effective_r,
            "rho": core + shell + recirculation,
            "core_remnant": core,
            "d_shell": shell,
            "recirculation": recirculation,
            "boundary_tail_integral": boundary_tail_integral,
            "lapse_synthetic": lapse_synthetic,
            "boundary_measure": "spatial_integral_of_shell_density_tail_not_conserved_flux",
            "model_status": "illustrative_recirculation_probe_not_field_equation_solution",
        }
