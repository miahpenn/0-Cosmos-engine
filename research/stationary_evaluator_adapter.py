"""Research-only coordinate-index evaluator for the stationary-surfer hypothesis.

This module is deliberately not imported by the production kernel. It evaluates
an illustrative profile on a caller-supplied spatial grid; it does not integrate
field equations, derive proper time, or establish that the profile is a physical
4D solution. tau_local and t_cosmos are independent query parameters.

A production BSSNMetricSlice-like object may be passed directly: only its
existing r and alpha arrays are read, and neither is modified.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class StationaryEvaluatorAdapter:
    """Evaluate toy core/wake profiles on an existing grid without time stepping."""

    core_scale: float = 5.0
    stream_speed: float = 0.8
    shell_width: float = 3.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.core_scale) or self.core_scale <= 0.0:
            raise ValueError("core_scale must be finite and positive")
        if not np.isfinite(self.shell_width) or self.shell_width <= 0.0:
            raise ValueError("shell_width must be finite and positive")
        if not np.isfinite(self.stream_speed):
            raise ValueError("stream_speed must be finite")

    @staticmethod
    def _finite_scalar(name: str, value: float) -> float:
        result = float(value)
        if not np.isfinite(result):
            raise ValueError(f"{name} must be finite")
        return result

    @staticmethod
    def _validate_grid(r: Any) -> np.ndarray:
        grid = np.asarray(r, dtype=float)
        if grid.ndim != 1 or grid.size < 2:
            raise ValueError("r must be a one-dimensional grid with at least 2 points")
        if not np.all(np.isfinite(grid)):
            raise ValueError("r must contain only finite values")
        if np.any(np.diff(grid) <= 0.0):
            raise ValueError("r must be strictly increasing")
        return grid

    def evaluate_slice(
        self,
        r: Any,
        *,
        tau_local: float,
        t_cosmos: float,
        lapse: Any | None = None,
    ) -> dict[str, Any]:
        """Evaluate one independent query on the supplied spatial grid.

        lapse is optional and, when supplied, is passed through unchanged in
        value from the caller's metric state. No synthetic lapse is substituted.
        The toy density profiles are illustrative only.
        """
        tau = self._finite_scalar("tau_local", tau_local)
        t = self._finite_scalar("t_cosmos", t_cosmos)
        grid = self._validate_grid(r)

        scale = 1.0 + 0.05 * tau
        if scale <= 0.0:
            raise ValueError("tau_local produces a non-positive profile scale")

        width_factor = 1.0 + 0.02 * t
        if width_factor <= 0.0:
            raise ValueError("t_cosmos produces a non-positive shell width")
        shell_amp = 0.5 + 0.1 * t
        if shell_amp < 0.0:
            raise ValueError("t_cosmos produces a negative shell amplitude")

        effective_r = grid / scale
        core = 10.0 * np.exp(-((effective_r / self.core_scale) ** 2))
        shell_center = 15.0 + self.stream_speed * t
        shell_width = self.shell_width / width_factor
        shell = shell_amp * np.exp(
            -0.5 * ((effective_r - shell_center) / shell_width) ** 2
        )

        lapse_values = None
        if lapse is not None:
            lapse_values = np.asarray(lapse, dtype=float)
            if lapse_values.shape != grid.shape:
                raise ValueError("lapse must have the same shape as r")
            if not np.all(np.isfinite(lapse_values)):
                raise ValueError("lapse must contain only finite values")
            lapse_values = lapse_values.copy()

        return {
            "r": grid.copy(),
            "tau_local": tau,
            "t_cosmos": t,
            "effective_r": effective_r,
            "core_remnant": core,
            "d_shell": shell,
            "rho": core + shell,
            "shell_center_effective": shell_center,
            "shell_width_effective": shell_width,
            "lapse": lapse_values,
            "model_status": "illustrative_profile_not_field_equation_solution",
        }

    def evaluate_metric_slice(
        self,
        metric_state: Any,
        *,
        tau_local: float,
        t_cosmos: float,
    ) -> dict[str, Any]:
        """Adapt a metric-state object exposing r and alpha arrays.

        This is a read-only bridge: it does not alter or evolve the supplied
        metric state and does not imply that the toy density is coupled to it.
        """
        if not hasattr(metric_state, "r") or not hasattr(metric_state, "alpha"):
            raise TypeError("metric_state must expose existing r and alpha arrays")
        return self.evaluate_slice(
            metric_state.r,
            tau_local=tau_local,
            t_cosmos=t_cosmos,
            lapse=metric_state.alpha,
        )
