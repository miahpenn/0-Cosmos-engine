"""Diagnostic-only helper for the R=160 alpha history at r=79.75.

The 49 source slices are spaced about 0.25 in time, not at the PIRK stage
interval dt=0.00375. Linear interpolation is not a stage-level error bound.
Held-out interpolation max absolute error on 0.5-spacing was about 1.5e-4
(RMS about 5.6e-5). This is index 159 (r=79.75), not the R=160 outermost cell.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

_DEFAULT_SERIES = Path(__file__).resolve().parent / "r160_alpha_boundary_r79_75.json"


class BoundaryAlphaSeries:
    """Linear interpolator for the R=160 alpha(r=79.75) reference history."""

    def __init__(self, path: str | Path | None = None):
        path = Path(path) if path is not None else _DEFAULT_SERIES
        with open(path) as f:
            data = json.load(f)
        if data.get("schema") != "r160_boundary_alpha_r79_75_v1":
            raise ValueError(f"unexpected schema in {path}")
        samples = data["samples"]
        self.t = np.asarray([s["t"] for s in samples], dtype=float)
        self.alpha = np.asarray([s["alpha"] for s in samples], dtype=float)
        self.source_run = data.get("source_run")
        self.source_index = data.get("source_index")
        self.source_radius = data.get("source_radius")
        self.path = str(path)
        if self.t.size < 2:
            raise ValueError("alpha series must contain at least two samples")
        if not np.all(np.isfinite(self.t)) or not np.all(np.isfinite(self.alpha)):
            raise ValueError("alpha series contains non-finite values")
        if not np.all(np.diff(self.t) > 0):
            raise ValueError("alpha sample times must be strictly increasing")

    def __call__(self, t: float) -> float:
        """Linear interpolation; clamps to endpoints outside [t0, tN]."""
        t = float(t)
        if t <= self.t[0]:
            return float(self.alpha[0])
        if t >= self.t[-1]:
            return float(self.alpha[-1])
        return float(np.interp(t, self.t, self.alpha))

    def describe(self) -> dict:
        return {
            "path": self.path,
            "n_samples": int(self.t.size),
            "t_range": [float(self.t[0]), float(self.t[-1])],
            "alpha_range": [float(self.alpha.min()), float(self.alpha.max())],
            "source_run": self.source_run,
            "source_index": self.source_index,
            "source_radius": self.source_radius,
            "interpolation": "linear_in_t_clamp_endpoints",
            "limitation": "49-slice history (~0.25 spacing) does not bound stage-level interpolation error at dt=0.00375; held-out max absolute error ~1.5e-4 on 0.5 spacing",
        }
