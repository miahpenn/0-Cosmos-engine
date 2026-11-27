"""Engine-owned V1 homogeneous research state contract; no evolution entry point.

This is a formulation-specific contract, not a ProductionKernel implementation.
The frozen native_v1 package supplies spatial mathematics. State construction
checks representation only; it does not establish physical admissibility.
"""
from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

SCHEMA = "native-v1-homogeneous-state/v1"
GEOMETRY_ORDER = ("p", "k", "ell", "J", "RV", "L0", "K", "kh", "eh", "incoming", "w")
SPECIES = ("dark_matter", "baryons", "radiation")
SCALAR_ORDER = ("S", "P_S", "D", "P_D", "phi", "P_phi")
TICK_DT = 0.000125
FINAL_TICK = 8


def _integer(value, label):
    if type(value) is not int:
        raise ValueError(f"{label}: Python integer required")


def _array(value, shape, label):
    a = np.asarray(value)
    if a.shape != shape or a.dtype != np.dtype("float64") or not np.all(np.isfinite(a)):
        raise ValueError(f"{label}: finite float64 with shape {shape} required")
    # Backing bytes are immutable, including against setflags(write=True).
    return np.frombuffer(a.tobytes(order="C"), dtype=np.float64).reshape(shape)


@dataclass(frozen=True, slots=True)
class V1RunSpec:
    resolution: int
    stride: int

    def __post_init__(self):
        _integer(self.resolution, "resolution")
        _integer(self.stride, "stride")
        if self.resolution not in (8, 16) or self.stride not in (1, 2):
            raise ValueError("outside frozen homogeneous research scope")

    @property
    def radius(self):
        return 1.0

    @property
    def dt(self):
        return self.stride * TICK_DT

    @property
    def geometry_shapes(self):
        n = self.resolution
        return {k: ((3*n+4,) if k in ("p", "k", "ell") else
                    (3*n+5,) if k in ("J", "RV") else
                    (n,) if k == "incoming" else ()) for k in GEOMETRY_ORDER}

    @property
    def gauge_shape(self):
        return (2*self.resolution+4,)


@dataclass(frozen=True, slots=True, eq=False)
class V1AcceptedState:
    spec: V1RunSpec
    tick: int
    geometry: Mapping
    matter: Mapping
    scalars: tuple
    packed_gauge: np.ndarray

    def __post_init__(self):
        if type(self.spec) is not V1RunSpec:
            raise ValueError("explicit V1RunSpec required")
        _integer(self.tick, "accepted tick")
        if not 0 <= self.tick <= FINAL_TICK or self.tick % self.spec.stride:
            raise ValueError("tick is not an accepted endpoint of the frozen schedule")
        if tuple(self.geometry) != GEOMETRY_ORDER or tuple(self.matter) != SPECIES:
            raise ValueError("frozen geometry/species layout or ordering mismatch")
        n = self.spec.resolution
        geometry = {k: _array(self.geometry[k], shape, k) for k, shape in self.spec.geometry_shapes.items()}
        matter = {}
        for k in SPECIES:
            if len(self.matter[k]) != 3:
                raise ValueError(f"{k}: three raw moments required")
            matter[k] = tuple(_array(v, (n,), f"{k}[{j}]") for j, v in enumerate(self.matter[k]))
        if len(self.scalars) != 6:
            raise ValueError("six scalar coefficient arrays required")
        scalars = tuple(_array(v, (n,), k) for k, v in zip(SCALAR_ORDER, self.scalars))
        gauge = _array(self.packed_gauge, self.spec.gauge_shape, "packed gauge")
        if gauge[-1] != 1.0:
            raise ValueError("augmented gauge must end in one")
        object.__setattr__(self, "geometry", MappingProxyType(geometry))
        object.__setattr__(self, "matter", MappingProxyType(matter))
        object.__setattr__(self, "scalars", scalars)
        object.__setattr__(self, "packed_gauge", gauge)

    @property
    def t(self):
        return self.tick * TICK_DT

    def kernel_input(self):
        """Detached, ordered physical input and incoming gauge; never alias accepted state."""
        geometry = {k: (float(v) if v.ndim == 0 else v.copy()) for k, v in self.geometry.items()}
        matter = {k: tuple(v.copy() for v in vs) for k, vs in self.matter.items()}
        scalars = [v.copy() for v in self.scalars]
        return (geometry, matter, scalars), self.packed_gauge.copy()
