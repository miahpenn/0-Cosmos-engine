"""Production-engine contracts.

The final runner is intentionally blocked from silently falling back to the
old CMC/reference driver. A production kernel must provide reference-metric
BSSN/PIRK state evolution and the frozen V5.5 matter source projections.
"""
from dataclasses import dataclass
from typing import Protocol


@dataclass(frozen=True)
class KernelCapabilities:
    reference_metric_bssn: bool
    pirk2: bool
    moving_gauge: bool
    dm_momentum: bool
    baryon_momentum: bool
    radiation: bool
    invariant_trapping: bool
    misner_sharp_current: bool
    # This is intentionally optional: current V5.5 has no Maxwell state/RHS.
    radial_maxwell: bool = False


class ProductionKernel(Protocol):
    capabilities: KernelCapabilities

    def initialize(self, resolution: int, r_max: float): ...
    def step(self, state, dt): ...
    def diagnostics(self, state): ...


REQUIRED = KernelCapabilities(
    reference_metric_bssn=True,
    pirk2=True,
    moving_gauge=True,
    dm_momentum=True,
    baryon_momentum=True,
    radiation=True,
    invariant_trapping=True,
    misner_sharp_current=True,
    radial_maxwell=False,
)


def require_production_capabilities(
    kernel: ProductionKernel,
    *,
    require_radial_maxwell: bool = False,
) -> None:
    """Check baseline runtime support and any explicitly requested feature.

    The baseline engine does not implement radial Maxwell dynamics. A caller
    that requires that feature must opt in here and fail closed unless the
    kernel advertises it; radiation support is not a proxy for Maxwell support.
    """
    missing = [
        name for name in REQUIRED.__dataclass_fields__
        if getattr(REQUIRED, name) and not getattr(kernel.capabilities, name, False)
    ]
    if require_radial_maxwell and not getattr(
        kernel.capabilities, "radial_maxwell", False
    ):
        missing.append("radial_maxwell")
    if missing:
        raise RuntimeError(
            "Production kernel is incomplete; missing: " + ", ".join(missing)
        )
