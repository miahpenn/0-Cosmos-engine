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
    # Runtime implementation availability only; this does not mean the model
    # has passed scientific admission. Older/external kernels default to false.
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
    """Require the baseline V5.5 runtime and any explicitly requested option.

    `radial_maxwell` means the kernel can execute the explicit opt-in research
    path; it is not a claim of physical normalization or scientific admission.
    The option remains optional for ordinary runs, while a caller that requires
    it fails closed if the selected kernel cannot provide it.
    """
    missing = [
        name
        for name, required in REQUIRED.__dict__.items()
        if required and not getattr(kernel.capabilities, name, False)
    ]
    if require_radial_maxwell and not getattr(
        kernel.capabilities, "radial_maxwell", False
    ):
        missing.append("radial_maxwell")
    if missing:
        raise RuntimeError(
            "Production kernel is incomplete; missing: " + ", ".join(missing)
        )
