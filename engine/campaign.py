"""Final numerical campaign orchestration contract.

This module can host the complete campaign once the production kernel is
promoted. It deliberately refuses to run against the incomplete CMC/reference
driver, preventing an accidental acceptance run during construction.
"""
from dataclasses import dataclass
from pathlib import Path
import json
import math

from .production_contract import require_production_capabilities


@dataclass(frozen=True)
class CampaignConfig:
    # Preserve the original dr values while moving the finite-radius
    # worldtube outward so the strong outgoing D geometry remains interior.
    resolutions: tuple = (80, 120, 160)
    r_max: float = 80.0
    cfl: float = 0.03
    # Coordinate-time safety cap; evolution may stop earlier at final_proper_time.
    final_time: float = 160.0
    final_proper_time: float | None = None
    checkpoint_interval: float = 2.0
    output_dir: str = "runs"


def campaign_endpoint_reached(t: float, tau: float, config: CampaignConfig) -> bool:
    """Whether a campaign endpoint was reached, without treating a cap as physics."""
    return (
        float(t) >= float(config.final_time)
        or (
            config.final_proper_time is not None
            and float(tau) >= float(config.final_proper_time)
        )
    )


def write_checkpoint(path: Path, metadata: dict, trajectory: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"metadata": metadata, "trajectory": trajectory}
    path.write_text(json.dumps(payload, indent=2, default=float))


def prepare_campaign(kernel, config: CampaignConfig = CampaignConfig()):
    """Validate the production kernel and return an immutable campaign plan."""
    require_production_capabilities(kernel)
    if not math.isfinite(float(config.final_time)) or config.final_time <= 0.0:
        raise ValueError("final_time must be finite and positive")
    if config.final_proper_time is not None and (
        not math.isfinite(float(config.final_proper_time))
        or config.final_proper_time <= 0.0
    ):
        raise ValueError("final_proper_time must be finite and positive when set")
    if config.r_max <= 0.0:
        raise ValueError("r_max must be positive")
    if not config.resolutions:
        raise ValueError("at least one resolution is required")
    return {
        "resolutions": tuple(config.resolutions),
        "r_max": config.r_max,
        "cfl": config.cfl,
        "final_time": config.final_time,
        "final_proper_time": config.final_proper_time,
        "checkpoint_interval": config.checkpoint_interval,
        "output_dir": config.output_dir,
        "status": "prepared_not_executed",
    }
