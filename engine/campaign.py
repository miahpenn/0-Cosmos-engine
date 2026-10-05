"""Final numerical campaign orchestration contract.

This module can host the complete campaign once the production kernel is
promoted. It deliberately refuses to run against the incomplete CMC/reference
driver, preventing an accidental acceptance run during construction.
"""
from dataclasses import dataclass
from pathlib import Path
import json

from .production_contract import require_production_capabilities


@dataclass(frozen=True)
class CampaignConfig:
    resolutions: tuple = (40, 60, 80)
    r_max: float = 40.0
    cfl: float = 0.03
    final_time: float = 160.0
    checkpoint_interval: float = 2.0
    output_dir: str = "runs"


def write_checkpoint(path: Path, metadata: dict, trajectory: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"metadata": metadata, "trajectory": trajectory}
    path.write_text(json.dumps(payload, indent=2, default=float))


def prepare_campaign(kernel, config: CampaignConfig = CampaignConfig()):
    """Validate the production kernel and return an immutable campaign plan."""
    require_production_capabilities(kernel)
    if config.final_time <= 0.0:
        raise ValueError("final_time must be positive")
    if config.r_max <= 0.0:
        raise ValueError("r_max must be positive")
    if not config.resolutions:
        raise ValueError("at least one resolution is required")
    return {
        "resolutions": tuple(config.resolutions),
        "r_max": config.r_max,
        "cfl": config.cfl,
        "final_time": config.final_time,
        "checkpoint_interval": config.checkpoint_interval,
        "output_dir": config.output_dir,
        "status": "prepared_not_executed",
    }
