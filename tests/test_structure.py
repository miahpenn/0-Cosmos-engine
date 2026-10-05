from engine.production_contract import (
    KernelCapabilities,
    REQUIRED,
    require_production_capabilities,
)
from engine.spherical_scope import SphericalMatterScope
from engine.stress_energy import TotalStressEnergy


def test_required_capability_contract_is_complete():
    assert all(REQUIRED.__dict__.values())


def test_spherical_scope_rejects_shear():
    try:
        SphericalMatterScope(include_shear=True).validate()
    except ValueError:
        return
    raise AssertionError("shear must be rejected in exact spherical scope")


def test_total_stress_energy_shape_contract():
    assert TotalStressEnergy.__dataclass_fields__.keys()
