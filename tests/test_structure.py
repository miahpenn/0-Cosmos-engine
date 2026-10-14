from engine.production_contract import (
    KernelCapabilities,
    REQUIRED,
    require_production_capabilities,
)
from engine.spherical_scope import SphericalMatterScope
from engine.stress_energy import TotalStressEnergy


def test_required_capability_contract_is_complete():
    assert all(
        getattr(REQUIRED, name)
        for name in REQUIRED.__dataclass_fields__
        if name != "radial_maxwell"
    )
    # This field is optional for the baseline contract, not a required feature.
    assert REQUIRED.radial_maxwell is False


def test_radial_maxwell_capability_is_opt_in_and_fails_closed():
    class KernelStub:
        def __init__(self, capabilities):
            self.capabilities = capabilities

    unsupported = KernelStub(REQUIRED)
    # Historical/default production runs do not require the research option.
    require_production_capabilities(unsupported)

    try:
        require_production_capabilities(
            unsupported, require_radial_maxwell=True
        )
    except RuntimeError as exc:
        assert "radial_maxwell" in str(exc)
    else:
        raise AssertionError("an EM-required call must fail closed on an unsupported kernel")

    from engine.production_kernel import V55ProductionKernel

    assert V55ProductionKernel.capabilities.radial_maxwell is True
    require_production_capabilities(
        V55ProductionKernel(), require_radial_maxwell=True
    )


def test_spherical_scope_rejects_shear():
    try:
        SphericalMatterScope(include_shear=True).validate()
    except ValueError:
        return
    raise AssertionError("shear must be rejected in exact spherical scope")


def test_total_stress_energy_shape_contract():
    assert TotalStressEnergy.__dataclass_fields__.keys()
