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
    # The active kernel must not imply that radial Maxwell dynamics exist.
    assert REQUIRED.radial_maxwell is False


def test_radial_maxwell_capability_is_opt_in_and_fails_closed():
    class KernelStub:
        def __init__(self, capabilities):
            self.capabilities = capabilities

    baseline = KernelStub(REQUIRED)
    # The existing non-EM production path remains valid.
    require_production_capabilities(baseline)

    try:
        require_production_capabilities(
            baseline, require_radial_maxwell=True
        )
    except RuntimeError as exc:
        assert "radial_maxwell" in str(exc)
    else:
        raise AssertionError(
            "an EM-required run must fail closed on the non-EM kernel"
        )

    advertised = KernelStub(
        KernelCapabilities(**{**REQUIRED.__dict__, "radial_maxwell": True})
    )
    require_production_capabilities(
        advertised, require_radial_maxwell=True
    )


def test_spherical_scope_rejects_shear():
    try:
        SphericalMatterScope(include_shear=True).validate()
    except ValueError:
        return
    raise AssertionError("shear must be rejected in exact spherical scope")


def test_total_stress_energy_shape_contract():
    assert TotalStressEnergy.__dataclass_fields__.keys()
