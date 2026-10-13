import numpy as np

from engine.production_kernel import V55ProductionKernel


def test_diagnostics_do_not_fail_when_domain_does_not_reach_fixed_cutoff():
    kernel = V55ProductionKernel()
    # r_max=4 is accepted by CampaignConfig and predates the fixed-r=5
    # comparison mask. The diagnostic should disclose missing exterior
    # coverage rather than turning this report-only limitation into a failure.
    state = kernel.initialize(
        resolution=40,
        r_max=4.0,
        include_radiation=True,
    )

    diagnostics = kernel.diagnostics(state)

    assert diagnostics["constraint_exclusion_radius"] == 5.0
    assert diagnostics["constraint_exterior_coverage_available"] is False
    assert diagnostics["constraint_inner_coverage_available"] is True
    assert (
        diagnostics["constraint_exterior_norm_region"]
        == "full_domain_fallback_no_r_ge_cutoff"
    )
    # The scalar aliases retain a defined value through the explicit full-domain
    # fallback; the coverage marker prevents misreading it as a fixed-r norm.
    assert diagnostics["hamiltonian_max"] == diagnostics["hamiltonian_full_absmax"]
    assert diagnostics["momentum_max"] == diagnostics["momentum_full_absmax"]
    assert diagnostics["connection_max"] == diagnostics["connection_full_absmax"]
    for key in (
        "hamiltonian_max",
        "momentum_max",
        "connection_max",
        "hamiltonian_full_absmax",
        "momentum_full_absmax",
        "connection_full_absmax",
    ):
        assert np.isfinite(diagnostics[key])
