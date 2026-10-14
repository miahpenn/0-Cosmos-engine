# Fixed-R residual discriminators (research-only)

This module adds named JSON witnesses for the charged-minus-neutral fixed-areal-radius Misner–Sharp residual on PR #184.

- \`frozen_geometry_stress_balance(...)\` evaluates the factor-interpolated stress-flux plus radial-pressure-work expression on the supplied snapshot profiles. It does not re-evolve geometry. The interpolated lapse must be finite and strictly positive; invalid values fail closed.
- \`pure_sector_rates(...)\` reports signed stress-flux and pressure-work for scalar, fluid, and electromagnetic sectors, with a 64-epsilon additive consistency check.

## Integration and interpretation

Both results are recorded by \`_worldtube_sample\` and copied into actual start/predictor-stage snapshots. The pure-sector helper is cross-checked against the pre-existing \`tube_sector_rates\` implementation at roundoff scale; the existing field is retained for continuity.

The audit separately records factor-interpolated and composite-interpolated flux/work. These are different discrete operators and must not be silently substituted for one another. The frozen-geometry result is a continuum-expression evaluation on supplied snapshot fields, not an independent geometry evolution and not proof that geometry terms have been eliminated.

No production equations, defaults, or source pins are changed. Disposition remains **draft / research-only / NOT ADMITTED** until the discrete identity, interpolation conventions, and actual evolution-stage pairing are independently reconciled.


## Bracket-secant moving-interpolant rate witness

Each captured predictor stage now records \`tube_piecewise_linear_fixed_R_rate_witness\`. It uses the exact active \`np.interp\` bracket and theta, with the piecewise-linear secant \`(M[right]-M[left])/(R[right]-R[left])\` and sweep \`-interp(Rdot)*secant\`. It reconstructs the directly differentiated moving-interpolant tangent without substituting the fourth-order \`D_r M\) operator. The factor/composite stress-flux and pressure-work convention is explicitly supplied and recorded. Manufactured-profile tests cover exact bracket pairing, the secant operator, and fail-closed theta/nonfinite/operator mismatch. This remains a research diagnostic; it does not assert the continuum balance closes or alter production physics.

