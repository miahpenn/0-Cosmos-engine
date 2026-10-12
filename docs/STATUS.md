# Build Status — 2026-10-11

## Current state

The canonical repository contains the unified production-facing state graph:

1. corrected archive COSMOS background equations;
2. local GEAR-03 S/D sector with GEAR-63/65 barrier;
3. scalar S/D/COSMOS evolution and stress projections;
4. conservative spherical dark-matter, baryon, and radiation matter sectors;
5. locked beta DM/scalar four-force exchange;
6. one total Einstein stress-energy assembly;
7. pinned reference-metric BSSN/PIRK numerical geometry;
8. synchronized geometry/scalar/matter evolution state;
9. dynamically derived proper time, H_eff, e-folds, trapping, and Misner-Sharp observables;
10. cycle/handoff ledgers without a synthetic lambda_cycle;
11. checkpoint and long-horizon multi-resolution campaign infrastructure.

## Explicit exclusions

- Bianchi-I shear remains outside the exact spherical local source.
- No D-to-matter identification law.
- No shell EOS/action.
- No phenomenological interface source.
- No fitted coefficient.
- No imposed clock conversion.
- No bounce/reset/branch flip/lapse floor/ejection threshold/physical stop.

## Validation state

An independent repository review reported the following checks on 2026-10-11:

- Full test suite: 79/79 passed, including tests in `tests/` and the vendor BSSN/PIRK operator suite.
- Production smoke campaign: `resolution=40`, `r_max=40`, `final_time=1.0`; initialization and a kernel step completed without NaNs or assertion failures.
- Default-path kernel initialization and stepping succeeded.

These results establish basic test and smoke-run functionality. They do **not** close the full production campaign gate.

Current classification:

**IMPLEMENTED / TESTED AT UNIT AND SMOKE-RUN LEVEL / NOT YET FULL-CAMPAIGN-VALIDATED**

## Active numerical audits

- Draft PR #2 tracks the N=320 reconstruction audit and host/runtime-dependent outer Hamiltonian L2 discrepancies. Frozen tolerances and the missing preregistration source documents must be resolved transparently; do not relax tolerances solely to make runners pass.
- Draft PR #1 contains experimental strong-field representation and central-clock diagnostics. It is not production validation and must remain separate from `main`.
- The N=640 signed-pole audit through t=12 is running on `diag/early-signed-pole-splice-audit`; its numerical result is pending until its artifact is available and admitted.
- The bond-sweep diagnostic artifact completed successfully. Its coupling-dependent mode differences are diagnostic evidence, not proof that the bond explains the full evolving constraint discrepancy.

## Final engineering gate

The remaining work is the complete regression and long-horizon multi-resolution campaign, followed by analysis of:

- Hamiltonian, momentum, connection, and determinant witnesses;
- invariant trapped-surface roots;
- Misner-Sharp current balance;
- proper-time and H_eff trajectories;
- dark-sector exchange closure;
- resolution convergence and host sensitivity;
- turnaround and re-expansion events;
- genuinely emergent multi-cycle behavior, if the solved equations produce it.

A numerical failure is classified as a numerical/gauge result unless invariant, resolution-independent evidence supports a physical interpretation. No production physics change is authorized by a passing unit suite or short smoke campaign alone.

## Numerical reproducibility note

Host/runtime-dependent reconstruction discrepancies have been reported in the N=320 audit. Preserve the preregistered tolerance and report failing diagnostics, runtime details, and source hashes. Do not silently loosen tolerances or omit failed fields.

Campaign implementation note: the radiation inversion bracket uses the adjacent representable value only when endpoint roundoff reverses the analytically positive residual; the physical EOS and bound remain unchanged.
