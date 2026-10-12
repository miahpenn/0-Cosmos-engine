# High-Resolution Radiation Cone: Actual-Step Check After RHS Contract Repair

**Date:** 2026-10-12  
**Disposition:** bounded one-step actual-stage admissibility passed at N=240 and N=320. This invalidates the *old one-step artifact* as evidence of a genuine high-resolution finite-step radiation failure. It does not establish long-run convergence or full production readiness.

## Provenance and admission

- Candidate base: PR #10 branch head `6f6fd9d81af739cacce16fcd428bdfeac996e343`.
- Diagnostic branch: `diag/high-resolution-radiation-real-step-after-rhs-contract-2026-10-12`.
- Exact tested head: `a5ca89e57599048503cf259d0fb10f6bfd518700`.
- Pinned vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Workflow [38318958640](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318958640): **success**.
- Existing actual-state radiation guard tests: **4 passed**.
- Full repository suite: **26 passed**.
- A/B-sized diagnostic artifact: [11776889362](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318958640/artifacts/11776889362).
- Artifact ZIP SHA-256: `fa4c915b1bb24daca3c7454319278c78fc7b65a9e04608da9a2c04721b23e301`.
- Embedded `report.json` SHA-256: `21742fcefebb3155b61aab4963086a0e31875756fe17c1af03cce09146f07453`.
- Archived comparison artifact used for source provenance: [11772240963](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304517568/artifacts/11772240963), report.json SHA-256 verified as `0a938cfca543e838cf0415ecd37c450375c34f4c103faab9ba28f212b862d90f`.

The experiment performed one normal `V55ProductionKernel.step` for each resolution from the same unadvanced initialization recipe, `r_max=80), `dt=0.03,dr`. It observes—not modifies—the production validator at `initial`, `explicit-predictor`, `primary-predictor`, and `final`. No retry, projection, floor, clamp, damping, source edit, geometry edit, boundary change, or second time step was introduced.

## Result

| Resolution | (dt) | Minimum cone margin, initial | Minimum margin, explicit predictor | Minimum margin, final | Maximum (|S|/E), final | Status |
|---:|---:|---:|---:|---:|---:|---|
| 240 | 0.0100 | (2.46389496\times10^{-10}) | (2.43876101\times10^{-10}) | (2.44025071\times10^{-10}) | 0.00991926171 | All four stage checks pass |
| 320 | 0.0075 | (1.38593035\times10^{-10}) | (1.37181308\times10^{-10}) | (1.37265093\times10^{-10}) | 0.00991880272 | All four stage checks pass |

The smallest margin in each profile is at centre cell (i=0), as expected for this metric/spherical volume weight. At the explicit predictor, the maximum full-domain ratios are 0.0101439889997 (N=240) and 0.0101434458555 (N=320). At the final state, they are 0.00991926171186 and 0.00991880271985. All cells have finite conservative values and valid metrics, and all satisfy (E\ge0) and (|S|\le E) in these four checked states.

At N=240, the final first-cell values are (U_E=2.46469870253229\times10^{-10}), (U_r=-2.4449822703864208\times10^{-12}), (E=8.871229659084301\times10^{-9}), (|S|=8.799604869450436\times10^{-11}), and (|S|/E=0.00991926171186).

At N=320, the final first-cell values are (U_E=1.3864023783397872\times10^{-10}), (U_r=-1.3752179394719555\times10^{-12}), (E=8.871739986627669\times10^{-9}), (|S|=8.799703870916744\times10^{-11}), and (|S|/E=0.00991880271985).

## Why the earlier N=240/320 rejection must not be read as a real step failure

The combined CMC row A/B artifact above was generated from candidate source commit `9f7ebf1aab2076fe8925737591366b416427305e`. The exact source at that SHA has:

1. `engine/matter_rhs.py::species_rhs` call `evolve_species(..., dt=1.0)` to derive a rate via subtraction.
2. An unconditional radiation realizability check inside `evolve_species` on the resulting output state (U+L(U)).

The archived N=240 exception says `radiation RHS created inadmissible state at cell i=0`, with (|S|/E=1.0202778); the N=320 exception reports 1.3602899. Those values were calculated from the surrogate `dt=1` output, not the actual production stage at the requested (dt=0.01) / (0.0075).

On the corrected PR #10 source, `species_rhs` disables only that finite-time admissibility check for the derivative-only operation. A caller explicitly using `evolve_species` for a real finite-time update retains the direct validator by default. The production kernel validates each actual stage. This is why the new experiment checks the kernel's named stages directly rather than estimating a dt-scaled state from the rejected unit-step output.

**Correction to prior interpretation:** the old combined CMC A/B does not establish that the production initial step fails at N=240 or N=320. Under the corrected API contract, both initial-state one-step tests pass. The previous rejection remains useful evidence of the API-contract defect and of the guard's overly early placement, but not as radiation-physics evidence.

## What remains unresolved

This experiment is deliberately one step. It does not prove multi-step high-resolution stability, constraint propagation, Misner–Sharp balance, or a physical radiation-coupling law. The separate saved N=80 long-run trace still records an actual outer-cell radiation-cone crossing around (t=46.98), when the explicit predictor metric drives (|S|/E) just above unity. Its geometry/source budget is documented in [RADIATION_PREDICTOR_FAILURE_GEOMETRY_SOURCE_ATTRIBUTION_2026-10-12.md](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/RADIATION_PREDICTOR_FAILURE_GEOMETRY_SOURCE_ATTRIBUTION_2026-10-12.md); it is distinct from the old `dt=1) centre-cell artifact.

Next work should use the admitted outer-cell failing-step witness and existing boundary/CMC analysis; no repeated high-resolution RHS derivation or new long trajectory is justified solely by the superseded one-step error. Keep source repair, CMC row candidates, radiation rest-current representation, and physical admission separate. `main` and the GEAR archive remain unchanged.
