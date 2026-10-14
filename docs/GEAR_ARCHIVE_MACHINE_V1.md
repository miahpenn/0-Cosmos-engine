# GEAR Archive Machine v1 — independent-lane build

**Status:** research-only scaffold; not a production or scientific-admission claim  
**Source branch:** \`research/gear-archive-machine-v1-2026-10-13\`  
**Preserved references:** \`CYCLIC-COSMOS ENGINE-.docx\`, \`GEAR_FULL_CYCLIC_COSMOS_ENGINE_ARCHIVE_SNAPSHOT_2026-10-02.md\`, and \`GEAR_MASTER_ARCHIVE_HOUSEKEEPING_INDEX_2026-10-03_v2.md\`.

## Purpose and scope

This starts the older, archive-derived two-lane machine as a versioned runner over the existing source-mapped sector RHSs. It fills the missing orchestration/test seam only. The canonical archived file and the current production kernel are not rewritten, and the local and homogeneous COSMOS sectors remain independent.

Run the small deterministic harness with:

\`\`\`bash
pytest -q tests/test_gear_archive_machine_v1.py
\`\`\`

The Python entry point is \`research/gear_archive_machine_v1.py:run_machine\`. Its local lane uses \`tau_local\`; its cosmology lane uses \`t_cosmos\`. The time axes may be advanced with different final times and step sizes. The runner does not align them, compare equal numerical times as if they were physically simultaneous, or use either clock to rescale the other.

## Source-to-kernel contract

| Sector / diagnostic | Governing/source map | Current v1 behavior | Interpretation limit |
|---|---|---|---|
| Local stable mode S | GEAR-03: \(\omega_S^2=\kappa+g\), \(\ddot S=-\omega_S^2 S\) on this independent local-time lane | Uses \`engine.local.local_rhs\` with \`H=0\`; verifies positive \(\kappa+g\) | This is the reduced mode sector, not the full radial Einstein–scalar PDE |
| Local D mode and geometric barrier | GEAR-03 and GEAR-63/65: \(\alpha_D^2=g-\kappa\), cubic barrier in \`engine/local.py\` | Uses the existing local RHS; checks \(\alpha_D^2=g-\kappa\) instead of accepting inconsistent inputs | Does not imply D energy has become ordinary matter |
| Local energy/activity | Archive ledgers distinguish raw D energy from positive kinetic/excursion activity; \(F_{\rm proxy}\sim |D\dot D|\) | Emits \(E_S\), raw \(E_D\), kinetic/excursion proxies, a generated-activity proxy, and a transport-activity proxy | Proxies are not particle abundance, ejection proof, a physical flux, or an exchange source |
| Homogeneous COSMOS | \(\dot a=Ha\); \(\dot H=-\tfrac12(\rho+p)\); scalar/DM exchange and baryon, radiation, shear dilution in \`engine/cosmos.py\` | Runs the existing signed-H RHS from archive-reconstructed initial data; records component densities, normalized budgets, and \(H^2-\rho/3\) | This is homogeneous background evolution, not a global solution of the exact spherical geometry |
| Expansion/contraction | Archive permits distinct branches; current RHS evolves H as a state variable | Branch choice changes only the *initial* sign of \(H\); the subsequent sign and trajectory come from the RHS | No forced turnaround, branch splice, or bounce |
| Local phase | Archive oscillator phase \(\theta=\operatorname{unwrap}(\operatorname{atan2}(-\dot S/\omega_S,S))\) | Reports local \(\theta\) from the sampled S/Sdot history | This phase is not a mapping to cosmic time or a global cycle coordinate |
| Cosmic scalar/cycle data | Archive does not establish a complete global cycle-coordinate map | Reports raw \(\phi(t)\) in the COSMOS state | No inferred \(\lambda_{\rm cycle}\) map |

## Explicitly unpaid gates

These gates are returned in every run, not converted to numerical assumptions:

1. Local proper time \(\leftrightarrow\) cosmic time.
2. D-mode energy \(\rightarrow\) localized matter source.
3. Local source \(\rightarrow\) global cosmological exchange \(Q\).
4. Full-GR oscillaton endpoint.
5. First-principles bounce.
6. Global anisotropic three-curvature / shear completeness beyond the spherical lane.

There is no cross-lane RHS term. The runner does not call \`exchange_from_flux\` to synthesize a reciprocal physical source, and does not invent a matter-conversion efficiency, clock scale, friction law, ejection threshold, or bounce prescription.

## Known source discrepancy to audit before calling this an exact archive reproduction

The historic DOCX local-mode listing includes an \(-\eta_G D^2\dot D\) term. The current \`engine/local.py\` RHS does **not** use its \`eta_G\` field and does not include that nonlinear damping term. v1 deliberately follows the current source module instead of silently reviving a historic term; this is a source mismatch to resolve through a separate equation/provenance review, not a claim that the two formulations are identical.

Likewise, the old DOCX and later corrected COSMOS lane have different potential-normalization histories. This runner uses the current, normalized \`engine/cosmos.py\` definitions, including its archive-reconstructed \(V_0\), rather than re-running an obsolete source snapshot under a new name.

## Validation boundary and next gate

The new tests cover RK4 plumbing, independent clock axes, lane-ledger shape/provenance, initial expansion/contraction sign, the GEAR-03 parameter identity, invalid controls, and open-gate reporting. They are software contracts for short deterministic runs only.

They do **not** establish long-horizon stability, the 47.203-e-fold contraction target, the physical validity of the proxy ledgers, conservation across the two lanes, or scientific admission. Before extending the harness, the next narrow task is a source-by-source comparison of the preserved archive equations against the currently canonical local/COSMOS RHS modules, with discrepancies recorded before choosing any physical branch.
