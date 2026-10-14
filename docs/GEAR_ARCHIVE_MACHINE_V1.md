# GEAR Archive Machine v1 — independent-lane build

**Status:** research-only scaffold; not a production or scientific-admission claim  
**Source branch:** \`research/gear-archive-machine-v1-2026-10-13\`  
**Preserved references:** \`CYCLIC-COSMOS ENGINE-.docx\`, \`GEAR_FULL_CYCLIC_COSMOS_ENGINE_ARCHIVE_SNAPSHOT_2026-10-02.md\`, and \`GEAR_MASTER_ARCHIVE_HOUSEKEEPING_INDEX_2026-10-03_v2.md\`.

## Purpose and scope

This starts the older, archive-derived two-lane machine as a versioned runner. The default local mode equation is transcribed from the preserved archive DOCX; a named comparison option invokes the current \`engine/local.py\` RHS. The homogeneous COSMOS lane uses the corrected, normalized \`engine/cosmos.py\` equations. This is a deliberate source-crosswalk harness, not a byte-for-byte reproduction of the historical DOCX program.

Run the small deterministic harness with:

\`\`\`bash
pytest -q tests/test_gear_archive_machine_v1.py
\`\`\`

The Python entry point is \`research/gear_archive_machine_v1.py:run_machine\`. The local lane uses \`tau_local\`; the homogeneous lane uses \`t_cosmos\`. They can use different final times and step sizes. The runner does not align them, compare equal numerical times as if physically simultaneous, or use either clock to rescale the other.

## Source-to-kernel contract

| Sector / diagnostic | Archive equation or source map | Current v1 behavior | Interpretation limit |
|---|---|---|---|
| Local stable mode S | GEAR-03: \(\omega_S^2=\kappa+g\), \(\ddot S=-\omega_S^2 S\) | Default \`archive_docx\` RHS transcribes the archived reduced mode equation; \`current_source\` calls \`engine.local.local_rhs\` with \(H=0\) | Reduced mode-space system, not the full radial Einstein–scalar PDE |
| Local D mode and geometric barrier | \(\alpha_D^2=g-\kappa\), \(\ddot D=\alpha_D^2D-\eta_GD^2\dot D-\frac{2\pi}{3}G\alpha_D^2D^3\) in the historical DOCX reduced RHS | Default \`archive_docx\` mode includes the archive's damping term; \`current_source\` intentionally omits it. The runner checks \(\alpha_D^2=g-\kappa\). | A source-version distinction, not evidence for D-to-matter conversion |
| Local energy/activity | Raw saddle energy \(E_D=\dot D^2/2-\alpha_D^2D^2/2\); positive kinetic/excursion activity are separate diagnostics; \(F_{\rm proxy}\sim|D\dot D|\) | Emits \(E_S\), raw \(E_D\), kinetic/excursion proxies, generated-activity proxy, and transport-activity proxy | Proxies are not particle abundance, ejection proof, a physical flux, or an exchange source |
| Homogeneous COSMOS | Scalar + coupled DM, baryons, radiation and shear background equations | Uses the corrected normalized \`engine.cosmos\` state and signed-\(H\) RHS; records component density/budget and \(H^2-\rho/3\) witnesses | Not the same as the historical branch-explicit DOCX implementation or a global solution of the spherical spacetime |
| Expansion/contraction | Historical lane ran expansion and contraction branches separately | Branch selector sets only the initial sign of \(H\); the current signed-\(H\) ODE then evolves the state | No forced turnaround, branch splice, or bounce |
| Local phase | \(\theta=\operatorname{unwrap}(\operatorname{atan2}(-\dot S/\omega_S,S))\) | Reports local \(\theta\) from sampled S/\(\dot S\) history | Not a local-to-cosmic time map or global cycle coordinate |

## Explicitly unpaid gates

These gates are returned in every run, not converted to numerical assumptions:

1. Local proper time \(\leftrightarrow\) cosmic time.
2. D-mode energy \(\rightarrow\) localized matter source.
3. Local source \(\rightarrow\) global cosmological exchange \(Q\).
4. Full-GR oscillaton endpoint.
5. First-principles bounce.
6. Global anisotropic three-curvature / shear completeness beyond the spherical lane.

There is no cross-lane RHS term. The runner does not call \`exchange_from_flux\` to synthesize a reciprocal physical source and does not invent a matter-conversion efficiency, clock scale, friction law, ejection threshold, or bounce prescription.

## Verified differences retained rather than silently merged

**Local damping:** the preserved DOCX reduced \`local_rhs\` contains \(-\eta_GD^2\dot D\). Current \`engine/local.py\` accepts \`eta_G\` but does not use it, and instead offers an optional \(-3H\dot D\) term. In the runner's comparison lane, \(H=0\), so a direct same-state test isolates the missing historical damping term. This update makes the archived equation the default *in this research harness*; it does not edit or declare the production module wrong.

**COSMOS potential normalization and state formulation:** the older DOCX potential and branch-explicit \(\pm\sqrt{\rho/3}\) integration differ from the later corrected normalized \`engine/cosmos.py\` potential (including \(V_0\)) and its evolved signed-\(H\) equation. v1 explicitly chooses the corrected current COSMOS module and exposes this as a source-version boundary. It does not claim the old and corrected COSMOS formulations are identical.

## Validation boundary and next gate

Tests now cover RK4 plumbing, independent clock axes, ledger shape/provenance, initial expansion/contraction sign, the GEAR-03 parameter identity, signed damping contribution, parity with the current local RHS when \(\eta_G=0\), comparison-mode selection, invalid controls and open-gate reporting. These are small deterministic software contracts only.

They do **not** establish long-horizon stability, the 47.203-e-fold contraction target, the physical validity of proxy ledgers, conservation across lanes, the full archived run, or scientific admission. Next gate: independent equation review of the archived local damping/barrier and corrected COSMOS normalization/branch choices, then a bounded low-cost run of each lane separately. No coupled campaign should be scheduled before the source/formulation choices are accepted.
