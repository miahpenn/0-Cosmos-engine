# External literature crosswalk under the 0* admission gate — 25 October 2026

**Status:** research-only, documentation-only. This note adds a focused crosswalk to the existing [scientific literature map](../SCIENTIFIC_LITERATURE_MAP.md) and [source-evidence crosswalk](../RESEARCH_SOURCE_CROSSWALK_2026-10-15.md). It does not replace them, reopen already reviewed claims, alter production code, or change scientific-admission status.

**Active engineering target:** [Issue #4 — momentum-constraint propagation](https://github.com/miahpenn/0-Cosmos-engine/issues/4), using the SHA-verified accepted states and sector-resolved diagnostics already recorded in [PR #420](https://github.com/miahpenn/0-Cosmos-engine/pull/420). Do not repeat accepted-state captures or generic interpolation experiments without a new source-pinned discriminator.

## Admission rule for borrowed results

Every borrowed result is a *candidate*, not an input to production by default. Its evidence record must keep four gates distinct:

1. **Published-source validity.** Verify the paper/version, exact equation and assumptions, errata/supplement, and whether the cited statement is a derivation, theorem, numerical observation, or conjecture.
2. **Convention and scope mapping.** Map variables, signs, index positions, lapse/shift, metric/reference-connection convention, stress-energy definition, units, branch, and domain of validity to the pinned archive and production sources. Similar notation is not a mapping.
3. **0* eligibility.** Apply the current archived admission criterion exactly as written: \(0^*_{\mathrm{eff}}\iff t\land b\land(w\geq1)\). Preserve the authoritative archive definitions of \(t\) and \(b\); this note does not redefine them. The archive's current interpretation distinguishes \(w\geq1\) as the stated threshold for avoiding BKL chaos from \(w>1\) as the stronger active-smoothing regime. Record whether each condition holds pointwise or only asymptotically, on which branch, and whether it fails at/after turnaround. Do not silently substitute \(w>1\) for \(w\geq1\), or infer the full criterion from one passing condition.
4. **Independent engine validation.** Test the mapped claim against the current pinned implementation and appropriate analytic/manufactured controls. A passing algebraic identity, unit suite, local flux balance, or frozen-state response is evidence only for that scope; it does not by itself prove a physical 0* state, full-system conservation, or scientific admission.

A numerical method paper need not itself posit a 0* cosmology. In that case, the **proposed application** of the method to a 0*-effective state must still pass the archive eligibility and convention-mapping gates. A method cannot be admitted as a way to bypass a failed physical condition. Conversely, the 0* gate does not replace proof that the numerical method is correct.

## Priority crosswalk to current Issue #4

### A. Constraint propagation: Sarbach & Tiglio (2012)

**Source:** Olivier Sarbach and Manuel Tiglio, “Continuum and Discrete Initial-Boundary Value Problems and Einstein’s Field Equations,” *Living Reviews in Relativity* 15, 9 (2012), especially Sec. 4.3.2, Eqs. (4.73)–(4.76). [Open full review](https://link.springer.com/article/10.12942/lrr-2012-9), [DOI](https://doi.org/10.12942/lrr-2012-9).

**Potential use:** The review writes the BSSN Hamiltonian, momentum, and conformal-connection constraints as a coupled propagation system. The momentum-constraint evolution contains a Hamiltonian-gradient contribution, algebraic momentum-constraint terms, and conformal-connection-constraint coupling. Its derivation assumes the total stress-energy tensor is divergence-free.

**Required mapping before use:**
- Translate the review's \(\mathcal H\), \(\mathcal M_i\), and \(C^i\) into the engine's exact variable definitions, including whether the compared momentum is covariant \(M_r\) or contravariant \(M^r\).
- Carry the engine's lapse/shift, conformal normalization, reference-metric connection terms, and any formulation parameters through the derivation; do not copy Eq. (4.75) by symbol resemblance alone.
- Establish how the archive-locked matter-sector exchange terms combine into the *total* divergence condition. Individual sectors may exchange energy/momentum; their exchange terms must not be assumed to vanish separately.
- Compare the resulting independently derived prediction with the measured tangent on the existing accepted states, retaining signed component budgets and a negative control.

**Limit:** This review establishes a continuum propagation structure for the cited BSSN formulation under its stated conservation assumption. It does not establish the correct source signs/normalizations for this engine, its exact discrete propagation identity, or 0* eligibility of any physical trajectory. The latter must be evaluated under the archive gate independently.

### B. Covariant/reference-metric conventions: Brown (2009)

**Source:** J. David Brown, “Covariant formulations of Baumgarte, Shapiro, Shibata, and Nakamura and the standard gauge,” *Physical Review D* 79, 104029 (2009). [Journal page](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.79.104029), [DOI](https://doi.org/10.1103/PhysRevD.79.104029).

**Potential use:** Convention cross-check for coordinate-covariant BSSN variables and the distinction between determinant-evolution variants.

**Required mapping before use:** Identify which variant is represented by the exact pinned production source; document tensor/density weights, reference metric and connection-difference definition, determinant condition, and momentum-constraint additions to the connection evolution. Keep this separate from whether the engine's current discrete right-hand side matches the selected formulation.

**Limit:** Brown's formulation is not a blanket authorization to replace the current vendor/native operator or connection convention. It does not establish matter closure, a scalar-sector face flux, or 0* eligibility.

### C. Spherical origin and conservative fluid fluxes: Montero & Cordero-Carrión (2012)

**Source:** Pedro J. Montero and Isabel Cordero-Carrión, “BSSN equations in spherical coordinates without regularization: Vacuum and nonvacuum spherically symmetric spacetimes,” *Physical Review D* 85, 124037 (2012). [Journal page](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.85.124037), [DOI](https://doi.org/10.1103/PhysRevD.85.124037).

**Potential use:** Origin regularity, spherical BSSN/PIRK methodology, and the separation between flux-conservative relativistic hydrodynamics and geometric evolution.

**Limit:** It is a different evolution scheme and implementation. It does not determine how this engine discretizes the scalar \(S\)-sector spatial term, nor may its origin treatment be transplanted without matching the pinned operator and testing the same accepted states.

### D. Coupled finite-volume source balance: Gaburro, Castro & Dumbser (2021)

**Source:** Elena Gaburro, Manuel J. Castro, and Michael Dumbser, “A Well Balanced Finite Volume Scheme for General Relativity,” *SIAM Journal on Scientific Computing* 43(6), B1226–B1251 (2021). [Journal page](https://epubs.siam.org/doi/10.1137/21M1399154), [DOI](https://doi.org/10.1137/21M1399154).

**Potential use:** Design principles for balancing fluxes, nonconservative products, algebraic source terms, and equilibrium reconstruction in a fully coupled Einstein–matter scheme.

**Limit:** The paper's FO-CCZ4 + GRMHD method is not the engine's present spherical BSSN/PIRK + HLL fluid update. Its well-balanced construction is a comparator for source accounting, not evidence that HLL defines a scalar \(S\)-flux or that a lower residual identifies the physically correct operator.

## Source-level mapping already established (read-only review, pinned production source)

The following is a *source inspection*, not a new evolution and not a closure result. Engine source was checked at `4aec88d0e3ff028e0a65a8ac36b30948625df792`; the vendor gitlink is `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

- **Raw geometry constraint:** the pinned vendor's [`vacuum_bssn.py`](https://github.com/RonBibb/bb-palatini-unified-r0/blob/d6052d605673ce9d82cdc99b0df79855fa2ea215/vacuum_bssn.py) computes
  \[
  M_r^{\rm geom}=D_r A_a-\tfrac23D_rK+6A_a\chi'
  +(A_a-A_b)\left(\tfrac2r+\frac{D_rb}{b}\right),
  \quad A_b=-\tfrac12A_a,\quad \chi'=-\tfrac12\frac{D_rX}{X}.
  \]
  This returned `momentum` is geometry-only. The same vendor's `constraints()` returns raw geometric Hamiltonian, momentum, connection, and determinant witnesses.

- **Matter-corrected engine witness:** the pinned [tools/audit_coupled_lambda_constraint_propagation.py](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/tools/audit_coupled_lambda_constraint_propagation.py) constructs \(H_{\mathrm{BSSN}}=H_{\mathrm{raw}}-16\pi\rho_{\mathrm{total}}\) and \(M_r=M_r^{\mathrm{geom}}-8\pi j_{\mathrm{total}}\); its independent metric-only Hamiltonian is retained as a distinct witness. Keep this factor/sign convention pinned; do not adjust it from analogy alone.

- **Hamiltonian normalization and operator caveat:** Sarbach–Tiglio define their geometric Hamiltonian as \(H_{\mathrm{ST}}=\tfrac12(R-\tilde A_{ij}\tilde A^{ij}+\tfrac23K^2)\), and \(\mathcal H_{\mathrm{ST}}=H_{\mathrm{ST}}-8\pi G_N\rho\). If the same Ricci operand is used, the engine's raw Hamiltonian is twice \(H_{\mathrm{ST}}\), so \(H_{\mathrm{raw}}-16\pi G_N\rho_{\mathrm{total}}\) is the correspondingly doubled residual; the factors \(8\pi\) and \(16\pi\) are not in themselves inconsistent. Expanding Sarbach–Tiglio's normalized Hamiltonian-gradient term then gives
  \[
  \frac{\alpha}{6}D_r H_{\mathrm{raw}}-\frac{D_r\alpha}{3}H_{\mathrm{raw}},
  \]
  matching the coefficients written in the existing Brown scaffold. **This does not settle the Hamiltonian operand:** the Brown observer deliberately uses a metric-only Ricci reconstruction, whereas the pinned vendor's BSSN Hamiltonian depends on the independently evolved connection variable. Preserve both witnesses until the exact formulation mapping determines which enters each term; do not substitute one because its residual is smaller. Likewise, map Sarbach–Tiglio's shift-adapted time operator and explicit shift/index terms consistently before comparing to a coordinate-time tangent.

- **Connection-evolution coupling:** vendor \`lambda_l2_rhs\` includes \(+\lambda_m\alpha M_r^{\mathrm{geom}}/a\). The local [engine/v55_pirk_adapter.py](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/v55_pirk_adapter.py) adds the matter contribution \(-8\pi\lambda_m\alpha j_{\mathrm{total}}/a\) in \`lambda_l3_with_matter\`. Taken together, those explicit terms form
  \[
  \frac{\lambda_m\alpha}{a}\left(M_r^{\mathrm{geom}}-8\pi j_{\mathrm{total}}\right).
  \]
  This is a useful, exact source-level correspondence to the matter-corrected \(M_r\) witness, **not proof** that the full nonlinear constraint-propagation tangent closes.

- **Remaining derivative ledger:** the pinned vendor matter-coupled geometry RHS also contributes the anisotropic-stress source \(-16\pi\alpha(p_r-p_t)/3\) to \(\dot A_a\) and \(+4\pi\alpha(\rho+p_r+2p_t)\) to \(\dot K\). A complete \(\dot M_r\) comparison must propagate these through the discrete derivative and the remaining metric/connection coefficients, differentiate the \(-8\pi j_{\mathrm{total}}\) term, include the Hamiltonian-gradient and connection-constraint terms, and use the actual shift/lapse and reference-metric convention. It must pair these with the **total** sector exchange/divergence budget, not the fluid flux budget alone.

- **Scope limit of the existing Brown observer:** [tools/audit_brown_vacuum_momentum_propagation.py](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/tools/audit_brown_vacuum_momentum_propagation.py) labels the Brown Eq. (11b) prediction as a lower-index \(M_r\) candidate and uses a metric-only Hamiltonian in the \(H,H_r\) slots. Its \`vacuum_rhs_rates\` are built from vendor vacuum RHS blocks, without the local matter L3 source insertion. It therefore cannot be treated as the full matter-coupled acceptance oracle requested by Issue #4. Preserve it as a useful geometric scaffold; derive/extend the matter-coupled identity rather than upgrading its status.

- **Exchange caution:** pinned [engine/v55_matter.py](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/v55_matter.py) contains the archive-locked scalar-\(\phi\)/DM exchange pair. Algebraic equal-and-opposite pair entries alone do not prove that the implemented stage RHS and stress projections satisfy \(\nabla_\mu T^{\mu\nu}_{\mathrm{total}}=0\) with the required normalization, index positions, and metric-work terms. That must be checked against the existing signed accepted-state budgets.

**Best next bounded discriminator:** reassemble the matter-corrected \(M_r\) directional tangent from the already pinned accepted-state arrays, separating (i) derivatives of the \(\dot A_a,\dot K\) rates, (ii) time variation of metric coefficients in \(M_r^{\mathrm{geom}}\), (iii) \(-8\pi\dot j_{\mathrm{total}}\), (iv) the independent Hamiltonian-gradient and connection-constraint prediction, and (v) sector-exchange/divergence terms. Compare signed terms on the same saved states and report centre cells separately from fixed-radius/annular witnesses. No new trajectory is needed unless this reassembly exposes a specific missing operand.


## Existing accepted-state artifact cross-check — no rerun

Read-only analysis of the three existing artifacts from [workflow run 40207821733](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40207821733) was performed with report and profile hashes checked against each archive manifest. The t=0.60 accepted-state snapshots were reused unchanged.

| Resolution | Profile artifact SHA-256 | Direct covariant \(M_r\) candidate residual at fixed \(r=2\) | At fixed \(r=8\) |
|---:|---|---:|---:|
| 60 | \`6aade0ced9ce41e472a491f620782b3569926c1193f0b2d17220bcd60414dd4f\` | \(+8.9752\times10^{-6}\) | \(-1.1941\times10^{-6}\) |
| 120 | \`84a0f7c8c37a6fbc57e8106b16aeba583813133649f61f400238d579a245b84f\` | \(+6.1218\times10^{-7}\) | \(-7.3282\times10^{-8}\) |
| 240 | \`b0ddbacbc0803c3c07bc28d3d14e3a45a9e23275224f0ba8d5913023c088b4be\` | \(+4.0296\times10^{-8}\) | \(-4.5908\times10^{-9}\) |

These are fixed **coordinate** radii. The absolute residual decreases at approximately fourth order: observed orders \(3.87,3.93\) at \(r=2\), and \(4.03,4.00\) at \(r=8\). This is consistent with a high-order discrete-truncation trend, but is not proof that the continuum source/index mapping is correct. The N=240 relative RMS residual remains \(0.657\) over the full domain and \(0.818\) for \(r\geq5\), so the candidate is not accepted as closed.

### What the large source-mapping residual does—and does not—say

The accepted-state decomposition reports the “actual matter source composite” as the explicit matter injection into the geometric momentum tangent minus \(8\pi\dot j_{\mathrm{total}}\). That object is **not automatically the same isolated term** as the provisional \(-8\pi\alpha Q_r\) divergence contribution: the complete propagation identity also contains matter-corrected Hamiltonian/momentum operands, connection terms, index/time-operator conversion, and the corresponding convention-specific coefficients.

At \(r=2\), the geometric Brown vacuum candidate and actual matter composite nearly cancel:
- N=60: \(-5.7350\times10^{-5} + 5.8301\times10^{-5} \approx 9.51\times10^{-7}\).
- N=120: \(-5.7744\times10^{-5} + 5.7812\times10^{-5} \approx 6.80\times10^{-8}\).
- N=240: \(-5.8347\times10^{-5} + 5.8351\times10^{-5} \approx 4.92\times10^{-9}\).

Meanwhile the independently measured-minus-predicted residual remains larger than those sums (values in the table above). So the right conclusion is **not** “the large source-mapping residual proves conservation is broken.” It identifies a mapping that has not yet been justified as an independent physical law; the full signed propagation equation and its operand conventions still need to be derived and compared as a whole.

For cross-checking bookkeeping only, the exact discrete \(M_r^{\mathrm{geom}}\) chain-rule maximum error is \(1.87\times10^{-13}\), \(1.12\times10^{-12}\), and \(7.13\times10^{-13}\) at N=60/120/240. The decomposition of the total tangent as \(\dot M_r^{\mathrm{geom}}-8\pi\dot j_{\mathrm{total}}\) closes to \(1.55\times10^{-16}\), \(3.84\times10^{-16}\), and \(5.22\times10^{-16}\). Those are strong implementation-bookkeeping checks, not continuum conservation or 0* admission.


## Immediate research disposition

The highest-value next step is **not** another broad paper search or a new evolution. It is an equation-by-equation mapping of Sarbach–Tiglio Eqs. (4.73)–(4.76) onto the pinned engine formulation, with Brown's reference-metric conventions used to resolve notation. Then compare that prediction against the existing accepted-state signed tangent and the already-recorded scalar/matter source budgets.

For any imported physical premise, first record its archive source and test the full \(t\land b\land(w\geq1)\) condition in the actual regime claimed. A pointwise/asymptotic mismatch or failure beyond turnaround must remain visible; it cannot be averaged away by a successful numerical-method test.

## Current non-claims / fail-closed status

- This crosswalk does not prove that a paper is novel to the project, that its full method is compatible, or that another team is solving the exact same open problem.
- The accepted-state discrete scalar-current tangent closing near roundoff is implementation bookkeeping; it is not total stress-energy conservation.
- The pointwise and finite-volume \(S\)-sector stress-divergence comparisons are different discrete operators and must not be selected by whichever gives the smaller residual.
- Do not infer individual-sector conservation from conservation of the total tensor, or the reverse.
- Do not change production equations, source signs, normalizations, flux choices, or archive assumptions from this literature note.

**Status remains:** \(\texttt{full\_matter\_coupled\_propagation\_identity\_closed=false}\), \(\texttt{total\_stress\_energy\_sector\_conservation\_closed=false}\), \(\texttt{scientific\_admission=false}\).

This document is a scoped literature crosswalk, not a renewed scientific review of already admitted archive claims.