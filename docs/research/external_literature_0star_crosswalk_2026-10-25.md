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

- **Matter-corrected engine witness:** the pinned [`tools/audit_coupled_lambda_constraint_propagation.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/tools/audit_coupled_lambda_constraint_propagation.py) constructs (H_{\rm BSSN}=H_{\rm raw}-16\pi\rho_{\rm total}) and (M_r=M_r^{\rm geom}-8\pi j_{\rm total}); its independent metric-only Hamiltonian is retained as a distinct witness. Keep this factor/sign convention pinned; do not adjust it from analogy alone.


- **Hamiltonian normalization and operator caveat:** Sarbach–Tiglio define their geometric Hamiltonian (H_{ST}) with an overall factor (1/2), and (ℒ_{ST}=H_{ST}-8\pi\rho). The engine's raw Hamiltonian expression is twice that geometric quantity, so the source-pinned corrected witness (H_{raw}-16\pi\rho_{total}) is the correspondingly doubled residual; the factors (8\pi) and (16\pi) are not in themselves inconsistent. In the expanded momentum-propagation RHS, this scaling yields the (+(\alpha/6)D_r H_{raw}-(D_r\alpha/3)H_{raw}) coefficients used in the existing Brown scaffold. **But** its Hamiltonian operand must match the formulation being mapped: the existing Brown observer deliberately uses a metric-only Ricci reconstruction, whereas the BSSN Hamiltonian used by the pinned vendor's geometric constraints includes the independently evolved connection variable. Preserve both witnesses until the source/formulation mapping establishes which enters each term; never substitute one because its residual is smaller. Sarbach–Tiglio's shift-adapted time operator and the explicit shift/index terms must likewise be mapped consistently when comparing with a coordinate-time tangent.


- **Connection-evolution coupling:** vendor `lambda_l2_rhs` includes (+\lambda_m\alpha M_r^{\rm geom}/a). The local [`engine/v55_pirk_adapter.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/v55_pirk_adapter.py) adds the matter contribution (-8\pi\lambda_m\alpha j_{\rm total}/a) in `lambda_l3_with_matter`. Taken together, those explicit terms form (\lambda_m\alpha(M_r^{\rm geom}-8\pi j_{\rm total})/a). This is a useful, exact source-level correspondence to the matter-corrected (M_r) witness, **not proof** that the full nonlinear constraint-propagation tangent closes.

- **Remaining derivative ledger:** the pinned vendor matter-coupled geometry RHS also contributes the anisotropic-stress source (-16\pi\alpha(p_r-p_t)/3) to (\dot A_a) and (+4\pi\alpha(\rho+p_r+2p_t)) to (\dot K). A complete (\dot M_r) comparison must propagate these through the discrete derivative and the remaining metric/connection coefficients, differentiate the (-8\pi j_{\rm total}) term, include the Hamiltonian-gradient and connection-constraint terms, and use the actual shift/lapse and reference-metric convention. It must pair these with the **total** sector exchange/divergence budget, not the fluid flux budget alone.

- **Scope limit of the existing Brown observer:** [`tools/audit_brown_vacuum_momentum_propagation.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/tools/audit_brown_vacuum_momentum_propagation.py) labels the Brown Eq. (11b) prediction as a lower-index (M_r) candidate and uses a metric-only Hamiltonian in the (H,H_r) slots. Its `vacuum_rhs_rates` are built from vendor vacuum RHS blocks, without the local matter (L3) source insertion. It therefore cannot be treated as the full matter-coupled acceptance oracle requested by Issue #4. Preserve it as a useful geometric scaffold; derive/extend the matter-coupled identity rather than upgrading its status.

- **Exchange caution:** pinned [`engine/v55_matter.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/4aec88d0e3ff028e0a65a8ac36b30948625df792/engine/v55_matter.py) contains the archive-locked scalar-(\phi)/DM exchange pair. Algebraic equal-and-opposite pair entries alone do not prove that the implemented stage RHS and stress projections satisfy (\nabla_\mu T^{\mu\nu}_{\rm total}=0) with the required normalization, index positions, and metric-work terms. That must be checked against the existing signed accepted-state budgets.

**Best next bounded discriminator:** reassemble the matter-corrected (M_r) directional tangent from the already pinned accepted-state arrays, separating (i) derivatives of the (A_a,K) rates, (ii) time variation of metric coefficients in (M_r^{\rm geom}), (iii) (-8\pi\dot j_{\rm total}), (iv) the independent Hamiltonian-gradient and connection-constraint prediction, and (v) sector-exchange/divergence terms. Compare signed terms on the same saved states and report centre cells separately from fixed-radius/annular witnesses. No new trajectory is needed unless this reassembly exposes a specific missing operand.

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