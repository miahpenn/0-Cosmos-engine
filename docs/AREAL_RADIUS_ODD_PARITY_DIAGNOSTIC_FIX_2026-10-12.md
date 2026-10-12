# Areal-radius regular-centre parity correction — diagnostic report

Date: 2026-10-12  
Branch: \`fix/areal-radius-odd-parity-2026-10-12\`  
Strict-admission base: \`d5a5b97eadcd5c565722e61262a22195baf3c1d5\`

## Defect found

The areal radius in the spherical metric is
\[
R(r) = \frac{r\sqrt{b}}{X}.
\]
For a regular centre, \(b\) and \(X\) extend as even functions of signed radius, so \(R\) is odd and \(dR/dr\) is even. Two production diagnostic paths incorrectly differentiated \(R\) using the even-function ghost rule:
- \`engine/worldtube.py::invariant_chi\`;
- \`engine/production_kernel.py::diagnostics\`.

For the flat-space control \(a=b=X=1\), this means \(R=r\) exactly, but the old parity choice yields \(R'_0=7/12\) and \(R'_1=13/12\), rather than \(R'=1\). It therefore makes the first-cell trapping indicator \(\chi\) about \(49/144\approx0.3403\) in flat, static space, where the exact value is 1. That closely matches the repeated \(0.356\)-scale minimum in the older ledgers and identifies it as a diagnostic artifact, not evidence of physical trapping.

## Change

A shared helper, \`areal_radius_derivative(grid, geometry)\`, now differentiates \(R\) with odd parity. Both invariant-\(\chi\) and production-kernel trapping/Misner–Sharp diagnostics use this helper. The fix changes neither the evolved state nor the production source terms.

Regression tests in \`tests/test_areal_radius_odd_parity.py\` assert that:
- \(R=r\) in flat space;
- the derivative is 1 at every cell, including the first two cells;
- flat, static \(\chi\) is 1 across the entire domain;
- the old even-parity stencil fails the first-cell controls.

## Validation

Frozen production checkpoints are being audited by [workflow 38283346865](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38283346865). The audit compares the old and corrected \(R'\), \(\chi\), and Misner–Sharp mass for N=40/60/80 at t=2, 6, and 12. It verifies exact equality of the old and new \(\chi\) profiles for cells \(i\ge2\), where the derivative stencil no longer uses centre ghost points.

This is an invariant-diagnostic repair. It does not address the Hamiltonian/momentum constraint drift: those equations are evaluated independently of \(R'\). Any previous trapping minima/root claims that depend on the first two cells must be recalculated using the corrected parity before interpretation.

## Research consistency

This follows the regular-centre parity treatment required by spherical BSSN/PIRK discretizations. The published spherical-coordinate formulations emphasize that origin/axis handling is integral to accurate numerical relativity, even when a method avoids an explicit centre regularization step: Montero & Cordero-Carrión, [Phys. Rev. D 85, 124037 (2012)](https://doi.org/10.1103/PhysRevD.85.124037); Baumgarte et al., [Phys. Rev. D 87, 044026 (2013)](https://doi.org/10.1103/PhysRevD.87.044026). These references motivate the parity discipline; they do not validate the rest of the production kernel.

## Decision

Keep this diagnostic fix isolated until the regression and frozen-checkpoint report pass. Do not change \`main\` or infer that the production constraint problem is solved. Continue the central Hamiltonian/momentum investigation separately, with the corrected invariant diagnostics as the new basis for any trapping or Misner–Sharp interpretation.


## Frozen-checkpoint result — completed

Workflow: [38283346865](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38283346865)  
Artifact: \`areal-radius-odd-parity-audit\`

The full repository test suite and the new dedicated parity regression both passed. The saved-state audit covered all nine checkpoints (N=40/60/80 at t=2,6,12) from the original admitted campaign. For every state, the old/new \(\chi\) values were bitwise identical for all cells \(i\ge2\); the change is confined exactly to the first two centre cells.

At t=12:
- N=40: the reported \(\min\chi\) changes from \(0.35604\) at r=1 to \(0.78863\) at the outer cell r=79.
- N=60: it changes from \(0.35635\) at r=2/3 to \(0.78684\) at r=79.333.
- N=80: it changes from \(0.35619\) at r=0.5 to \(0.78595\) at r=79.5.

The first-cell \(dR/dr\) correction is large and consistent with the parity derivation:
- N=80 at t=12: old \(R'_0=0.62607\), corrected \(R'_0=1.05840\); old \(\chi_0=0.35619\), corrected \(\chi_0=1.01797\).
- At the adjacent cell, old \(R'_1=1.15666\) becomes \(1.06937\), and old \(\chi_1=1.14655\) becomes \(0.98001\).

The old first-cell Misner–Sharp mass is also strongly affected: N=80 at t=12 changes from \(0.16859\) to \(-0.00471\). This does not alter the reported worldtube mass around r=10, since that derivative stencil does not use the centre ghosts. Neither old nor corrected profile has a sign-crossing marginal root in these nine saved states.

Thus the previously repeated \(\min\chi\approx0.34\)–0.36 was a centre-stencil artifact. The corrected \(\min\chi\) now occurs near the outer edge of the domain; it is a separate witness that needs interpretation in light of the finite-radius boundary, not evidence for a new trapped surface. No trajectory or constraint equation was changed.

## Updated decision

Use the corrected areal-radius parity for subsequent trapping and full-domain Misner–Sharp profiles. Mark old first-two-cell \(\chi\) minima and masses as superseded diagnostics. Keep the Hamiltonian/momentum root-cause investigation separate: this fix does not change H/M/C profiles at any cell, so it is not their repair.
