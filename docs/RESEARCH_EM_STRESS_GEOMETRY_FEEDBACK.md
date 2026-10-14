# Research note: electromagnetic stress feedback into geometry

**Status:** diagnostic-only, opt-in, research branch. This note describes an audit, not a new evolution law or a claim of global energy conservation.

## Hypothesis to test

In the engine's radial-electric sector, the Maxwell field contributes to the stress-energy tensor used by the gravitational source. The useful numerical discriminator is therefore two-part:

1. Does the direct electromagnetic stress contribution to the actual PIRK primary L3 curvature source reconstruct from the same-stage projected stress, with geometry and non-electromagnetic state held fixed?
2. How does that verified source sit alongside the metric contribution to the electromagnetic field-energy tangent at the actual start and predictor stages?

The second quantity is not expected to equal either curvature-source component algebraically. The source increment acts on the geometry variables; the field-energy metric tangent is evaluated from the explicit \(a,b,X\) geometry RHS. Their time/stage relationship must be measured, not imposed.

## Source identities under the current model conventions

For a pure radial electric field, the active projection is

\[
\rho_{\rm EM}=\frac{E_{\hat r}^2}{8\pi},\qquad
p_{r,\rm EM}=-\rho_{\rm EM},\qquad
p_{t,\rm EM}=+\rho_{\rm EM},\qquad j_{\rm EM}=0.
\]

The current engine's primary L3 matter-source conventions give

\[
\dot K_{\rm EM}=4\pi\alpha(\rho_{\rm EM}+p_{r,\rm EM}+2p_{t,\rm EM})
=8\pi\alpha\rho_{\rm EM},
\]

\[
\dot A_{a,\rm EM}=-\frac{16\pi}{3}\alpha
(p_{r,\rm EM}-p_{t,\rm EM})
=\frac{32\pi}{3}\alpha\rho_{\rm EM}.
\]

The audit obtains the direct source increment by taking the production primary L3 result and subtracting a same-stage matter-only primary L3 evaluation with the same geometry, lapse, scalar fields and non-EM matter, but without the EM stress projection. It then independently reconstructs the increments above from the captured EM \(\rho,p_r,p_t\) profiles. The source-split software identity fails closed at an operand-scaled float64 roundoff bound. This frozen-stage subtraction is an operator-level counterfactual, not a second evolved spacetime.

The independent EM field-energy functional has metric tangent

\[
\left.\dot U_{\rm EM}\right|_{\rm metric}
=\frac{\Delta r}{2}\sum_i w_i Q_i^2
\left(\frac{\dot X_i}{X_i}
+\frac{1}{2}\frac{\dot a_i}{a_i}
-\frac{\dot b_i}{b_i}\right).
\]

Keep this distinct from the face-charge/current tangent and the cell-centred Lorentz power. The latter pair still has the known face/cell projection mismatch; the HLL dissipative jump is far too small to explain that mismatch in the current witnesses.

## Interpretation limits

A successful source split shows that the EM stress is included in the active discrete curvature source with the expected signs and conventions. It does **not** by itself establish that curvature feedback explains the measured energy remainder. The metric power is a downstream term in the evolving geometry; charged and neutral arms also differ dynamically, and the global scalar/gravity/Killing-energy ledger is not yet closed.

All physical residuals remain signed and unthresholded. No production equations, defaults, carrier assumptions, normalization choices or vendor files are changed by this diagnostic.


## Cross-check against the local GEAR archive (14 October 2026)

Sources reviewed in the local research library include `GEAR_DYNAMIC_INTERFACE_BLUEPRINT_AUDIT_2026-10-12.md`, `GEAR_V5_5_MASTER_CHECKPOINT_2026-10-05.md`, `GEAR-66-Force-Q-Bridge-Evidentiary-Sep29-2026.pdf`, `KERNEL_RESEARCH_REVIEW_2026-10-09.md`, and the provenance note `SOLVER-SOURCE-4-bianchi9-py-recovered-from-the-record-and-six-still-missing-Sep7-2026`. These are design/archive evidence, not replacements for exact-head source and runtime measurements.

### What transfers to this diagnostic

- **Keep the V5.5 ownership model intact.** The checkpoint specifies one evolving spherical metric and one total stress-energy tensor, with no phenomenological interface source or fitted feedback coefficient. The blueprint audit distinguishes this unified-spacetime architecture from the older two-lane local/COSMOS experiments. A separate homogeneous `CosmosState` or a guessed exchange `Q` should not be wired into the production step as a way to force the EM remainder to close.
- **GEAR-66 is an accounting discipline, not an EM law.** Its evidentiary chain is local force → dissipated power → event energy → event/ejection rate → cosmological exchange `Q`. The event population and production history, as well as normalization against the independently obtained target, remained unresolved in that archive. This is a useful precedent for deriving every intermediate term and forbidding a fitted closure; it does not supply the radial-Maxwell carrier, current law, normalization, or a value for an EM source.
- **Geometry and the interface ledger need their own witnesses.** The blueprint audit records the older V5.4 moving-worldtube prototype as diagnostic-only: it was resolution-sensitive and did not satisfy the required two-sided Israel matching. A fixed outer numerical face is not automatically a physical interface. In the present EM ledger, the outer summand of the discrete charge tangent must therefore remain distinct from physical Poynting flux unless a separate finite-domain derivation establishes that identification. The charged-dust boundary flux, field metric tangent, Lorentz power, and Misner–Sharp/gravity terms also remain separately signed.
- **Do not transfer a diagnostic coordinate derivation into an exchange law.** The archived proper-radius derivation can be translated to the current engine's radial metric as (\gamma_{rr}=a/X^2) and (d\ell/dr=\sqrt{a}/|X|); the corresponding profile is a geometric observable (tracked separately in [PR #195](https://github.com/miahpenn/0-Cosmos-engine/pull/195)), not a derivation of the charge-current pairing or of global energy exchange.
- **Respect source provenance.** The saved GEAR-41/Bianchi-IX solver source labels itself as transcribed from a conversation record and not rerun in the present recovery. It covers that canonical shear test problem only. It must not be treated as executable evidence for the live V5.5 EM path. As with the current audit, a mathematical/archive claim, source-level implementation, algebraic reconstruction, and live stage result are distinct evidence classes.

### Consequence for the next review

The actionable question remains whether the already-present EM stress/current operators and geometry RHS form the intended compatible *discrete* exchange at matched stages. The current same-stage (K/A_a) source reconstruction checks the direct stress-to-curvature contribution; it does not prove that this contribution explains the energy remainder or closes the global ledger. Retain the signed, unthresholded physical residual, use the saved raw operands for independent reconstruction, and require a separate scientific review before any physical interpretation is admitted. No production source, evolution equation, tolerance, default, carrier assumption, or vendor pin is changed by this archive cross-check.
