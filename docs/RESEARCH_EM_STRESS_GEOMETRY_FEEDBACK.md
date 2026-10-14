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
