# Conditional canonical D-mode trace audit

**Status:** research-only conditional derivation; not an archival action, not a production change, and not scientific admission.  
**Date:** 15 October 2026 UTC.

## Question

Can the GEAR-63/67E curvature candidates be recovered from the signed inverted-potential D field by treating it as a minimally coupled canonical scalar? This note computes the stress-energy tensor and its trace under that explicit assumption. It does not claim that the missing primary GEAR action has been found.

## Assumptions and conventions

Take metric signature \((- + + +)\), the Einstein equation
\[
G_{\mu\nu}=\kappa T_{\mu\nu},\qquad \kappa=8\pi G,
\]
and the standard Ricci convention \(G_{\mu\nu}=R_{\mu\nu}-\tfrac12 g_{\mu\nu}R\). Consider only a minimally coupled canonical real scalar D with
\[
S_D=\int d^4x\sqrt{-g}\left[-\frac12 g^{\mu\nu}\nabla_\mu D\nabla_\nu D-V_D(D)\right],
\qquad
V_D(D)=-\frac12\mu^2D^2.
\]
Here \(\mu\) is the growth rate (the archive's \(\alpha\)); current V5.5 code uses the dimensionless choice \(\mu=1\) in its signed potential \(V_D=-D^2/2\).

These are declared assumptions for a conditional calculation, not evidence that this action, signature, or curvature-coupling choice is the original archive's intended model. A cosmological constant, other matter, a surface layer, or non-minimal terms would change the complete source/trace.

## Metric variation and trace

Metric variation of the stated minimal action gives
\[
T_{\mu\nu}=\nabla_\mu D\nabla_\nu D
-g_{\mu\nu}\left[\frac12(\nabla D)^2+V_D\right].
\]
Taking the four-dimensional trace,
\[
T\equiv g^{\mu\nu}T_{\mu\nu}=-(\nabla D)^2-4V_D.
\]
For a spatially homogeneous field, \((\nabla D)^2=-\dot D^2\), hence
\[
\rho_D=\frac12\dot D^2+V_D,\qquad
p_D=\frac12\dot D^2-V_D,\qquad
T_D=-\rho_D+3p_D=\dot D^2-4V_D.
\]
Evaluate this stress tensor on the archive's homogeneous growing-branch relation \(\dot D=\mu D\), the signed inverted potential implies
\[
\rho_D=0,\qquad p_D=\mu^2D^2,\qquad T_D=3\mu^2D^2.
\]
This is an algebraic evaluation on that branch relation, not a claim that \\(\\dot D=\\mu D\\) solves the full curved/expanding production scalar equation when extrinsic-curvature, lapse-gradient, or other-background terms are present.\n\nThus zero canonical energy density does **not** mean zero stress-energy or zero curvature: the pressure and trace are nonzero. If D is the only source and there is no cosmological constant, taking the trace of Einstein's equation gives
\[
-R=\kappa T_D
\quad\Longrightarrow\quad
R=-3\kappa\mu^2D^2=-24\pi G\mu^2D^2.
\]
This result is conditional on the stated Einstein/Ricci sign convention and the minimal canonical action.

## Comparison with the archived curvature candidates

The archive notes report distinct hypotheses, which should not be conflated:

| Route recorded in the research archive | Reported curvature and reduced force | Status |
|---|---|---|
| GEAR-67E's converted/dust branch | \(R=4\pi G\alpha^2d^2\), so \((R/6)d=(2\pi/3)G\alpha^2d^3\) | Algebraic closure conditional on the converted/dust source |
| GEAR-67E's unconverted/inverted branch | \(R=-16\pi G\alpha^2d^2\), so \((R/6)d=-(8\pi/3)G\alpha^2d^3\) | Different branch/source assumption; physical origin remains open |
| GEAR-68's stated test datum | \(\rho=\tfrac12\alpha^2d^2\), \(R\approx8\pi G\alpha^2d^2\), giving \(F_{\rm geom}\approx(4\pi/3)G\alpha^2d^3\) | The audit records this as a supplied estimate/hypothesis, not a controlled derivation |
| Minimal canonical inverted scalar, derived above | \(R=-24\pi G\mu^2D^2\) on its homogeneous growing branch | Conditional trace result; not equal to the above routes |

The minimal-canonical result does **not** pick the archive's \(2\pi/3\) or \(4\pi/3\) cubic coefficient. It shows why an energy-density-only shortcut is insufficient: Einstein's scalar-curvature trace depends on \(T=-\rho+3p\), not on \(\rho\) alone. A conversion to pressureless dust changes the stress-energy model and must be independently derived; it cannot be silently substituted for the canonical scalar.

## What is still missing

The minimal action above has no explicit \(\xi R D^2\) term, so its Euler–Lagrange equation does not contain a separately postulated \(R D/6\) force. Such a term requires an explicit non-minimal curvature coupling in the action, or a distinct effective closure justified from the continuum model. The algebraic identity \(F=(R/6)D\) does not establish that coupling or fix its coefficient.

This matches the stated limits in **GEAR-67E — Rd/6 Algebraic Shortcut Finding**, **GEAR-61 Amendment — Retarded Conversion and Separated Reaction**, and **GEAR-68 Audit — Exact Local Self-Gravity and Retarded Memory**: the physical origin of the curvature closure and the localized remnant/source projection remain open. The current production scalar projection uses \(V_D=-D^2/2\), but code convention is not proof of the intended historical action.

## Disposition / next gate

- **Established:** under the explicitly declared minimal canonical action, the homogeneous growing D branch has \(\rho_D=0\), \(p_D=\mu^2D^2\), and \(R=-24\pi G\mu^2D^2\) under the stated Einstein/Ricci convention.
- **Not established:** that this is the source action intended by GEAR; a dust conversion; a non-minimal \(R D^2\) coupling; the physical \(F_{\rm geom}\) coefficient; or COSMOS/local density provenance.
- **Required next evidence:** an authoritative action-bearing archive source/revision with equation/page, metric and curvature sign conventions, field normalization and any \(\xi R D^2\) term; then an independent variation and source projection. If no such source exists, an independently reviewed scope decision must explicitly leave the archive COSMOS lane as a normalized reference and S/D normalization unresolved.

No production equation, source coefficient, coupling, vendor code, or campaign setting is changed. Issue #210 remains open and scientific admission remains **NOT ADMITTED**.
