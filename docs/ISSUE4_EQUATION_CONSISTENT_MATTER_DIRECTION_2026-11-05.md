# Issue #4 equation-consistent matter-direction control

**Date:** 2026-11-05  
**Class:** `equation_consistent_matter_direction`  
**Status:** analytic construction specification and constructor/test implementation; not yet projected, archived, loaded, or evaluated  
**Scope:** research-only, zero-step, frozen flat geometry; all admission and closure flags remain false.

## Purpose

Create a second, independent analytic matter-direction control. The existing
`manufactured_regular_family` archive remains immutable and is not reused as
a source of tangent values. This new direction aligns scalar conjugate
momenta and fluid velocities with the instantaneous matter equations at a
zero-velocity base state.

## Equation-consistent matter direction

The base profiles retain the declared Gaussian scalar and density amplitudes
and widths of the arbitrary-family control. At the base, all scalar conjugate
momenta and fluid velocities vanish.

For (f(r)=A e^{-r^2/w^2}), the independently evaluated continuum operators are

[
f'(r)=-\frac{2r}{w^2}f(r),\qquad
\Delta f=f''+\frac{2}{r}f'
=\left(\frac{4r^2}{w^4}-\frac{6}{w^2}\right)f.
]

The pinned scalar equations on the exact flat slice ((\alpha=1,\beta=K=0))
give

[
\dot S=PS=0,\quad \dot{PS}=\Delta S-S,
\]
[
\dot D=PD=0,\quad \dot{PD}=\Delta D+D,
\]
[
\dot\phi=\Pi=0,\quad
\dot\Pi=\Delta\phi-V'_{\rm cosmos}(\phi)+\beta_{DM}\rho_{DM}.
]

Here (V'_{\rm cosmos}=-V_0 e^{-\phi}), (V_0=8.242522415500654\times10^{-5}),
and (\beta_{DM}=-0.04), all from the pinned scalar source.

For each fluid at rest, spherical continuity gives (\dot\rho=0) at the
base. The radial Euler/Valencia momentum equation is

[
(\rho+p)\dot v^r=-\partial_r p+Q_r.
]

Thus pressureless baryons have (\dot v^r=0); dark matter has
(\dot v^r=\beta_{DM}\partial_r\phi), using the pinned exchange covector
(Q_{DM,r}=\beta_{DM}\rho_{DM}\partial_r\phi); and radiation has
(\dot v^r=-\partial_r(\rho/3)/(4\rho/3)).

The family is the instantaneous first-order direction
(U(s)=U_0+s\dot U_{matter}), with symmetric plus/minus states for the
locked epsilons. Conservative fluid arrays are reconstructed from these
primitives using the pinned flat-areal Valencia conventions, including
\(\sqrt{\gamma}=r^2\) and zero radiation rest-current.

## Deliberate geometry boundary

**This is not a full Einstein-matter tangent.** Geometry, lapse, shift, and
connection are held exactly flat/fixed. The direction is equation-consistent for the pinned scalar spatial
operators and the fluid continuity/Valencia momentum equations on that
background, but it does not supply the metric/extrinsic-curvature tangent
required for a coupled Einstein-matter solution. Consequently, any residual from the full Bianchi
predictor must be interpreted as a frozen-geometry matter-direction control,
not as a physical closure verdict or evidence of a missing term.

The channel runner injects the pinned vendor grid's fourth-order first- and
second-derivative operators into the analytic-rate constructor. Thus the
scalar accelerations and the DM/radiation radial gradients use the same pinned
spatial operators as the scalar/Valencia equations, without calling either
production RHS. The constructor's no-operator default remains an exact
continuum-Gaussian reference used by its isolated unit tests. The archived
control records the pinned-operator rates as explicit analytic-rate arrays.

## Guardrails

- Do not edit or regenerate the existing arbitrary-family archive.
- Do not call production initializer, scalar/matter RHS, kernel step, vendor
  step, or any time integrator.
- Do not fit coefficients, tune the predictor, or select terms by residual.
- Preserve the diagnostic source, helper baseline, vendor, tangent-helper,
  total-stress-helper, and reviewed-code blob pins.
- Use a separate class label, archive, report, and fail-closed loader path.
- Keep all scientific admission, acceptance-oracle, physical-exchange, and
  conservation-closure flags false.
- The next gates are constructor tests, pinned six-sector projections, both
  independent divergence channels, analytic-family tangent, archive/loader
  validation, then same-predictor evaluation. Compare with the arbitrary
  family and accepted-state class only after these gates pass.
