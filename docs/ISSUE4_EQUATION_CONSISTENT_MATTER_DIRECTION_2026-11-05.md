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


## Bounded execution and three-class comparison

The equation-consistent constructor and full observer completed under the
locked pins in workflow run
[42082807688](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42082807688).
The full three-class comparison and complete 28-test suite passed in workflow run
[42084531386](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42084531386).
The latter run used the exact accepted-state N=60 archive whose SHA-256 is
`41d646090e05c08684806c5e39d57851f119cf0f38aa90f1d26ab5f72f87d1ce`,
and the immutable arbitrary-family archive SHA-256
`ac4a12d5c45c37b505c4f20db83a56b399fb4232cfc0ca7102371cf1717355b7`.

The same current full-Bianchi evaluator was run on all three classes:

| Class | Pointwise max residual, r<=10 | FV max residual, r<=10 |
|---|---:|---:|
| Arbitrary manufactured family | 1.148320316987254e-2 | 1.143576348708114e-2 |
| Equation-consistent matter direction, frozen geometry | 1.148320316987254e-2 | 1.147475977637545e-2 |
| Accepted t=0.60, full metric | 1.840585461201662e-5 | 9.640437042539719e-6 |

For the two manufactured controls, the measured tangent changes by a maximum
`2.122036622890433e-2`. The pointwise predictor changes by the same amount,
leaving the pointwise residual arrays equal to within
`1.65e-17`. This confirms that the pointwise residual under frozen geometry
does not discriminate between these matter directions: the time-derivative
piece moves in both the measured constraint tangent and the stress-divergence
source. The finite-volume residual differs by up to
`1.452304462221208e-3`, reflecting its distinct spatial operator.

The accepted full-metric state is about 624 times smaller in pointwise peak
residual and about 1,189 times smaller in finite-volume peak residual than
the equation-consistent frozen-geometry control. This is a comparison of
different geometric contracts, not proof that the matter control has found a
missing term. It demonstrates why a full Einstein-matter tangent is needed
before interpreting a Bianchi propagation residual physically.

The accepted-state reference regression passed its archived signed-witness
orientation, pointwise/FV residual, and epsilon-sweep checks. All three
evaluations ran zero production steps, did not advance state, and left the
zero-step tripwire armed and untriggered. Scientific admission,
acceptance-oracle, physical-exchange, and conservation-closure flags remain
false.

The pointwise time-density component-accounting flag remains false, with an
unexplained reconstruction remainder at approximately `8.4e-18`; the raw
six-sector sum reconstructs the total. The finite-volume corrected accounting
passes. Neither remainder is remotely large enough to explain the
approximately `1.1e-2` frozen-geometry residual.

The final full three-class artifact is GitHub Actions artifact
`13442329240` (outer ZIP SHA-256
`acfc123ff1748c8f7adcbf19e1b94f6fafd24771fd0787daa5542b17ee9f8e1b`).
It contains the equation-consistent fixture, the accepted fixture and report,
the immutable arbitrary-family fixture and report, and
`three-class-comparison.json`.
