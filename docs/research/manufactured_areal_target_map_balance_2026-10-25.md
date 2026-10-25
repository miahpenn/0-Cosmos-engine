# Exact manufactured mass/current balance on the target areal maps — 2026-10-25

**Status:** research-only calibration of fixed-\(R\) sampling. The manufactured continuum identity is exact; the observer measures only the interpolation/query remainder on the actual target map trajectory. It does not assert or test the physical matter source law.

**PR:** [#420](https://github.com/miahpenn/0-Cosmos-engine/pull/420)  
**All-resolution workflow:** [run 40175812454](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40175812454), all three matrix jobs passed. Artifacts: [N=60](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40175812454/artifacts/12591729334), [N=120](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40175812454/artifacts/12592370465), [N=240](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40175812454/artifacts/12591584633).  
**Accepted source:** SHA-verified full continuous-source profile artifact from [run 40156226086](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40156226086).  
**Pinned commits:** source `4aec88d0e3ff028e0a65a8ac36b30948625df792`; vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. All jobs passed the protected-source/vendor guard; accepted and copied shadow-state digests stayed unchanged.

## Method

For each N=60/120/240, the observer loads the immutable accepted state at \(t=0.60\), advances **copies** through two pinned-kernel shadow steps at \(\Delta t=0.03\Delta r\) and \(\Delta t/2\), and captures the actual areal maps \(R(t,r)=r\sqrt{b}/X\). It does not use the evolved target stress tensor or target mass in the manufactured identity.

On each captured map it evaluates:
\[
m_*(t,R)=(M_0+\dot M t)f(R),\quad
\chi_*(t,R)=1-\frac{2m_*}{R},\quad
F_*(t,R)=\dot M f(R),
\]
where \(M_0=0.30\), \(\dot M=0.12\). Therefore \(-R\partial_t\chi_*|_R/2=F_*\) exactly. Three smooth, predeclared profiles are tested:
- rational: \(f=R^3/(7^3+R^3)\)
- Gaussian: \(f=e^{-(R/6)^2}\)
- broad rational: \(f=R^2/(12^2+R^2)\)

The balance is queried at common fixed radii over \(2\le R\le10\), using linear and four-point cubic interpolation. This tests the actual target grid motion without invoking a physical source relation.

## Actual areal-map motion in the annulus

The domain-wide map-speed RMS is about 0.262, but it is dominated by large-radius motion. In the measured annulus, the map moves much less:

| N | Annulus cells | RMS \(|\partial_t R|\) | Max \(|\partial_t R|\) | RMS one-step shift / cell | Max shift / cell |
|---:|---:|---:|---:|---:|---:|
| 60 | 6 | 0.03255 | 0.04869 | 0.000976 | 0.001461 |
| 120 | 12 | 0.03610 | 0.05442 | 0.001083 | 0.001633 |
| 240 | 24 | 0.03611 | 0.05538 | 0.001083 | 0.001661 |

These cell-fraction shifts are for \(\Delta t\); at \(\Delta t/2\), each RMS shift is halved (to approximately 0.000488, 0.000542, and 0.000542). Thus the fixed-radius query stencil moves by only about one-tenth of one percent of a cell per full step in this annulus. A generic moving-map control with faster local motion should not be assumed to represent this target annulus.

## Exact-identity balance residuals

Centered mass-rate minus midpoint-current RMS over \(2\le R\le10\), using cubic interpolation and \(\Delta t/2\):

| Analytic profile | N=60 | N=120 | N=240 | Observed spatial orders |
|---|---:|---:|---:|---:|
| Rational | \(4.3845\times10^{-5}\) | \(2.9584\times10^{-6}\) | \(2.0827\times10^{-7}\) | 3.89, 3.83 |
| Gaussian | \(3.0732\times10^{-3}\) | \(2.0867\times10^{-4}\) | \(9.6274\times10^{-6}\) | 3.88, 4.44 |
| Broad rational | \(2.5797\times10^{-6}\) | \(1.6987\times10^{-7}\) | \(1.1167\times10^{-8}\) | 3.92, 3.93 |

The Gaussian fine-interval order 4.44 is treated as a descriptive superconvergent/cancellation result, not a general order claim. With linear interpolation the measured orders are roughly 1.7–2.0, depending on profile.

At N=240, the rational cubic centered witness is \(2.0830\times10^{-7}\) for \(\Delta t\) and \(2.0827\times10^{-7}\) for \(\Delta t/2\) (about 0.015% change). The corresponding linear witness is \(3.4860\times10^{-5}\) and \(3.4861\times10^{-5}\). This confirms that halving the time step has negligible effect on these interpolation-dominated controls.

## Interpretation and disposition

1. With the actual target maps, the exact manufactured identity has an approximately fourth-order cubic-sampling remainder over N=60/120/240 for the rational and broad-rational profiles. The target annulus’s very small per-step cell shift explains why the generic, more rapidly moving map’s third-order result does not carry over.
2. The error amplitude is profile-dependent: at N=240 the cubic centered RMS is \(2.08\times10^{-7}\) for the rational profile, \(9.63\times10^{-6}\) for the Gaussian, and \(1.12\times10^{-8}\) for the broad rational. These are calibration controls, not interchangeable acceptance bounds. The target evolved-state witness near \(1.03\times10^{-7}\) cannot be declared closed or bounded from an arbitrary manufactured profile.
3. This answers the sampling question for the target map trajectory: the apparent fourth-order trend is compatible with fourth-order sampling in that annulus, and generic moving-map third-order behavior is not decisive here. It does **not** establish the physical continuum current identity, close the parity-aware face-shell FV residual, or complete the full matter-coupled momentum-propagation identity.
4. Next effort should return to the source-pinned discrete constraint/equation and integrated sector-balance crosswalk, keeping the annulus, inner shell, native center, and face-shell FV witnesses separate. Do not add a production correction chosen by residual minimization.

**No production equation, initializer, engine, vendor source, or dependency changed.** Full matter-coupled propagation, total stress-energy sector conservation, and scientific admission remain open/false.

Observer: `tools/audit_manufactured_areal_target_map_balance.py`. Workflow: `.github/workflows/manufactured-areal-target-map-balance.yml`.
