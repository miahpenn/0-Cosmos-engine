# GEAR 69 / 70B independent verification — 15 October 2026

**Review disposition: PARTIAL PASS — perturbative structure verified; stated GEAR-69 coefficient target not reproduced from the displayed equations.**

This report verifies the equations and numerical claims in the preserved GEAR-69 and GEAR-70B materials against an independent implementation. It is a research-only calculation. It does not change production physics or verify a global/cosmological arrest mechanism.

## Provenance

- Repository `miahpenn/0-Cosmos-engine`; initial audit preflight `main`: `f754492d6f652c4b00c6511c65a5475d5f3382ee`. The final PR branch was refreshed onto current `main` `3a64d0b1f989595e4669740ed2ba4be85c51bf78`. The intervening main commits add the already-claimed momentum-tangent diagnostic/workflow only; they do not change this standalone weak-field equation check.
- Preserved derivation: `GEAR_69_First_Order_Einstein_Scalar_Backreaction.md` (Library copy).
- Preserved numerical report: `GEAR-70B_radial_feedback_result_2026-10-02.md` (Library copy).
- Preserved calculation script: `gear70B_direct_radial_einstein_scalar_2026-10-02.py`, SHA-256 `49bf809e0813d68ce874ee2c3510c9eff10c56ca27ba0fb7b379635f7a807704`.
- New standalone verifier: [`research/verify_gear69_70b.py`](../research/verify_gear69_70b.py). It imports no production module and derives the background/source directly from the equations printed below.
- Equation-convention cross-check: Liebling & Palenzuela, “Dynamical boson stars,” *Living Reviews in Relativity* 20, 5 (2017; updated version available), [DOI: 10.1007/s41114-017-0007-y](https://doi.org/10.1007/s41114-017-0007-y), §2.3–2.4, especially polar-areal equations (35)–(40). Their scalar potential/model is not the GEAR inverted potential; this source is used to cross-check the standard metric constraint form only.

## 1. Derivation audit

Assumptions being verified—not inferred beyond their stated scope:

- Polar-areal metric `ds^2=-alpha^2 dt^2 + a^2 dr^2 + r^2 dOmega^2`.
- Canonical real scalar with inverted quadratic potential `V(phi)=-mu^2 phi^2/2`.
- Dirichlet cavity at `r=R`, regular centre, and declared exterior-clock condition `alpha(R)=1/a(R)`.
- Einstein constraints and field equations as written in GEAR 69; no extra matter, phenomenological feedback, or cosmological coupling.

The leading mode `u=sin(k r)/(k r)`, `k=pi/R`, satisfies
`u'' + (2/r)u' + k^2 u = 0`. Therefore the flat growing rate is
`gamma^2=mu^2-k^2`. With `R=10, mu=G=1`, the independently evaluated rate is

```
gamma = 0.9493702944526474
```

The second-order scalar energy density is
`rho2=(gamma^2*u^2 + u'^2 - mu^2*u^2)/2`; the mass constraint is
`m2'=4*pi*G*r^2*rho2`, with `a2=m2/r`. From the radial-pressure Einstein equation the lapse constraint expands to

```
alpha2' = a2/r + 2*pi*G*r*(gamma^2*u^2 + u'^2 + mu^2*u^2)
```

and the boundary condition gives `alpha2(R)=-a2(R)`. Define `q=alpha2-a2` and `s=alpha2+a2`. Expanding
`phi=A*u*exp(gamma*t) + A^3*w*exp(3*gamma*t) + O(A^5)`
in the stated first-order momentum system, including `phi_tt = D_t*Pi + D*Pi_t` for `D=alpha/a`, gives

```
w'' + 2*w'/r + (mu^2 - 9*gamma^2)*w = S3(r)
S3 = -(1/r^2)*d/dr(r^2*q*u') - mu^2*s*u - 3*gamma^2*q*u
w'(0) = 0,  w(R) = 0
```

The radial derivative expands as
`(1/r^2)*d/dr(r^2*q*u') = 2*q*u'/r + q'*u' + q*u''`.
Thus the correction to include `q*u''` in 70B is mathematically correct. The cubic mode correction determines the central coordinate-time growth rate:
`Gamma_meas = gamma*(1 + 2*w(0)*phi_c^2 + O(phi_c^4))`.
With `X=G*R^2*phi_c^2`, this yields the stated coefficient identity
`1-Gamma_meas/gamma = C*X + O(X^2)`, `C=-2*w(0)/(G*R^2)`.
**These algebraic steps pass, under the assumptions above.**

## 2. Replay of preserved GEAR-70B code

Running the preserved script unchanged reproduces its recorded output:

| Quantity | Preserved output |
|---|---:|
| `gamma` | 0.949370294453 |
| `w(0)` | -72.160002567356 |
| `C` | 1.443200051347 |

However, instrumenting the returned `solve_bvp` object shows that its solve is **not marked converged** at the default parameters: `success=False`, status 1, message “The maximum number of mesh nodes is exceeded,” 31,655 nodes, and maximum reported RMS residual about `2.00e-5`. The script reads `sol.sol(1e-5)` and prints a value without checking `sol.success`. Grid refinements of that original implementation continue to print values, but the solver failure flag remains false-success; therefore those prints alone are not a convergence certificate.

Original-script source-grid sweep, retaining its piecewise-linear `np.interp` source and original BVP routine:

| Background points | Printed `C` | `solve_bvp.success` |
|---:|---:|---|
| 1,000 | 1.443198910425 | false |
| 2,000 | 1.443200051347 | false |
| 4,000 | 1.443200345680 | false |
| 8,000 | 1.443200426468 | false |
| 16,000 | 1.443200443289 | false |

This sequence trends to approximately 1.44320045 but is not, by itself, a converged solve because the BVP integrator does not report success.

## 3. Independent, converged re-solve

The new `research/verify_gear69_70b.py` rebuilds the second-order mass/lapse profile using cumulative Simpson quadrature, uses a smooth cubic-spline source, enforces the regular centre in the BVP, and requires `solve_bvp.success=True`. It then cross-checks `w(0)` using a separate regular-center shooting construction. No values are fitted to the GEAR target.

Default run: `python research/verify_gear69_70b.py`.

| Diagnostic | Independent result |
|---|---:|
| `gamma` | 0.9493702944526474 |
| `m2(R)` | `4.55e-16` |
| `alpha2(R)` | `-4.55e-17` |
| `w(0)`, BVP | -72.160022487208 |
| `C`, BVP | 1.443200449744 |
| `w(0)`, separate shooting | -72.160022487091 |
| `C`, separate shooting | 1.443200449742 |
| Absolute BVP–shooting difference in `C` | `2.33e-12` |
| BVP status | converged; 1,015 adaptive nodes |

The independently computed `C` remained `1.443200449744` when the background grid was increased from 5,001 to 100,001 points (the first 12 decimals shown did not change). This is stable numerical evidence for the BVP defined by the equations above.

## 4. Material inconsistency requiring correction or source reconciliation

GEAR 69 reports `C_target=1.443222`. That target implies
`w(0)=-C*G*R^2/2=-72.1611` at `G=1,R=10`. The independently solved BVP from the displayed equations instead gives `w(0)=-72.160022487208` and `C=1.443200449744`.

- Independent `C` minus reported GEAR-69 target: `-2.1550256e-5`.
- Relative discrepancy against the target: `-0.0014932%`.
- Independent `C` minus the reported GEAR-70B number: `+3.9840e-7`.

So the original GEAR-70B output is very close to the independently converged answer, but **the displayed equations do not reproduce the stated GEAR-69 target of 1.443222**. The small discrepancy is reproducible, not a matter of rounding the printed digits. The prior direct nonlinear results (about 1.44135 at their largest listed resolution) and GEAR-70C.1 time-domain values (about 1.4458–1.4475) are sub-percent checks, but they do not resolve this (2.16\times10^{-5}) coefficient-level discrepancy.

### Required next gate

1. Review the exact original calculation/source that produced the `1.443222` target and reconcile it term by term against `rho2`, the lapse constraint, `S3`, and the boundary conditions.
2. Either correct the target to the independently verified solution of the displayed equations or identify and document the additional convention/term in the original derivation that changes the coefficient. Do not tune the source to match the target.
3. Update 70B to check `sol.success` explicitly and preserve the BVP residual/status in its evidence report; use a smooth source representation or another solver whose convergence is demonstrated.
4. Obtain independent review before marking GEAR 69/70B fully verified.

## Scope boundary

**Verified:** the leading mode/rate, the second-to-third-order feedback structure, the `q*u''` contribution, the central growth-rate expansion, and a stable numerical solution of the displayed radial BVP.

**Not verified:** the specific GEAR-69 target value `1.443222`; a unique original derivation yielding that value; strong-field arrest, matter formation, local-to-cosmological coupling, or a global cyclic solution.

No production equation, physical coefficient, default, threshold, or campaign configuration was changed. The density/source normalization and local S/D-to-COSMOS provenance issues remain separate.
