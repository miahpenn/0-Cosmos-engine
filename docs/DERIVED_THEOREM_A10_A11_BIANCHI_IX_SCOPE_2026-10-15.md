# GEAR A10–A11 — Bianchi IX identity, wall encounter, and quiescence scope

**Status:** Research-only equation and inference audit; this note does not amend a locked GEAR or REV16.  
**Date:** 15 October 2026  
**Owner provenance:** DERIVED AND DEDUCED REV4 attributes the claims to the author's independent GEAR derivations. This audit preserves that provenance separately from mathematical scope and literature priority.  
**Source anchor:** COSMOS Engine `main` at `1411c060b9b9b61c22e4840647dd9026a253b616`.  
**Primary archive evidence:** GEAR 41 (numerical Bianchi IX controls), GEAR 42 (Kasner-exponent identity), GEAR 43 (Heron/S identity and its own “honest weight”), GEAR 47 (matter-trace correction and incomplete \(\gamma_*>2\) scan), REV16 (current cyclic-criterion scope), and GEAR 55 (flatness removes curvature walls in the current model).

## 1. Executive finding

The archive supports three different claims, which must **not** be collapsed into one:

1. **Exact geometric identity:** for the stated diagonal Bianchi IX equations and time gauge, the derivative of \(S=\alpha_1+\alpha_2+\alpha_3\) is given by a Heron-discriminant term plus a matter-trace correction.
2. **Wall-encounter statement:** the sign of the Heron term fixes the sign of \(dS/d\tau\), but a dynamically significant wall encounter also requires a non-negligible magnitude. \(T<0\) alone is not proof that a wall is currently being hit.
3. **Quiescence/finite-bounce statement:** a massless scalar can make the Kasner-exponent quadratic sum decrease when \(S\) falls, but the sign of this change does not bound the decrement per encounter or by itself prove that only finitely many encounters occur.

The first is an algebraic identity under its hypotheses. The second needs a size/encounter diagnostic. The third retains BKL/model assumptions; GEAR 47 does **not** establish a universal attractor when the matter trace moves above the stiff boundary.

## 2. The exact Bianchi IX identity

Use the notation in GEAR 43/47:

- Positive scale factors \(a,b,c\), with \(A=a^2\), \(B=b^2\), \(C=c^2\).
- Logarithmic variables \(\alpha_1=\ln a\), \(\alpha_2=\ln b\), \(\alpha_3=\ln c\).
- \(S=\alpha_1+\alpha_2+\alpha_3\), volume \(V=abc\), and time gauge \(dt=V\,d\tau\).
- Directional expansion rates \(H_i\), volume expansion \(\theta=H_1+H_2+H_3\), density \(\rho\), and directional pressures \(P_i\).

For the archive's diagonal Bianchi equations,

\[
\dot H_i+\theta H_i
=\tfrac12\left(\rho-\sum_j P_j+2P_i\right),
\]

GEAR 47's conversion to \(\tau\)-time gives the summed equation

\[
\boxed{2\frac{dS}{d\tau}
=-T+V^2\left(3\rho-\sum_i P_i\right),}
\]

where

\[
T\equiv 2AB+2BC+2CA-A^2-B^2-C^2.
\]

This expression is a Heron discriminant: with \(A=a^2\), etc.,

\[
T=2a^2b^2+2b^2c^2+2c^2a^2-a^4-b^4-c^4
=(a+b+c)(-a+b+c)(a-b+c)(a+b-c).
\]

So \(T>0\) precisely when the positive lengths satisfy the strict triangle inequalities; \(T<0\) when one exceeds the sum of the other two. The pure geometric part of the equation is therefore sign-exact:

- \(T>0\): \(dS/d\tau<0\) in vacuum or when the matter-trace term vanishes.
- \(T<0\): \(dS/d\tau>0\) in vacuum or when the matter-trace term vanishes.

In the archive's singularity direction, \(\tau\) decreases. Consequently \(dS/d\tau>0\) means that \(S\) decreases along the trajectory toward that singularity.

### Scope of the matter term

The Heron-only form \(2S'=-T\) is recovered when

\[
3\rho-\sum_i P_i=0.
\]

For isotropic stiff matter \(P_i=P=\rho\), this is automatic. It also holds for an anisotropic source whose **pressure trace** obeys \(\sum_i P_i=3\rho\); the directional pressures need not individually be equal.

For GEAR 47's parameterization \(P_i=(\gamma_i-1)\mu\), with \(\rho=\mu\) and \(\gamma_*=(\gamma_1+\gamma_2+\gamma_3)/3\),

\[
3\rho-\sum_iP_i=3\mu(2-\gamma_*).
\]

Thus:

- \(\gamma_*=2\): the correction vanishes; the Heron identity is exact even for anisotropic pressure.
- \(\gamma_*>2\): the correction is negative and opposes the sign-direction associated with the wall's reduction of \(S\) toward the singularity.
- \(\gamma_*<2\): the correction is positive and aids that reduction.

**Correction to broad wording:** “stiff on average” here must mean the equality \(\gamma_*=2\) (average \(w=1\)) for the *unmodified* Heron identity. For ultra-stiff \(\gamma_*>2\), the trace correction is not zero and must be retained.

## 3. Why the sign identity is not an exact wall-activity criterion

GEAR 43's own honest-weight section explicitly notes that \(T<0\) gives the sign of \(S'\), not the magnitude of an encounter. Deep in a Kasner epoch, all polynomial terms in \(T\) can be tiny, so \(S'\) can be effectively zero despite the triangle inequality being violated. A substantial change is associated with \(|T|\) becoming large; GEAR 43 used a top-1% cut in the sampled trajectory to identify encounters.

Therefore, the following statements have different evidentiary status:

- “\(T<0\) iff the triangle inequality fails” — exact geometric statement.
- “\(T<0\) fixes the sign of the vacuum/stiff \(S'\)” — exact algebraic consequence.
- “A curvature wall is dynamically active iff \(T<0\)” — **too strong** without a magnitude/dynamical definition.
- “Every encounter reduces \(S\) by enough to guarantee a finite number of encounters” — **not proved** by the identity.

If an operational wall detector is built, its magnitude threshold must be explicitly named as a diagnostic convention and studied for sensitivity. It must not be reported as if it were derived by Heron's formula.

## 4. Quiescence: what the sign argument proves, and what it does not

For the massless-scalar Bianchi IX regime used in GEAR 42–43, the archive records the Kasner-exponent identity

\[
K_2\equiv\sum_i p_i^2=1-\frac{\phi'^2}{S^2}.
\]

Here \(p_i\) denotes the directional Kasner exponents, not pressure. If \(\phi'\ne0\) is conserved in the stated reduced system, a decrease in positive \(S\) increases \(\phi'^2/S^2\), so \(K_2\) decreases. The archive's positivity boundary is \(K_2=1/2\); with \(\sum_i p_i=1\), \(K_2<1/2\) excludes a negative Kasner exponent. The isotropic floor \(S=|\phi'|\sqrt{3/2}\) corresponds to \(K_2=1/3\).

The limiting control is essential:

- With \(\phi'=0\), \(K_2=1\) identically. \(S\) can fall substantially while this quadratic sum does not move, so the scalar-mediated quiescence mechanism is absent in vacuum.
- With \(\phi'\ne0\), each sufficiently substantial encounter that lowers \(S\) lowers \(K_2\), under the stated reduced equations.

This establishes the **direction** of the scalar mechanism. It does not give a positive lower bound on the decrement in \(S\) or \(K_2\) per encounter, so it does not alone prove a finite number of encounters. The finite-bounce result remains dependent on the relevant BKL/Andersson–Rendall hypotheses; it is not proven for every matter model or every cyclic trajectory by this identity.

### How to read GEAR 41 and GEAR 47 numerical evidence

GEAR 41 reports that its vacuum control reproduced multiple bounces while scalar cases in its tested family became quiescent after zero or one encounter. These are useful model-specific checks and controls, not a universal proof.

GEAR 47 derives the trace correction and reports:
- the \(\gamma_*=2\) anisotropic-pressure case remains quiescent in the tested run, consistent with the vanishing trace correction;
- one valid \(\gamma_*=3\), \((1,3,5)\) row leaves the quiescent region (minimum exponent \(-0.1904\), \(\sum p_i^2=0.7904\));
- three other tested \(\gamma_*>2\) rows lost constraint validity and are not meaningful outcomes.

GEAR 47 itself therefore labels the above-boundary failure **indicated, not established**. The single valid row supports the predicted direction but does not establish a sharp general boundary, a universal attractor failure, or the Barrow–Ganguly Type IX result (whose model includes a negative-energy ghost bounce that the archive's matter-lane argument separately challenges).

## 5. Verdicts for A10 and A11

### A10 — Bianchi IX / Heron identity

**Verdict: exact conditional identity; inference overstated in REV4's prose.**

The identity and the sign equivalence follow from the displayed equations. The matter-trace term must be included except when \(\sum_i P_i=3\rho\). The statement “the wall is active precisely when the geometry cannot close” confuses a sign condition with encounter magnitude, and the identity gives no bound on the size of the \(S\) change.

### A11 — quiescence as an attractor

**Verdict: conditional monotonic mechanism; finite-bounce/global-attractor claim not proven by this derivation alone.**

With a nonzero conserved massless-scalar term, falling \(S\) drives \(\sum_i p_i^2\) down in the reduced Bianchi IX system. Vacuum is a negative control in which that mechanism vanishes. The finite encounter count requires an independent theorem or further proof bounding the accumulated decrease and is not obtained merely from a derivative sign.

The current COSMOS model is spatially flat and has no Bianchi IX curvature walls (REV16/GEAR 55); this mechanism must not be silently inserted as an attractor or production rule in the existing flat engine.

## 6. Suggested research-machine use

An **offline, separate Bianchi IX diagnostic** could report:

1. \(T\), the factorized Heron form, and the signed residual
   \[
   R_S=2S'+T-V^2\left(3\rho-\sum_iP_i\right);
   \]
2. \(|T|\) separately from its sign, with any wall-encounter cutoff labelled as an empirical diagnostic choice;
3. the matter-trace variable \(3\rho-\sum_iP_i\) and \(\gamma_*\) where defined;
4. \(K_2=\sum_i p_i^2\), \(\phi'\), \(S\), and constraint residuals, separating monotonic direction from bounce count;
5. terminal events and constraint-validity failures, including failed runs rather than silently dropping them.

Before any new expensive integration, verify the exact source equations and existing artifacts; predeclare the range of wall thresholds and admissible constraint residuals. The existing GEAR 47 run should not be repeated unchanged: its unresolved question is threshold/parameter robustness of the **single valid above-boundary result**, not whether the recorded failed rows can produce numbers.

## Provenance and review checklist

- [x] Heron factorization and trace correction written explicitly.
- [x] Exact sign relation separated from wall magnitude/encounter.
- [x] Quiescence sign mechanism separated from finite-bounce/global-attractor claim.
- [x] GEAR 47's one-valid-row/three-invalid-row limitation retained.
- [x] Current flat COSMOS model kept distinct from a Bianchi IX study.
- [ ] Independent review of time gauge, matter normalization, and Kasner identity against the original equations.
- [ ] If a new scan is justified after review, use predeclared stop thresholds and constraint gates.

**Admission status:** NOT ADMITTED. This document is research-only; it changes no production equation, parameter, GEAR archive, or scientific-admission decision.
