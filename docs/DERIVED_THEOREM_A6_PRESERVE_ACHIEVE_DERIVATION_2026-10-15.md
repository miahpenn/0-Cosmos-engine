# GEAR A6 — Preserve versus achieve: conditional evolution identity

**Status:** Research derivation note; not an amendment to a locked GEAR or to REV16.  
**Date:** 15 October 2026  
**Owner provenance:** The archive lists “preserve versus achieve” as an independently derived owner result (DERIVED AND DEDUCED REV4, A6). This note reconstructs the identity from the controlling equations; it does not make a priority claim.  
**Source anchor:** COSMOS Engine `main` at `1411c060b9b9b61c22e4840647dd9026a253b616`.  
**Relevant archive notes:** GEAR 35 (free shear scaling), GEAR 39 (definition and 1/4 criterion), GEAR 45 (anisotropic-stress source), GEAR 58 (radiation/turnaround caveat), and *The 0* Cycle — Authoritative, REV16* (asymptotic scope and (H=0) exclusion).

## 1. Result in one sentence

For a homogeneous Bianchi-I model with **no anisotropic-stress source**, a separately conserved isotropic non-shear sector of density \(\rho>0\) and effective equation of state \(w=p/\rho\), and away from \(H=0\),

\[
\boxed{\frac{d\Omega_{\sigma}}{d\ln a}
=3(w-1)\,\Omega_{\sigma}(1-\Omega_{\sigma}) .}
\]

This is an exact identity under those assumptions, not a universal law for all GEAR channels. In contraction, \(H<0\) and \(a\) decreases with increasing proper time. Thus a **positive** derivative with respect to \(\ln a\) means \(\Omega_\sigma\) falls as contraction proceeds.

## 2. Definitions and assumptions

Use the archive convention

\[
3H^2=\rho+\sigma^2,\qquad
\sigma^2\equiv\tfrac12\sum_i\sigma_i^2,\qquad
\Omega_\sigma\equiv\frac{\sigma^2}{3H^2}.
\]

This agrees with GEAR 39's \(\sum_i\sigma_i^2=6\Omega_\sigma H^2\). Define the unnormalised shear-to-fluid ratio

\[
r\equiv\frac{\sigma^2}{\rho}.
\]

The derivation assumes:
1. Bianchi-I homogeneous expansion/contraction with the above sign and normalization conventions.
2. Isotropic pressure, so the shear equation is unsourced: \(\dot\sigma_{ab}+3H\sigma_{ab}=0\).
3. The non-shear density is separately conserved: \(\dot\rho+3H(\rho+p)=0\). If several isotropic components are present, \(w\) can be their effective \(p/\rho\), provided their sum has no net exchange with another sector.
4. \(\rho>0\), \(\sigma^2>0\), and \(H\ne0\) on the interval where logarithmic derivatives and \(\Omega_\sigma\) are used.

Do not silently substitute a solver variable called `w_total`: confirm that it denotes the equation of state of the non-shear density \(\rho\), not a definition that already includes the effective shear contribution.

## 3. Derivation

From the unsourced shear equation,

\[
\dot\sigma_{ab}=-3H\sigma_{ab}
\quad\Longrightarrow\quad
\dot{\sigma^2}=-6H\sigma^2,
\]

and therefore, for \(H\ne0\),

\[
\frac{d\ln\sigma^2}{d\ln a}
=\frac{\dot{\ln\sigma^2}}{\dot{\ln a}}=-6.
\]

Fluid conservation gives

\[
\dot\rho=-3H(1+w)\rho
\quad\Longrightarrow\quad
\frac{d\ln\rho}{d\ln a}=-3(1+w).
\]

Subtracting the logarithmic slopes yields

\[
\frac{d\ln r}{d\ln a}
=\frac{d\ln\sigma^2}{d\ln a}
-\frac{d\ln\rho}{d\ln a}
=3(w-1),
\]

or \(dr/d\ln a=3(w-1)r\).

Now use the Friedmann constraint to write the Hubble-normalised fraction in terms of \(r\):

\[
\Omega_\sigma
=\frac{\sigma^2}{\rho+\sigma^2}
=\frac{r}{1+r}.
\]

Thus \(d\Omega_\sigma/dr=(1+r)^{-2}\), so

\[
\frac{d\Omega_\sigma}{d\ln a}
=\frac{1}{(1+r)^2}\,3(w-1)r
=3(w-1)\Omega_\sigma(1-\Omega_\sigma).
\]

The same result holds for time-dependent effective \(w(a)\) under the same conservation assumptions, with its instantaneous value \(w(a)\) on the right-hand side.

## 4. Meaning on a contracting branch

Define a contraction-increasing e-fold variable \(N_c\equiv\ln(a_{\rm ref}/a)\), so \(dN_c/dt=-H>0\) while \(H<0\). Then

\[
\boxed{\frac{d\Omega_\sigma}{dN_c}
=-3(w-1)\Omega_\sigma(1-\Omega_\sigma).}
\]

| Non-shear equation of state | Slope in \(\ln a\) | Evolution as contraction proceeds | Interpretation |
|---|---:|---|---|
| \(w>1\) | positive | \(\Omega_\sigma\) decreases | **Achieve**: dynamics drives the shear fraction down |
| \(w=1\) | zero | \(\Omega_\sigma\) is constant | **Preserve**: the entering shear fraction is retained |
| \(w<1\) | negative | \(\Omega_\sigma\) increases | shear fraction grows; the quiescence margin can be lost |

This resolves the sign ambiguity: the sign of \(d\Omega_\sigma/d\ln a\) is a useful label only if the contracting-branch orientation is stated. It is not the same sign as the derivative along increasing physical time during contraction.

## 5. Limiting and consistency checks

- **Marginal case:** if \(w=1\), then \(r\) and \(\Omega_\sigma\) are exactly constant throughout any interval satisfying the assumptions.
- **Sub-stiff/super-stiff direction:** for \(0<\Omega_\sigma<1\), the prefactor \(\Omega_\sigma(1-\Omega_\sigma)\) is positive, so the sign is exactly the sign of \(w-1\) in \(d/d\ln a\).
- **Small-shear limit:** as \(\Omega_\sigma\to0^+\), the slope tends to \(3(w-1)\Omega_\sigma\); a zero-shear solution remains zero in the unsourced model.
- **Normalized-shear endpoint:** the formula's slope tends to zero as \(\Omega_\sigma\to1^-\), but the fluid-density fraction also tends to zero; the value \(w=p/\rho\) may cease to be defined at the endpoint. Do not infer an independent fluid evolution there.
- **Radiation:** a separately conserved radiation component has \(w=1/3\) and ratio \(\rho_r/\rho_{\rm kin}\propto a^2\) against kination, so it becomes less important as \(a\) falls. A sourced radiation component requires the source term in the generalized equation below; it is not covered by the passive-fluid case.

## 6. Where the identity stops applying

### The turnaround or any \(H=0\) crossing

\(\Omega_\sigma=\sigma^2/(3H^2)\) is Hubble-normalized. For nonzero shear it becomes singular as \(H\to0\), and \(d/d\ln a=(1/H)d/dt\) is undefined at the crossing. REV16 says the \(\Omega_\sigma<1/4\) condition is asymptotic; GEAR 58 also records the divergence of normalized fractions at turnaround. Therefore this identity must **not** be used as a pointwise test through \(H=0\), and it does not invalidate the archive's finite turnaround-exposure warning. Use a regular proper-time diagnostic and restart the logarithmic-\(a\) analysis only on a branch where \(H\ne0\) and its direction is known.

### Anisotropic stress or energy exchange

If a source changes the shear magnitude, or the non-shear sector exchanges energy, define

\[
\dot{\sigma^2}=-6H\sigma^2+\mathcal S_\sigma,\qquad
\dot\rho=-3H(1+w)\rho+Q.
\]

Where the logarithms are defined,

\[
\frac{d\Omega_\sigma}{d\ln a}
=\Omega_\sigma(1-\Omega_\sigma)
\left[3(w-1)+\frac{\mathcal S_\sigma}{H\sigma^2}
-\frac{Q}{H\rho}\right].
\]

The quantities \(\mathcal S_\sigma\) and \(Q\) must be derived from the selected model and conventions; they are not free fitted terms. In particular, GEAR 45 records that \(\dot\sigma_{ab}+3H\sigma_{ab}=\mu P_{ab}\) when anisotropic stress is present. Then the free \(a^{-6}\) shear law and this unsourced preserve/achieve conclusion do not follow without evaluating that source. Similarly, if other interactions change \(\rho\), the \(Q\) term can change the sign.

## 7. Verdict for claim 6

**Verdict: conditionally proven after sharpening the statement.**

A defensible version of A6 is:

> In the asymptotic, positive-density, homogeneous Bianchi-I contraction with isotropic pressure and no shear/energy-exchange source, the sign of \(d\Omega_\sigma/d\ln a=3(w-1)\Omega_\sigma(1-\Omega_\sigma)\), read with \(\ln a\) decreasing in physical time, distinguishes preservation at \(w=1\) from dynamical suppression for \(w>1\) and growth for \(w<1\). The result is not valid as a normalized diagnostic at \(H=0\), and does not cover the anisotropic-stress channel.

This does **not** prove finite BKL bounce count, global cyclic viability, or production readiness. The finite-bounce/quiescence statement has its own literature/model hypotheses and is not a corollary of the sign identity alone.

## 8. Suggested machine use (research-only)

Implement an offline diagnostic that reports, separately:
- \(\Omega_\sigma\), \(r=\sigma^2/\rho\), \(w_{\rm eff}=p/\rho\), \(H\), and the branch sign;
- the analytic slope \(3(w_{\rm eff}-1)\Omega_\sigma(1-\Omega_\sigma)\);
- a finite-difference slope computed only where \(|H|\) is safely nonzero, with resolution/tolerance declared;
- the proper-time trend \(\dot\Omega_\sigma=H\,d\Omega_\sigma/d\ln a\);
- a scope/status flag excluding turnaround samples and identifying any nonzero anisotropic-stress or exchange source.

Use exact limiting controls \(w=1\), \(w>1\), \(w<1\), zero shear, and passive radiation. Do not infer \(\Omega_\sigma\) from an initialization fraction such as GEAR 35's \(f_i\); record the realized value, as the GEAR 38 correction requires. No production change is authorized by this note.

## Provenance and review checklist

- [x] Reconstructed from REV16/GEAR 39 definitions and GEAR 35/45/58 source-scope notes.
- [x] Explicit branch direction, units/conventions and excluded \(H=0\) layer.
- [x] Conditional source/energy-exchange generalization written.
- [x] Limiting cases stated algebraically.
- [ ] Independent reviewer to verify every normalization against the original GEAR equations and confirm the proposed offline diagnostics' variable mapping.
- [ ] If implemented, exact-head checks and source-pinned tests; no long campaign until a specific research question requires one.

**Admission status:** NOT ADMITTED. This file is research documentation only; it changes no production equation, model parameter, GEAR archive, or scientific-admission decision.
