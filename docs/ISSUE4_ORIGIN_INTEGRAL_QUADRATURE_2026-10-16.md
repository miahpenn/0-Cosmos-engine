# Issue #4 first-cell integral quadrature discriminator

**Status:** bounded manufactured diagnostic; scientific admission remains **NOT ADMITTED**.

## Purpose

The analytic regular-metric test has shown that the existing proper-radial and vendor BSSN curvature paths can reproduce the closed-form centre Ricci limit. The remaining leading hypothesis is a nonuniform error in the *radial-source-to-metric integration*, not a generic Ricci tensor formula failure.

The manufactured source is

\[
s(r)=1+s_2r^2,\qquad s_2=-4\pi A^2,\quad A=0.01
\]

over the entire centre stencil. Its exact metric reconstruction is

\[
B(r)=\frac1r\int_0^r s(x)\,dx=1+\frac{s_2}{3}r^2,
\]

for which the exact trace-free Ricci scalar is zero:

\[
R_{\rm TF}=\frac23\left[\frac{B-1}{r^2}-\frac{B'}{2r}\right]=0.
\]

The source is cut off smoothly nowhere only at a remote, declared radius \(r=8\) to keep the full-grid metric positive; all reported centre stencils lie well inside that cutoff. This outer continuation is not part of the centre acceptance witness.

## Three integration arms

At \(N=160,240,320\), compare:

1. Production-style point-sampled first half-cell \(I_0=r_0s(r_0)\), followed by centre-to-centre trapezoids.
2. Analytic first half-cell \(I_0=\int_0^{r_0}s\,dr\), followed by the same centre-to-centre trapezoids.
3. A parity-aware composite Simpson candidate: use the known regular value at the origin for the first half-cell, then integrate the even extension with Simpson 1/3 panels and a Simpson 3/8 end panel.
4. Exact cumulative polynomial integration at every centre.

The signed identities provide a priori controls, rather than a tolerance fitted to an observed curvature mismatch. For the point-sampled first cell,

\[
B^{\rm point}_0-B^{\rm exact}_0=\frac23s_2r_0^2.
\]

Even after correcting the first half-cell, the first centre-to-centre trapezoid has signed integral error \(s_2\Delta r^3/6\), therefore

\[
B^{\rm half+trap}_1-B^{\rm exact}_1=\frac{s_2\Delta r^2}{9}.
\]

This second \(O(\Delta r^2)\) metric error is concentrated in the first few cell-index rows and can become a finite centre-curvature residual when a second-derivative stencil divides it by \(O(\Delta r^2)\). The diagnostic tests whether this predicted mechanism quantitatively reproduces the observed pattern, whether a practical parity-aware composite-Simpson candidate removes it, and whether full analytic integration removes it.

For each arm, report the signed errors in B at cells 0–2; vendor BSSN and proper-radial \(R_{\rm TF}\) against the analytic zero; and the signed difference between the two curvature routes. Retain source/vendor pins, metric mapping checks, state digests, CSV/NPZ arrays, and SHA-256 manifest. No production step is advanced and no engine/vendor/dependency path or tolerance is changed.

## Interpretation boundary

A match to the manufactured quadrature identities would support a specific discretization-mechanism hypothesis, not alone authorize a production patch. Any production change must be reviewed against the actual initializer quadrature contract and regression controls. The independent same-stage total-stress-energy balance and Lumen's Brown/Bianchi propagation oracle remain separate scientific-admission gates. Keep Issue #4 open and scientific admission false until those independent requirements are satisfied.
