# Issue #4 analytic-s2 origin half-cell discriminator

**Status:** diagnostic-only candidate. Do not use to justify a production initializer patch; scientific admission remains **NOT ADMITTED**.

## Question

Does enforcing the analytically known regular-origin coefficients \(s_0=1\) and \(s_2=-4\pi A^2\), while fitting only \(s_4,s_6,s_8\), reduce the model-form ambiguity in the copied initializer's first half-cell integral?

For the standard uncharged initializer, the radial source is

\[
s(r)=1-8\pi r^2\rho(r)+2r^2 K_rK_t+r^2K_t^2.
\]

Given \(K_t=-H_0\), \(K_r=K_t+O(r^2)\), \(3H_0^2=\rho_{\rm archive}\), and \(8\pi\rho(0)=\rho_{\rm archive}+4\pi A^2\), the archive terms cancel and \(s_2=-4\pi A^2\). This scope does not include the optional charged-dust/radial-Maxwell initializer.

For cell-centre \(r_0=\Delta r/2\), production uses \(I_{\rm point}=r_0s(r_0)\), whereas a regular even-source integral is \(I_0=\int_0^{r_0}s(r)\,dr\). The leading correction is \(8\pi A^2r_0^3/3+O(r_0^5)\).

## Design

- Pin diagnostic branch to live main SHA 27da39e1e4d3a37ce0f7ef70b228009323183048, all eleven production-source blobs, and vendor d6052d605673ce9d82cdc99b0df79855fa2ea215.
- At \(N=160,240,320\), reconstruct the original six-pass point-sampled initializer and require it to reproduce the production \(X^6\) profile.
- Solve four independent copied-source/metric fixed points: existing unconstrained degree-3 \(r\le2.5\), degree-4 \(r\le3\), degree-5 \(r\le3.5\), and the new constrained \(s_0/s_2\) arm fitting only \(s_4,s_6,s_8\) over fixed physical radius \(r\le3.5\).
- Preserve inter-centre trapezoids. Report fit residuals and leave-one-out witnesses, half-cell moment, signed first-two-cell and fixed-annulus Ricci differences, constraints, metric identities, fixed-point closure, state digests, arrays, and hashes.
- No production step is advanced; no production/vendor/dependency path or tolerance is changed.

## Decision gate

Treat a successful workflow as validation of the copied diagnostic only. Compare the constrained arm with each unconstrained fit family. If the source fit remains sensitive, if a fixed point fails to converge, if input state changes, or if any metric identity exceeds its declared tolerance, retain Issue #4 as open and do not propose a production change based solely on this result.
