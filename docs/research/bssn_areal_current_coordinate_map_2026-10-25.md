# BSSN isotropic-to-areal coordinate map for the integrated-current gate — 2026-10-25

**Purpose:** pin the coordinate and measure transformation required before comparing the archive's areal-radius Misner–Sharp current identity with the current BSSN finite-volume divergence observer. This is an analytic coordinate-map note; it is not a new conservation pass.

## 1. Exact metric map in the current adapter

The pinned \`engine/valencia.py:spherical_metric_from_bssn\` constructs
\[
\gamma_{rr}=\frac{a}{X^2},\qquad
\gamma_{\theta\theta}=\frac{b r^2}{X^2},\qquad
\sqrt{\gamma}=\frac{\sqrt{a}\,b r^2}{X^3}.
\]
Therefore the areal radius is
\[
R(t,r)=\sqrt{\gamma_{\theta\theta}}=\frac{r\sqrt b}{X},
\]
and
\[
\sqrt{-g}_{(t,r)}=\alpha\sqrt{\gamma}
=\frac{\alpha\sqrt a\,b r^2}{X^3}.
\]

The radial derivatives are
\[
R_r=R\left(\frac1r+\frac12\frac{b_r}{b}-\frac{X_r}{X}\right),\qquad
R_t=R\left(\frac12\frac{b_t}{b}-\frac{X_t}{X}\right).
\]
These require the pinned derivative and actual same-stage metric-rate convention; they must not be replaced by a guessed conformal-factor power.

## 2. Geometric mass and the moving areal coordinate

For the current line element
\[
ds^2=-\alpha^2dt^2+\gamma_{rr}(dr+\beta\,dt)^2+R^2d\Omega^2,
\]
the invariant radial norm is
\[
\chi=g^{ab}\nabla_aR\nabla_bR
=\frac{X^2}{a}R_r^2-\frac{(R_t-\beta R_r)^2}{\alpha^2},
\]
so in \(G=1\) units the Misner–Sharp mass is
\[
m_{\rm MS}=\frac R2(1-\chi).
\]

The archive's areal-coordinate statement \(\partial_t m_{\rm MS}|_R=4\pi R^2 T_t{}^R\) is a derivative at **fixed areal radius**, not fixed isotropic coordinate \(r\). Since \(R=R(t,r)\),
\[
\left.\partial_t\right|_R
=\left.\partial_t\right|_r-\frac{R_t}{R_r}\left.\partial_r\right|_t.
\]
The stress component must also be transformed as a mixed tensor under \(y=(t,R(t,r))\), not copied from \(T_t{}^r\) in the isotropic chart. With \(x=(t,r)\), the required component is
\[
T_{t}{}^{R}\big|_{(t,R)}
=\frac{\partial x^\mu}{\partial y^t}
 \frac{\partial y^R}{\partial x^\nu}T_\mu{}^\nu,
\]
where \(\partial x^\mu/\partial y^t=(1,-R_t/R_r)\) and \(\partial y^R/\partial x^\nu=(R_t,R_r)\).

The volume density transforms with the coordinate Jacobian:
\[
\sqrt{-g}_{(t,R)}=\frac{\sqrt{-g}_{(t,r)}}{R_r}
\]
for a monotone outward map \(R_r>0\). This is separate from the spherical face-area factor \(4\pi R^2\).

## 3. Why this matters for the current results

The existing pointwise observer computes a local coordinate-r divergence in the isotropic chart. Its FV companion computes a cell-averaged balance using the chart's \(\sqrt{-g}\) measure and reconstructed spherical face fluxes. Both converge, and the FV residual is much smaller on the full-profile accepted states, but neither quantity is automatically identical to the archive's fixed-areal-radius Misner–Sharp current without the transformation above.

The archive record's areal-coordinate formula is therefore a valid source for the conservation target, but it cannot be applied by simply comparing \(\partial_t m_{\rm MS}|_r\) with the untransformed isotropic \(T_t{}^r\). That would mix both the time derivative and the tensor component.

## 4. Next executable check

Using the already SHA-verified accepted \(t=0.60\) states and the exact same-stage tangent endpoints:

1. compute \(R,R_r,R_t,\chi,m_{\rm MS}\) from the pinned \(a,b,X,\alpha,\beta\) fields and geometry rates;
2. transform the total stress tensor to the \((t,R)\) chart with the Jacobian above;
3. compare the fixed-\(R\) mass-rate identity against an independently reconstructed integrated shell balance;
4. repeat for N=60/120/240 and both directional epsilons, preserving source/vendor pins.

Until that comparison passes, retain pointwise and FV divergence as distinct convention-labelled diagnostics. No production or initializer change is implied; the full matter-coupled momentum-propagation identity and total sector conservation remain open, and \`scientific_admission=false\`.
