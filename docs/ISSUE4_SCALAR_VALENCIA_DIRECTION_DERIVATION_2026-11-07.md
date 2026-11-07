# Issue #4 scalar and Valencia direction derivation — pinned source transcription

**Date:** 2026-11-07  
**Status:** derivation note; no fixture implementation, archive, production import, or production change  
**Purpose:** define the equation-derived scalar rates and the pinned semi-discrete fluid rates needed by the proposed `coupled_analytic_matter_direction` control.

## 1. Provenance and scope

The transcription below is tied to the Issue #4 helper baseline `7f47669a1c0edfafef9fc5f27b2db80840d3f035` and these exact source blobs:

| Source | Blob SHA-1 |
|---|---|
| `engine/scalar_system.py` | `ffc119b1a6f2d495d972323901f76a85783b0a1b` |
| `engine/matter_system.py` | `e52a117690253ddad7b24d116145a65474a351d2` |
| `engine/matter_rhs.py` | `f41180efbcff2c72cf36477a765318ab539fac80` |
| `engine/valencia.py` | `03ed86398f424c03d9d6a3e646da9222a34fac5c` |
| `engine/v55_matter.py` | `0ee3c51fac7dc46c63feba2747a750542bcee3d0` |
| `engine/production_kernel.py` | `9685e76bb5a6f26d6bf383e240de24f37523f265` |

The exact `matter_system.py` / `matter_rhs.py` path is the production semi-discrete finite-volume operator to transcribe. Do not substitute a continuum Euler formula, the alternate mixed-tensor convenience methods in `engine/valencia.py`, or a generic textbook Valencia flux for the actual stored-array and source conventions in `matter_system.py`.

The work remains research-only. No production module is to be imported by the eventual independent constructor; these source blobs are references for a separately written transcription. No equation coefficient, sign, branch, stencil, limiter, tolerance, or source normalization may be fitted.

## 2. Scalar field rates

Use the same fixed parity-aware linear operators (D_1,D_2) as the pinned scalar source. Define
[
q=\frac{X^2}{a},\qquad
C=\frac{D_1a}{2a}-\frac{D_1b}{b}+\frac{D_1X}{X}-\frac2r,
\qquad
\mathcal L_f=q\left(D_2f-(D_1f)C\right).
]
For each (f\in\{S,D,\phi\}), with conjugate momentum (p_f\in\{PS,PD,\Pi\}), the pinned source map is
[
\dot f=\alpha_{out}p_f+\beta D_1f,
\]
[
\dot p_f=\beta D_1p_f+\alpha_{out}\mathcal L_f
+\alpha_{out}Kp_f+q(D_1\alpha_{out})(D_1f)
-\alpha_{out}V_f'(f)+s_f.
]
For the zero-shift CMC fixture, (eta=\dot\beta=0) is an explicit gauge condition. The scalar evolution uses the refreshed **output** lapse (alpha_{out}) and its direction. The incoming lapse (alpha_{in}) is a separate input to the CMC refresh map and must not be silently substituted for (alpha_{out}).

The pinned potential derivatives are:
[
V'_S=S,\quad V''_S=1;\qquad
V'_D=-D,\quad V''_D=-1;
]
[
V'_\phi=-V_0e^{-\phi},\quad V''_\phi=V_0e^{-\phi},
\quad V_0=8.242522415500654\times10^{-5}.
]
The scalar exchange source is (s_\phi=\beta_{DM}\rho_{DM}), (\beta_{DM}=-0.04); (s_S=s_D=0). The DM density here is the recovered archive-normalized primitive density, not a freely chosen density profile. Its direction includes the conservative-state and metric response derived in §4.

### 2.1 Analytic directional derivative

Let a dot on an operand denote the independent instantaneous direction, and let (D_1,D_2) remain fixed linear maps. Then
[
\frac{\dot q}{q}=2\frac{\dot X}{X}-\frac{\dot a}{a},
]
[
\dot C=
\frac12\left(\frac{D_1\dot a}{a}
-\frac{(D_1a)\dot a}{a^2}\right)
-\left(\frac{D_1\dot b}{b}
-\frac{(D_1b)\dot b}{b^2}\right)
+\left(\frac{D_1\dot X}{X}
-\frac{(D_1X)\dot X}{X^2}\right),
]
[
\dot{\mathcal L}_f=
\dot q\left(D_2f-(D_1f)C\right)
+q\left(D_2\dot f-(D_1\dot f)C-(D_1f)\dot C\right).
]
The field-rate direction is
[
\ddot f=\dot\alpha_{out}p_f+\alpha_{out}\dot p_f
+\dot\beta D_1f+\beta D_1\dot f.
]
The momentum-rate direction is
[
\begin{aligned}
\ddot p_f={}&\dot\beta D_1p_f+\beta D_1\dot p_f
+\dot\alpha_{out}\mathcal L_f+\alpha_{out}\dot{\mathcal L}_f
+\dot\alpha_{out}Kp_f+\alpha_{out}\dot Kp_f
+\alpha_{out}K\dot p_f\\
&+\dot q(D_1\alpha_{out})(D_1f)
+q(D_1\dot\alpha_{out})(D_1f)
+q(D_1\alpha_{out})(D_1\dot f)\\
&-\dot\alpha_{out}V'_f(f)
-\alpha_{out}V''_f(f)\dot f+\dot s_f.
\end{aligned}
]
For fixed coupling, (dot s_\phi=\beta_{DM}\dot\rho_{DM}). All terms, including the lapse-gradient and spatial-operator responses, remain present even when a particular base profile makes one term small.

The independent constructor must compute these rates without calling `scalar_rhs_arrays`; compare each reconstructed scalar rate to this transcription and to symmetric (+\epsilon/-\epsilon) reconstructions at both locked epsilons.

## 3. Actual semi-discrete Valencia state and RHS

For the pinned `matter_system.py` path, the evolved arrays are represented by the following Eulerian Valencia quantities:
[
D=\rho W,\quad E=(\rho+p)W^2-p,\quad
S_r=(\rho+p)W^2\gamma_{rr}v^r,
\]
[
U=(\sqrt\gamma D,\sqrt\gamma E,\sqrt\gamma S_r).
]
Here (W=(1-v^2)^{-1/2}), (v^2=\gamma_{rr}(v^r)^2), (\gamma_{rr}=a/X^2), and
[
\sqrt\gamma=\frac{r^2\sqrt a\,b}{X^3}.
]
The names `rest`, `energy_t`, and `momentum_r` in the stored structure are array/API names; use the above actual `matter_system.py` arithmetic and do not infer a different sign or lapse factor from a similarly named tensor component in another helper.

For each face, `matter_system.py` linearly averages the *primary* (a,b,X,\alpha,\beta) to the face radius and constructs the face metric from those averaged fields. It reconstructs cell primitive ((\rho,p,v^r)) with the pinned minmod piecewise-linear rule, using the face (\gamma_{rr}), then evaluates the species-specific HLL flux. The outer face uses the last-cell primitive without a reconstructed slope and the outer-extrapolated primary metric; the inner face flux is exactly zero.

The Valencia flux vector for dust and baryons is
[
F=\sqrt\gamma\left(
D(\alpha v^r-\beta),\;
\alpha S^r-\beta E,\;
\alpha\big((\rho+p)W^2v^2+p\big)-\beta S_r
\right),
\qquad S^r=\gamma^{rr}S_r.
]
Radiation has an identically zero first (rest-current) component; its energy and momentum fluxes use the same second and third expressions. The HLL wave speeds are those in the pinned helper, with (c_s=0) for dust and (c_s=1/\sqrt3) for radiation.

At each cell, define
[
S_{rr}=(\rho+p)W^2v^2\gamma^{rr}+p\gamma^{rr},
\qquad S_{\theta\theta}=p\gamma_{\theta\theta}^{-1},
]
[
K_{rr}=\gamma_{rr}(K/3+A_a),\qquad
K_{\theta\theta}=\gamma_{\theta\theta}(K/3-A_a/2).
]
The exact pinned geometric source pair is
[
S_E=\sqrt\gamma\left[
\alpha(K_{rr}S_{rr}+2K_{\theta\theta}S_{\theta\theta})
-S^rD_1\alpha
\right],
]
[
S_S=\sqrt\gamma\left[
-E D_1\alpha+S_rD_1\beta
+\frac{\alpha}{2}\left(S_{rr}D_1\gamma_{rr}
+2S_{\theta\theta}D_1\gamma_{\theta\theta}\right)
\right].
]
For dark matter only, with the locked (Q_{DM,\nu}=\beta_{DM}\rho_{DM}\partial_\nu\phi),
[
Q_E=-\sqrt\gamma(Q_{DM,t}-\beta Q_{DM,r}),\qquad
Q_S=\sqrt\gamma Q_{DM,r},
]
where (Q_{DM,t}=\beta_{DM}\rho_{DM}\dot\phi) and
(Q_{DM,r}=\beta_{DM}\rho_{DM}D_1\phi). Baryons and radiation have zero exchange arrays. The scalar (phi) source has the opposite sector-exchange sign, as pinned above.

The cell-centered semi-discrete rates returned by the pinned `matter_rhs.species_rhs` arithmetic are
[
\dot U_i=-\frac{F_{i+1/2}-F_{i-1/2}}{\Delta r}
+(0,S_{E,i}+Q_{E,i},S_{S,i}+Q_{S,i}).
]
For radiation, (dot U_{D}=0) identically. Preserve the helper's operation ordering and boundary-face rules. The plus/minus cross-check must reconstruct the complete face fluxes and sources independently; it must not use a finite timestep or invoke the production stepper.

## 4. Conservative-to-primitive directions

The primitive tangent must include the metric response of (\sqrt\gamma) and (\gamma^{rr}):
[
\frac{\dot{\sqrt\gamma}}{\sqrt\gamma}
=\frac12\frac{\dot a}{a}+\frac{\dot b}{b}-3\frac{\dot X}{X},
\qquad
\frac{\dot\gamma^{rr}}{\gamma^{rr}}
=2\frac{\dot X}{X}-\frac{\dot a}{a}.
]

### 4.1 Dust (dark matter and baryons)

The pinned recovery uses
[
D=U_D/\sqrt\gamma,\quad S_r=U_S/\sqrt\gamma,\quad
s^2=\gamma^{rr}S_r^2,\quad Z=\sqrt{D^2+s^2},
]
[
W=Z/D,\qquad \rho=D^2/Z,\qquad
v^r=\gamma^{rr}S_r/Z,qquad p=0.
]
Thus
[
\dot D=\dot U_D/\sqrt\gamma-D\,\dot{\sqrt\gamma}/\sqrt\gamma,
\quad
\dot S_r=\dot U_S/\sqrt\gamma-S_r\,\dot{\sqrt\gamma}/\sqrt\gamma,
]
[
\dot{s^2}=\dot\gamma^{rr}S_r^2+2\gamma^{rr}S_r\dot S_r,
\quad
\dot Z=(D\dot D+\tfrac12\dot{s^2})/Z,
]
[
\dot W=\dot Z/D-W\dot D/D,
\quad
\dot\rho=2D\dot D/Z-D^2\dot Z/Z^2,
\]
[
\dot v^r=\dot\gamma^{rr}S_r/Z
+\gamma^{rr}\dot S_r/Z-\gamma^{rr}S_r\dot Z/Z^2.
]
The stored energy array remains an independently evolved conservative channel and must be retained in its own RHS and residual ledger; the pinned dust primitive recovery does not use it to recover ((\rho,v^r)).

### 4.2 Radiation

The pinned radiation recovery uses
[
E=U_E/\sqrt\gamma,\quad S_r=U_S/\sqrt\gamma,\quad
s^2=\gamma^{rr}S_r^2.
]
Its pressure (p) is the root of
[
R(p;E,s^2)=3p-E+\frac{4p s^2}{(E+p)^2-s^2}=0,
\quad 0\le p\le E/3.
]
Then
[
\rho=3p,\quad v^2=\frac{s^2}{(E+p)^2},\quad
v^r=\operatorname{sgn}(S_r)\sqrt{v^2/\gamma_{rr}}.
]
For a regular interior root, put (T=(E+p)^2-s^2). The implicit derivative is
[
R_p=3+\frac{4s^2}{T}
-\frac{8ps^2(E+p)}{T^2},\quad
R_E=-1-\frac{8ps^2(E+p)}{T^2},\quad
R_{s^2}=\frac{4p(E+p)^2}{T^2},
]
[
\dot p=-\frac{R_E\dot E+R_{s^2}\dot{s^2}}{R_p},
\qquad \dot\rho=3\dot p,
]
[
\dot{v^2}=\frac{\dot{s^2}}{(E+p)^2}
-\frac{2s^2(\dot E+\dot p)}{(E+p)^3}.
]
The radial velocity derivative follows by differentiating its signed expression on a fixed nonzero-(S_r) branch. The (s^2=0) branch is separate: do not differentiate through the sign/square-root representation as if it were smooth when a perturbation changes the momentum branch. A base at or near a recovery branch boundary is inadmissible for this tangent control.

The radiation pressure root is differentiated implicitly, not by differentiating the bisection iteration count. The actual conservative energy and momentum rates remain primary archived outputs.

## 5. Piecewise numerical branches are part of the derivative contract

The pinned finite-volume RHS is not globally differentiable across all limiter and wave-speed branches. The independent tangent must explicitly record the active branches and strict margins for:

- minmod slope selection (including zero-slope/tie cases);
- primitive reconstruction's nonnegative density/pressure guards and velocity bound;
- HLL left/right signal-speed extrema and the zero-speed branch;
- radiation admissibility and the implicit pressure-root interior;
- the dust/radiation primitive recovery branches.

For a valid centered-difference witness, the base and both (+\epsilon/-\epsilon) members must remain on the same differentiable branch, with nonzero branch margins. If a tie, clamp, limiter switch, or wave-speed branch switch occurs, report that witness as non-differentiable/inconclusive; do not smooth, clip, freeze an invalid branch, or relax the gate to manufacture agreement. A separate deliberately branch-crossing control may test rejection behavior.

The analytic tangent must differentiate the *selected active formula* for face reconstruction, face metric, HLL flux, geometric source, and exchange source. The independent (+\epsilon/-\epsilon) reconstruction must execute the actual piecewise algorithm afresh on each member so branch stability is demonstrated rather than imposed.

## 6. Continuum identities versus the discrete Valencia contract

For each pressureless species, the continuum continuity witness in spherical coordinates is
[
\partial_t(\sqrt\gamma\rho W)
+\partial_r\left[\sqrt\gamma\rho W(\alpha v^r-\beta)\right]=0.
]
The continuum Euler/acceleration witness is the declared covariant equation
[
(\rho+p)u^\mu\nabla_\mu u_\nu
+(\delta_\nu{}^\mu+u_\nu u^\mu)\nabla_\mu p
=(\delta_\nu{}^\mu+u_\nu u^\mu)Q_\mu.
]
For baryons (Q_\mu=0); for dark matter (Q_{DM,\mu}=+\beta_{DM}\rho_{DM}\partial_\mu\phi); radiation has (p=\rho/3), no independently evolved rest-current channel, and its pressure-gradient acceleration is retained.

These continuum identities are separate equation witnesses. The actual archived directional arrays are generated by the pinned finite-volume semi-discrete RHS in §3. Do not demand that a finite-volume truncation residual vanish pointwise as though the numerical flux divergence were the continuum derivative. Report (a) continuum continuity/acceleration residuals using their declared differential operators and (b) exact reconstruction of the pinned discrete conservative RHS as separate channels. The continuum witnesses are not a substitute for the finite-volume directional contract, and the discrete RHS pass is not by itself a continuum conservation proof.

## 7. Required next checks before coupled fixture construction

1. Implement an independent scalar transcription and its analytic tangent; compare all six scalar rates at both locked epsilons.
2. Implement independent dust primitive/conservative tangent controls and radiation implicit-recovery tangent controls, with explicit branch and admissibility margins.
3. Implement the pinned face metric, minmod reconstruction, HLL flux, geometric source, and DM exchange transcription plus analytic directional derivative. Keep raw face/source arrays in the diagnostic output.
4. Compare independently reconstructed conservative RHS arrays and their tangents against centered (+\epsilon/-\epsilon) reconstructions at both locked epsilons.
5. Separately compute continuum continuity and acceleration witnesses from the recovered primitive histories; label their operator and truncation scope.
6. Only after those pass, combine these rates with the already specified geometric/primary/connection and stage-aware lapse directions, then derive six-sector projection tangents.

No archive, fixture schema, observer run, production call, or admission flag change is authorized by this derivation note.

## 8. Explicit unresolved details

- The regularity-complete base profile is still required; the Gaussian incoming-lapse centre gap recorded in the main specification is not waived.
- The exact parity of every fluid primitive and conservative array must be declared (radial momentum is odd; scalar densities and energy channels are even), with a regular centre construction and parity-aware derivative/ghost handling.
- The scalar DM exchange density direction must be derived through the dust recovery tangent, including metric volume and inverse-metric response.
- A radiation base must remain strictly inside the recovery domain (E>|S|) wherever (S\ne0), and its pressure root must remain interior with (R_p\ne0).
- No claim is made here that the pinned continuum acceleration witness, full total stress-energy divergence, six-sector closure, or physical exchange admission has passed.

**Disposition:** source equations and a derivative implementation contract are now transcribed for the scalar and fluid lanes. This is a derivation milestone, not a passed scalar/fluid test. The next implementation remains research-only and must reproduce the exact pinned formulas independently before any coupled fixture is assembled.
