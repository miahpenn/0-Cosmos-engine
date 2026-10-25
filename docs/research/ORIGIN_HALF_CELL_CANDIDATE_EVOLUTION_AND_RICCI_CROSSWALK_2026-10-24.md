# Origin-half-cell candidate: evolved source and same-state Ricci crosswalk — 2026-10-24

**Disposition:** research/diagnostic only. No production equation, initializer, coefficient, gauge, boundary, vendor pin, threshold, or scientific-admission flag is changed. The baseline and candidate trajectories below are both already-completed, source-pinned artifacts. The Ricci crosswalk is offline and does not advance either state.

## Provenance and reuse

- Production source anchor: `4aec88d0e3ff028e0a65a8ac36b30948625df792`.
- Vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Baseline trajectory: [workflow #40127716287](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40127716287), artifacts [N=120 #12573419382](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40127716287/artifacts/12573419382) and [N=240 #12572764318](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40127716287/artifacts/12572764318).
- Origin-half-cell candidate trajectory: [workflow #40128647367](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40128647367), artifacts [N=120 #12574010545](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40128647367/artifacts/12574010545) and [N=240 #12574710341](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40128647367/artifacts/12574710341).
- Both artifacts at each resolution use the same production source blobs, vendor pin, `r_max=80`, CFL (=0.03), time steps, and accepted output times. Internal `SHA256SUMS.json` manifests are verified by the source-pinned observer and the new offline crosswalk.
- The candidate is the historical PR [#84](https://github.com/miahpenn/0-Cosmos-engine/pull/84)/[#106](https://github.com/miahpenn/0-Cosmos-engine/pull/106) origin-only half-cell quadrature: (I_0=dr\,\mathrm{dot}([401/720,-31/480,11/1440],S_{0:3})). It preserves the six-pass source/metric recurrence and every inter-centre trapezoid. Nothing in `engine/` or the vendor is changed in this PR.

## The evolved (A_a) L2 source

The source channel is computed on each saved accepted state from the pinned production profile. The initial centre source is strongly suppressed at both resolutions, but the long-enough bounded run is not resolution-independent.

| Resolution | Baseline (t=0) | Candidate (t=0) | Initial absolute reduction | Baseline (t=0.6) | Candidate (t=0.6) | Final candidate/baseline absolute ratio |
|---:|---:|---:|---:|---:|---:|---:|
| N=120 | (+1.0883\times10^{-3}) | (-8.6144\times10^{-5}) | 92.1% | (-1.2060\times10^{-3}) | (+9.9420\times10^{-5}) | 8.24% |
| N=240 | (+1.0875\times10^{-3}) | (-9.3882\times10^{-5}) | 91.4% | (+1.0250\times10^{-3}) | (-4.4358\times10^{-4}) | 43.3% |

At initialization, the candidate changes mainly the trace-free Ricci contribution; the signed lapse-Hessian contribution is almost unchanged. The initial metric quotient (q_0=(1-a_0/b_0)/r_0^2) moves much closer to the continuation fitted from cells 1–3, but crosses the fit rather than matching it exactly:

| Resolution | Baseline (q_0) | Candidate (q_0) | Candidate fitted continuation |
|---:|---:|---:|---:|
| N=120 | (-1.2514\times10^{-3}) | (-4.1945\times10^{-4}) | (-4.9921\times10^{-4}) |
| N=240 | (-1.2553\times10^{-3}) | (-4.1872\times10^{-4}) | (-5.0481\times10^{-4}) |

The evolved centre (A_a) changes sign across resolutions. The centre (K)-versus-even-fit residual is essentially unchanged at (t=0.6): approximately (-5.6551\times10^{-4}) for both arms at N=120 and (-5.6325\times10^{-4}) at N=240. So the candidate suppresses a particular initial/evolved source component; it does not close the full centre regularity problem.

## Same-state proper-radial Ricci crosswalk

Using the stored metric profiles ((a,b,X)), the audit recomputes the geometric trace-free Ricci through the areal-radius/proper-radial representation:

[
R=r\frac{\sqrt b}{X},\quad f=\frac{\sqrt a}{X},\quad
R_\ell=\frac{D_1^{\rm odd}R}{f},\quad
R_{\ell\ell}=\frac{D_2^{\rm odd}R-R_rD_1^{\rm even}f/f}{f^2}.
]

Then (R^\ell{}_\ell=-2R_{\ell\ell}/R), (R^\theta{}_\theta=(1-R_\ell^2-RR_{\ell\ell})/R^2), and (R_{\rm TF}=\frac23(R^\ell{}_\ell-R^\theta{}_\theta)). The vendor trace-free Ricci is recovered from the saved signed source `Aa_tracefree_ricci / alpha`. Both calculations use the same pinned cell-centred D1/D2 operators, so this is an independent geometric representation—not an independent discretization.

Key measure: candidate-to-baseline ratio of max absolute (|R_{\rm TF}^{\rm vendor}-R_{\rm TF}^{\rm proper}|) over the first eight cell centres. RMS ratios are over those same eight cells.

| N | Time | Signed gap at cell 0, baseline → candidate | First-8-cell max ratio | First-8-cell RMS ratio |
|---:|---:|---:|---:|---:|
| 120 | 0.0 | (+1.5944\times10^{-3}\to-1.7838\times10^{-4}) | 0.112 | 0.124 |
| 120 | 0.1 | (+1.3334\times10^{-3}\to-1.2681\times10^{-4}) | 0.095 | 0.106 |
| 120 | 0.3 | (-2.6120\times10^{-4}\to+1.7450\times10^{-4}) | 0.668 | 0.576 |
| 120 | 0.6 | (-1.4924\times10^{-3}\to+2.9165\times10^{-4}) | 0.195 | 0.295 |
| 240 | 0.0 | (+1.5900\times10^{-3}\to-1.9294\times10^{-4}) | 0.121 | 0.133 |
| 240 | 0.1 | (+6.3872\times10^{-4}\to+4.8264\times10^{-6}) | 0.016 | 0.022 |
| 240 | 0.3 | (-1.4241\times10^{-3}\to+3.6513\times10^{-4}) | 0.256 | 0.385 |
| 240 | 0.6 | (+3.6115\times10^{-3}\to+6.0611\times10^{-4}) | 0.168 | 0.169 |

The candidate initially reduces the first-eight-cell Ricci representation gap by about 9x at both resolutions. It remains materially smaller at (t=0.6), but its evolution is non-monotonic and resolution-sensitive. The proper-radial reconstruction itself remains irregular at the final high-resolution state: its cell-0 TF Ricci is about (-2.579\times10^{-3}) baseline versus (-1.044\times10^{-3}) candidate at N=240. A smaller representation gap is not the same as a regular or converged geometry.

The first few cells also show why reporting only cell 0 is inadequate: at (t=0) the candidate changes signed gaps in cells 0–3 and may reduce an apparent cell-0 spike while leaving alternating signed residuals nearby. The new offline tool reports those signed profiles, first-4/first-8/r≤10 maxima and RMS values for the vendor Ricci, proper-radial Ricci, their difference, and the stored (A_a) source channels.

## Interpretation and next gate

**Supported:** the initial origin half-cell moment is a significant contributor to the centre (A_a) source and the near-origin Ricci representation gap. The candidate is not merely algebraically cosmetic: the change persists in the stored evolution and improves a multi-cell geometric comparison.

**Not supported:** treating this quadrature candidate as a production fix; claiming convergence from two resolutions; concluding the current Ricci formulation alone is wrong; or treating the unchanged (K) residual as solved. The archived candidate is exact only for even source polynomials through (r^4), not for the fully coupled source without a truncation error.

The next work should use the already-saved stage profiles to identify the earliest captured PIRK operand at which the Ricci gap and (A_a) source diverge again. Do not launch another trajectory until that offline stage join has been exhausted. Scientific admission remains false.

## Reproducible offline tool

`tools/audit_origin_half_cell_saved_profile_ricci_crosswalk.py` verifies artifact manifests and source/vendor pins, then regenerates the table from the existing baseline/candidate artifact ZIPs or extracted directories. The dedicated [offline workflow](../.github/workflows/origin-half-cell-saved-profile-ricci-crosswalk.yml) downloads the already-produced artifacts and does not run the production integrator. Focused tests include a flat-space Ricci control and signed-profile norm accounting.
