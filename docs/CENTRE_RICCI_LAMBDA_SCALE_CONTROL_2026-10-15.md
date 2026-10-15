# Centre-Ricci Λ radial-term length-scale control — 15 October 2026

**Disposition:** research-only, manufactured profiles, no evolution. Scientific admission remains **NOT ADMITTED**.

## Purpose and non-overlap

This is a complementary scale-sensitivity experiment for the derivative-representation question in [PR #254](https://github.com/miahpenn/0-Cosmos-engine/pull/254). It is deliberately distinct from the regular-centre C5/Taylor-series work: it does not fit or assume a truncated (q_2,q_4,q_6) series. Instead it uses a smooth analytic profile at three length scales and independently differentiates its closed form. It requires no saved machine artifact and performs no engine initialization, source capture, or time evolution.

## Manufactured geometry and observable

For amplitude (A=0.1), length (L\in\{4,8,16\}), and fixed domain (0<r<32), define

\[
q(r)=1+A\left(1-e^{-(r/L)^2}\right),\qquad
a=q^{-4/3},\quad b=q^{2/3},\quad X=q^{-2/3}.
\]

The parametrization is determinant-consistent (\(ab^2=1\)); it gives \(d\ell=(\sqrt{a}/X)dr=dr\) and areal radius \(R=r\sqrt{b}/X=rq\). The singular-looking radial component is

\[
\Lambda_r=\frac{2}{r}\left(\frac{1}{b}-\frac{1}{a}\right)
=\frac{2}{r}\left(q^{-2/3}-q^{4/3}\right).
\]

The independent continuum derivative is evaluated from the analytic profile:

\[
q'=\frac{2Ar}{L^2}e^{-(r/L)^2},\quad
\delta=q^{-2/3}-q^{4/3},\quad
\delta'=\left(-\frac{2}{3}q^{-5/3}-\frac{4}{3}q^{1/3}\right)q',
\]
\[
\Lambda_r'=\frac{2\delta'}{r}-\frac{2\delta}{r^2}.
\]

This makes the error comparison independent of both discrete differentiation routes. All sample points are cell centres (r_i=(i+1/2)\Delta r), so no quotient is numerically evaluated at (r=0).

## Compared channels

1. **Sampled-odd route:** apply the pinned fourth-order first derivative with odd parity to the sampled (Lambda_r).
2. **D1 chain route:** assemble (-2\delta/r^2 + 2[-D_1(b)/b^2+D_1(a)/a^2]/r) with the pinned fourth-order even-field (D_1).
3. **Independent reference:** use the closed-form (Lambda_r') above.

We report the signed error for each discrete route at cells 0 and 1, each error divided by (r_i^2), and the maximum/signed/rms error on the **same fixed physical annulus** (1\le r\le4). We also report the D4-minus-chain route gap and verify the existing three-term Ricci response recombination and isolated radial-term identity. Refinement order is computed only from fixed-annulus maxima; the centre cells are not treated as fixed-coordinate convergence samples.

## Results

Both discrete representations converge against the independent analytic derivative at approximately fourth order on the fixed annulus for all three length scales. The D4-minus-chain-rule gap also converges at fourth order: it is a finite-stencil representation difference on these smooth controlled profiles rather than a resolution-independent offset.

| Length scale L | Fixed-annulus error | N=80 | N=160 | N=320 | Observed order 80→160 / 160→320 |
|---:|---|---:|---:|---:|---:|
| 4 | Odd-sampled D4 minus exact | 1.1940e-6 | 6.7787e-8 | 4.4801e-9 | 4.139 / 3.919 |
| 4 | D1 chain minus exact | 1.2099e-5 | 7.3333e-7 | 4.7176e-8 | 4.044 / 3.958 |
| 4 | D4 minus chain | 1.0905e-5 | 6.6554e-7 | 4.2695e-8 | 4.034 / 3.962 |
| 8 | Odd-sampled D4 minus exact | 2.6160e-8 | 1.6045e-9 | 1.0146e-10 | 4.027 / 3.983 |
| 8 | D1 chain minus exact | 2.2793e-7 | 1.4147e-8 | 8.9063e-10 | 4.010 / 3.989 |
| 8 | D4 minus chain | 2.0177e-7 | 1.2542e-8 | 7.8917e-10 | 4.008 / 3.990 |
| 16 | Odd-sampled D4 minus exact | 4.4169e-10 | 2.7484e-11 | 1.7252e-12 | 4.006 / 3.994 |
| 16 | D1 chain minus exact | 3.7331e-9 | 2.3293e-10 | 1.4587e-11 | 4.002 / 3.997 |
| 16 | D4 minus chain | 3.2915e-9 | 2.0544e-10 | 1.2861e-11 | 4.002 / 3.998 |

For example, at (L=4), the signed first-cell D4-minus-chain gap is (-1.3190\times10^{-5}), (-8.4708\times10^{-7}), and (-5.3307\times10^{-8}) at (N=80,160,320); dividing by (r_0^2) gives (-3.2975\times10^{-4}), (-8.4708\times10^{-5}), and (-2.1323\times10^{-5}). These centre-scaled witnesses are retained separately from the fixed-annulus order fit. Full three-term recombination and the isolated radial-term identity are checked to (2\times10^{-12}) in the test.

## Interpretation and limits

The controlled family shows that the sampled-odd D4 and D1 chain-rule representations each approximate the same known continuum derivative at fourth order over the tested scale/domain matrix, and their difference decreases at fourth order too. This is useful calibration for expected finite-stencil error and units/sign handling.

It does **not** explain the roughly (10^{-3}) first-cell response on the saved production predictor states, prove the production initializer lies on this analytic family, choose one production representation, identify a unique defective term, or certify evolved centre regularity. The exact production gap remains a separate diagnostic question. No production, vendor, dependency, physical equation, or campaign was changed by this experiment.
