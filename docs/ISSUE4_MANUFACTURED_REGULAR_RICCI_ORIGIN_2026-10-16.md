# Issue #4 manufactured regular-centre Ricci oracle

**Status:** test-only diagnostic; scientific admission remains **NOT ADMITTED**.

## Purpose

Classify the persistent first-cell trace-free Ricci discrepancy without using the production source fit or a time-evolved state. The comparison uses an exact continuum oracle for a known regular metric, the current proper-radial formula, and the vendor BSSN Ricci route with both the production-discrete metric connection and a continuum metric-derived connection.

## Manufactured metric

Use the areal-radius form

\[
dl^2 = \frac{dr^2}{B(r)} + r^2d\Omega^2,\qquad
B(r)=1+b_2r^2+b_4r^4+b_6r^6,
\]

with \(b_2=2\times10^{-5}\), \(b_4=2\times10^{-7}\), and \(b_6=10^{-12}\). The engine metric representation is set by \(a=B^{-2/3}\), \(b=B^{1/3}\), and \(X=B^{1/6}\), so the identities \(ab^2=1\), \(b=X^2\), \(a=X^{-4}\), \(X^6=B\), and areal radius \(r\) are testable independently.

The warped-product identity gives

\[
R^l{}_l=-B'/r,\qquad
R^\theta{}_\theta=(1-B-rB'/2)/r^2
\]

and hence

\[
R_{\rm TF}=R^r{}_r-\frac13R
=\frac23\left[\frac{B-1}{r^2}-\frac{B'}{2r}\right].
\]

For the even sextic, the \(b_2\) terms cancel exactly:

\[
R_{\rm TF}=-\frac23b_4r^2-\frac43b_6r^4,\qquad
\lim_{r\to0}\frac{R_{\rm TF}}{r^2}=-\frac23b_4.
\]

This analytic identity is the oracle. The proper-radial and vendor BSSN results are measured against it; no tolerance is chosen from the observed production gap.

## Negative control and reporting

A second profile adds \(c_1r\), with \(c_1=10^{-5}\), intentionally breaking even regularity. Its \(R_{\rm TF}/r^2\) has a leading \(r^{-3}\) signature, so the magnitude should grow by approximately eight when the first-cell radius halves from \(N=160\) to \(N=320\). This control confirms that the analytic witness is sensitive to the missing origin parity.

At \(N=160,240,320\), retain the signed cell-0/cell-1 errors, fixed-annulus maximum and RMS, per-cell \(R_{\rm TF}/r^2\), both Lambda constructions, all metric identities, state digests, profile arrays, and SHA-256 manifest. The vendor/proper-radial errors are reported diagnostically rather than gated against a tolerance inferred after inspection.

## Scope boundary

The test uses copied geometry only: zero production steps, no changes to engine equations, initializer, vendor, dependency, campaign configuration, or scientific-acceptance tolerances. Passing validates the manufactured-oracle harness and its controls; it does **not** prove total stress-energy conservation, the matter-coupled Brown/Bianchi constraint-propagation identity, or long-time scientific admission. Any production repair still requires an isolated, source-derived discrete identity failure and the independent matter/propagation gates.
