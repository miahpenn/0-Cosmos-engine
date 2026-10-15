# Analytic half-cell injection into the self-consistent copied initializer — 15 October 2026

**Status: research-only; no production change; scientific admission NOT ADMITTED.** This experiment follows the continuous integral/enclosure audit and the source-fit-family copied-initializer diagnostics. It asks a more specific question: does substituting the independently computed continuous first-half-cell integral into the self-consistent discrete outer-grid recurrence produce a reproducible signed centre-Ricci and constraint response, and how much can that response vary across the analytic integral enclosure?

## Method

The default uncharged initializer admits the continuous-source form

\[
s(r)=A(r)-q(r)B(r), \qquad I(r)=rB(r), \qquad
I'(r)+\frac{q(r)}{r}I(r)=A(r), \quad I(0)=0.
\]

The dedicated continuous-source oracle provides (i) a Richardson reference for I(r0), where r0=dr/2, and (ii) an analytic model-conditional enclosure I_A-Q_max <= I(r0) <= I_A. This experiment reruns the copied initializer at each of N=160,240,320, using three first-half-cell inputs: the lower enclosure endpoint, Richardson nominal, and upper endpoint.

Only the first half-cell moment changes. The source on all sampled centres is recomputed from each candidate B, and the usual inter-centre trapezoidal recurrence is preserved. Each arm reports the six-iteration analogue and the converged fixed point separately. The selected polynomial fit is still calculated as a diagnostic comparator but is **not applied** in these analytic-integral arms.

On each candidate geometry, the diagnostic records signed vendor-vs-independent-proper-radial trace-free Ricci gaps and native-grid Hamiltonian, momentum, connection, and determinant witnesses. The Lambda derivative comparison is reported as a distinct offline arm: it compares the pinned vendor D4(Lambda) contribution with the derivative of the metric-to-Lambda identity on the exact same candidate metric. It is not folded into the primary Ricci gap and is not proposed as a production operator substitution.

## Source and safety fence

- Engine source anchor: 285a2ee0c1a3b97cddeb4519d8e749109b115aab.
- Vendor commit: d6052d605673ce9d82cdc99b0df79855fa2ea215.
- No edits to production engine/, vendor, dependency files, physical equations, or campaign controls are permitted in this branch.
- Every observer checks initialized state hashes before and after; the audit takes zero evolution steps.
- The copied geometry retains the initialized scalar/matter fields and unchanged lapse/shift/extrinsic-curvature fields except the spatial-metric variables required by the radial B reconstruction. The report explicitly labels this a hybrid continuous-half-cell/discrete-outer-grid counterfactual.
- The analytic enclosure is conditional on the default uncharged Gaussian source, constant phi/Pi/Kt, radiation enabled, Kt<0. Charged dust/Maxwell are outside scope. Endpoints are float64-evaluated, not directed-rounding interval-certified.

## Discriminating outcomes

The primary question is not whether a candidate gap gets numerically smaller at one cell. It is whether the sign and magnitude of the Ricci response are stable across resolutions and across the lower/nominal/upper analytic integral inputs, while metric identities and signed constraints remain explicitly reported.

A useful result would show the three integral inputs yield nearly identical metric/Ricci/constraint witnesses compared with the larger production-to-counterfactual change. A material endpoint sensitivity would instead indicate that the analytic error enclosure is too wide for this coupled discrete map, even if it encloses the local continuum integral.

This experiment can narrow the initialization hypothesis, but it cannot by itself prove full-horizon convergence or conservation. Keep scientific admission NOT ADMITTED pending independent numerical-method review, signed response interpretation, multi-resolution/full-horizon checks, and the repository's separate physics/release gates.

## Reproduction

The dedicated workflow compiles the diagnostic, checks source/vendor pins, runs focused tests, performs all nine counterfactuals (three integral inputs at three resolutions, each with six-step and fixed-point witnesses), and uploads deterministic JSON/CSV/NPZ files with SHA-256 checksums.
