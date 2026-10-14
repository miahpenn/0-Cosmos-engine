# Applying the local GEAR derivation to V5.5 diagnostics

**Status:** research diagnostic only; no evolution equations or physical-source assumptions changed.

## Derived relation and metric translation

The local-library derivation dated 2026-10-13 supports the sequence

\[
(S,D)\ \to\ T_{\mu\nu}\ \to\ \text{solved metric/lapse}
\ \to\ (\ell,\tau)\ \to\ \text{measured flux}.
\]

In the archived polar-areal local model, the spatial line element is
\(d\ell=a_{\rm PA}(t,r)\,dr\), so
\(\ell(t,r)=\int_0^r a_{\rm PA}(t,r')\,dr'\). The archived cell-centred code uses midpoint quadrature.

The current V5.5 engine is written in different metric variables. Its solved radial spatial metric is

\[
\gamma_{rr}=\frac{X^2}{a},
\qquad
d\ell=\sqrt{\gamma_{rr}}\,dr.
\]

Thus the corresponding main-engine diagnostic at uniform cell centres is

\[
\ell_j=\Delta r\left[
\sum_{i=0}^{j}\sqrt{\gamma_{rr,i}}
-\frac12\sqrt{\gamma_{rr,j}}
\right].
\]

This map is implemented as engine.worldtube.proper_radial_coordinate. It returns a diagnostic coordinate computed from the current solved metric; it is not an evolved field.

## Where it is applied

V55ProductionKernel.diagnostics now emits proper_radius_sigma alongside the existing worldtube surface observables. With profiles=True, it also emits a proper_radius array corresponding to the same radial cells as areal_radius, chi, and misner_sharp. Existing tau_rate = alpha(0) and trapezoidal proper-time accumulation are unchanged. Existing total-stress coordinate_energy_flux, flux_T, and work_pR are unchanged.

The proper_radius and areal_radius outputs are intentionally both retained: areal radius is the invariant sphere-area label R, while proper radius measures distance along the spatial slice. They are not interchangeable in a curved slice.

## Validation and interpretation boundary

Tests check that the map gives ell=r on the flat production metric, gives the expected linear scaling on a manufactured constant radial metric, rejects invalid metric data, and is exposed by the production profile output.

This applies the part of the local derivation that can be translated directly into the current engine's geometry. It deliberately does not copy the archived polar-areal scalar-mode flux formula into the BSSN kernel: the metric variables, canonical momenta, and source normalization must first be mapped and independently verified. It also does not alter the energy/current ledger, infer a D-to-matter conversion, add a guessed COSMOS source, or claim gravitational balance closure.

The archive's reported Misner–Sharp-type residual at r=10,T=10 decreases by about a factor of two under each tested resolution doubling, but remains nonzero at finite resolution. Separately, the saved proper-radial flux-peak CSV has an unresolved source-provenance caveat: the older experiment helper has a centre-divergence form inconsistent with the later corrected face-consistent source. The peak track is therefore not admitted here as a corrected-source regression witness. A full gravitational balance/continuum gate remains separate from the proper-coordinate diagnostic.
