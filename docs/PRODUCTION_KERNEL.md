# Production kernel provenance

The repository pins RonBibb/bb-palatini-unified-r0 as a git submodule at
commit d6052d605673ce9d82cdc99b0df79855fa2ea215.

The pinned source supplies the reference-metric spherical BSSN/PIRK numerical
operators and moving-puncture gauge machinery. Its Palatini/complex-scalar
matter equations are not imported as V5.5 physics.

The local adapter adds only the frozen V5.5 stress-energy projections:
- local S scalar with V=+S^2/2;
- local inverted D scalar with V=-D^2/2;
- corrected COSMOS scalar normalization;
- conservative beta-coupled dark matter;
- baryonic dust;
- radiation p=rho/3.

The final production engine is not declared complete until the conservative
matter states are advanced in lockstep with the BSSN stages and the invariant
witnesses are recorded throughout the strong-field continuation.
