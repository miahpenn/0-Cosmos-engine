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

The production path now advances geometry, scalar fields, and conservative
matter in lockstep through `engine/production_kernel.py`, with capability checks
in `engine/production_contract.py` and invariant diagnostics in
`engine/invariant_diagnostics.py`.

That closes the software-integration requirement; it does **not** by itself
establish scientific admission. Long-horizon radiation continuation, centre
regularity and constraint propagation, multi-resolution convergence, and
fixed-areal-radius Misner–Sharp balance remain independently gated. See
[docs/STATUS.md](STATUS.md) for the current source-pinned status.
