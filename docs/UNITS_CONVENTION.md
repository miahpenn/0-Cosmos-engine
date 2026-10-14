# Gravitational source units and Friedmann consistency

This note records the density normalization used by the COSMOS background and the local V5.5 production source path. The convention is software-level bookkeeping; it does not admit the model's physics or identify the local S/D sectors with COSMOS components.

## Convention

Use geometrized time and length units (`c=1`) and fix `G=1`. Define the archive/background density

\[
\rho_E \equiv \kappa\rho_{\rm physical},\qquad
\kappa=8\pi G=8\pi.
\]

The COSMOS lane evolves `rho_E` and its matching pressure `p_E`; therefore its background equations are

\[
H^2=\frac{\rho_E}{3},\qquad
\dot H=-\frac{\rho_E+p_E}{2}.
\]

The COSMOS density numbers (for example, `RHO_DM_PRESENT`) are in this Einstein-normalized convention. The `G` field on `CosmosParams` is retained as compatibility metadata for archived parameter records, but it is constrained to `1`: changing it without re-normalizing every density and potential is unsupported.

## Local production conversion

The Valencia conserved states and archive initializer carry the same normalized fluid-density values `rho_E`. The local matter projection converts them once to physical stress-energy:

\[
\rho_{\rm local}=\frac{\rho_E}{\kappa}.
\]

The production Hamiltonian diagnostic then uses the standard local Einstein source `16*pi*rho_local`, and the local geometry's flat-FLRW constraint is

\[
3H^2=\kappa\rho_{\rm local}
\quad\Longleftrightarrow\quad
H^2=\frac{\rho_E}{3}.
\]

Thus the two lanes agree only when this conversion is applied exactly once. Do not multiply a projected local `total.rho` by `kappa` a second time inside the BSSN source, and do not interpret `total.rho` itself as the unprojected archive-normalized density.

The local scalar sector is not one uniform normalization: the COSMOS `phi` projection is divided by `kappa`, while the archive-locked local S/D projections retain their own normalization. This note makes no claim that the independent S/D sector has a homogeneous COSMOS counterpart, and it makes no change to their equations. The density conversion establishes the source scale used by the Friedmann/Hamiltonian identity; it does not by itself prove the local scalar–DM exchange identity, derive the value or sign of beta, or construct the missing local-to-COSMOS reciprocal source. Those remain separate gates.

## Regression scope and limits

`tests/test_friedmann_unit_convention.py` checks (1) the COSMOS present-state Friedmann relation and total continuity including cancellation of the declared scalar–DM exchange, and (2) a synthetic homogeneous, flat local state through the actual production diagnostic path. In the local test, the Hamiltonian residual must close with `16*pi*total.rho` and the measured `H_eff` must satisfy `H_eff**2 = kappa*total.rho/3`.

These are normalization/consistency regressions, not evidence for a reciprocal local↔COSMOS coupling, the outer-boundary time identification, a closed radiation trajectory, a physically admitted turnaround, or full scientific validity. The local test intentionally omits local S/D excitations and shear so that it isolates the shared conventional-fluid density scale.
