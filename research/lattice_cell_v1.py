"""RESEARCH ONLY. The lattice-cell wall (GEAR-107 and AMENDMENT-GEAR-107; HANDOFF to projects Nov 14 2026 (4), Task 1).

The spherical machine holds one cell of a lattice of identical pockets. The outer boundary r = r_max is a symmetry
surface between identical cells. No production file is changed: this module subclasses the grid and the kernel.

What the wall is, in the machine's variables (cell-centred grid, the wall is the last cell FACE):
  W1  Even fields. Every field that is even at the regular centre (a, b, X, K, Aa, alpha, S, PS, D, PD, phi, Pi, and
      every scalar built from them) gets the mirror ghost f(r_max + s) = f(r_max - s); every field that is odd at the
      centre (radial-vector components: Lambda, beta, B) gets f(r_max + s) = -f(r_max - s). All derivative, face-gradient
      and face-interpolation operators then use their ordinary centred stencils up to the wall. No one-sided closure.
  W1b Trace-free tensor variables. Aa, a and b cannot be even at the wall: the momentum constraint there (K uniform,
      zero flux, X even) reads  d_r(r^3 b^(3/2) Aa / X^3) = 0,  so Aa' = -3 Aa / r_max when Aa(wall) != 0, and Aa(wall) is
      fixed by the interior (the r^3-moment of the momentum density inside the cell: GEAR-107 E2). The first version
      of this wall took Aa even; the momentum residual at the wall was then exactly (3/r_max) Aa(wall) (4e-5 at t = 4).
      Since d_t ln a = -2 alpha Aa and d_t ln b = alpha Aa, the same holds for ln a and ln b to first order. Wall
      ghosts: r^3 b^(3/2) Aa / X^3 mirrored, and r^3 ln a, r^3 ln b mirrored (exact for Aa; first order in ln a for
      a, b). This is the image of the centre's regularity rule (there Aa / r^2 is even). Lambda stays odd: with these
      rules Lambda(wall) vanishes to first order in ln a.
  W2  No flux. The fluid flux through the wall face is the HLL flux between the last cell's state and its mirror image
      (v_r -> -v_r): rest-mass and energy fluxes vanish identically, the momentum flux is the wall pressure. The scalar
      fluxes vanish by W1 (zero radial gradient at the wall).
  W3  No boundary data. The vendor's incoming-light reconstruction of Aa at the last cell is NOT applied: Aa there is
      evolved by the same equations as inside (GEAR-107 E2: the wall's geometry is whatever the interior equations
      and constraints give; nothing is prescribed). Shift stays zero (production's CMC branch).
  W4  The clock. The CMC lapse solves  D^2 alpha - Q alpha = -Kdot,  Q = K_ij K^ij + 4 pi (rho + S),  with zero lapse
      gradient at the wall (and at the centre). With both ends Neumann the solution is alpha = Kdot u, where u solves
      D^2 u - Q u = -1: the equation fixes the SHAPE of alpha and leaves one number, Kdot (the rate of the time label).
        closure "ii"  Kdot = <Q>, the proper-volume average of Q over the cell (the cell-averaged acceleration
                      equation). Integrating the lapse equation over the closed cell gives <alpha Q> = Kdot, so this is
                      the normalisation <alpha Q>/<Q> = 1: the Q-weighted mean clock keeps the label.
        closure "i"   alpha = 1 at the wall: Kdot = 1/u(wall). The wall shell's proper time is the label.
      DERIVED consequence: the two closures give the SAME slices (the same u) and differ only by a time-dependent
      relabelling dt_i = alpha_ii(wall) dt_ii. They are one foliation with two clocks; they agree identically when
      alpha_ii(wall) = 1, which is the uniform cell. (GEAR-107 §3 states (ii) as "Kdot fixed by the volume average";
      the reading Kdot = <Q> is this module's, flagged for chat.)
  W5  No sponge, no layer, no friction (E6).

Initial data for a cell (`cell_initial`, `conformal_cell_slice`): fields as research/dm_one_budget_v1.initial, but the
slice is conformally flat (a = b = 1, X = psi^-2, Lambda = 0, uniform K, Aa = 0) instead of areal, so that every metric
function can be even at the wall. (In the areal slice of PR #475 a/b = 1/B at the wall differs from 1 by the cell's
excess mass, 3 m_ex/r_max = 0.16 at r_max = 24, and no choice of K makes all of a, b, X, Lambda consistent with W1:
tried first, Hamiltonian residual 0.13 at the wall, rejected.) K is not taken from a far-field budget: with psi' = 0
at both ends the Hamiltonian constraint is solvable only for  K^2 = 24 pi <rho psi^5>/<psi^5>,  the cell average.
For a uniform cell this is the Friedmann equation with psi = 1. The pocket's profile is laid down in the coordinate
r, which is proper radius at the centre (psi = 1 there; the first version normalised psi at the wall, which made the
pocket's proper size depend on the cell: 1.4% between r_max = 24 and 36, and the ring transient with it).

FORCED / imposed (AUDIT §D): the spherical cell stands in for a space-filling one (Lindquist-Wheeler approximation);
all neighbours identical and in step; zero shift; the pocket and its start are placed by hand.
Global side effect: constructing a LatticeCellKernel replaces engine.matter_rhs._species_rhs_components by the
wall-flux wrapper (`install_wall_flux`); `uninstall_wall_flux` restores it.
"""
import dataclasses
import math
import sys

sys.path.insert(0, ".")

import numpy as np  # noqa: E402
from scipy.linalg import solve_banded  # noqa: E402

import engine.matter_rhs as matter_rhs  # noqa: E402
import engine.matter_system as msys  # noqa: E402
from engine import v55_pirk_adapter as adapter  # noqa: E402
from engine.cmc_gauge import _proper_volume_weights  # noqa: E402
from engine.native_geometry import NativeSphericalCellGrid, native_cmc_explicit_metric_rhs  # noqa: E402
from engine.production_kernel import LAMBDA_M, V55ProductionKernel  # noqa: E402
from engine.v55_matter import total_matter_projection  # noqa: E402
from engine.valencia import FluidPrimitive, spherical_metric_from_bssn  # noqa: E402

grid_ops, _vacuum, _ = adapter.vendor_modules()


import os  # noqa: E402

_REG = {"a": None, "b": None, "Aa": None, "Lambda": None}
OPT = {"ab": os.environ.get("CELL_AB", "exact"), "lam": os.environ.get("CELL_LAM", "geom"),
       "ghost": os.environ.get("CELL_GHOST", "poly"), "lapse": os.environ.get("CELL_LAPSE", "fourth")}


def _poly_weights(g, value_constraint):
    """Ghost weights at s = +h/2, +3h/2 (units of h) from the quartic through the last four cells (s = -1/2, -3/2,
    -5/2, -7/2) with zero slope at the wall (value_constraint False) or a prescribed wall value (True; the fifth
    weight multiplies the wall value)."""
    xs = np.array([-0.5, -1.5, -2.5, -3.5])
    A = np.vstack([np.vstack([xs**p for p in range(5)]).T,
                   np.array([1.0, 0, 0, 0, 0]) if value_constraint else np.array([0, 1.0, 0, 0, 0])])
    inv = np.linalg.inv(A)
    sg = np.array([0.5, 1.5][:g])
    P = np.vstack([sg**p for p in range(5)]).T
    return P @ inv


_PW = {(g, v): _poly_weights(g, v) for g in (1, 2) for v in (False, True)}


def register(geometry):
    """Tell the wall operators which arrays are the trace-free tensor variables of the current stage (W1b)."""
    _REG["a"], _REG["b"], _REG["Aa"] = (np.asarray(geometry.a, float), np.asarray(geometry.b, float),
                                        np.asarray(geometry.Aa, float))
    _REG["Lambda"] = np.asarray(geometry.Lambda, float)


def _kind(u):
    for name in ("Aa", "a", "b"):
        t = _REG[name]
        if t is not None and (u is t or (u.shape == t.shape and np.array_equal(u, t))):
            return name
    return None


class _WallOperators:
    """Centred stencils with ghosts at both ends (replaces the one-sided outer closures)."""

    def _wall_ghosts(self, f, g, wall_value=None):
        """g ghost values of f beyond the wall. "mirror": f(R + s) = f(R - s) (or the odd image about wall_value).
        "poly": the quartic through the last four cells with zero slope at the wall (or f(R) = wall_value). The mirror
        assumes every odd derivative vanishes at the wall, which a spherical wall does not give (the lapse there has
        third derivative -(2/R) times its second derivative, plus source terms); the quartic imposes only the one
        condition the wall supplies."""
        last = f[:-5:-1]                                   # f(-1/2), f(-3/2), f(-5/2), f(-7/2)
        if OPT["ghost"] == "mirror":
            return f[:-g - 1:-1] if wall_value is None else 2.0 * wall_value - f[:-g - 1:-1]
        W = _PW[(g, wall_value is not None)]
        return W[:, :4] @ last + (W[:, 4] * wall_value if wall_value is not None else 0.0)

    def _ext(self, values, parity, g):
        """Extend by g ghosts at the centre (parity mirror) and at the wall (W1, or W1b for a, b, Aa, Lambda)."""
        u = np.asarray(values, dtype=float)
        left = parity * u[g - 1::-1] if g > 1 else parity * u[:1]
        R, h = self.r_max, self.dr
        kind = _kind(u) if parity == 1 else None
        tl = _REG["Lambda"]
        if parity == 1 and kind is None:
            right = self._wall_ghosts(u, g)
        elif parity == -1:
            lw = 0.0
            if (OPT["lam"] == "geom" and tl is not None and np.any(u != 0.0)
                    and (u is tl or (u.shape == tl.shape and np.array_equal(u, tl)))):
                aw = 1.5 * _REG["a"][-1] - 0.5 * _REG["a"][-2]
                lw = 2.0 / (aw * R) * (aw**0.75 - 1.0) ** 2   # Lambda's wall value, second order in ln a
            right = self._wall_ghosts(u, g, wall_value=lw)
        else:
            rc = R - (np.arange(len(u))[::-1] + 0.5) * h    # cell radii, grid order
            rg = R + (np.arange(g) + 0.5) * h
            if kind in ("a", "b"):
                pw = -0.75 if kind == "a" else 1.5          # r^3 (1 - a^(-3/4)) = r^3 (1 - b^(3/2)) has zero slope
                if OPT["ab"] == "exact":
                    right = (1.0 - self._wall_ghosts(rc**3 * (1.0 - u**pw), g) / rg**3) ** (1.0 / pw)
                else:
                    right = np.exp(self._wall_ghosts(rc**3 * np.log(u), g) / rg**3)
            else:
                b = _REG["b"]
                bg = (1.0 - self._wall_ghosts(rc**3 * (1.0 - b**1.5), g) / rg**3) ** (1.0 / 1.5)
                right = self._wall_ghosts(rc**3 * b**1.5 * u, g) / (rg**3 * bg**1.5)   # X even: its factor cancels
        return np.concatenate([left, u, right])

    def cell_derivative(self, values, *, parity):
        e = self._ext(values, parity, 1)
        return (e[2:] - e[:-2]) / (2.0 * self.dr)

    def cell_second_derivative(self, values, *, parity):
        e = self._ext(values, parity, 1)
        return (e[2:] - 2.0 * e[1:-1] + e[:-2]) / self.dr**2

    def _cell_derivative_fourth(self, values, *, parity, order):
        e = self._ext(values, parity, 2)
        if order == 1:
            return (e[:-4] - 8.0 * e[1:-3] + 8.0 * e[3:-1] - e[4:]) / (12.0 * self.dr)
        return (-e[:-4] + 16.0 * e[1:-3] - 30.0 * e[2:-2] + 16.0 * e[3:-1] - e[4:]) / (12.0 * self.dr**2)

    def cell_derivative_fourth(self, values, *, parity):
        return self._cell_derivative_fourth(values, parity=parity, order=1)

    def cell_second_derivative_fourth(self, values, *, parity):
        return self._cell_derivative_fourth(values, parity=parity, order=2)

    def gradient_faces_fourth(self, values, *, parity, inner_flux_gradient=0.0, outer_flux_gradient=0.0):
        e = self._ext(values, parity, 2)                 # e[i + 2] = u[i]
        res = np.empty(self.n + 1)
        res[1:] = (e[1:self.n + 1] - 27.0 * e[2:self.n + 2] + 27.0 * e[3:self.n + 3] - e[4:self.n + 4]) / (24.0 * self.dr)
        res[0] = inner_flux_gradient
        if outer_flux_gradient != 0.0:
            res[-1] = outer_flux_gradient
        return res

    def interpolate_faces_fourth(self, values, *, parity):
        e = self._ext(values, parity, 2)
        res = np.empty(self.n + 1)
        res[:] = (-e[0:self.n + 1] + 9.0 * e[1:self.n + 2] + 9.0 * e[2:self.n + 3] - e[3:self.n + 4]) / 16.0
        return res


class LatticeCellGrid(_WallOperators, grid_ops.SphericalCellGrid):
    """The vendor grid with the lattice-cell wall at r_max."""


class LatticeNativeGrid(_WallOperators, NativeSphericalCellGrid):
    """The native-geometry grid with the lattice-cell wall at r_max."""


# ------------------------------------------------------------------------------------------------ W2: wall flux
_ORIG_COMPONENTS = msys._species_rhs_components


def _wall_components(metric, metric_derivatives, state, species, **kw):
    out = _ORIG_COMPONENTS(metric, metric_derivatives, state, species, **kw)
    face_flux = out[0]
    n = len(metric.r)
    dr = float(metric.r[1] - metric.r[0])
    wall = spherical_metric_from_bssn(float(metric.r[-1]) + 0.5 * dr, float(metric.a[-1]), float(metric.b[-1]),
                                      float(metric.X[-1]), float(metric.alpha[-1]), 0.0)
    last = msys._metric_arrays(type(metric)(*(np.asarray(getattr(metric, f.name))[-1:] for f in dataclasses.fields(metric))))
    q = msys.primitives(last, msys.ConservedSpecies(state.rest[-1:], state.energy_t[-1:], state.momentum_r[-1:]), species)[0]
    # face states: the last cell's value extrapolated to the wall with its mirror-limited slope (v_r -> 0), and
    # the mirror image. Mass and energy fluxes vanish identically; the momentum flux is the wall pressure.
    ql = FluidPrimitive(rho=max(0.0, q.rho), pressure=max(0.0, q.pressure), v_r=0.0, gamma_rr=wall.gamma_rr)
    qr = ql
    face_flux[n] = msys._hll_flux(wall, ql, qr, species)
    return out


_ORIG_RECON = msys._reconstructed_primitive


def _wall_reconstructed_primitive(prim, i, side, gamma_rr):
    """As the production reconstruction, but the last cell uses the wall's mirror ghost (rho, p even; v_r odd)
    for its limited slope instead of falling back to piecewise constant. (Without this the last cell's face value of
    v_r is its centre value, half the true face velocity for a flow that vanishes at the wall, and the dust density
    develops a +-0.5% two-cell sawtooth at the wall: found at t = 8, r_max = 12.)"""
    n = len(prim)
    if i != n - 1 or n < 2:
        return _ORIG_RECON(prim, i, side, gamma_rr)
    attrs = {}
    for name, par in (("rho", 1.0), ("pressure", 1.0), ("v_r", -1.0)):
        c, lft = float(getattr(prim[i], name)), float(getattr(prim[i - 1], name))
        slope = msys._minmod(c - lft, par * c - c)
        attrs[name] = c - 0.5 * slope if side == "left" else c + 0.5 * slope
    vmax = (1.0 - 1.0e-12) / math.sqrt(gamma_rr)
    return FluidPrimitive(rho=max(0.0, attrs["rho"]), pressure=max(0.0, attrs["pressure"]),
                          v_r=max(-vmax, min(vmax, attrs["v_r"])), gamma_rr=gamma_rr)


def install_wall_flux():
    matter_rhs._species_rhs_components = _wall_components
    msys._reconstructed_primitive = _wall_reconstructed_primitive


def uninstall_wall_flux():
    matter_rhs._species_rhs_components = _ORIG_COMPONENTS
    msys._reconstructed_primitive = _ORIG_RECON


# ------------------------------------------------------------------------------------------------ W1b: registration
_VENDOR_ENTRY = ("geometry_terms", "primary_l2_rhs", "primary_l3_rhs", "momentum_constraint", "lambda_l2_rhs",
                 "lambda_l3_rhs", "constraints")
_VENDOR_ORIG = {}


def install_registration():
    """Wrap the vendor entry points so each registers the geometry it is evaluated on (only for lattice grids)."""
    for name in _VENDOR_ENTRY:
        if name in _VENDOR_ORIG:
            continue
        orig = getattr(_vacuum, name)
        _VENDOR_ORIG[name] = orig

        def wrapper(grid, state, *args, _orig=orig, **kw):
            if isinstance(grid, _WallOperators):
                register(state)
            return _orig(grid, state, *args, **kw)
        setattr(_vacuum, name, wrapper)


def uninstall_registration():
    for name, orig in _VENDOR_ORIG.items():
        setattr(_vacuum, name, orig)
    _VENDOR_ORIG.clear()


# ------------------------------------------------------------------------------------------------ W4: the clock
def solve_cell_lapse(grid, geometry, scalars, matter, closure="ii"):
    """CMC lapse with zero gradient at the centre and at the wall. Returns (alpha, Kdot)."""
    register(geometry)
    weights, total_weight = _proper_volume_weights(grid, geometry)
    total = total_matter_projection(grid, geometry, scalars, matter)
    a, b, X = np.asarray(geometry.a), np.asarray(geometry.b), np.asarray(geometry.X)
    r = np.asarray(grid.centers)
    inv = X * X / a
    d1 = lambda f: grid.cell_derivative_fourth(f, parity=1)  # noqa: E731
    c = -0.5 * d1(a) / a + d1(b) / b - d1(X) / X + 2.0 / r
    Q = (1.5 * geometry.Aa * geometry.Aa + geometry.K * geometry.K / 3.0
         + 4.0 * math.pi * (total["rho"] + total["pr"] + 2.0 * total["pt"]))
    m = Q / inv
    n, h = grid.n, grid.dr
    band = np.zeros((3, n))
    ci = c / (2.0 * h)
    band[1, :] = -2.0 / h**2 - m                       # diagonal
    band[0, 1:] = (1.0 / h**2 + ci)[:-1]               # super-diagonal
    band[2, :-1] = (1.0 / h**2 - ci)[1:]               # sub-diagonal
    band[1, 0] += 1.0 / h**2 - ci[0]                   # even ghost at the centre
    band[1, -1] += 1.0 / h**2 + ci[-1]                 # even ghost at the wall (zero lapse gradient)
    u = solve_banded((1, 1), band, -1.0 / inv, check_finite=False)
    if OPT["lapse"] == "fourth":
        # Deferred correction to the operator the K equation itself uses (fourth-order stencils with the wall's
        # ghosts), so that the K right-hand side equals Kdot at every cell including the last one.
        d2 = lambda f: grid.cell_second_derivative_fourth(f, parity=1)  # noqa: E731
        for _ in range(4):
            res = d2(u) + c * d1(u) - m * u + 1.0 / inv
            u = u - solve_banded((1, 1), band, res, check_finite=False)
    if not np.all(np.isfinite(u)) or float(np.min(u)) <= 0.0:
        raise ValueError(f"cell lapse solve has no positive solution: min u = {float(np.min(u)):.6e}")
    if closure == "ii":
        kdot = float(np.sum(weights * Q) / total_weight)
    elif closure == "i":
        kdot = 1.0 / (0.5 * (u[-1] + u[-1]))           # even ghost: the wall-face value is the last cell's
    else:
        raise ValueError("closure must be 'i' or 'ii'")
    return kdot * u, kdot


class LatticeCellKernel(V55ProductionKernel):
    """Production kernel with the lattice-cell wall (W1-W5). Research setting; nothing in engine/ is modified."""

    def __init__(self, closure="ii"):
        self.closure = closure
        install_wall_flux()
        install_registration()

    def _rhs(self, state):
        register(state.geometry)
        return super()._rhs(state)

    def _solve_lapse(self, grid, geometry, scalars, matter, *, electromagnetic_projection=None):
        if electromagnetic_projection is not None:
            raise NotImplementedError("lattice cell: no electromagnetic sector")
        return solve_cell_lapse(grid, geometry, scalars, matter, self.closure)

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        register(geometry)
        native_grid = LatticeNativeGrid(grid.n, grid.r_max)
        explicit_native = native_cmc_explicit_metric_rhs(native_grid, geometry)
        return adapter.geometry_stage_terms(grid, geometry, scalars, matter, lambda_m=LAMBDA_M,
                                            explicit_rhs_override=explicit_native)

    @staticmethod
    def _apply_outer_light_boundary(grid, geometry, scalars, matter, *, electromagnetic_projection=None):
        return None                                    # W3: nothing is prescribed at the wall


def as_cell(state):
    """The same state on the lattice-cell grid."""
    g = state.grid
    return dataclasses.replace(state, grid=LatticeCellGrid(g.n, g.r_max))


# ------------------------------------------------------------------------------------------------ initial data
def conformal_cell_slice(r, dr, grad2, Vm, tol=1e-13, itmax=60):
    """Conformally flat CMC initial slice in a cell: gamma = psi^4 (dr^2 + r^2 dOmega^2), K uniform, Aa = 0.
    Hamiltonian constraint  psi'' + (2/r) psi' = (K^2/12 - 2 pi rho) psi^5,  rho = grad2 / psi^4 + Vm,
    with psi' = 0 at the centre and at the wall (W1) and psi = 1 at the centre (r is proper radius there, so the
    same coordinate profile is the same physical pocket in every cell: psi(r) - psi(0) depends only on what lies
    inside r and on K). Integrating over the
    cell with both ends Neumann gives the solvability condition  K^2 = 24 pi <rho psi^5>/<psi^5>: K is fixed by the
    cell average, not chosen. Second-order finite differences, Newton iteration on (psi, K^2)."""
    n = len(r)
    psi = np.ones(n); K2 = 24 * math.pi * float(np.sum(Vm * r * r) / np.sum(r * r))
    lo = 1.0 / dr**2 - 1.0 / (r * dr); hi = 1.0 / dr**2 + 1.0 / (r * dr)
    for _ in range(itmax):
        rho = grad2 / psi**4 + Vm
        F = np.empty(n + 1)
        lap = -2.0 / dr**2 * psi
        lap[1:] += lo[1:] * psi[:-1]; lap[:-1] += hi[:-1] * psi[1:]
        lap[0] += lo[0] * psi[0]; lap[-1] += hi[-1] * psi[-1]
        F[:n] = lap - (K2 / 12 - 2 * math.pi * rho) * psi**5
        F[n] = psi[0] - 1.0
        J = np.zeros((n + 1, n + 1))
        idx = np.arange(n)
        # d/dpsi of -(K2/12 - 2 pi (grad2 psi^-4 + Vm)) psi^5 = -(5 K2/12) psi^4 + 2 pi (grad2 + 5 Vm psi^4)
        J[idx, idx] = -2.0 / dr**2 - 5 * K2 / 12 * psi**4 + 2 * math.pi * (grad2 + 5 * Vm * psi**4)
        J[idx[1:], idx[:-1]] = lo[1:]; J[idx[:-1], idx[1:]] = hi[:-1]
        J[0, 0] += lo[0]; J[n - 1, n - 1] += hi[-1]
        J[:n, n] = -psi**5 / 12
        J[n, 0] = 1.0
        dx = np.linalg.solve(J, -F)
        psi = psi + dx[:n]; K2 = K2 + dx[n]
        if np.max(np.abs(dx[:n])) < tol and abs(dx[n]) < tol * K2:
            break
    if K2 <= 0 or np.any(psi <= 0) or not np.all(np.isfinite(psi)):
        raise FloatingPointError("no expanding conformally flat cell slice")
    return psi, K2


def cell_initial(case, x_c, beta_D, N, RMAX, closure="ii", start="cold", R_pocket=None):
    """One-budget initial data in a lattice cell. case: full | noD | nomatter; start: cold | hot (D released)."""
    import research.dm_one_budget_v1 as dm
    import research.gear91_coupling_v1 as g91
    from engine import scalar_system as ss, v55_initial
    from engine.matter_system import _metric_arrays, initialize_dust, initialize_radiation
    from engine.v55_matter import V55MatterState, metric_slice_from_q
    if R_pocket is not None:
        dm.R_POCKET = R_pocket
    _, st, r, _ = dm.initial(case, x_c, beta_D, N, RMAX)
    st = as_cell(st)
    grid = st.grid; dr = float(grid.dr); g = st.geometry
    S, D, phi, Pi = (np.array(getattr(st.scalars, k), dtype=float) for k in ("S", "D", "phi", "Pi"))
    LAM = dm.LAM
    if start == "hot" and case == "full":
        x = LAM * S * S
        eD = math.exp(-2 * beta_D * dm.PHI_TODAY)
        D = np.sqrt(np.clip((2 / LAM) * (eD - 0.75 * x), 0, None)) + 1e-10 * np.exp(-(r / (0.5 * dm.R_POCKET)) ** 2)
    dm0, b0, rr0 = dm.C2 * dm.vm.RHO_DM_PRESENT, dm.C2 * dm.vm.RHO_B_PRESENT, dm.C2 * dm.vm.RHO_R_PRESENT
    d1 = lambda f: grid.cell_derivative_fourth(f, parity=1)  # noqa: E731
    Sp, Dp = d1(S), d1(D)
    Vm = (g91.potential(S, D, phi, 0.0, beta_D, 1.0, LAM)
          + (0.5 * Pi**2 + ss.cosmos_potential(phi)) / ss.KAPPA + (dm0 + b0 + rr0) / ss.KAPPA)
    psi, K2 = conformal_cell_slice(r, dr, 0.5 * (Sp**2 + Dp**2), Vm)
    H0 = math.sqrt(K2) / 3.0
    g.a[:] = 1.0; g.b[:] = 1.0; g.X[:] = psi ** (-2.0)
    g.K[:] = -math.sqrt(K2); g.Aa[:] = 0; g.beta[:] = 0; g.B[:] = 0; g.Lambda[:] = 0
    mets = _metric_arrays(metric_slice_from_q(grid, g))
    z = np.zeros_like(r)
    st = dataclasses.replace(
        st, geometry=g,
        scalars=dataclasses.replace(st.scalars, S=S, PS=z.copy(), D=D, PD=z.copy(), phi=phi, Pi=Pi),
        matter=V55MatterState(initialize_dust(mets, np.full_like(r, dm0)), initialize_dust(mets, np.full_like(r, b0)),
                              initialize_radiation(mets, np.full_like(r, rr0)), None))
    k = LatticeCellKernel(closure)
    st.geometry.alpha[:] = k._solve_lapse_for_state(grid, st.geometry, st.scalars, st.matter, None)[0]
    return k, st, r, H0
