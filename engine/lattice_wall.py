"""The zero-flux (lattice-wall) outer closure for EXPANDING runs.

Label: "identical neighbours, in step: one-stone simplification" (HANDOFF Nov 15 2026 (3), decision 2,
confirmed by Jeremiah). The spherical machine holds one cell of a lattice of identical cells that all evolve in step;
the outer boundary r = r_max is the symmetry surface between neighbours. Nothing crosses it and nothing is
prescribed there. The wall REFLECTS outgoing waves: correct for identical neighbours (each neighbour sends its own
ringing in), WRONG for an isolated set-up. Static and isolated runs keep production's open edge
(engine.production_kernel, the pinned incoming-light reconstruction and the projected-K-dot lapse).

Promoted from research/lattice_cell_v1.py (W1-W5 there; results/zero-flux-closure-2026-11-16). The research version
installed its wall flux and its geometry registration module-wide; here nothing global is touched. The closure is a
property of the GRID the state carries (``LatticeWallGrid``), and the kernel threads it explicitly:
  - operators: ``bind(grid, geometry)`` returns a copy of the wall grid that knows the geometry it is evaluating
    (for the W1b ghosts of a, b, Aa and Lambda); the open-edge grid is returned unchanged;
  - fluid flux: ``wall_species_rhs`` (W2) replaces ``engine.matter_rhs.species_rhs`` for wall grids only;
  - lapse: ``solve_wall_lapse`` (W4) replaces the open edge's CMC solve for wall grids only;
  - boundary data: none (W3), the incoming-light reconstruction is skipped for wall grids only.

What the wall is (cell-centred grid; the wall is the last cell FACE):
  W1  Even fields get a ghost with zero slope at the wall (the quartic through the last four cells with f'(R) = 0);
      odd fields (Lambda, beta, B) the quartic with a prescribed wall value (0, or Lambda's geometric value below).
      All derivative, face-gradient and face-interpolation operators use their centred stencils up to the wall.
  W1b Trace-free tensor variables. The momentum constraint at the wall (K uniform, zero flux, X even) reads
      d_r(r^3 b^(3/2) Aa / X^3) = 0, so Aa is not even there. Ghosts: r^3 b^(3/2) Aa (X even, its factor cancels),
      r^3 (1 - a^(-3/4)) and r^3 (1 - b^(3/2)) have zero slope at the wall. Lambda's wall value is
      2 (a^(3/4) - 1)^2 / (a R), second order in ln a. The image of the centre's regularity rule.
  W2  No flux. The fluid flux through the wall face is the HLL flux between the last cell's state at rest and its
      mirror image: rest-mass and energy fluxes vanish identically; the momentum flux is the wall pressure. The last
      cell reconstructs its face value of v_r with the mirror ghost (v_r odd). Scalar fluxes vanish by W1.
  W3  No boundary data. Aa in the last cell is evolved by the interior equations. Shift stays zero (CMC branch).
  W4  The clock. The CMC lapse solves D^2 alpha - Q alpha = -Kdot with zero lapse gradient at the centre AND the
      wall, Q = K_ij K^ij + 4 pi (rho + S). Its solution is alpha = Kdot u with D^2 u - Q u = -1: the equation fixes
      the lapse's SHAPE; Kdot is one number fixed by the normalisation, never by the stored lapse. So no step's lapse
      feeds the next and the open edge's undamped lapse mode (results/floor-swap-2026-11-15/KERNEL_DRIFT.md) cannot
      exist. Clocks (one foliation, two labels; Derived: the same u, dt_i = alpha_ii(wall) dt_ii):
        "i"   alpha = 1 at the wall: the time label is the proper time of the wall shell. PRODUCTION DEFAULT.
        "ii"  Kdot = <Q>, the proper-volume mean of Q over the cell (GEAR-107 sec. 3 read as Kdot = <Q>; flagged for chat).
  W5  No sponge, no layer, no friction.

Initial data. Production's background start is exactly flat (psi = 1) and needs no change. A non-flat start (the
labelled pocket fixtures) needs a conformally flat slice so every metric function can obey W1 at the wall:
``conformal_wall_slice`` (K^2 = 24 pi <rho psi^5>/<psi^5> is fixed by the cell average, not chosen).

FORCED / imposed: the spherical cell stands in for a space-filling one (Lindquist-Wheeler); all neighbours identical
and in step; zero shift. Unpaid: the one-stone simplification itself; clock "i" (the wall shell's proper time).
"""
from __future__ import annotations

import contextlib
from dataclasses import dataclass, field, replace
import math

import numpy as np
from scipy.linalg import solve_banded

from .cmc_gauge import _proper_volume_weights
from .matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    Species,
    _hll_flux,
    _metric_arrays,
    _minmod,
    _reconstructed_primitive,
    _species_rhs_components,
    primitives,
)
from .native_geometry import NativeSphericalCellGrid
from .v55_matter import total_matter_projection
from .valencia import FluidPrimitive, spherical_metric_from_bssn

WALL_LABEL = "identical neighbours, in step: one-stone simplification"
OPEN_EDGE = "open_edge"
LATTICE_WALL = "lattice_wall"
OUTER_EDGES = (OPEN_EDGE, LATTICE_WALL)
WALL_CLOCKS = ("i", "ii")
DEFAULT_WALL_CLOCK = "i"


# The edge production's default start (today, background alone: an expanding run) gets when none is asked for.
_EXPANDING_DEFAULT = {"edge": LATTICE_WALL}


def expanding_default_edge() -> str:
    return _EXPANDING_DEFAULT["edge"]


@contextlib.contextmanager
def open_edge_reference():
    """LABELLED REFERENCE: inside this block production's default start gets the open edge, as it did before the
    wall was added (bit-identical to that code). For contracts written against the open edge's implementation seams
    (tests/conftest.py::OPEN_EDGE_REFERENCE_MODULES); not for expanding production runs."""
    saved = _EXPANDING_DEFAULT["edge"]
    _EXPANDING_DEFAULT["edge"] = OPEN_EDGE
    try:
        yield
    finally:
        _EXPANDING_DEFAULT["edge"] = saved


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


@dataclass(frozen=True)
class WallTensorFields:
    """The geometry a wall grid is evaluating: the trace-free tensor variables that take the W1b ghosts."""

    a: np.ndarray
    b: np.ndarray
    Aa: np.ndarray
    Lambda: np.ndarray

    @classmethod
    def of(cls, geometry):
        return cls(np.asarray(geometry.a, dtype=float), np.asarray(geometry.b, dtype=float),
                   np.asarray(geometry.Aa, dtype=float), np.asarray(geometry.Lambda, dtype=float))


class _WallOperators:
    """Centred stencils with ghosts at both ends (replaces the open edge's one-sided outer closures).

    The W1b ghosts need to know which array is a, b, Aa or Lambda of the geometry being evaluated; the kernel binds
    the grid to that geometry (``bind``). An unbound wall grid treats every field as W1 (even or odd)."""

    def _kind(self, u):
        t = self.tensor_fields
        if t is None:
            return None
        for name in ("Aa", "a", "b"):
            ref = getattr(t, name)
            if u is ref or (u.shape == ref.shape and np.array_equal(u, ref)):
                return name
        return None

    def _wall_ghosts(self, f, g, wall_value=None):
        """g ghost values of f beyond the wall: the quartic through the last four cells with zero slope at the wall
        (or f(R) = wall_value). It imposes only the one condition the wall supplies; a plain mirror would also force
        every odd derivative to vanish there, which a spherical wall does not give."""
        last = f[:-5:-1]                                   # f(-1/2), f(-3/2), f(-5/2), f(-7/2)
        W = _PW[(g, wall_value is not None)]
        return W[:, :4] @ last + (W[:, 4] * wall_value if wall_value is not None else 0.0)

    def _ext(self, values, parity, g):
        """Extend by g ghosts at the centre (parity mirror) and at the wall (W1, or W1b for a, b, Aa, Lambda)."""
        u = np.asarray(values, dtype=float)
        left = parity * u[g - 1::-1] if g > 1 else parity * u[:1]
        R, h = self.r_max, self.dr
        kind = self._kind(u) if parity == 1 else None
        t = self.tensor_fields
        tl = None if t is None else t.Lambda
        if parity == 1 and kind is None:
            right = self._wall_ghosts(u, g)
        elif parity == -1:
            lw = 0.0
            if (tl is not None and np.any(u != 0.0)
                    and (u is tl or (u.shape == tl.shape and np.array_equal(u, tl)))):
                aw = 1.5 * t.a[-1] - 0.5 * t.a[-2]
                lw = 2.0 / (aw * R) * (aw**0.75 - 1.0) ** 2   # Lambda's wall value, second order in ln a
            right = self._wall_ghosts(u, g, wall_value=lw)
        else:
            rc = R - (np.arange(len(u))[::-1] + 0.5) * h    # cell radii, grid order
            rg = R + (np.arange(g) + 0.5) * h
            if kind in ("a", "b"):
                pw = -0.75 if kind == "a" else 1.5          # r^3 (1 - a^(-3/4)) = r^3 (1 - b^(3/2)) has zero slope
                right = (1.0 - self._wall_ghosts(rc**3 * (1.0 - u**pw), g) / rg**3) ** (1.0 / pw)
            else:
                b = t.b
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


def _check_clock(clock):
    if clock not in WALL_CLOCKS:
        raise ValueError(f"wall clock must be one of {WALL_CLOCKS}, not {clock!r}")


@dataclass(frozen=True)
class LatticeNativeWallGrid(_WallOperators, NativeSphericalCellGrid):
    """The native-geometry grid with the lattice wall at r_max (for the native CMC a/b/X right-hand side)."""

    clock: str = DEFAULT_WALL_CLOCK
    tensor_fields: WallTensorFields | None = field(default=None, compare=False, repr=False)


_VENDOR_WALL_GRID = {}


def _vendor_wall_grid_class():
    """The vendor grid with the lattice wall (built on first use, so importing this module needs no vendor)."""
    cls = _VENDOR_WALL_GRID.get("cls")
    if cls is None:
        from . import v55_pirk_adapter as adapter

        grid_ops, _, _ = adapter.vendor_modules()

        @dataclass(frozen=True)
        class LatticeWallGrid(_WallOperators, grid_ops.SphericalCellGrid):
            """The production (vendor-operator) grid with the lattice wall at r_max. Expanding runs only;
            label: "identical neighbours, in step: one-stone simplification"."""

            clock: str = DEFAULT_WALL_CLOCK
            tensor_fields: WallTensorFields | None = field(default=None, compare=False, repr=False)

            def __post_init__(self):
                super().__post_init__()
                _check_clock(self.clock)

        LatticeWallGrid.__module__ = __name__
        LatticeWallGrid.__qualname__ = "LatticeWallGrid"     # picklable as engine.lattice_wall.LatticeWallGrid
        cls = _VENDOR_WALL_GRID["cls"] = LatticeWallGrid
    return cls


def __getattr__(name):
    if name == "LatticeWallGrid":
        return _vendor_wall_grid_class()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def lattice_wall_grid(n, r_max, clock=DEFAULT_WALL_CLOCK):
    """A production grid whose outer boundary is the lattice wall (clock "i" by default)."""
    _check_clock(clock)
    return _vendor_wall_grid_class()(n, r_max, clock=clock)


def as_lattice_wall(grid, clock=DEFAULT_WALL_CLOCK):
    """The same cells with the lattice wall at r_max."""
    return lattice_wall_grid(grid.n, grid.r_max, clock)


def is_lattice_wall(grid) -> bool:
    return isinstance(grid, _WallOperators)


def outer_edge_of(grid) -> str:
    return LATTICE_WALL if is_lattice_wall(grid) else OPEN_EDGE


def bind(grid, geometry):
    """The grid that evaluates operators on ``geometry``: a wall grid bound to its tensor variables (W1b), or the
    open-edge grid itself, unchanged."""
    if not is_lattice_wall(grid):
        return grid
    return replace(grid, tensor_fields=WallTensorFields.of(geometry))


def native_wall_grid(grid):
    """The native-geometry twin of a (bound) wall grid."""
    return LatticeNativeWallGrid(grid.n, grid.r_max, clock=grid.clock, tensor_fields=grid.tensor_fields)


# ------------------------------------------------------------------------------------------------ W2: no flux
def _wall_reconstructed_primitive(prim, i, side, gamma_rr):
    """As the production reconstruction, but the last cell uses the wall's mirror ghost (rho, p even; v_r odd) for
    its limited slope instead of falling back to piecewise constant. Without it the last cell's face value of v_r is
    its centre value, half the true face velocity of a flow that vanishes at the wall, and the dust density develops a
    +-0.5% two-cell sawtooth at the wall (research: t = 8, r_max = 12)."""
    n = len(prim)
    if i != n - 1 or n < 2:
        return _reconstructed_primitive(prim, i, side, gamma_rr)
    attrs = {}
    for name, par in (("rho", 1.0), ("pressure", 1.0), ("v_r", -1.0)):
        c, lft = float(getattr(prim[i], name)), float(getattr(prim[i - 1], name))
        slope = _minmod(c - lft, par * c - c)
        attrs[name] = c - 0.5 * slope if side == "left" else c + 0.5 * slope
    vmax = (1.0 - 1.0e-12) / math.sqrt(gamma_rr)
    return FluidPrimitive(rho=max(0.0, attrs["rho"]), pressure=max(0.0, attrs["pressure"]),
                          v_r=max(-vmax, min(vmax, attrs["v_r"])), gamma_rr=gamma_rr)


def _tail(metric: BSSNMetricSlice, k: int) -> BSSNMetricSlice:
    return BSSNMetricSlice(*(np.asarray(getattr(metric, name))[-k:] for name in BSSNMetricSlice.__dataclass_fields__))


def wall_face_fluxes(metric, state, species, face_flux):
    """Overwrite the last two face fluxes of the open edge's flux array with the wall's (W2), in place.

    Face n - 1 (between the last two cells): the same HLL flux, with the last cell reconstructed against its mirror.
    Face n (the wall): HLL between the last cell's state at rest and its mirror image."""
    n = len(metric.r)
    if n < 4:
        raise ValueError("the lattice wall needs at least four cells")
    tail = _tail(metric, 3)
    prim = primitives(_metric_arrays(tail), ConservedSpecies(np.asarray(state.rest)[-3:],
                                                             np.asarray(state.energy_t)[-3:],
                                                             np.asarray(state.momentum_r)[-3:]), species)
    i = n - 2
    mf = spherical_metric_from_bssn(
        0.5 * (float(metric.r[i]) + float(metric.r[i + 1])),
        0.5 * (float(metric.a[i]) + float(metric.a[i + 1])),
        0.5 * (float(metric.b[i]) + float(metric.b[i + 1])),
        0.5 * (float(metric.X[i]) + float(metric.X[i + 1])),
        0.5 * (float(metric.alpha[i]) + float(metric.alpha[i + 1])),
        0.5 * (float(metric.beta[i]) + float(metric.beta[i + 1])),
    )
    ql = _reconstructed_primitive(prim, 1, "right", mf.gamma_rr)
    qr = _wall_reconstructed_primitive(prim, 2, "left", mf.gamma_rr)
    face_flux[n - 1] = _hll_flux(mf, ql, qr, species)

    dr = float(metric.r[1] - metric.r[0])
    wall = spherical_metric_from_bssn(float(metric.r[-1]) + 0.5 * dr, float(metric.a[-1]), float(metric.b[-1]),
                                      float(metric.X[-1]), float(metric.alpha[-1]), 0.0)
    q = prim[2]
    at_rest = FluidPrimitive(rho=max(0.0, q.rho), pressure=max(0.0, q.pressure), v_r=0.0, gamma_rr=wall.gamma_rr)
    face_flux[n] = _hll_flux(wall, at_rest, at_rest, species)
    return face_flux


def wall_species_rhs(metric, metric_derivatives, state, species, dphi_t=None, dphi_r=None, beta_dm=-0.04):
    """engine.matter_rhs.species_rhs with the wall's face fluxes (W2). Same sources; same arithmetic otherwise."""
    face_flux, source_e, source_s, q_e, q_s, q_rest = _species_rhs_components(
        metric, metric_derivatives, state, species, dphi_t=dphi_t, dphi_r=dphi_r, beta_dm=beta_dm,
    )
    wall_face_fluxes(metric, state, species, face_flux)
    inv_dr = 1.0 / (metric.r[1] - metric.r[0])
    return ConservedSpecies(
        rest=-inv_dr * (face_flux[1:, 0] - face_flux[:-1, 0]) + q_rest,
        energy_t=(-inv_dr * (face_flux[1:, 1] - face_flux[:-1, 1]) + (source_e + q_e)),
        momentum_r=(-inv_dr * (face_flux[1:, 2] - face_flux[:-1, 2]) + (source_s + q_s)),
    )


# ------------------------------------------------------------------------------------------------ W4: the clock
def wall_lapse_shape(grid, geometry, scalars, matter):
    """u with D^2 u - Q u = -1, zero gradient at the centre and at the wall; returns (u, Q, proper-volume weights,
    total weight). ``grid`` must be a wall grid bound to ``geometry``."""
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
    # Deferred correction to the operator the K equation itself uses (fourth-order stencils with the wall's ghosts),
    # so that the K right-hand side equals Kdot at every cell including the last one.
    d2 = lambda f: grid.cell_second_derivative_fourth(f, parity=1)  # noqa: E731
    for _ in range(4):
        res = d2(u) + c * d1(u) - m * u + 1.0 / inv
        u = u - solve_banded((1, 1), band, res, check_finite=False)
    if not np.all(np.isfinite(u)) or float(np.min(u)) <= 0.0:
        raise ValueError(f"lattice-wall lapse solve has no positive solution: min u = {float(np.min(u)):.6e}")
    return u, Q, weights, total_weight


def solve_wall_lapse(grid, geometry, scalars, matter):
    """CMC lapse with zero gradient at the centre and at the wall (W4). Returns (alpha, Kdot); Kdot from the grid's
    clock ("i": alpha = 1 at the wall; "ii": Kdot = <Q>). The lapse already in ``geometry`` is never read."""
    if not is_lattice_wall(grid):
        raise TypeError("solve_wall_lapse needs a lattice-wall grid")
    grid = bind(grid, geometry)
    u, Q, weights, total_weight = wall_lapse_shape(grid, geometry, scalars, matter)
    if grid.clock == "ii":
        kdot = float(np.sum(weights * Q) / total_weight)
    else:
        kdot = 1.0 / (0.5 * (u[-1] + u[-1]))           # even ghost: the wall-face value is the last cell's
    return kdot * u, kdot


# ------------------------------------------------------------------------------------------------ initial data
def conformal_wall_slice(r, dr, grad2, Vm, tol=1e-13, itmax=60):
    """Conformally flat CMC initial slice in a cell: gamma = psi^4 (dr^2 + r^2 dOmega^2), K uniform, Aa = 0.

    Hamiltonian constraint  psi'' + (2/r) psi' = (K^2/12 - 2 pi rho) psi^5,  rho = grad2 / psi^4 + Vm,
    with psi' = 0 at the centre and at the wall (W1) and psi = 1 at the centre (r is proper radius there, so the same
    coordinate profile is the same physical pocket in every cell). Integrating over the cell with both ends Neumann
    gives the solvability condition  K^2 = 24 pi <rho psi^5>/<psi^5>: K is fixed by the cell average, not chosen.
    Second-order finite differences, Newton iteration on (psi, K^2). Returns (psi, K^2)."""
    n = len(r)
    psi = np.ones(n)
    K2 = 24 * math.pi * float(np.sum(Vm * r * r) / np.sum(r * r))
    lo = 1.0 / dr**2 - 1.0 / (r * dr)
    hi = 1.0 / dr**2 + 1.0 / (r * dr)
    for _ in range(itmax):
        rho = grad2 / psi**4 + Vm
        F = np.empty(n + 1)
        lap = -2.0 / dr**2 * psi
        lap[1:] += lo[1:] * psi[:-1]
        lap[:-1] += hi[:-1] * psi[1:]
        lap[0] += lo[0] * psi[0]
        lap[-1] += hi[-1] * psi[-1]
        F[:n] = lap - (K2 / 12 - 2 * math.pi * rho) * psi**5
        F[n] = psi[0] - 1.0
        J = np.zeros((n + 1, n + 1))
        idx = np.arange(n)
        J[idx, idx] = -2.0 / dr**2 - 5 * K2 / 12 * psi**4 + 2 * math.pi * (grad2 + 5 * Vm * psi**4)
        J[idx[1:], idx[:-1]] = lo[1:]
        J[idx[:-1], idx[1:]] = hi[:-1]
        J[0, 0] += lo[0]
        J[n - 1, n - 1] += hi[-1]
        J[:n, n] = -psi**5 / 12
        J[n, 0] = 1.0
        dx = np.linalg.solve(J, -F)
        psi = psi + dx[:n]
        K2 = K2 + dx[n]
        if np.max(np.abs(dx[:n])) < tol and abs(dx[n]) < tol * K2:
            break
    if K2 <= 0 or np.any(psi <= 0) or not np.all(np.isfinite(psi)):
        raise FloatingPointError("no expanding conformally flat cell slice")
    return psi, K2


def conformal_wall_initial_data(grid, geometry, scalars, include_radiation=True):
    """The labelled pocket fixtures on the wall: the same scalar fields and fluid densities on the conformally flat
    slice (``conformal_wall_slice``) instead of the open edge's areal slice. Returns (geometry, matter, H0).

    The slice's source is the engine's own stress projection on the flat metric, split into the gradient energy
    (which scales as psi^-4) and the rest. The fluids start at rest with production's densities today."""
    from .scalar_system import KAPPA
    from .v55_matter import initialize_from_archive

    if not is_lattice_wall(grid):
        raise TypeError("conformal_wall_initial_data needs a lattice-wall grid")
    g = geometry.copy()
    g.a[:] = 1.0
    g.b[:] = 1.0
    g.X[:] = 1.0
    g.Aa[:] = 0.0
    g.Lambda[:] = 0.0
    g.beta[:] = 0.0
    g.B[:] = 0.0
    g.alpha[:] = 1.0
    flat = bind(grid, g)
    r = np.asarray(grid.centers, dtype=float)
    dr = float(grid.dr)
    d1 = lambda f: flat.cell_derivative_fourth(f, parity=1)  # noqa: E731
    grad2 = 0.5 * (d1(scalars.S) ** 2 + d1(scalars.D) ** 2) + 0.5 * d1(scalars.phi) ** 2 / KAPPA
    flat_matter = initialize_from_archive(flat, g, include_radiation=include_radiation)
    rho_flat = total_matter_projection(flat, g, scalars, flat_matter)["rho"]
    Vm = rho_flat - grad2
    psi, K2 = conformal_wall_slice(r, dr, grad2, Vm)
    g.X[:] = psi ** (-2.0)
    g.K[:] = -math.sqrt(K2)
    matter = initialize_from_archive(grid, g, include_radiation=include_radiation)
    return g, matter, math.sqrt(K2) / 3.0
