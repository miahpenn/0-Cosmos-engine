"""Corrected V5.5 production initial data.

The production initializer uses the archive-corrected COSMOS operating point:

    a=33.8983
    phi=0.179055
    phidot=0.00341913
    rho_DM=2.5857e-5
    rho_b =4.0306e-6

Radiation is enabled by default using rho_r = Omega_r*rho_total.
The archived Bianchi-I shear is deliberately excluded from the exact
spherical geometry.

The older reference_pirk_unified.make_initial() is retained unchanged as an
audit/reference artifact; it contains a historical variable-mapping mismatch
and is therefore not used here.

Initial structure (hand-floor removal step 2, Nov 16 2026; decision of Nov 15):
production's default start is TODAY, the background alone, labelled
"no structure placed": S = 0 and D at the sea's rest D_s(phi_today)
everywhere, so the floor is carried by D (engine/model_constants.py).
phi_today is a chart (borrowed). Settled stones enter only as the S-D species
the lane carries; a reset start is the goal once the floor question resolves.
The old S pocket survives only as labelled fixtures:

  "none"                  production default: background alone, no structure placed.
  "rescaled_pocket"       TEST FIXTURE (a): the old pocket's shape with its amplitude
                          rescaled so that x = lambda A^2 equals pocket_x, CHOSEN BY HAND.
  "inherited_pocket"      TEST FIXTURE (b): the old pocket's amplitude (0.01) carried
                          over unchanged to the new lambda, so x = lambda A^2 ~ 157
                          (INHERITED, not chosen: an artefact of the old lambda = 0 data).
  "hand_floor_reference"  REFERENCE ONLY: the pre-step-2 default (pocket of amplitude
                          0.01, D seed 1e-10, PD = D), built with the old arithmetic.
                          The default under engine.scalar_system.hand_floor_reference().

In both fixtures D is relaxed with friction against the time-averaged S in
flat space (the cold start of PR #470), from its adiabatic rest
D^2 = (2/lambda)(e^(-2 beta_D phi) - (3/4) lambda S(r)^2) clipped at 0, plus the
D_amplitude seed. S is placed at a turning point (PS = 0), not settled: a
fixture, not a settled stone.

The slice built here is the AREAL slice (the open edge's data). Since the
pocket-fixtures step (Nov 16 2026) the production kernel runs the RESCALED
fixture on the lattice wall by default, where its geometry is rebuilt on the
conformally flat cell slice from the same S and D profiles (LABELLED: the
fixture's data changed; engine/lattice_wall.py, POCKET_FIXTURE_WALL_LABEL). The
INHERITED fixture keeps the open edge and this areal slice by default, LABELLED
"wall behaviour unresolved" (lattice_wall.INHERITED_POCKET_WALL_LABEL). The
areal slice stays each fixture's open-edge reference (outer_edge="open_edge").
"""
from dataclasses import dataclass
import math
import numpy as np

from .scalar_system import cosmos_potential, sd_potential
from . import scalar_system as _scalar_system
from .v55_matter import (
    RHO_TOTAL_PRESENT,
    RHO_DM_PRESENT,
    RHO_B_PRESENT,
    RHO_R_PRESENT,
    V55MatterState,
    initialize_from_archive,
)


@dataclass(frozen=True)
class CorrectedInitialData:
    geometry: object
    scalars: object
    matter: V55MatterState
    H0: float
    include_radiation: bool


def _d1(grid, x, parity=1):
    if hasattr(grid, "cell_derivative_fourth"):
        return grid.cell_derivative_fourth(x, parity=parity)
    return grid.d1(x, parity)


def _cumulative_origin_integral(grid, r, dr, f):
    """Fourth-order cumulative integral F(r_i) = int_0^{r_i} f(s) ds of an even f.

    Used for the areal reconstruction (r B)' = src. Any constant error in F is a
    homogeneous solution, B ~ eps / r: smooth in the continuum but odd and
    singular at the origin, so the even-parity fourth-order stencils turn
    eps ~ h^3 into a resolution-independent Hamiltonian-constraint violation at
    the centre cells. The former rule (a rectangle over [0, r0] and plain
    trapezoids) left exactly such an O(h^3) constant.

    Over [0, r0]: the exact integral of the even quartic c0 + c2 s^2 + c4 s^4
    through cells 0-2. Then trapezoids with the Euler-Maclaurin end correction
    -(h^2/12) [f'(r_i) - f'(r0)], f' from the grid's even-parity derivative.
    """
    A = np.vstack([np.ones(3), r[:3] ** 2, r[:3] ** 4]).T
    c = np.linalg.solve(A, f[:3])
    first = c[0] * r[0] + c[1] * r[0] ** 3 / 3.0 + c[2] * r[0] ** 5 / 5.0
    fp = _d1(grid, f, 1)
    integ = np.cumsum(np.r_[first, 0.5 * dr * (f[1:] + f[:-1])])
    return integ - (dr * dr / 12.0) * (fp - fp[0])


INITIAL_STRUCTURES = ("none", "rescaled_pocket", "inherited_pocket", "hand_floor_reference")


def _pocket_D(S, phi0, D_amplitude, r, dr, width, T=200.0, gamma=0.3):
    """D for the pocket fixtures: relaxed with friction against the time-averaged S, in flat space, then the
    D_amplitude seed added (the cold start of PR #470; research/dm_one_budget_v1.cold_D, ported). Fixtures only.

    Starts from D's adiabatic rest against the time-averaged S, D^2 = (2/lambda)(e^(-2 beta_D phi) - (3/4) lambda S^2),
    clipped at 0, and relaxes the gradient energy away (damped wave equation for w = r D, odd ghost at the centre,
    Neumann at r_max) so the profile is smooth. Not a settled stone: S itself is placed, not settled.
    """
    lam = _scalar_system.SD_QUARTIC_LAMBDA
    eD = math.exp(-2.0 * _scalar_system.BETA_D * phi0)
    x_avg = 0.75 * lam * S * S                      # (3/2) lambda <S^2>, <S^2> = S^2/2 at a turning point
    D = np.sqrt(np.clip((2.0 / lam) * (eD - x_avg), 0.0, None))
    w = r * D
    wp = w.copy()
    stiff = float(np.max(np.abs(eD - x_avg))) + 2.0 * eD
    dt = min(0.4 * dr, 0.5 / math.sqrt(stiff))
    c = 0.5 * gamma * dt
    for _ in range(int(math.ceil(T / dt))):
        lap = np.empty_like(w)
        lap[1:-1] = (w[2:] - 2.0 * w[1:-1] + w[:-2]) / dr**2
        lap[0] = (w[1] - 3.0 * w[0]) / dr**2
        lap[-1] = (w[-2] - 2.0 * w[-1] + w[-1] * (r[-1] + dr) / r[-1]) / dr**2
        Dv = w / r
        aw = lap + w * (eD - 0.5 * lam * Dv * Dv - x_avg)
        wn = (2.0 * w - (1.0 - c) * wp + dt * dt * aw) / (1.0 + c)
        wp, w = w, wn
    return w / r + D_amplitude * np.exp(-(r / width) ** 2)


def build_initial_data(
    grid,
    vacuum_state_factory,
    amplitude=None,
    width=None,
    D_amplitude=None,
    include_radiation=True,
    *,
    initial_structure=None,
    pocket_x=None,
    charged_dust_density=None,
    radial_maxwell_charge_to_mass=None,
):
    r = np.asarray(grid.centers)
    n = len(r)
    dr = float(grid.dr)
    rho_r = RHO_R_PRESENT if include_radiation else 0.0

    # Optional conventional radial-Maxwell experiment: a separate charged
    # dust sector starts at rest. Its density contributes to the initializer's
    # matter budget, and its Gaussian electric energy is included in the same
    # iterative radial Hamiltonian reconstruction.
    if charged_dust_density is None and radial_maxwell_charge_to_mass is None:
        rho_charged = np.zeros(n, dtype=float)
        charge_to_mass = None
    elif charged_dust_density is None or radial_maxwell_charge_to_mass is None:
        raise ValueError(
            "charged_dust_density and radial_maxwell_charge_to_mass must be supplied together"
        )
    else:
        raw_density = np.asarray(charged_dust_density)
        if raw_density.dtype.kind not in "iuf" or raw_density.shape != (n,):
            raise ValueError("charged_dust_density must be a real profile matching the grid")
        rho_charged = np.asarray(raw_density, dtype=float)
        if not np.all(np.isfinite(rho_charged)) or np.any(rho_charged < 0.0):
            raise ValueError("charged_dust_density must be finite and nonnegative")
        if isinstance(radial_maxwell_charge_to_mass, (bool, np.bool_)):
            raise ValueError("radial_maxwell_charge_to_mass must be finite and real")
        try:
            charge_to_mass = float(radial_maxwell_charge_to_mass)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("radial_maxwell_charge_to_mass must be finite and real") from exc
        if not math.isfinite(charge_to_mass):
            raise ValueError("radial_maxwell_charge_to_mass must be finite and real")

    structure = (
        _scalar_system.DEFAULT_INITIAL_STRUCTURE if initial_structure is None else initial_structure
    )
    if structure not in INITIAL_STRUCTURES:
        raise ValueError(f"initial_structure must be one of {INITIAL_STRUCTURES}, not {structure!r}")
    lam = _scalar_system.SD_QUARTIC_LAMBDA
    if structure != "hand_floor_reference" and lam <= 0.0:
        raise ValueError(
            "D has no floor at lambda = 0 (D is unbounded); only the labelled reference "
            "initial_structure='hand_floor_reference' runs there"
        )
    # The background-only default places no structure, so it refuses pocket arguments instead of silently
    # dropping them (Jeremiah's review of PR #503). The pockets and the reference keep the old defaults.
    pocket_args = {"amplitude": amplitude, "width": width, "D_amplitude": D_amplitude}
    if structure == "none":
        given = sorted(k for k, v in pocket_args.items() if v is not None)
        if given:
            raise ValueError(
                f"the background-only default ('no structure placed') takes no pocket arguments; got {given}. "
                "Use initial_structure='rescaled_pocket' / 'inherited_pocket' (labelled fixtures) or the "
                "hand_floor_reference() path"
            )
    amplitude = 0.01 if amplitude is None else amplitude
    width = 7.0 if width is None else width
    D_amplitude = 1.0e-10 if D_amplitude is None else D_amplitude
    if structure == "rescaled_pocket":
        if pocket_x is None or not math.isfinite(float(pocket_x)) or float(pocket_x) <= 0.0:
            raise ValueError("the rescaled pocket needs pocket_x > 0, chosen by hand (label it so)")
    elif pocket_x is not None:
        raise ValueError("pocket_x belongs to initial_structure='rescaled_pocket' only")

    phi0 = 0.179055
    pi0 = 0.00341913
    dm0 = RHO_DM_PRESENT
    b0 = RHO_B_PRESENT

    # Archive-normalized spherical budget. The exact spherical branch
    # excludes anisotropic shear, so it is not silently folded into rho.
    rho0 = 0.5 * pi0**2 + float(cosmos_potential(np.array([phi0]))[0])
    rho0 += dm0 + b0 + rho_r
    if structure != "hand_floor_reference":
        # The sea's floor (D at rest, S = 0 far out), Einstein-normalised: KAPPA V_SD(0, D_s; phi_today).
        # At the today-matched lambda this equals the old hand floor, so H0 is unchanged.
        rho0 += float(_scalar_system.sea_floor_density(phi0))
    H0 = math.sqrt(max(rho0 / 3.0, 0.0))

    if structure == "hand_floor_reference":
        # REFERENCE ONLY: the pre-step-2 default initial data, unchanged.
        S = amplitude * np.exp(-(r / width) ** 2)
        D = D_amplitude * np.exp(-(r / width) ** 2)
        PD = D.copy()
    elif structure == "none":
        # Production default: background alone, "no structure placed".
        S = np.zeros_like(r)
        D = np.full(n, float(_scalar_system.sd_rest_value(phi0)))
        PD = np.zeros_like(r)
    else:
        if structure == "rescaled_pocket":
            A = math.sqrt(float(pocket_x) / lam)          # FIXTURE (a): x chosen by hand
        else:
            A = amplitude                                  # FIXTURE (b): inherited amplitude, x = lambda A^2
        S = A * np.exp(-(r / width) ** 2)
        D = _pocket_D(S, phi0, D_amplitude, r, dr, width)
        PD = np.zeros_like(r)
    PS = np.zeros_like(r)
    phi = np.full(n, phi0)
    Pi = np.full(n, pi0)
    Kt = np.full(n, -H0)

    scalars = _scalar_system.ScalarFields(S, PS, D, PD, phi, Pi)
    Sp = _d1(grid, S, 1)
    Dp = _d1(grid, D, 1)
    phip = _d1(grid, phi, 1)
    Ktp = _d1(grid, Kt, 1)
    # Archive fluids start at rest, and the radial electric field carries no j_r.
    j = _scalar_system.scalar_momentum_projection(scalars, (Sp, Dp, phip))
    B = np.ones_like(r)

    # Same constrained areal initial-data reconstruction used by the archived
    # V5.5 gate, with corrected COSMOS variables and optional radiation.
    for _ in range(6):
        if _scalar_system.coupling_active():
            # GEAR-91 (step 2): the full V_SD(S, D; phi), counted once.
            rhoS = 0.5 * (PS**2 + B * Sp**2)
            rhoD = 0.5 * (PD**2 + B * Dp**2) + sd_potential(S, D, phi=phi)
        else:
            rhoS = 0.5 * (PS**2 + B * Sp**2) + 0.5 * S**2
            rhoD = 0.5 * (PD**2 + B * Dp**2) - 0.5 * D**2
            if _scalar_system.SD_QUARTIC_LAMBDA != 0.0:
                # GEAR-67 quartic without the coupling (uncoupled research setting).
                rhoD = rhoD + sd_potential(S, D) - (0.5 * S**2 - 0.5 * D**2)
        rhoPhi = (0.5 * Pi**2 + cosmos_potential(phi)) / (8.0 * math.pi)
        rhoMatter = (dm0 + b0 + rho_r + rho_charged) / (8.0 * math.pi)
        rho_em = np.zeros(n, dtype=float)
        if charge_to_mass is not None and np.any(rho_charged != 0.0):
            # For the initialized isotropic metric, R_areal=r and sqrt(gamma)
            # per unit solid angle is r^2*A, A=1/sqrt(B). The cumulative Q
            # follows the exact shell rest-mass budget; u_E=Ehat^2/(8*pi)
            # is then included in the initial Hamiltonian-source iteration.
            A_iter = 1.0 / np.sqrt(B)
            sqrt_gamma_iter = r * r * A_iter
            shell_charge = (
                4.0 * math.pi * charge_to_mass
                * rho_charged * sqrt_gamma_iter * dr
            )
            q_faces = np.r_[0.0, np.cumsum(shell_charge)]
            q_centres = 0.5 * (q_faces[:-1] + q_faces[1:])
            ehat = q_centres / (r * r)
            rho_em = ehat * ehat / (8.0 * math.pi)
        rho = rhoS + rhoD + rhoPhi + rhoMatter + rho_em
        Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
        src = (
            1.0 - 8.0 * math.pi * r * r * rho
            + 2.0 * r * r * Kr * Kt
            + r * r * Kt * Kt
        )
        integ = _cumulative_origin_integral(grid, r, dr, src)
        B = integ / r

    if np.any(~np.isfinite(B)) or np.any(B <= 0.0):
        raise FloatingPointError("invalid constrained initial metric")

    A = 1.0 / np.sqrt(B)
    Kr = Kt + r * Ktp + 4.0 * math.pi * r * j
    K = Kr + 2.0 * Kt
    Aa = Kr - K / 3.0
    a = A ** (4.0 / 3.0)
    b = A ** (-2.0 / 3.0)
    X = A ** (-1.0 / 3.0)
    Lambda = (
        _d1(grid, a, 1) / (2.0 * a * a)
        - _d1(grid, b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )

    state = vacuum_state_factory(grid)
    state.a[:] = a
    state.b[:] = b
    state.X[:] = X
    state.alpha[:] = 1.0
    state.beta[:] = 0.0
    state.Aa[:] = Aa
    state.K[:] = K
    state.Lambda[:] = Lambda
    state.B[:] = 0.75 * Lambda

    matter = initialize_from_archive(
        grid,
        state,
        include_radiation=include_radiation,
        charged_dust_density=(
            rho_charged if charge_to_mass is not None else None
        ),
    )
    return CorrectedInitialData(
        geometry=state,
        scalars=scalars,
        matter=matter,
        H0=H0,
        include_radiation=include_radiation,
    )
