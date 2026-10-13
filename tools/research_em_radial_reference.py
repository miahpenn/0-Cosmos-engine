"""Analytic reference oracle for a proposed spherical Maxwell/Coulomb sector.

Research-only flat-space fixture. This module is deliberately not imported by
or connected to the production evolution engine: it is an exact benchmark for
a future radial Gauss-law discretization, not a Maxwell solver or a GEAR
electromagnetism claim.

Conventions: Gaussian units, c = 1, spherical Minkowski background, positive
or negative finite central charge density rho0, and
    rho_e(r) = rho0 * (1 - (r/R0)**2)**2, 0 <= r < R0,
with rho_e = 0 outside. The orthonormal radial field is E_hat_r = Q(r)/r**2.
"""
from __future__ import annotations

from dataclasses import dataclass
import math


FOUR_PI = 4.0 * math.pi


def _finite_real(value: float, name: str) -> float:
    """Return a validated finite real; reject bools to keep the API unambiguous."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


@dataclass(frozen=True)
class RadialElectricStress:
    """Eulerian orthonormal projections for a static pure radial E field."""

    energy_density: float
    radial_pressure: float
    tangential_pressure: float
    radial_momentum_density: float


@dataclass(frozen=True)
class RadialCoulombReference:
    """Closed-form compact charge profile and its flat-space Coulomb field.

    The charge profile is a manufactured numerical fixture only. It does not
    assign charge to any existing GEAR field or matter species.
    """

    rho0: float
    source_radius: float

    def __post_init__(self) -> None:
        rho0 = _finite_real(self.rho0, "rho0")
        radius = _finite_real(self.source_radius, "source_radius")
        if radius <= 0.0:
            raise ValueError("source_radius must be positive")
        object.__setattr__(self, "rho0", rho0)
        object.__setattr__(self, "source_radius", radius)

    @staticmethod
    def _radius(radius: float) -> float:
        value = _finite_real(radius, "radius")
        if value < 0.0:
            raise ValueError("radius must be nonnegative")
        return value

    @property
    def total_charge(self) -> float:
        """Exact total charge: (32 pi / 105) rho0 R0**3."""
        return (32.0 * math.pi / 105.0) * self.rho0 * self.source_radius**3

    def charge_density(self, radius: float) -> float:
        """Return the compactly supported analytic charge density rho_e(r)."""
        r = self._radius(radius)
        if r >= self.source_radius:
            return 0.0
        x = r / self.source_radius
        return self.rho0 * (1.0 - x * x) ** 2

    def enclosed_charge(self, radius: float) -> float:
        """Return exact Q(r) = 4 pi integral_0^r rho_e(s) s**2 ds."""
        r = self._radius(radius)
        if r >= self.source_radius:
            return self.total_charge
        r2 = r * r
        r3 = r2 * r
        R2 = self.source_radius * self.source_radius
        return FOUR_PI * self.rho0 * (
            r3 / 3.0 - 2.0 * r3 * r2 / (5.0 * R2)
            + r3 * r2 * r2 / (7.0 * R2 * R2)
        )

    def electric_field(self, radius: float) -> float:
        """Return orthonormal E_hat_r = Q(r)/r**2, with regular E_hat_r(0)=0."""
        r = self._radius(radius)
        if r == 0.0:
            return 0.0
        return self.enclosed_charge(r) / (r * r)

    def stress_energy(self, radius: float) -> RadialElectricStress:
        """Return exact pure-radial-E projections in Gaussian units.

        For this static purely radial field there is no Poynting flux:
        rho_EM = E**2/(8 pi), p_r = -rho_EM, p_t = +rho_EM, j_r = 0.
        """
        field = self.electric_field(radius)
        energy_density = field * field / (8.0 * math.pi)
        return RadialElectricStress(
            energy_density=energy_density,
            radial_pressure=-energy_density,
            tangential_pressure=energy_density,
            radial_momentum_density=0.0,
        )

    def gauss_shell_residual(
        self,
        r_inner: float,
        r_outer: float,
        E_inner: float,
        E_outer: float,
    ) -> float:
        """Return flat-space finite-volume Gauss residual on one radial shell.

        The integral identity is
            [r**2 E_hat_r]_inner^outer
                - 4 pi integral_inner^outer rho_e(r) r**2 dr = 0.
        The source integral is evaluated exactly as Q(outer)-Q(inner). Field
        values are inputs so a future numerical field can be tested against
        this source oracle. This residual is a shell-integrated Gauss check,
        not a claim that a differential field evolution has been solved.
        """
        inner = self._radius(r_inner)
        outer = self._radius(r_outer)
        if outer < inner:
            raise ValueError("r_outer must be greater than or equal to r_inner")
        e_inner = _finite_real(E_inner, "E_inner")
        e_outer = _finite_real(E_outer, "E_outer")
        flux_difference = outer * outer * e_outer - inner * inner * e_inner
        source_charge = self.enclosed_charge(outer) - self.enclosed_charge(inner)
        return flux_difference - source_charge



@dataclass(frozen=True)
class CurvedRadialCoulombReference:
    """Exact Gauss-law oracle on a fixed, nontrivial spherical conformal metric.

    This manufactured spatial metric uses the production BSSN variable
    convention a=b=1, X=1/psi, psi(r)=1+eta*(r/R0)**2. Therefore
        gamma_rr = psi**2,
        gamma_thth = (r*psi)**2,
        sqrt_gamma = r**2 * psi**3,
        R_areal = r*psi.
    The compact charge profile is specified against coordinate r; the
    enclosed charge integrates rho_e against the proper spatial volume.

    This is a frozen-background metric/constraint fixture, not evolved EM or
    an Einstein-Maxwell trajectory. eta must be positive so the metric is
    nontrivial and regular on the whole domain.
    """

    rho0: float
    source_radius: float
    metric_strength: float = 0.4

    def __post_init__(self) -> None:
        rho0 = _finite_real(self.rho0, "rho0")
        radius = _finite_real(self.source_radius, "source_radius")
        strength = _finite_real(self.metric_strength, "metric_strength")
        if radius <= 0.0:
            raise ValueError("source_radius must be positive")
        if strength <= 0.0:
            raise ValueError("metric_strength must be positive")
        object.__setattr__(self, "rho0", rho0)
        object.__setattr__(self, "source_radius", radius)
        object.__setattr__(self, "metric_strength", strength)

    @staticmethod
    def _radius(radius: float) -> float:
        value = _finite_real(radius, "radius")
        if value < 0.0:
            raise ValueError("radius must be nonnegative")
        return value

    def conformal_factor(self, radius: float) -> float:
        """Return psi(r) for the manufactured frozen metric."""
        r = self._radius(radius)
        x = r / self.source_radius
        return 1.0 + self.metric_strength * x * x

    def areal_radius(self, radius: float) -> float:
        """Return R=r*sqrt(b)/X=r*psi for a=b=1, X=1/psi."""
        r = self._radius(radius)
        return r * self.conformal_factor(r)

    def sqrt_gamma(self, radius: float) -> float:
        """Return the spatial volume factor without sin(theta)."""
        r = self._radius(radius)
        return r * r * self.conformal_factor(r) ** 3

    def charge_density(self, radius: float) -> float:
        """Return rho_e(r); its compact support is defined in coordinate r."""
        r = self._radius(radius)
        if r >= self.source_radius:
            return 0.0
        x = r / self.source_radius
        return self.rho0 * (1.0 - x * x) ** 2

    def _charge_integral_to(self, radius: float) -> float:
        """Analytic 4*pi integral rho_e*sqrt(gamma) dr up to radius <= R0."""
        x = radius / self.source_radius
        # Expand psi**3=(1+eta*x**2)**3. Integrate each polynomial term
        # exactly, retaining the proper-volume factor in the source.
        integral = 0.0
        for k in range(4):
            coefficient = math.comb(3, k) * self.metric_strength**k
            integral += coefficient * (
                x ** (2 * k + 3) / (2 * k + 3)
                - 2.0 * x ** (2 * k + 5) / (2 * k + 5)
                + x ** (2 * k + 7) / (2 * k + 7)
            )
        return FOUR_PI * self.rho0 * self.source_radius**3 * integral

    @property
    def total_charge(self) -> float:
        """Exact total charge integrated over the curved proper spatial volume."""
        return self._charge_integral_to(self.source_radius)

    def enclosed_charge(self, radius: float) -> float:
        """Return Q(r)=4*pi integral_0^r rho_e*sqrt(gamma) dr."""
        r = self._radius(radius)
        if r >= self.source_radius:
            return self.total_charge
        return self._charge_integral_to(r)

    def electric_field(self, radius: float) -> float:
        """Return the orthonormal radial field E_hat_r=Q/R_areal**2."""
        r = self._radius(radius)
        if r == 0.0:
            return 0.0
        R = self.areal_radius(r)
        return self.enclosed_charge(r) / (R * R)

    def coordinate_electric_field(self, radius: float) -> float:
        """Return contravariant E^r, using E_hat_r=sqrt(gamma_rr)*E^r."""
        r = self._radius(radius)
        return self.electric_field(r) / self.conformal_factor(r)

    def stress_energy(self, radius: float) -> RadialElectricStress:
        """Return pure radial-E projections in the orthonormal Eulerian frame."""
        field = self.electric_field(radius)
        energy_density = field * field / (8.0 * math.pi)
        return RadialElectricStress(
            energy_density=energy_density,
            radial_pressure=-energy_density,
            tangential_pressure=energy_density,
            radial_momentum_density=0.0,
        )

    def gauss_shell_residual(
        self,
        r_inner: float,
        r_outer: float,
        E_inner_r: float,
        E_outer_r: float,
    ) -> float:
        """Return the curved-metric shell residual in coordinate E^r variables.

        The production convention is sqrt(gamma)=sqrt(a)*b*r**2/X**3 and
        R_areal=r*sqrt(b)/X. With E_hat_r=sqrt(gamma_rr)*E^r, the boundary
        flux is sqrt(gamma)*E^r=R_areal**2*E_hat_r. The exact proper-volume
        source is Q(outer)-Q(inner). This is an integral constraint check,
        not a solved Maxwell evolution.
        """
        inner = self._radius(r_inner)
        outer = self._radius(r_outer)
        if outer < inner:
            raise ValueError("r_outer must be greater than or equal to r_inner")
        e_inner = _finite_real(E_inner_r, "E_inner_r")
        e_outer = _finite_real(E_outer_r, "E_outer_r")
        flux_difference = (
            self.sqrt_gamma(outer) * e_outer
            - self.sqrt_gamma(inner) * e_inner
        )
        source_charge = self.enclosed_charge(outer) - self.enclosed_charge(inner)
        return flux_difference - source_charge
