'''
DARK SECTOR SOLVER - linear perturbations for coupled quintessence.

Conformal Newtonian gauge, conformal time tau, ' = d/dtau:
  ds^2 = a^2[ -(1+2Psi)dtau^2 + (1-2Phi)dx^2 ], Phi = Psi.
Units 8 pi G = 1. H below is the CONFORMAL Hubble rate a'/a.

MODEL. Dark matter mass depends on the field, m_dm(phi) = m_0 exp(-beta phi),
giving the Amendola coupling and the background continuity equation.

BACKGROUND
  a'    = a H, H = a sqrt(rho_tot/3)
  phi'' = -2H phi' - a^2 V_,phi + a^2 beta rho_dm
  rho_dm' = -3H rho_dm - beta rho_dm phi'
  rho_b' = -3H rho_b
  rho_phi = phi'^2/(2a^2) + V , p_phi = phi'^2/(2a^2) - V

PERTURBATIONS. Varying mass means number and energy perturbations differ:
  rho_dm = m(phi) n => delta_dm = delta_n/n - beta delta-phi

  delta_dm' = -theta_dm + 3 Phi' - beta delta-phi'            (1)
  theta_dm' = -H theta_dm + k^2 Phi - beta k^2 delta-phi
         + beta phi' theta_dm                       (2)
  delta_b' = -theta_b + 3 Phi'
   theta_b' = -H theta_b + k^2 Phi                           (3)
   delta-phi'' = -2H delta-phi' - (k^2 + a^2 V_,phiphi) delta-phi
            + 4 phi' Phi' - 2 a^2 (V_,phi - beta rho_dm) Phi
            + a^2 beta rho_dm delta_dm                        (4)
   Phi' = -H Phi - [k^2 Phi + (a^2/2) delta-rho_tot] / (3H)        (5)
with delta-rho_tot = rho_dm delta_dm + rho_b delta_b + delta-rho_phi,
delta-rho_phi = (phi' delta-phi' - phi'^2 Phi)/a^2 + V_,phi delta-phi.
(5) is evaluated FIRST and substituted into (1),(3),(4).

*** EQUATION (2) HAS HAD TWO SIGN ERRORS. BOTH TERMS ARE NOW CORRECT. ***

FIFTH FORCE, -beta k^2 delta-phi. A particle with m = m_0 exp(-beta phi) feels
force per unit mass -grad(Psi) + beta grad(delta-phi) = -grad(Psi - beta
delta-phi), and -grad(X) maps to +k^2 X in theta'. The opposite sign makes the
scalar force REPULSIVE and suppresses clustering. Caught by rung 3.

DRAG, +beta phi' theta. Free streaming conserves comoving momentum a m v, so
v ~ 1/(a m) and dln(theta)/dtau = -H - dln(m)/dtau = -H + beta phi'. The
opposite sign implies a v / m is conserved instead, and loses a factor of ~300
over a test span. Caught by rung 6, and independently confirmed against
Amendola 1999 eq (36), whose drag term reads +C phi_dot theta_c with his
rho_m = rho_m0 a^-3 exp(-C phi) matching ours exactly, so C <-> beta with no
sign flip.

(NOTE ON THAT PAPER: its eq (36) fifth-force term reads +C k^2 varphi, which
contradicts its OWN eq (1) coupling definition and its own radiation-sector
convention theta_gamma_dot = k^2 delta_gamma/4. Derived from his eq (1), the
term must be MINUS. Our sign is also favoured by measurement, since rung 3
returns the 1 + 2 beta^2 ENHANCEMENT and a flipped sign gives suppression.
This is the archive's one unpublished item.)

BOTH ERRORS WERE INVISIBLE wherever beta = 0 or the field is frozen. Any future
change to equation (2) is checked against a conservation law before it is
trusted.
'''
import numpy as np
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings("ignore")

OM_DM, OM_B, OM_DE = 0.26, 0.05, 0.69          # units H0 = 1, a0 = 1 => rho_crit0 = 3
class Model:
    def __init__(self, beta=0.0, lam=1.0, D=0.0, V0=None, w_const=None):
        self.beta, self.lam, self.D = beta, lam, D
        self.V0 = 3.0 * OM_DE if V0 is None else V0
        self.w_const = w_const

    def V(self, phi):
        with np.errstate(over='ignore'):
            return self.V0 * (np.exp(-self.lam * phi) - self.D)

    def dV(self, phi):
        with np.errstate(over='ignore'):
            return -self.lam * self.V0 * np.exp(-self.lam * phi)

    def d2V(self, phi):
        with np.errstate(over='ignore'):
            return self.lam**2 * self.V0 * np.exp(-self.lam * phi)



def background_rhs(tau, y, M):
  a, phi, pp, rdm, rb = y
  if a <= 0 or not np.isfinite(a):
      return [0, 0, 0, 0, 0]
  rphi = pp**2 / (2 * a**2) + M.V(phi)
  rtot = rdm + rb + rphi
  if rtot <= 0:
      return [0, 0, 0, 0, 0]
  H = a * np.sqrt(rtot / 3.0)
  return [a * H,
         pp,
         -2 * H * pp - a**2 * M.dV(phi) + a**2 * M.beta * rdm,
         -3 * H * rdm - M.beta * rdm * pp,
         -3 * H * rb]



def full_rhs(tau, y, M, k):
  a, phi, pp, rdm, rb, ddm, tdm, db, tb, dphi, dpp, Phi = y
  if a <= 0 or not np.isfinite(a):
      return [0] * 12
  rphi = pp**2 / (2 * a**2) + M.V(phi)
  rtot = rdm + rb + rphi
  if rtot <= 0:
      return [0] * 12
  H = a * np.sqrt(rtot / 3.0)

  drho_phi = (pp * dpp - pp**2 * Phi) / a**2 + M.dV(phi) * dphi
  drho = rdm * ddm + rb * db + drho_phi

  Phip = -H * Phi - (k * k * Phi + 0.5 * a**2 * drho) / (3.0 * H) # (5) first

  dphip = dpp
  ddpp = (-2 * H * dpp - (k * k + a**2 * M.d2V(phi)) * dphi
       + 4 * pp * Phip
       - 2 * a**2 * (M.dV(phi) - M.beta * rdm) * Phi
       + a**2 * M.beta * rdm * ddm)                  # (4)

  dddm = -tdm + 3 * Phip - M.beta * dpp                  # (1)
  # (2) BOTH SIGNS CORRECTED - see module docstring.
  dtdm = (-H * tdm + k * k * Phi - M.beta * k * k * dphi
        + M.beta * pp * tdm)                      # (2)
  ddb = -tb + 3 * Phip                           # (3)
  dtb = -H * tb + k * k * Phi

  return [a * H, pp,
       -2 * H * pp - a**2 * M.dV(phi) + a**2 * M.beta * rdm,
       -3 * H * rdm - M.beta * rdm * pp,
       -3 * H * rb,
       dddm, dtdm, ddb, dtb, dphip, ddpp, Phip]



def initial_conditions(M, a_i, Phi_i=1.0):
  """Adiabatic, super-horizon, matter dominated. delta_m = -2 Phi."""
  rdm = 3 * OM_DM / a_i**3
  rb = 3 * OM_B / a_i**3
  return [a_i, 0.0, 0.0, rdm, rb,
         -2 * Phi_i, 0.0, -2 * Phi_i, 0.0, 0.0, 0.0, Phi_i]



def evolve(M, k, a_i=1e-4, a_f=1.0, n=4000):
  y0 = initial_conditions(M, a_i)
  ev = lambda t, y, *A: y[0] - a_f
  ev.terminal, ev.direction = True, 1
  s = solve_ivp(full_rhs, (0, 1e6), y0, args=(M, k), method='LSODA',
             rtol=1e-10, atol=1e-14, events=ev, dense_output=True,
             max_step=np.inf)
  tt = np.linspace(s.t[0], s.t[-1], n)
  Y = s.sol(tt)
  return dict(a=Y[0], phi=Y[1], pp=Y[2], rdm=Y[3], rb=Y[4],
          ddm=Y[5], tdm=Y[6], db=Y[7], tb=Y[8],
          dphi=Y[9], dpp=Y[10], Phi=Y[11], tau=tt, sol=s)
