"""RESEARCH ONLY. Thomas-Fermi (large lambda) self-held real S: averaged EOS from S = A cos(w t), V = S^2/2 + lam S^4/8, then TOV.
x = lam A^2.  rho = x(1 + 9x/32)/(2 lam), P = 3 x^2/(64 lam).  Units G = c = m = 1; results scale with lam."""
import math, numpy as np
from scipy.integrate import solve_ivp
lam = 5000.0
def rho_of_x(x): return x*(1+9*x/32)/(2*lam)
def P_of_x(x): return 3*x*x/(64*lam)
def x_of_P(P): return math.sqrt(max(P,0)*64*lam/3)
def tov(xc):
    Pc = P_of_x(xc)
    def f(r, y):
        m, P, nu = y
        if P <= 0: return [0,0,0]
        rho = rho_of_x(x_of_P(P))
        dnu = 2*(m + 4*math.pi*r**3*P)/(r*(r-2*m))
        return [4*math.pi*r*r*rho, -(rho+P)*dnu/2, dnu]
    ev = lambda r, y: y[1] - 1e-12*Pc; ev.terminal = True
    r0 = 1e-6
    s = solve_ivp(f, [r0, 1e4], [0, Pc, 0], events=ev, rtol=1e-10, atol=1e-14, max_step=1.0)
    R = s.t[-1]; M = s.y[0, -1]
    # lapse: alpha_s = sqrt(1-2M/R); alpha_c = alpha_s*exp((nu_c - nu_s)/2)
    alpha_s = math.sqrt(1-2*M/R); alpha_c = alpha_s*math.exp(-(s.y[2,-1])/2)
    # radius where x drops to 4/3 (averaged hold) and 2/3 (peak hold)
    xs = np.array([x_of_P(p) for p in s.y[1]])
    r43 = s.t[xs >= 4/3].max() if (xs >= 4/3).any() else 0.0
    r23 = s.t[xs >= 2/3].max() if (xs >= 2/3).any() else 0.0
    return M, R, 2*M/R, alpha_c, r43, r23
print("x_c   M      R      2M/R   alpha_c  r(x>4/3) r(x>2/3)   A_c")
for xc in [0.3, 0.67, 1.0, 1.33, 2, 3, 4, 6, 8, 12, 20, 40]:
    M, R, C, ac, r43, r23 = tov(xc)
    print(f"{xc:5.2f} {M:6.3f} {R:6.2f} {C:6.3f} {ac:7.3f} {r43:8.2f} {r23:8.2f}  {math.sqrt(xc/lam):.4f}")
