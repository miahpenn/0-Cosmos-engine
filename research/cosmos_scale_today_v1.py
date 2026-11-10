"""RESEARCH ONLY. The COSMOS lane's scale freedom (V0 is a unit in the isolated lane) and markers along today's trajectory.
Usage (repo root): python research/cosmos_scale_today_v1.py. Backward integration from today leaves the forward-built
history within a few e-folds (scalar decaying modes grow backwards); markers far in the past are not GEAR 20's history."""
import sys; sys.path.insert(0, ".")
import numpy as np, dataclasses
from scipy.integrate import solve_ivp
from engine.cosmos import CosmosParams, construct_present_day, cosmos_rhs, cosmos_rho_p, potential
p = CosmosParams(); s = construct_present_day(p)
y0 = np.array([s.a, s.H, s.phi, s.pi_phi, s.rho_dm, s.rho_b, s.rho_r, s.rho_shear])
def run(p, y0, T, back=False):
    f = lambda t, y: np.array(cosmos_rhs(tuple(y), p))
    ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = -1
    return solve_ivp(f, (0, -T if back else T), y0, events=None if back else ev, rtol=1e-12, atol=1e-30, method="DOP853", dense_output=True)
def obs(p, y):
    a, H, phi, pi, dm, b, r, sh = y
    rho, pr = cosmos_rho_p(tuple(y), p); rde = 0.5*pi*pi + potential(phi, p)
    return dict(rho=rho, w=pr/rho, Ode=rde/rho, Om=(dm+b)/rho, Or=r/rho, wphi=(0.5*pi*pi-potential(phi,p))/rde)
# 1) scale invariance: rescale densities and V0 by c^2, H and pi by c, time by 1/c
A = run(p, y0, 5000); NT = np.log(A.y[0,-1]/y0[0]); tT = A.t[-1]
for c2 in (10.0, 1e-6):
    c = np.sqrt(c2); p2 = dataclasses.replace(p, V0=p.V0*c2)
    y2 = y0.copy(); y2[1]*=c; y2[3]*=c; y2[4:]*=c2
    B = run(p2, y2, 5000/c)
    print(f"V0 x {c2:g}: e-folds today->turnaround {np.log(B.y[0,-1]/y2[0]):.10f} (base {NT:.10f});  H0*t_turn {B.t[-1]*y2[1]:.8f} (base {tT*y0[1]:.8f});  phi_turn {B.y[2,-1]:.8f} (base {A.y[2,-1]:.8f})")
o = obs(p, y0); print(f"today: phi={y0[2]} Omega_de={o['Ode']:.4f} Omega_m={o['Om']:.4f} w_phi={o['wphi']:.4f} H0={y0[1]:.5e}; e-folds to turnaround {NT:.4f}; t_turn*H0 {tT*y0[1]:.3f}")
# 2) backward from today: markers along the trajectory
Bk = run(p, y0, 3000, back=True)
ts = np.linspace(0, Bk.t[-1], 200001); Y = Bk.sol(ts)
lnA = np.log(Y[0]/y0[0])
ob = [obs(p, Y[:,i]) for i in range(0, Y.shape[1], 50)]; idx = list(range(0, Y.shape[1], 50))
Ode = np.array([o['Ode'] for o in ob]); Om = np.array([o['Om'] for o in ob]); w = np.array([o['w'] for o in ob]); Or=np.array([o['Or'] for o in ob])
def mark(name, cond):
    k = np.where(cond)[0]
    if len(k)==0: print(f"  {name}: not reached (back to ln a = {lnA[idx[-1]]:.2f})"); return
    i = idx[k[0]]; print(f"  {name}: ln(a/a_today) = {lnA[i]:+.4f}, phi = {Y[2,i]:.5f}, t*H0 = {ts[i]*y0[1]:+.4f}")
print("backward from today (first point going back where...):")
mark("Omega_de < Omega_m (DE-matter equality)", Ode < Om)
mark("w_total > -1/3 (acceleration starts)", w > -1/3)
mark("Omega_r > Omega_m (radiation-matter equality)", Or > Om)
print(f"  earliest reached: ln a = {lnA[-1]:.2f}, phi = {Y[2,-1]:.5f}, Omega_de = {Ode[-1]:.3e}")
