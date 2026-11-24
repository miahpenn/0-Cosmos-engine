import numpy as np
from numpy.polynomial.legendre import leggauss
from research import b7_w2_origin_check_v1 as oc
def interval_reference(r, fn, n):
    # y=rB obeys y'=a-b*y. With z=y-r, z'=s-b*z and s=a-1-b*r.
    Kt = -oc.clone.background_H0(True, True)
    left = np.r_[0.,r[:-1]]
    h = r-left
    nodes, weights = leggauss(n)
    nodes, weights = (nodes+1)/2, weights/2
    x = left[:,None]+h[:,None]*nodes
    v = fn(x)
    def coef(x,v):
        return 4*np.pi*x*(v["Sp"]**2+v["Dp"]**2)
    b = coef(x,v)
    rho0 = oc._rho(v["S"],v["PS"],v["D"],v["PD"],v["Sp"],v["Dp"],v["phi"],v["Pi"],0.,oc.RHO_R_PRESENT)
    j = -v["PS"]*v["Sp"]-v["PD"]*v["Dp"]
    Kr = Kt+4*np.pi*x*j
    s = x*x*(-8*np.pi*rho0+2*Kr*Kt+Kt*Kt)-b*x
    # Within each interval evaluate the attenuation from each node to its right edge.
    tail = r[:,None]-x
    u = x[:,:,None]+tail[:,:,None]*nodes
    beta_tail = tail*(coef(u,fn(u))@weights)
    forcing = h*((s*np.exp(-beta_tail))@weights)
    decay = np.exp(-h*(b@weights))
    z = np.empty_like(r)
    last = 0.
    for i in range(len(z)):
        last = decay[i]*last+forcing[i]
        z[i] = last
    return 1+z/r

