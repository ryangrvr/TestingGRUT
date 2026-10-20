#!/usr/bin/env python3
"""Cross-check of cw_floor.py: finite-difference Hessian (no analytic second derivatives),
and minimization including D != 0."""
import numpy as np, json
from scipy.optimize import minimize
gX, lp, lX, mD0sq, d, lD, muR = 0.1, 0.01, 0.001, 0.01, 1e-4, 0.001, 0.1
def V0q(q):
    phi, x, y, D = q; R2 = x*x+y*y
    return lp/4*(phi**2-1)**2 + gX**2/2*phi**2*R2 + lX/4*R2**2 + 0.5*(mD0sq+d*phi**2)*D**2 + lD/4*D**4
def fd_hess(q, h=1e-4):
    q = np.asarray(q, float); H = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            e_i = np.eye(4)[i]*h; e_j = np.eye(4)[j]*h
            H[i,j] = (V0q(q+e_i+e_j)-V0q(q+e_i-e_j)-V0q(q-e_i+e_j)+V0q(q-e_i-e_j))/(4*h*h)
    return H
def V1q(q):
    ev = np.linalg.eigvalsh(fd_hess(q))
    return sum(e*e*(np.log(abs(e)/muR**2)-1.5) for e in ev if abs(e)>1e-14)/(64*np.pi**2)
vac = V1q([1,0,0,0])
Rs = 10**0.25
strict = 2*(V1q([0,Rs,0,0]) - vac)/Rs**2
f = lambda p: 2*(V0q([p[0],p[1],0,p[2]]) + V1q([p[0],p[1],0,p[2]]) - vac)/p[1]**2
r = minimize(f, [0.0, Rs, 0.05], method='Nelder-Mead', options={'xatol':1e-10,'fatol':1e-16,'maxiter':20000})
print(json.dumps({'fd_vacuum_V1': vac, 'fd_strict_delta_mstar2': strict,
                  'min_ratio_with_D_free': r.fun, 'minimizer_phi_R_D': [float(x) for x in r.x],
                  'relative_to_tree': r.fun/np.sqrt(lp*lX)-1}, indent=1))
