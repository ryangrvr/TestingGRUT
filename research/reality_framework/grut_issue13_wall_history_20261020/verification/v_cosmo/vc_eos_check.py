#!/usr/bin/env python3
"""Equation-of-state checks for homogeneous hidden-field configurations of the unchanged action
(Minkowski, period/time averages).  Units v = 1.  numpy/scipy only, deterministic.

w = <p>/<rho>, rho = K + U, p = K - U, K = sum fdot^2/2 (+|Xdot|^2 for complex X),
U measured from the true vacuum (U = 0 at phi = +-v, X = D = 0).
"""
import json
import math
import os

import numpy as np
from scipy.integrate import solve_ivp

lam, gX, lamX, hD, lamD = 0.01, 0.1, 0.001, 0.01, 0.001
v = 1.0
mD0 = gX * v


def U(phi, x1, x2, D):
    X2 = 0.5 * (x1**2 + x2**2)          # |X|^2 with X = (x1 + i x2)/sqrt(2)
    return (lam / 4) * (phi**2 - v**2) ** 2 + gX**2 * phi**2 * X2 + lamX * X2**2 \
        + 0.5 * (mD0**2 + hD**2 * phi**2) * D**2 + (lamD / 4) * D**4


def grad(phi, x1, x2, D):
    X2 = 0.5 * (x1**2 + x2**2)
    dphi = lam * phi * (phi**2 - v**2) + 2 * gX**2 * phi * X2 + hD**2 * phi * D**2
    dX2 = gX**2 * phi**2 + 2 * lamX * X2
    dx1 = dX2 * x1
    dx2 = dX2 * x2
    dD = (mD0**2 + hD**2 * phi**2) * D + lamD * D**3
    return dphi, dx1, dx2, dD


def rhs(t, y):
    phi, x1, x2, D, pphi, p1, p2, pD = y
    g = grad(phi, x1, x2, D)
    return [pphi, p1, p2, pD, -g[0], -g[1], -g[2], -g[3]]


def w_of(y0, T=4000.0):
    ts = np.linspace(0, T, 400001)
    sol = solve_ivp(rhs, [0, T], y0, method="DOP853", rtol=1e-10, atol=1e-12, t_eval=ts)
    phi, x1, x2, D, pphi, p1, p2, pD = sol.y
    K = 0.5 * (pphi**2 + p1**2 + p2**2 + pD**2)
    Uu = U(phi, x1, x2, D)
    rho = K + Uu
    p = K - Uu
    return float(np.mean(p) / np.mean(rho)), float(np.mean(rho)), float(rho.std() / rho.mean())


out = {}
# (a) phi oscillation inside one vacuum (amplitude 0.5 v about +v), and over the barrier
for start in (0.5, 1.1, 1.3, 1.45, 2.0, 5.0):
    y0 = [start, 0, 0, 0, 0, 0, 0, 0]
    w, rho, drift = w_of(y0)
    out[f"phi_start_{y0[0]:.2f}v"] = dict(w=w, w_minus_third=w - 1 / 3, rho_v4=rho,
                                          over_barrier=bool(U(y0[0], 0, 0, 0) > lam / 4 * v**4))
# analytic: 3p - rho = lam v^2 (<phi^2> - v^2)
# (b) rotating X (circular orbit) at phi = v (true vacuum), various radii; phi held by its own dynamics
for r in (0.1, 0.5, 1.0, 3.0, 10.0):
    # circular orbit in (x1,x2) at phi = v requires omega^2 = dU/d|X|^2 ; start phi at v (it will respond)
    X2 = 0.5 * r**2
    om = math.sqrt(gX**2 * v**2 + 2 * lamX * X2)
    y0 = [v, r, 0, 0, 0, 0, om * r, 0]
    w, rho, drift = w_of(y0, T=3000.0)
    out[f"X_rot_r{r}"] = dict(w=w, w_minus_third=w - 1 / 3, rho_v4=rho)
# (c) D oscillation
for amp in (0.5, 5.0, 50.0):
    y0 = [v, 0, 0, amp, 0, 0, 0, 0]
    w, rho, drift = w_of(y0, T=3000.0)
    out[f"D_amp{amp}"] = dict(w=w, w_minus_third=w - 1 / 3, rho_v4=rho)
out["note"] = ("w > 1/3 only for phi oscillations crossing both vacua with <phi^2> > v^2; "
               "requires rho > lam v^4/4 = 2.5e9 GeV^4 at v = 1 TeV")
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vc_eos_results.json")
with open(path, "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out, indent=1))
