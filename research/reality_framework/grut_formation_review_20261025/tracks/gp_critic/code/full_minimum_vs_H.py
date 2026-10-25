#!/usr/bin/env python3
"""gp_critic: does the full X, Y, Sigma inflationary minimum continue past the truncated-EFT existence edge?

Radial potential with phases fixed (psi1 = 3thX - thY aligned; psi2 = 2thY - thS aligned or at pi/2):
  V = -cH H^2 R_X^2 + lam_X R_X^4 + lam_Y R_Y^4 + lam_S (R_S^2 - h^2/2)^2 + M_Y^2 R_Y^2
      - 2 kap R_Y R_X^3 - 2 mu R_S R_Y^2 cos(psi2 - pi)     (R_i = |field|)
Own code; global minimum from a log grid of starts plus Nelder-Mead/BFGS polish.
Usage: python3 full_minimum_vs_H.py <out.json>
"""
import json
import sys
import warnings

import numpy as np
from scipy.optimize import minimize

warnings.filterwarnings('ignore')

M_Y, h, mu, kap, lam_X, lam_Y, lam_S = 1e12, 1e11, 5e11, 1e-3, 1e-3, 1.0, 1.0
Dm = M_Y**4 - 2 * (mu * h)**2
beta6 = M_Y**2 * kap**2 / Dm
X_edge_eft = np.sqrt(lam_X / (6 * beta6))


def V(p, H, cH, c2):
    RX, RY, RS = np.exp(p)
    return (-cH * H**2 * RX**2 + lam_X * RX**4 + lam_Y * RY**4 + lam_S * (RS**2 - h**2 / 2)**2 + M_Y**2 * RY**2
            - 2 * kap * RY * RX**3 - 2 * mu * RS * RY**2 * c2)


def global_min(H, cH, c2):
    best = None
    for lx in np.linspace(np.log(1e10), np.log(1e16), 7):
        for ly in np.linspace(np.log(1e6), np.log(1e15), 5):
            p0 = [lx, ly, np.log(h / np.sqrt(2))]
            s = 1.0 / H**4
            res = minimize(lambda p: V(p, H, cH, c2) * s, p0, method='Nelder-Mead',
                           options=dict(xatol=1e-12, fatol=1e-16, maxiter=20000))
            res = minimize(lambda p: V(p, H, cH, c2) * s, res.x, method='BFGS', options=dict(gtol=1e-12))
            if best is None or res.fun < best.fun:
                best = res
    RX, RY, RS = np.exp(best.x)
    return dict(RX=RX, RY=RY, RS=RS, V_over_H4=best.fun, kapX2_over_MY2=kap * RX**2 / M_Y**2,
                RX_over_X0=RX / np.sqrt(cH * H**2 / (2 * lam_X)))


out = dict(X_edge_eft_GeV=X_edge_eft, rows=[])
for cH in (1.0, 10.0):
    H_edge = X_edge_eft / np.sqrt(cH / (2 * lam_X))
    for H in (0.5 * H_edge, 0.9 * H_edge, H_edge, 1.5 * H_edge, 3 * H_edge, 1e13, 4.7e13):
        for label, c2 in (('psi2_aligned', 1.0), ('psi2_pi_over_2', 0.0)):
            m = global_min(H, cH, c2)
            out['rows'].append(dict(cH=cH, H=H, H_over_H_edge=H / H_edge, phase=label, **m))
            print(cH, f'{H:.3g}', f'{H / H_edge:.2f}', label,
                  {k: f'{v:.4g}' for k, v in m.items()})
json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
