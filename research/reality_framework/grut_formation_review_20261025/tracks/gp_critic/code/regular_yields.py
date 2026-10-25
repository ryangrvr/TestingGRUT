#!/usr/bin/env python3
"""gp_critic: yield and charge fraction in the regular (truncation-valid) regime versus H_inf, with the companion.

Uses run() from companion_ad.py (own code). Best supplied phase phi0 = pi/2, instant reheating, radiation release.
eta = |alpha6| X0^2/lam with X0 = sqrt(cH/(2 lam)) H; the Hubble-held minimum exists only for eta <= r/12.
Usage: python3 regular_yields.py <out.json>
"""
import json
import sys

import numpy as np

from companion_ad import run

lam, MRED, gstar, Y_Q = 1e-3, 2.435e18, 106.75, 4.8294009e-9
M_Y, h, mu, kap = 1e12, 1e11, 5e11, 1e-3
Dm = M_Y**4 - 2 * (mu * h)**2
alpha6 = mu * kap**2 * (h / np.sqrt(2)) / Dm
beta6 = M_Y**2 * kap**2 / Dm
r = 2 * alpha6 / beta6
D = 2 * np.pi**2 / 45 * gstar * (90 / (np.pi**2 * gstar))**0.75
out = dict(r=r, rows=[])
for cH in (1.0, 10.0, 100.0):
    X_edge = np.sqrt(lam / (6 * beta6))
    H_edge = X_edge / np.sqrt(cH / (2 * lam))
    for H in (2e10, 1e11, 0.5 * H_edge, 0.95 * H_edge):
        X0 = np.sqrt(cH / (2 * lam)) * H
        eta = alpha6 * X0**2 / lam
        if eta > r / 12:
            continue
        res = run(eta, r, np.pi / 2, cH=cH)
        Y = res['q'] * cH**1.5 * H**1.5 / (2 * lam * D * MRED**1.5)
        out['rows'].append(dict(cH=cH, H_GeV=H, H_edge_GeV=H_edge, eta=eta, eps=res['eps'], q=res['q'],
                                Y=Y, Y_over_YQ=Y / Y_Q))
json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
for row in out['rows']:
    print({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in row.items()})
