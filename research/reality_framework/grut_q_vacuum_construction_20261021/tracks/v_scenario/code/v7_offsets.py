"""v_scenario / v7: full-equation check (m/H_c = 100 toy, my constraint-solved formulation of v2) of the late
de Sitter vacuum energy for the two kinds of offset the scenario uses:
  type i  (heavy threshold at fixed mu):      eps -> eps + c        => rho_V,late = c [1 + (4/3) c/(K0 m^2) + ...]
  type ii (QV1-cal, mu = mu0 - c/q0):         rho_V(Y) += c (1+Y), R += 2c  => rho_V,late = c [1 + c/(K0 m^2) + ...]
(static-root coefficients from v1 checks H; the scenario quotes 4/3 for QV1-cal, NOTES §0 row 8 and §5.2).
Output: ../results/v7_offsets.json
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2_toy_full as v2  # noqa: E402  (my own toy, not the scenario's)


def late_rhoV(m, c, kind, N_end=5.0):
    toy = v2.Toy(m, 1.0, True, Ni=-1.0)
    s0, _ = v2.initial_state(toy, v2.qs_series_post())
    m2 = toy.m2

    def rhoV(y):
        r = 0.75 * m2 * y * y * (1 + 4 * y / 3 + y * y / 3)
        return r + (c * (1 + y) if kind == 'ii' else 0.0)

    def R(y):
        return 3 * m2 * y * (1 + 2 * y + 2 * y * y / 3) + (2 * c if kind == 'ii' else 0.0)

    ci = c if kind == 'i' else 0.0

    def rhs(t, s):
        y, u, N, rr, J = s
        dv = toy.dV(N)
        rd = toy.rho_d(N)
        rho = rr + rd + dv + ci
        S = rd + 4 * dv + 4 * ci
        A = 3 * (1 + y)
        H = (-3 * u + math.sqrt(9 * u * u + 4 * A * (rho + rhoV(y)))) / (2 * A)
        ud = (S + 4 * rhoV(y) - (1 + y) * R(y)) / 3 - 3 * H * u
        return [u, ud, H, H * (-4 * rr - toy.dV1(N)), H * H]
    ev = lambda t, s: s[2] - N_end
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 1e7), s0, method='DOP853', rtol=1e-12, atol=1e-30, events=ev, dense_output=True,
                    max_step=0.25 * 2 * math.pi / m)
    tf = sol.t[-1]
    ts = tf - 2 * math.pi / m * (np.arange(64) / 64)
    vals = []
    for t in ts:
        y, u, N, rr, J = sol.sol(t)
        rho = rr + toy.rho_d(N) + toy.dV(N)
        A = 3 * (1 + y)
        H = (-3 * u + math.sqrt(9 * u * u + 4 * A * (rho + ci + rhoV(y)))) / (2 * A)
        vals.append(3 * H * H * (1 + y) + 3 * H * u - rho)     # total gravitating vacuum (includes c for type i)
    rv = float(np.mean(vals))
    return dict(m=m, c=c, kind=kind, rhoV_late=rv, ratio_minus_1=rv / c - 1,
                four_thirds_c_over_m2=4 * c / (3 * m2), c_over_m2=c / m2, N_end=N_end)


out = [late_rhoV(100.0, 0.03, 'i'), late_rhoV(100.0, 0.03, 'ii')]
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, '..', 'results', 'v7_offsets.json'), 'w') as f:
    json.dump(out, f, indent=1)
for o in out:
    print(o)
