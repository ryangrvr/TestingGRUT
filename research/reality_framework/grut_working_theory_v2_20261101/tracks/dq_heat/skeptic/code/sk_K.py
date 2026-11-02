#!/usr/bin/env python3
"""dq_heat SKEPTIC (Claude, own code, written fresh): the late radiation constant K = a^2 H of the (m_I, Gamma_I) = (3, 0.5)
release from an inflaton-only background, independent of rptri_dq.c.

Quadratic inflaton phi (M_P = 1 units) released at rest at a = 1 with H = 1 (H_i units): phi_0 = sqrt(6) H/m;
  phi'' + (3H + Gamma) phi' + m^2 phi = 0,   rho_r' = -4 H rho_r + Gamma phi'^2,   3 H^2 = phi'^2/2 + m^2 phi^2/2 + rho_r.
The hidden sector (|V| ~ 3e4 H_i^4 against 3 M_P^2 H_i^2 ~ 1e15) is dropped.  dq_heat's K = 1.98581 (engine) and the
dependence of K on (m_I, Gamma_I) are compared.  Vector abundance carries K^{3/2}, R_vec carries K^2 (frozen part).
Output: ../results/sk_K.json.   Labels: [C].
"""
import json, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

RES = Path(__file__).resolve().parent.parent / 'results'


def K_of(m, G, tend=4000.0):
    phi0 = math.sqrt(6.0) / m

    def H_of(y):
        la, ph, pd, rr = y
        return math.sqrt(max((0.5 * pd * pd + 0.5 * m * m * ph * ph + rr) / 3.0, 0.0))

    def rhs(t, y):
        la, ph, pd, rr = y
        H = H_of(y)
        return [H, pd, -(3 * H + G) * pd - m * m * ph, -4 * H * rr + G * pd * pd]

    sol = solve_ivp(rhs, (0, tend), [0.0, phi0, 0.0, 0.0], method='DOP853', rtol=1e-12, atol=1e-16, max_step=0.05,
                    dense_output=True)
    out = {}
    for tt in (0.25, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, tend):
        y = sol.sol(tt); a = math.exp(y[0]); H = H_of(y)
        out[str(tt)] = dict(a=a, H=H, a2H=a * a * H)
    yend = sol.y[:, -1]
    a = math.exp(yend[0])
    return a * a * H_of(yend), out


def main():
    K, hist = K_of(3.0, 0.5)
    print('K(3, 0.5) = %.6f   (dq_heat engine 1.985805)' % K)
    for k, v in hist.items():
        print('  tau=%s a=%.6f H=%.6f a2H=%.6f' % (k, v['a'], v['H'], v['a2H']))
    scan = []
    for m, G in ((3.0, 0.25), (3.0, 1.0), (2.0, 0.5), (4.0, 0.5), (3.0, 3.0), (3.0, 10.0)):
        Kx, _ = K_of(m, G)
        scan.append(dict(m_I=m, Gamma_I=G, K=Kx, K15=Kx ** 1.5, K2=Kx ** 2))
        print('m_I=%.1f Gamma_I=%.2f K=%.4f K^1.5=%.3f K^2=%.3f' % (m, G, Kx, Kx ** 1.5, Kx ** 2))
    json.dump(dict(method=__doc__, K_3_0p5=K, history=hist, scan=scan, dq_heat_K=1.9858052754013398,
                   rel_diff=K / 1.9858052754013398 - 1), open(RES / 'sk_K.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
