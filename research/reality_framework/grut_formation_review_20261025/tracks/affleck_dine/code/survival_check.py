#!/usr/bin/env python3
"""Independent free-Green (Born) estimate of the X^d late core-charge loss (own code).

Derivation (independent of Codex): with X = (v/sqrt2) F(r) e^{-i w t} and
dV = lambda X^d + h.c., the linear response obeys (Helmholtz, free exterior)
  (nabla^2 + k^2) dX = d lambda^* (v/sqrt2)^(d-1) F^(d-1) e^{+i(d-1) w t}
(negative frequency: anti-X of energy (d-1)w).  Outgoing solution
  dX -> A e^{ikr}/r,  A = -(1/k) int sin(k r') r' S(r') dr'.
The number current of a complex scalar with L=|dX|^2 is 2k|A|^2/r^2, so
  dN/dt = 8 pi k |A|^2.  In psi units and rho = m_X r this is
  dN/dt = 4 pi m_X |eps_d|^2 K |A_unit|^2 / g^2,
  eps_d = d 2^(1-d/2) c_d v^(d-2) / (m_X^2 M^(d-4)).
Core charge loss per emitted anti-X = (d-1) (energy (d-1)w at chemical potential w).
The full response potential (Codex) is NOT included: this checks normalization and
order of magnitude only.
Usage: python3 -I survival_check.py <profile_csv> <Omega>
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline

G, V = 0.1, 0.924089637375362
MX = G * V
HBAR = 6.582119569e-25
M = 2.435e18
AGE = 13.8e9 * 365.25 * 86400
Q = 200000.


def main():
    path, Omega = sys.argv[1], float(sys.argv[2])
    arr = np.loadtxt(path, delimiter=',', skiprows=1)
    F = CubicSpline(arr[:, 0], arr[:, 1])
    rmax = arr[-1, 0]
    out = {}
    for d in (5, 6):
        K = np.sqrt(((d - 1) * Omega) ** 2 - 1.0)
        res = []
        for n in (20001, 80001):
            r = np.linspace(0, rmax, n)
            A = -simpson(np.sin(K * r) * r * F(r) ** (d - 1), x=r) / K
            res.append(A)
        A = res[-1]
        flux = K * A ** 2
        eps = d * 2 ** (1 - d / 2) * V ** (d - 2) / (MX ** 2 * M ** (d - 4))
        Ndot = 4 * np.pi * MX * eps ** 2 * flux / (G ** 2 * HBAR)
        frac = (d - 1) * Ndot / Q
        cmax = np.sqrt(0.01 / (frac * AGE))
        out[d] = dict(K=K, A_unit_free=A, A_grid_change=abs(res[0] - res[1]), unit_flux=flux,
                      eps_c1=eps, Ndot_anti_c1_per_s=Ndot, free_core_fraction_rate_c1=frac,
                      free_c_max_1pct=cmax, free_unit_c_scale_GeV=M / cmax ** (1 / (d - 4)),
                      profile_rmax=rmax)
        print(d, out[d])
    Path(sys.argv[3]).write_text(json.dumps(out, indent=1)) if len(sys.argv) > 3 else None


if __name__ == '__main__':
    main()
