#!/usr/bin/env python3
"""Track sp_vector (Claude, own code): cost check for a larger gauge coupling e.

Late release (Codex window a = 600-1000): the X condensate is a conformally periodic quartic orbit with comoving
E_X = 741.65 and N_X = -238.65 (H_i units). Its Higgs norm gives the vector a periodically modulated mass,
a^2 m_A^2(tau) = 2 e^2 |chi_X(tau)|^2 (+ the slowly growing Sigma part 36 e^2 a^2 h^2, omitted here, which only
detunes). Transverse vector modes (Maxwell is conformal) obey A'' + [k^2 + 2 e^2 |chi_X|^2] A = 0.
Floquet exponents mu(k) per unit conformal time (= per unit a in radiation) are computed from the monodromy matrix
over one radial period (Greene-Kofman-Linde-Starobinsky-type analysis [R]). Longitudinal modes are not computed;
the transverse result is a lower bound on resonant vector production.
Also the ratio m_A/omega_X of the vector mass to the orbit frequency during X domination.
Usage: python3 floquet_vector.py ../results/floquet_vector.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

LX = 1e-3
E = 741.64988
N = 238.64868


def orbit():
    # radial turning points of r'^2 = E - N^2/(4 r^2) - LX r^4
    f = lambda r: E - N**2 / (4 * r * r) - LX * r**4
    rmax = brentq(f, 5.0, 100.0)
    rmin = brentq(f, 0.5, rmax * 0.999) if f(0.5) < 0 else None
    chi0 = rmax + 0j
    th_dot = -N / (2 * rmax**2)
    y0 = [rmax, 0.0, 0.0, rmax * th_dot]

    def rhs(t, y):
        c = y[0] + 1j * y[1]
        acc = -2 * LX * abs(c)**2 * c
        return [y[2], y[3], acc.real, acc.imag]

    # find radial period: next maximum of r
    def ev(t, y):
        return y[0] * y[2] + y[1] * y[3]  # d(r^2)/2 = 0
    ev.direction = -1
    sol = solve_ivp(rhs, (0, 200), y0, method='DOP853', rtol=1e-12, atol=1e-12, events=ev, dense_output=True)
    T = [t for t in sol.t_events[0] if t > 1e-6][0]
    return sol, T, rmax, rmin


def floquet(sol, T, e, k):
    def rhs(t, y):
        c = sol.sol(t)
        m2 = 2 * e * e * (c[0]**2 + c[1]**2)
        w2 = k * k + m2
        return [y[1], -w2 * y[0], y[3], -w2 * y[2]]
    s = solve_ivp(rhs, (0, T), [1, 0, 0, 1], method='DOP853', rtol=1e-11, atol=1e-13)
    Mmat = np.array([[s.y[0, -1], s.y[2, -1]], [s.y[1, -1], s.y[3, -1]]])
    ev = np.linalg.eigvals(Mmat)
    return float(np.log(max(abs(ev))) / T)


def main(out):
    sol, T, rmax, rmin = orbit()
    ts = np.linspace(0, T, 2001)
    r2 = np.array([sol.sol(t)[0]**2 + sol.sol(t)[1]**2 for t in ts])
    omega = 2 * np.pi / T
    res = dict(E_X=E, N_X=N, r_max=rmax, r_min=float(np.sqrt(r2.min())), radial_period=T, omega_radial=omega,
               scope='transverse modes, X part of the Higgs norm only, conformal (exact for the quartic orbit)', scans={})
    for e in [1e-7, 1e-3, 0.01, 0.0287, 0.05, 0.1, 0.3, 1.0]:
        mA_max = e * np.sqrt(2 * r2.max()); mA_min = e * np.sqrt(2 * r2.min())
        ks = np.linspace(0.0, 3.0 * max(omega, mA_max), 241)
        mus = [floquet(sol, T, e, k) for k in ks]
        i = int(np.argmax(mus))
        res['scans'][str(e)] = dict(e=e, mA_conf_min=float(mA_min), mA_conf_max=float(mA_max),
                                     mA_over_omega=float(mA_max / omega), mu_max=float(mus[i]), k_at_max=float(ks[i]),
                                     growth_exponent_over_400_units=float(2 * mus[i] * 400))
        print(e, res['scans'][str(e)], flush=True)
    Path(out).write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
