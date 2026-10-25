#!/usr/bin/env python3
"""Independent ground-state profile solver (Claude review, track xcharge).

Model (quoted from the packets, re-derived in NOTES.md section 1):
  L = |dX|^2 + (1/2)(dphi)^2 - g^2 phi^2 |X|^2 - lam_X |X|^4 - (lam_phi/4)(phi^2-v^2)^2
Dimensionless: rho = m_X r, tau = m_X t, m_X = g v, X = (v/sqrt2) F(rho) exp(-i Omega tau), phi = v S.
  F'' + 2F'/rho = (S^2 + A F^2 - Omega^2) F,      A = lam_X/g^2 = 0.1
  S'' + 2S'/rho = (B (S^2-1) + F^2) S,           B = lam_phi/g^2 = 3
  g^2 Q = 4 pi Omega int rho^2 F^2 drho

Method (differs from Codex): solve at FIXED Omega with collocation, Robin
outer conditions, then root-find Omega for the target charge with Brent.
"""
import numpy as np
from scipy.integrate import solve_bvp, simpson
from scipy.optimize import brentq

G, VPHI = 0.1, 0.924089637375362
MX = G * VPHI
A, B = 0.1, 3.0


def solve_fixed_omega(Om, R=70.0, guess=None, tol=1e-10, n=3001, R0=4.5):
    kap = np.sqrt(1 - Om**2)
    r = np.linspace(0, R, n)

    def rhs(x, y):
        F, dF, S, dS = y
        return np.vstack([dF, (S * S + A * F * F - Om**2) * F, dS, (B * (S * S - 1) + F * F) * S])

    Smat = np.zeros((4, 4)); Smat[1, 1] = -2; Smat[3, 3] = -2

    def bc(a, b):
        return np.array([a[1], a[3], b[1] + kap * b[0], b[2] - 1.0])

    if guess is None:
        F0 = 2.6
        w = 0.5 * (1 - np.tanh((r - R0) / 0.9))
        y0 = np.vstack([F0 * w, np.gradient(F0 * w, r), 1 - w * (1 - 0.01), np.gradient(1 - w, r)])
    else:
        y0 = guess.sol(r)
    sol = solve_bvp(rhs, bc, r, y0, S=Smat, tol=tol, max_nodes=200000)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def charge(sol, Om, R=None):
    x = np.linspace(0, sol.x[-1], 40001)
    F = sol.sol(x)[0]
    return 4 * np.pi * Om * simpson(x * x * F * F, x=x)  # = g^2 Q


def ground_state(Q, Om_lo=0.85, Om_hi=0.875, R=70.0, tol=1e-10, R0=4.5):
    target = G * G * Q
    cache = {}

    def f(Om):
        prev = None
        if cache:
            k = min(cache, key=lambda o: abs(o - Om))
            prev = cache[k]
        s = None
        if prev is not None:
            try:
                s = solve_fixed_omega(Om, R=R, guess=prev, tol=tol)
            except RuntimeError:
                s = None
        if s is None or s.sol(0.0)[0] < 0.5:  # staged tolerance from a tanh guess
            s0 = solve_fixed_omega(Om, R=R, guess=None, tol=1e-6, R0=R0)
            s = solve_fixed_omega(Om, R=R, guess=s0, tol=tol)
        if s.sol(0.0)[0] < 0.5:
            raise RuntimeError('trivial branch')
        cache[Om] = s
        return charge(s, Om) - target

    Om = brentq(f, Om_lo, Om_hi, xtol=1e-13, rtol=1e-14)
    s = solve_fixed_omega(Om, R=R, guess=cache[min(cache, key=lambda o: abs(o - Om))], tol=tol)
    return Om, s


def energy(sol, Om):
    """Energy in units of v/g (E = (v/g) * e_bar)."""
    x = np.linspace(0, sol.x[-1], 40001)
    F, dF, S, dS = sol.sol(x)
    dens = 0.5 * (Om**2 * F**2 + dF**2) + 0.5 * dS**2 + 0.5 * S**2 * F**2 + 0.25 * A * F**4 + 0.25 * B * (S**2 - 1)**2
    return 4 * np.pi * simpson(x * x * dens, x=x)


if __name__ == '__main__':
    import json, sys
    Q = float(sys.argv[1]) if len(sys.argv) > 1 else 200000.
    Om, s = ground_state(Q)
    e = energy(s, Om)
    print(json.dumps(dict(Q=Q, Omega=Om, F0=float(s.sol(0)[0]), S0=float(s.sol(0)[2]),
                          E_over_QmX=float(e / (G * G * Q)), nodes=len(s.x)), indent=1))
