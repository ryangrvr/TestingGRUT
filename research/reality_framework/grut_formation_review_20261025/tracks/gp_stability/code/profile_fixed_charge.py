#!/usr/bin/env python3
"""Own stationary-profile solver for the reference Q = 200000 object (track gp_stability).

Method (differs from Codex's Background class and from track `conversion`): the
frequency Omega is an UNKNOWN PARAMETER of the boundary-value problem and the
charge constraint is imposed as an extra boundary condition on the running
charge q(rho) = 4 pi Omega int_0^rho u^2 drho', u = rho F.  No root finding.

Model (quoted, C1-HP-03 dimensionless form, a = lambda_X/g^2 = 0.1, b = lambda_phi/g^2 = 3):
    F'' + 2F'/rho = (S^2 + a F^2 - Omega^2) F
    S'' + 2S'/rho = (b (S^2 - 1) + F^2) S
X = (v/sqrt 2) F exp(-i Omega m_X t), phi = v S, rho = m_X r, m_X = g v.
Charge Q = g^-2 * 4 pi Omega int rho^2 F^2 ; energy E = (v/g) * eb with
eb = 4 pi int rho^2 [ (F'^2+S'^2)/2 + Omega^2 F^2/2 + S^2F^2/2 + aF^4/4 + b(S^2-1)^2/4 ].
"""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp

A, B, G = 0.1, 3.0, 0.1
V = 0.924089637375362          # GeV
MX = G * V
QBAR = G * G * 200000.0         # = 2000


def make(L):
    def rhs(r, y, p):
        om = p[0]
        u, up, z, zp, q = y
        rs = np.where(r > 0, r, 1.0)
        F = np.where(r > 0, u / rs, up)
        S = np.where(r > 0, 1 + z / rs, 1 + zp)
        return np.vstack([up, (S * S + A * F * F - om * om) * u,
                          zp, (B * (S * S - 1) + F * F) * S * r,
                          4 * np.pi * om * u * u])

    def bc(ya, yb, p):
        k = np.sqrt(max(1 - p[0] ** 2, 1e-12))
        return np.array([ya[0], ya[2], ya[4], yb[1] + k * yb[0], yb[2], yb[4] - QBAR])
    return rhs, bc


def solve(L=70.0, n=3001, om0=0.86, tol=1e-10):
    rhs, bc = make(L)
    r = np.linspace(0, L, n)
    t = 1 / (1 + np.exp((r - 4.3) / 0.9))
    F = 2.6 * t; S = 1 - 0.994 * t
    u = r * F; z = r * (S - 1)
    q = np.cumsum(4 * np.pi * om0 * u * u) * (r[1] - r[0])
    q *= QBAR / q[-1]
    y0 = np.vstack([u, np.gradient(u, r), z, np.gradient(z, r), q])
    sol = None
    for tl in [1e-5, 1e-7, 1e-9, tol]:   # tolerance continuation from a coarse solve
        new = solve_bvp(rhs, bc, r, y0, p=[om0] if sol is None else sol.p, tol=tl, max_nodes=400000)
        if not new.success:
            assert sol is not None and tl == tol, (tl, new.message)
            break                         # keep the last converged (1e-9) solution
        sol = new; sol.tol_reached = tl
        r = sol.x; y0 = sol.y
    return sol


def fields(sol, r):
    u, up, z, zp, q = sol.sol(r)
    rs = np.where(r > 0, r, 1.0)
    F = np.where(r > 0, u / rs, up)
    S = np.where(r > 0, 1 + z / rs, 1 + zp)
    Fp = np.where(r > 0, (up - F) / rs, 0.0)
    Sp = np.where(r > 0, (zp - (S - 1)) / rs, 0.0)
    return F, Fp, S, Sp, q


def gl_grid(L, panels=400, order=16):
    """Composite Gauss-Legendre nodes/weights on [0, L] (own quadrature)."""
    x, w = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(0, L, panels + 1)
    R = [];W = []
    for a, b in zip(edges[:-1], edges[1:]):
        R.append(0.5 * (b - a) * x + 0.5 * (a + b)); W.append(0.5 * (b - a) * w)
    return np.concatenate(R), np.concatenate(W)


def observables(sol, L):
    om = sol.p[0]
    r, w = gl_grid(L)
    F, Fp, S, Sp, q = fields(sol, r)
    dens = 0.5 * (Fp ** 2 + Sp ** 2) + 0.5 * om ** 2 * F ** 2 + 0.5 * S * S * F * F + A / 4 * F ** 4 + B / 4 * (S * S - 1) ** 2
    qn = 4 * np.pi * om * np.sum(w * r * r * F * F)
    eb = 4 * np.pi * np.sum(w * r * r * dens)
    I4 = 4 * np.pi * np.sum(w * r * r * F ** 4)           # int F^4 d^3rho
    G4 = 4 * np.pi * np.sum(w * r * r * (2 * F * Fp) ** 2)  # int |grad F^2|^2 d^3rho
    E = V / G * eb
    return dict(L=L, Omega=float(om), qbar=float(qn), Q=float(qn / G ** 2), ebar=float(eb), E_GeV=float(E),
                binding_margin_GeV=float(200000 * MX - E), F0=float(F[0]), S0=float(S[0]),
                int_F4_d3rho=float(I4), int_gradF2_sq_d3rho=float(G4), nodes=int(sol.x.size),
                rms_residual=float(np.max(sol.rms_residuals)), bvp_tol=float(sol.tol_reached))


def main():
    out = {}
    sols = {}
    for L, n in [(70.0, 3001), (90.0, 4001)]:
        s = solve(L, n)
        sols[L] = s
        out[f'L{int(L)}'] = observables(s, L)
    out['codex_quoted'] = dict(E_GeV=17064.537454786467, binding_margin_GeV=1417.2552927207762,
                               source='SCREENING_RESULTS.json (b5174e5) [Q]')
    path = Path(__file__).resolve().parents[1] / 'results' / 'profile_fixed_charge.json'
    path.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=1))
    return sols


if __name__ == '__main__':
    main()
