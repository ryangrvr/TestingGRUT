#!/usr/bin/env python3
"""Own ground-state profile for the Q=200000, lambda_phi=0.03 object (Claude review, track gp_walls).

Model (Candidate 1 light sector, quoted in NOTES.md section 1):
  V = (lam/4)(phi^2-v^2)^2 + g^2 phi^2 |X|^2 + lam_X |X|^4
Units: rho = m_X r, m_X = g v, X = (v/sqrt2) F e^{-i Omega m_X t}, phi = v S.
  F'' + 2F'/rho = (S^2 + a F^2 - Omega^2) F,   a = lam_X/g^2
  S'' + 2S'/rho = (b (S^2-1) + F^2) S,          b = lam/g^2
  g^2 Q = 4 pi Omega int rho^2 F^2 drho
Method (different from Codex's and from tracks/xcharge): one boundary-value problem in which
Omega is an unknown parameter and the charge integral is an extra state with two boundary
conditions, so the charge is imposed directly instead of by a root search over Omega.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp, simpson

G, VPHI, LAMX, LAMP = 0.1, 0.924089637375362, 0.001, 0.03
MX = G * VPHI
A, B = LAMX / G**2, LAMP / G**2


def solve(Q=200000.0, R=60.0, n=2401, tol=1e-9):
    target = G * G * Q
    r = np.linspace(0.0, R, n)

    def rhs(x, y, p):
        Om = p[0]
        F, dF, S, dS, I = y
        return np.vstack([dF, (S * S + A * F * F - Om**2) * F, dS, (B * (S * S - 1) + F * F) * S, x * x * F * F])

    Smat = np.zeros((5, 5)); Smat[1, 1] = -2.0; Smat[3, 3] = -2.0

    def bc(ya, yb, p):
        Om = p[0]
        kap = np.sqrt(max(1 - Om * Om, 1e-12))
        return np.array([ya[1], ya[3], ya[4], yb[1] + (kap + 1.0 / R) * yb[0], yb[2] - 1.0,
                         4 * np.pi * Om * yb[4] - target])

    w = 0.5 * (1 - np.tanh((r - 4.6) / 0.8))
    F0 = 2.6 * w
    S0 = 1 - 0.99 * w
    I0 = np.concatenate([[0.0], np.cumsum(0.5 * (r[1:]**2 * F0[1:]**2 + r[:-1]**2 * F0[:-1]**2) * np.diff(r))])
    y0 = np.vstack([F0, np.gradient(F0, r), S0, np.gradient(S0, r), I0])
    sol = solve_bvp(rhs, bc, r, y0, p=[0.86], S=Smat, tol=1e-5, max_nodes=300000)
    sol = solve_bvp(rhs, bc, sol.x, sol.y, p=sol.p, S=Smat, tol=tol, max_nodes=300000)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def energy_parts(sol):
    x = np.linspace(0, sol.x[-1], 60001)
    F, dF, S, dS, _ = sol.sol(x)
    Om = sol.p[0]
    dens = 0.5 * (Om**2 * F**2 + dF**2) + 0.5 * dS**2 + 0.5 * S**2 * F**2 + 0.25 * A * F**4 + 0.25 * B * (S**2 - 1)**2
    E = (VPHI / G) * 4 * np.pi * simpson(x * x * dens, x=x)
    q = 4 * np.pi * Om * simpson(x * x * F * F, x=x) / G**2
    return x, F, dF, S, dS, Om, E, q


def main(out):
    res = {}
    for R, n in [(60.0, 2401), (80.0, 3201)]:
        sol = solve(R=R, n=n)
        x, F, dF, S, dS, Om, E, q = energy_parts(sol)
        res[f'R{int(R)}'] = dict(Omega=float(Om), F0=float(F[0]), S0=float(S[0]), E_GeV=float(E), Q=float(q),
                                 E_over_mXQ=float(E / (MX * q)), nodes=int(len(sol.x)))
    sol = solve(R=60.0)
    x, F, dF, S, dS, Om, E, q = energy_parts(sol)
    # Bias energy on the unperturbed core: V_bias = Delta/4 (S-1)^2 (S+2).
    DELTA = 1e-22
    dE_bias = 4 * np.pi / MX**3 * DELTA / 4 * simpson(x * x * (S - 1)**2 * (S + 2), x=x)
    rhalf = float(np.interp(0.5, S, x))
    r09 = float(np.interp(0.9, S, x))
    np.savez(sys.argv[3] if len(sys.argv) > 3 else Path(out).with_suffix('.npz'), x=x, F=F, dF=dF, S=S, dS=dS, Om=Om)
    res.update(dict(bias_energy_on_core_GeV=float(dE_bias),
                    rho_S_half=rhalf, r_S_half_GeVinv=rhalf / MX,
                    rho_S_0p9=r09, r_S_0p9_GeVinv=r09 / MX,
                    binding_margin_GeV=float(MX * q - E), object_mass_GeV=float(E)))
    # Reproduction-type comparison with Codex's supplied profile file (input, not Codex code).
    csv = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if csv and csv.exists():
        d = np.loadtxt(csv, delimiter=',', skiprows=1)
        Sc = np.interp(d[:, 0], x, S); Fc = np.interp(d[:, 0], x, F)
        m = d[:, 0] <= 40
        res['max_abs_diff_vs_codex_profile'] = dict(F=float(np.max(np.abs(Fc[m] - d[m, 1]))),
                                                    S=float(np.max(np.abs(Sc[m] - d[m, 3]))))
    Path(out).write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
