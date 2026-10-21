"""Stationary Q-ball family of the unchanged Candidate 1 action (D = 0), benchmark couplings.

Track `reach` of the a69bd86 audit. Units follow Codex / charge_phase: S = phi/v, F^2 = 2|X|^2/v^2, length
ell = 1/(sqrt(lam_phi) v), a = gX^2/lam_phi = 1, b = lam_X/lam_phi = 0.1. Physical charge Q = Qbar/lam_phi = 100 Qbar;
E/(m_X Q) = Ebar/(sqrt(a) Qbar); R90 = radius enclosing 90 % of the energy (archived convention, qball.py line 124).

Continuation in Omega starting from the archived profile (grut_candidate_1_20261018/code/results_profile.csv, read only).
Check: Omega = 0.8 must reproduce Qbar = 399.44 (Q = 39,944), E/(m_X Q) = 0.91108652 and R90 = 5.1034.
"""
import json, os, sys
import numpy as np
from scipy.integrate import solve_bvp, trapezoid

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..', '..', '..'))
PROFILE = os.path.join(ROOT, 'research/reality_framework/grut_candidate_1_20261018/code/results_profile.csv')
OUT = os.path.join(HERE, '..', 'results', 'qball_family.json')

a, b = 1.0, 0.1
LAM_PHI = 0.01


def rhs(r, y, p):
    Om = p[0]
    F, Fp, S, Sp = y
    rr = np.where(r > 0, r, 1.0)
    Fpp = np.where(r > 0, -2 * Fp / rr, 0.0) + (a * S**2 - Om**2) * F + b * F**3
    Spp = np.where(r > 0, -2 * Sp / rr, 0.0) + S * (S**2 - 1) + a * S * F**2
    # at r=0 the regular limit is f'' = (1/3) * source
    Fpp = np.where(r > 0, Fpp, ((a * S**2 - Om**2) * F + b * F**3) / 3.0)
    Spp = np.where(r > 0, Spp, (S * (S**2 - 1) + a * S * F**2) / 3.0)
    return np.vstack([Fp, Fpp, Sp, Spp])


def make_bc(F0):
    def bc(ya, yb, p):
        return np.array([ya[1], ya[3], yb[0], yb[2] - 1.0, ya[0] - F0])
    return bc


def observables(r, F, Fp, S, Sp, Om):
    e = 0.5 * Sp**2 + 0.5 * Fp**2 + 0.5 * Om**2 * F**2 + 0.25 * (S**2 - 1)**2 + 0.5 * a * S**2 * F**2 + 0.25 * b * F**4
    q = Om * F**2
    w = 4 * np.pi * r**2
    E = trapezoid(w * e, r)
    Q = trapezoid(w * q, r)
    Ecum = np.concatenate([[0], np.cumsum(0.5 * (w[1:] * e[1:] + w[:-1] * e[:-1]) * np.diff(r))])
    Qcum = np.concatenate([[0], np.cumsum(0.5 * (w[1:] * q[1:] + w[:-1] * q[:-1]) * np.diff(r))])
    r90E = float(np.interp(0.9 * Ecum[-1], Ecum, r))
    r90Q = float(np.interp(0.9 * Qcum[-1], Qcum, r))
    rhalfQ = float(np.interp(0.5 * Qcum[-1], Qcum, r))
    return dict(Omega=float(Om), Qbar=float(Q), Ebar=float(E), Q_phys=float(Q / LAM_PHI),
                E_over_mXQ=float(E / (np.sqrt(a) * Q)), R90_energy_ell=r90E, R90_charge_ell=r90Q,
                Rhalf_charge_ell=rhalfQ, F0sq=float(F[0]**2), S0=float(S[0]))


def rhs_fixed(Om):
    def f(r, y):
        return rhs(r, y, [Om])
    return f


def bc_fixed(ya, yb):
    return np.array([ya[1], ya[3], yb[0], yb[2] - 1.0])


def solve_fixed_omega(Om, guess_r, guess_y, L, n=3000):
    r = np.linspace(0, L, n)
    y = np.vstack([np.interp(r, guess_r, guess_y[i], right=(0.0 if i in (0, 1, 3) else 1.0)) for i in range(4)])
    return solve_bvp(rhs_fixed(Om), bc_fixed, r, y, tol=1e-7, max_nodes=400000)


def solve_at(F0, guess_r, guess_y, Om_guess, L):
    r = np.linspace(0, L, 1500)
    y = np.vstack([np.interp(r, guess_r, guess_y[i], right=(0.0 if i in (0, 1, 3) else 1.0)) for i in range(4)])
    sol = solve_bvp(rhs, make_bc(F0), r, y, p=[Om_guess], tol=1e-7, max_nodes=200000)
    return sol


def main():
    d = np.loadtxt(PROFILE, delimiter=',', skiprows=1)
    r0, F, Fp, S, Sp = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4]
    gy = np.vstack([F, Fp, S, Sp])
    F0_arch = F[0]
    sol = solve_at(F0_arch, r0, gy, 0.8, 30.0)
    assert sol.status == 0, sol.message
    rr = np.linspace(0, sol.x[-1], 20001)
    yy = sol.sol(rr)
    ref = observables(rr, *yy, sol.p[0])
    ref['status'] = int(sol.status)
    print('archived check', ref)
    family = [ref]
    # (1) smaller charges: continuation in central amplitude F0 downward (Omega rises toward m_X)
    gr, gyy, gOm = rr, yy, sol.p[0]
    for F0 in np.arange(F0_arch - 0.025, 1.40, -0.025):
        s = solve_at(F0, gr, gyy, gOm, 40.0)
        if s.status != 0:
            print('stop small-Q branch at F0', F0); break
        rgrid = np.linspace(0, s.x[-1], 40001); ygrid = s.sol(rgrid)
        ob = observables(rgrid, *ygrid, s.p[0]); ob.update(status=0, L=float(s.x[-1]), branch='F0_continuation')
        family.append(ob); gr, gyy, gOm = rgrid, ygrid, s.p[0]
        print('smallQ', round(F0, 4), ob['Omega'], ob['Q_phys'], ob['E_over_mXQ'], ob['R90_energy_ell'])
    # (2) larger charges: continuation in fixed Omega downward toward m* (thin-wall limit)
    gr, gyy = rr, yy
    L = 30.0
    for Om in [0.79, 0.78, 0.76, 0.74, 0.72, 0.70, 0.68, 0.66, 0.65, 0.64, 0.63, 0.62, 0.615, 0.61, 0.605, 0.60]:
        ok = False
        for Ltry in (L, 1.5 * L, 2.2 * L):
            s = solve_fixed_omega(Om, gr, gyy, Ltry)
            if s.status == 0:
                rgrid = np.linspace(0, s.x[-1], 60001); ygrid = s.sol(rgrid)
                ob = observables(rgrid, *ygrid, Om)
                if ob['R90_energy_ell'] < 0.6 * s.x[-1] and ob['Q_phys'] > family[0]['Q_phys'] * 0.99:
                    ok = True; break
        if not ok:
            print('stop large-Q branch at Omega', Om, s.status, s.message); break
        ob.update(status=0, L=float(s.x[-1]), branch='Omega_continuation')
        family.append(ob); gr, gyy = rgrid, ygrid
        L = max(30.0, 2.5 * ob['R90_energy_ell'])
        print('largeQ', Om, ob['Q_phys'], ob['E_over_mXQ'], ob['R90_energy_ell'], ob['F0sq'], ob['S0'])
    family.sort(key=lambda o: o['Q_phys'])
    # thin-wall comparison: E/Q -> m* + 3 sigma_Q/(n_* R); n_* = 1.778279 (Qbar units per ell^3), m* = 0.5623413
    mstar = (b) ** 0.25
    nstar = 1.0 / np.sqrt(b) * mstar  # nbar_* = Omega_* x_* = m* / sqrt(b) = 1.778
    for o in family:
        Rtw = (3 * o['Qbar'] / (4 * np.pi * nstar)) ** (1 / 3)
        o['R_thinwall_ell'] = float(Rtw)
    out = dict(
        description='Stationary Q-ball family, unchanged Candidate 1, benchmark couplings, D=0, flat space',
        units='Q_phys = Qbar/lam_phi; mass M = (E/(m_X Q)) * gX v * Q_phys; radius in ell = 1/(sqrt(lam_phi) v) = 1/(100 GeV)(TeV/v)',
        archived_check=dict(target_Qphys=39944, target_E_over_mXQ=0.91108652, target_R90=5.1034322,
                            got_Qphys=ref['Q_phys'], got_E_over_mXQ=ref['E_over_mXQ'], got_R90=ref['R90_energy_ell']),
        mstar_over_mX=float(mstar), nstar_bar=float(nstar),
        family=family)
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
