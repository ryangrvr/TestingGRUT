#!/usr/bin/env python3
"""kin_thermo (Claude, own code): finite-Q energy E(Q), omega(Q) = dE/dQ and geometry of the light sector's
two-field Q-balls, tabulated along the stationary branch by continuation in the internal frequency Omega.

Model (quoted, bb211d1 MATHEMATICS.md sec. 2 / our conversion and xcharge tracks): rho = m_X r, X = (v/sqrt2) F(rho)
e^{-i Omega tau}, phi = v S(rho), a = lambda_X/g^2 = 0.1, b = lambda_phi/g^2 = 3:
    F'' + 2F'/rho = (S^2 + a F^2 - Omega^2) F,      S'' + 2S'/rho = (b (S^2 - 1) + F^2) S,
    g^2 Q = 4 pi Omega int rho^2 F^2 drho,   E = (v/g) 4 pi int rho^2 [ (F'^2+S'^2)/2 + Omega^2 F^2/2 + S^2F^2/2
                                                                      + a F^4/4 + b (S^2-1)^2/4 ] drho.
Method (own; differs from conversion/profile.py (u = rho F, z = rho(S-1), z(L) = 0) and xcharge/profile.py
(F, S with a singular-term matrix)): the unknowns are (F, F', S, S') on [rho0, L] with rho0 = 1e-4, the inner
condition is the regular power series through O(rho^2) at rho0, and the outer conditions carry the 1/rho of the
Yukawa tails, F' + (kappa + 1/L) F = 0 with kappa = sqrt(1 - Omega^2) and S' + (sqrt6 + 1/L)(S - 1) = 0.  L is
enlarged until the charge outside L/2 is below 1e-9 of Q.  scipy solve_bvp (collocation), tol 1e-9, continuation
in Omega from 0.86 towards the thin-wall (Omega -> 0.74) and the thick-wall (Omega -> 1) limits.  Observables by
Simpson quadrature on 60001 points.  Checks: dE/dQ = Omega m_X along the table (exact identity for stationary
Q-balls), the Derrick/virial identity E_grad + 3(E_pot - E_kin) = 0, and the archived Q = 200000 reference.
Output: results/objects_EQ.json.  Runtime about 10-20 min.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp, simpson
from scipy.optimize import brentq

A, B = 0.1, 3.0
G, VPHI = 0.1, 0.924089637375362
MX = G * VPHI
RHO0 = 1e-4
RES = Path(__file__).resolve().parents[1] / 'results'


def rhs(r, y, om):
    F, Fp, S, Sp = y
    return np.vstack([Fp, (S * S + A * F * F - om * om) * F - 2 * Fp / r,
                      Sp, (B * (S * S - 1) + F * F) * S - 2 * Sp / r])


def bc_factory(om, L):
    kap = np.sqrt(1 - om * om)

    def bc(ya, yb):
        F0, Fp0, S0, Sp0 = ya
        cF = S0 * S0 + A * F0 * F0 - om * om
        cS = B * (S0 * S0 - 1) + F0 * F0
        return np.array([Fp0 - RHO0 * F0 * cF / 3, Sp0 - RHO0 * S0 * cS / 3,
                         yb[1] + (kap + 1 / L) * yb[0], yb[3] + (np.sqrt(6) + 1 / L) * (yb[2] - 1)])
    return bc


def guess_fn(F0=2.62, S0=0.006, R=4.3, w=0.9):
    def g(r):
        t = 1 / (1 + np.exp((r - R) / w))
        F = F0 * t
        S = 1 - (1 - S0) * t
        return np.vstack([F, np.gradient(F, r), S, np.gradient(S, r)])
    return g


def solve(om, L, init=None, tol=1e-9, n=3001):
    r = np.linspace(RHO0, L, n)
    y0 = guess_fn()(r) if init is None else init(r)
    sol = solve_bvp(lambda x, y: rhs(x, y, om), bc_factory(om, L), r, y0, tol=tol, max_nodes=600000)
    if not sol.success:
        raise RuntimeError('solve_bvp: ' + sol.message)
    return sol


def observables(sol, om, L, n=60001):
    r = np.linspace(RHO0, L, n)
    F, Fp, S, Sp = sol.sol(r)
    kin = 0.5 * om * om * F * F
    grad = 0.5 * (Fp * Fp + Sp * Sp)
    pot = 0.5 * S * S * F * F + 0.25 * A * F**4 + 0.25 * B * (S * S - 1)**2
    w = r * r
    qn = 4 * np.pi * om * simpson(w * F * F, x=r)
    e_kin, e_grad, e_pot = (4 * np.pi * simpson(w * z, x=r) for z in (kin, grad, pot))
    eb = e_kin + e_grad + e_pot
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] * F[1:]**2 + w[:-1] * F[:-1]**2) * np.diff(r))])
    cum /= cum[-1]
    r50, r90, r99 = (float(np.interp(f, cum, r)) for f in (0.5, 0.9, 0.99))
    half = int(n // 2)
    outside_half = 1 - cum[half]
    Q = qn / G**2
    E = VPHI / G * eb
    return dict(Omega=float(om), Q=float(Q), E_GeV=float(E), E_over_QmX=float(eb / qn),
                E_minus_Q_omega_GeV=float(E - Q * om * MX), F0=float(F[0]), S0=float(S[0]),
                r50_over_mX=r50, r90_over_mX=r90, r99_over_mX=r99, L=float(L),
                virial=float((e_grad + 3 * (e_pot - e_kin)) / eb), charge_outside_L_half=float(outside_half),
                E_kin_frac=float(e_kin / eb), E_grad_frac=float(e_grad / eb), E_pot_frac=float(e_pot / eb),
                nodes=int(sol.x.size), rms_residual=float(np.max(sol.rms_residuals)))


class Branch:
    """Continuation with adaptive box size."""

    def __init__(self, om0=0.86, L0=60.0):
        self.last = None
        self.L = L0
        self.om = om0

    def __call__(self, om, tol=1e-9):
        init = None if self.last is None else (lambda r, s=self.last: s.sol(np.clip(r, RHO0, s.x[-1])))
        for _ in range(6):
            sol = solve(om, self.L, init, tol=tol)
            ob = observables(sol, om, self.L)
            if ob['charge_outside_L_half'] < 1e-9 and ob['r99_over_mX'] < 0.45 * self.L:
                break
            self.L *= 1.6                       # enlarge the box and resolve from the current solution
            init = (lambda r, s=sol: s.sol(np.clip(r, RHO0, s.x[-1])))
        if ob['F0'] < 0.2:
            raise RuntimeError('trivial branch reached at Omega=%r' % om)
        self.last = sol
        self.om = om
        return sol, ob


def omega_for_charge(branch, Q, lo, hi):
    return brentq(lambda om: branch(om)[1]['Q'] - Q, lo, hi, xtol=1e-13, rtol=1e-13)


def main():
    t0 = time.time()
    out = dict(method=__doc__, rows=[])
    br = Branch()
    sol, base = br(0.86)
    print('seed', {k: base[k] for k in ('Omega', 'Q', 'E_over_QmX', 'F0', 'S0', 'virial')}, flush=True)
    # reference charge
    om_ref = omega_for_charge(br, 200000.0, 0.855, 0.867)
    sol_ref, ref = br(om_ref)
    print('Q=200000', ref, flush=True)
    out['reference_Q200000'] = ref
    # save my own profile for the kernels (scratch, not results)
    scr = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/kin_thermo')
    scr.mkdir(parents=True, exist_ok=True)
    r = np.linspace(RHO0, br.L, 120001)
    F, Fp, S, Sp = sol_ref.sol(r)
    np.savez(scr / 'profile_Q200000.npz', r=r, F=F, S=S, Fp=Fp, Sp=Sp, Omega=om_ref)
    # finite-difference dE/dQ at the reference, two step sizes
    fd = []
    for h in (2e-4, 5e-5):
        _, op = br(om_ref + h)
        _, om_ = br(om_ref - h)
        fd.append(dict(dOmega=h, dEdQ_over_mX=(op['E_GeV'] - om_['E_GeV']) / ((op['Q'] - om_['Q']) * MX),
                       dQ_dOmega=(op['Q'] - om_['Q']) / (2 * h)))
    out['dEdQ_finite_difference_at_reference'] = fd
    br(om_ref)
    # table: thin-wall direction (down) then thick-wall direction (up)
    rows = {}
    down = list(np.arange(0.86, 0.745, -0.005))
    up = list(np.arange(0.86, 0.9961, 0.005)) + [0.998]
    for seq in (down, up):
        br.last = sol_ref
        br.L = max(60.0, ref['L'])
        for om in seq:
            try:
                _, ob = br(float(om))
            except Exception as exc:             # recorded, not hidden; the branch ends here
                print('stop at Omega=%.4f: %r' % (om, exc), flush=True)
                out.setdefault('branch_ends', []).append(dict(Omega=float(om), error=repr(exc)))
                break
            rows[round(float(om), 6)] = ob
            print('%.4f Q=%.6g E/Qm=%.7f F0=%.4f S0=%.5f r90=%.3f L=%.0f vir=%.1e nodes=%d  %.0fs' % (
                om, ob['Q'], ob['E_over_QmX'], ob['F0'], ob['S0'], ob['r90_over_mX'], ob['L'], ob['virial'],
                ob['nodes'], time.time() - t0), flush=True)
            out['rows'] = [rows[k] for k in sorted(rows)]
            (RES / 'objects_EQ.json').write_text(json.dumps(out, indent=1) + '\n')
    tab = out['rows']
    # dE/dQ along the table (central differences in the table index), the branch label and Q_min
    Om = np.array([x['Omega'] for x in tab]); Q = np.array([x['Q'] for x in tab]); E = np.array([x['E_GeV'] for x in tab])
    dEdQ = np.gradient(E, Om) / np.gradient(Q, Om) / MX
    dQdOm = np.gradient(Q, Om)
    for x, d, dq in zip(tab, dEdQ, dQdOm):
        x['dEdQ_over_mX_table'] = float(d)
        x['dQ_dOmega'] = float(dq)
        x['branch'] = 'stable (dQ/dOmega<0)' if dq < 0 else 'thick-wall (dQ/dOmega>0, classically unstable)'
    i = int(np.argmin(Q))
    out['Q_min'] = dict(Omega=float(Om[i]), Q=float(Q[i]), E_over_QmX=float(tab[i]['E_over_QmX']))
    out['max_abs_dEdQ_minus_Omega_stable_interior'] = float(np.max(np.abs(dEdQ[2:-2] - Om[2:-2])))
    out['seconds'] = time.time() - t0
    (RES / 'objects_EQ.json').write_text(json.dumps(out, indent=1) + '\n')
    print('done', out['Q_min'], out['max_abs_dEdQ_minus_Omega_stable_interior'], out['seconds'])


if __name__ == '__main__':
    main()
