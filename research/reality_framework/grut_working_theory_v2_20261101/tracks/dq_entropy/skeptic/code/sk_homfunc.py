#!/usr/bin/env python3
"""dq_entropy SKEPTIC (Claude, own code, written fresh): the T = 0 homogeneous M1 functional behind dq_entropy item 3
(spinodal latent heat), recomputed by direct minimisation at fixed charge density.

Units m_X = g v = 1 (g = 0.1, v = 10), lambda_X = 1e-3, lambda_phi = 0.03, eps = 0.002 (M1 tilt).  Complex X with
kinetic |dX|^2, charge density n = 2 omega |X|^2; energy density of a homogeneous rotating condensate
    e(n) = min_{F = |X|^2, phi} [ n^2/(4F) + g^2 phi^2 F + lambda_X F^2 + V_phi(phi) ],
    V_phi = lambda_phi/4 (phi^2 - v^2)^2 - eps v (phi - v)^3   (zero in the vacuum phi = v).
Bulk Q-matter: Omega_b = min_n e(n)/n (zero pressure); spinodal: minimum of omega(n) = de/dn on the homogeneous branch.
Different method from heat_onset.hom_state (which solves the phi equation at fixed F and differentiates along F): here
a 2-D minimisation at fixed n (Nelder-Mead + polish), then finite differences in n.
Latent heat of the supercooled condensate decomposing at n_sp into bulk lumps + vacuum: e(n_sp)/n_sp - Omega_b.
Usage: sk_homfunc.py OUT.json   Labels [C].
"""
import json, math, sys
import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq

G, VV, LX, LP, EPS = 0.1, 10.0, 1e-3, 0.03, 0.002


def Vphi(p):
    return LP / 4 * (p * p - VV * VV) ** 2 - EPS * VV * (p - VV) ** 3


def efun(n, F, p):
    return n * n / (4 * F) + G * G * p * p * F + LX * F * F + Vphi(p)


def e_of_n(n, guess=None):
    """minimise over (ln F, phi); phi in [-v, 1.2 v] by a 1-D inner minimisation at each F (the phi root is unique
    for the tilted model) and an outer 1-D minimisation in ln F."""
    def inner(F):
        r = minimize_scalar(lambda p: efun(n, F, p), bounds=(-0.2 * VV, 1.2 * VV), method='bounded', options=dict(xatol=1e-12))
        return r.fun, r.x
    # F at fixed phi=0 would be (n^2/(8 LX))^{1/3}; bracket widely
    F0 = (n * n / (8 * LX)) ** (1 / 3)
    r = minimize_scalar(lambda lf: inner(math.exp(lf))[0], bounds=(math.log(F0) - 6, math.log(F0) + 1), method='bounded',
                        options=dict(xatol=1e-12))
    F = math.exp(r.x); ev, p = inner(F)
    # polish jointly
    m = minimize(lambda z: efun(n, math.exp(z[0]), z[1]), [r.x, p], method='Nelder-Mead',
                 options=dict(xatol=1e-13, fatol=1e-15, maxiter=4000))
    if m.fun < ev:
        F, p, ev = math.exp(m.x[0]), m.x[1], m.fun
    return ev, F, p


def omega(n, h=None):
    h = h or 1e-4 * n
    return (e_of_n(n + h)[0] - e_of_n(n - h)[0]) / (2 * h)


def main():
    out = {}
    ns = np.geomspace(60, 1500, 80)
    rows = []
    for n in ns:
        e, F, p = e_of_n(n)
        rows.append(dict(n=float(n), e=e, EQ=e / n, F=F, phi_over_v=p / VV, omega_direct=n / (2 * F)))
    out['scan'] = rows
    EQ = np.array([r['EQ'] for r in rows])
    i = int(np.argmin(EQ))
    rb = minimize_scalar(lambda ln: e_of_n(math.exp(ln))[0] / math.exp(ln), bounds=(math.log(ns[i - 2]), math.log(ns[i + 2])),
                         method='bounded', options=dict(xatol=1e-12))
    nb = math.exp(rb.x); eb, Fb, pb = e_of_n(nb)
    out['bulk'] = dict(n_b=nb, Omega_b=eb / nb, omega_at_nb=nb / (2 * Fb), phi_over_v=pb / VV)
    # spinodal: min of omega(n) = n/(2F) (the rotation frequency of the minimising state equals de/dn)
    om = np.array([r['omega_direct'] for r in rows])
    j = int(np.argmin(om))
    rs = minimize_scalar(lambda ln: (lambda n: n / (2 * e_of_n(n)[1]))(math.exp(ln)), bounds=(math.log(ns[j - 2]), math.log(ns[j + 2])),
                         method='bounded', options=dict(xatol=1e-12))
    nsp = math.exp(rs.x); esp, Fsp, psp = e_of_n(nsp)
    out['spinodal'] = dict(n_sp=nsp, EQ_hom=esp / nsp, omega_sp=nsp / (2 * Fsp), omega_fd=omega(nsp), phi_over_v=psp / VV)
    out['latent_heat_T0'] = dict(dE_per_charge=esp / nsp - eb / nb,
                                 dq_entropy=dict(EQ_hom_sp=0.82727, Omega_b=0.77586, dE=0.05140, n_sp=233.22, n_b=449.31, phi_sp=0.34))
    print(json.dumps({k: v for k, v in out.items() if k != 'scan'}, indent=1))
    json.dump(out, open(sys.argv[1], 'w'), indent=1)


if __name__ == '__main__':
    main()
