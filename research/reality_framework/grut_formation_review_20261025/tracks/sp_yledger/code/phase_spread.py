#!/usr/bin/env python3
"""sp_yledger (Claude, own code): side check for the leftover verdict.  How large is the inflationary spread
of the source phase psi2 = 2 theta_Y - theta_Sigma at the leading point, and what net charge fraction does
a Gaussian source-phase spread of that size leave?

Stochastic equilibrium of the light angular mode (mass^2 m^2 < 9/4 H^2): <phi^2> = 3 H^4/(8 pi^2 m^2) [R:
Starobinsky-Yokoyama]; heavy and gauge directions excluded.  The ensemble perturbs only the Sigma phase about
delta = pi/2 (Codex's ensemble convention); radial fluctuations are not added.  Sampling SEM only.
Usage: python3 phase_spread.py <out.json>     (~3 min, 4 cores)
"""
import json, sys
from multiprocessing import Pool
import numpy as np
from ymodel import minimum, params, initial, run, KAP

H = 1.1211167304108865e13


def psi2_rms():
    m = minimum(H, 1.0); x, y, s = m['x'], m['y'], m['s']; p = params(H)
    v1 = np.array([3., -1., 0.]); v2 = np.array([0., 2., -1.])
    Mth = 2 * KAP * y * x**3 * np.outer(v1, v1) + 2 * p['mu'] * s * y * y * np.outer(v2, v2)
    D = np.diag(1 / (np.sqrt(2) * np.array([x, y, s])))
    w, U = np.linalg.eigh(D @ Mth @ D)
    var = sum((v2 @ D @ U[:, k])**2 * 3 / (8 * np.pi**2 * w[k]) for k in range(3) if 1e-8 < w[k] < 2.25)
    return float(np.sqrt(var)), w.tolist()


def member(args):
    d, amp = args
    x, y, s = amp
    z0 = np.array([x, -y, -s * np.exp(1j * d)])
    r = run(H, z0, t_end=0.5 * 1000**2, rtol=1e-10, atol=1e-13, nsamp=4)
    w = r['windows']['600-1000']
    return dict(delta=d, eps=w['eps'], NX=w['NX'], Ncirc=w['Ncirc'])


def main():
    sig, masses = psi2_rms()
    m = minimum(H, 1.0); amp = (m['x'], m['y'], m['s'])
    rng = np.random.default_rng(2026102501)
    out = dict(psi2_rms_equilibrium=sig, angular_mass2_over_H2=masses, ensembles=[])
    for s in (0.03, sig):
        ds = np.pi / 2 + s * rng.standard_normal(128)
        with Pool(4) as p:
            rows = p.map(member, [(d, amp) for d in ds])
        e = np.array([r['eps'] for r in rows]); nx = np.array([r['NX'] for r in rows]); nc = np.array([r['Ncirc'] for r in rows])
        out['ensembles'].append(dict(sigma=s, N=len(rows), mean_eps=float(e.mean()), sem=float(e.std(ddof=1) / np.sqrt(len(e))),
                                     net_eps=float(nx.mean() / nc.mean()), frac_negative=float((e < 0).mean()),
                                     rows=rows))
        print(s, e.mean(), e.std(ddof=1) / np.sqrt(len(e)), nx.mean() / nc.mean(), flush=True)
    json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
