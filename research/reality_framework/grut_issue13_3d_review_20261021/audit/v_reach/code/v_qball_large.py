"""Verifier v_reach: large-Q points of the constrained-minimisation family (Q = 6.19e5, 1e6, 1.52e6), h = 0.08, 0.04.
Usage: OMP_NUM_THREADS=1 python3 v_qball_large.py ../results/v_qball_large.json   (~4 min)"""
import json, sys
import numpy as np
from v_qball_min import Ball
res_all = {}
for h in (0.08, 0.04):
    prev = None
    for Qb in (6185.3, 10000.0, 15199.1):
        L = max(30.0, 3.0 * (3 * Qb / (4 * np.pi * 1.778)) ** (1 / 3) + 18.0)
        B = Ball(h, L)
        z0 = B.guess(Qb) if prev is None else B.interp(*prev)
        res, E, parts = B.solve(Qb, z0)
        ob = B.observables(res.x, Qb, E, parts)
        prev = (B.r, res.x, B.N)
        print(h, ob['Q_phys'], ob['e'], ob['Omega'], ob['R90E'], res.nit, flush=True)
        res_all[f'{h}_{Qb}'] = dict(e=float(ob['e']), Omega=float(ob['Omega']), R90=ob['R90E'])
json.dump(res_all, open(sys.argv[1], 'w'), indent=1)
