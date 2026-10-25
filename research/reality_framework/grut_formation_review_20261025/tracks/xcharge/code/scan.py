#!/usr/bin/env python3
"""Charge-loss rates versus object charge / density along the ground-state family (track xcharge).

Continuation in Omega from the reference Q = 200000 solution. For each ground state: Q, F0 (central
amplitude), E/(Q m_X), half-amplitude radius, the coherent single-quantum rate (forced.py), and the
pair estimate (pairs.py), for d = 5 and 6 with unit coefficient at the reduced Planck scale.
"""
import json, sys
import numpy as np

from profile import G, solve_fixed_omega, charge, energy, ground_state
from forced import Bg, solve, rates
from pairs import pair_rate


def main(out):
    Om0, sol0 = ground_state(200000.)
    grid_up = [0.875, 0.89, 0.91, 0.93, 0.95]
    grid_dn = [0.85, 0.84, 0.83, 0.825]
    rows = []
    for direction in (grid_dn, grid_up):
        prev = sol0
        for Om in direction:
            try:
                s = solve_fixed_omega(Om, R=90.0, guess=prev, tol=1e-9)
            except RuntimeError as e:
                print('fail', Om, e); break
            if s.sol(0)[0] < 0.3:
                print('trivial', Om); break
            prev = s
            q = charge(s, Om); Q = q / G**2
            bg = Bg(s, Om)
            x = np.linspace(0, s.x[-1], 20001); F = s.sol(x)[0]
            rhalf = float(x[np.argmax(F < 0.5 * F[0])])
            row = dict(Omega=Om, Q=Q, F0=float(F[0]), S0=float(s.sol(0)[2]), E_over_QmX=float(energy(s, Om) / q), r_half=rhalf)
            for d in (5, 6):
                r = solve(bg, d, R=60.0)
                rr = rates(r, Q)
                pr = pair_rate(bg, d, Q)
                row[f'd{d}_single_frac'] = rr['frac_rate_per_s']
                row[f'd{d}_single_flux'] = r['flux'].tolist()
                row[f'd{d}_pair_frac'] = pr.get('frac_core_rate_free_born', 0.0)
                row[f'd{d}_rel_work'] = r['rel_work']
            rows.append(row)
            print(json.dumps(row), flush=True)
    rows.sort(key=lambda r: r['Omega'])
    json.dump(dict(rows=rows), open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
