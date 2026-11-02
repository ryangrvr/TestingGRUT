#!/usr/bin/env python3
"""dq_critic item 4 (Claude, own driver, 2 Nov 2026): Higgs survival of xi_EW = 0.1 at the recalibrated H_i.

Uses the higgs_reheat SKEPTIC's normalisation-fixed copy of the track's mode code (skeptic/code/modes_canon.py),
imported read-only and unchanged: k-resolved Higgs modes through inflation's end, the oscillating-R era and reheating
for (3, 0.5), gravitational only (no plasma), danger exponent D60 of the scale-resolved criterion [E criterion, C runs].
New here: D60 at the H_i of this pass (free X 1.2923e11, Q-matter 1.4908e11, free X with M_P portals 1.0724e11) and
the top mass M_t,crit(H_i) at which D60 = 0 (bisection), i.e. the largest M_t for which xi_EW = 0.1 survives
gravitationally under that criterion.
usage: higgs_check.py points | mtcrit <H> [Mt_lo Mt_hi]   -> appends to ../results/higgs_check.json
"""
import json, os, sys, time
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
SK = os.path.abspath(os.path.join(HERE, '..', '..', 'higgs_reheat', 'skeptic', 'code'))
sys.path.insert(0, SK)
import modes_canon as M

OUT = os.path.join(HERE, '..', 'results', 'higgs_check.json')


def D60(Mt, H, cHs=()):
    lam = M.Lam(Mt)
    xb = float(lam.sm.xi(H, 0.1))
    r, _ = M.run(xb, Mt=Mt, H_GeV=H, lam=lam, nsave=800, cHs=cHs)
    d = dict(Mt=Mt, H=H, xi_H=xb, D60=r['danger']['best'][60]['D'])
    for k, v in r.get('danger_cH', {}).items():
        d['D60_' + k] = v[60]['D']
    return d


def save(key, val):
    data = json.load(open(OUT)) if os.path.exists(OUT) else dict(method=__doc__)
    data[key] = val
    json.dump(data, open(OUT, 'w'), indent=1, default=float)


def main():
    mode = sys.argv[1]
    t0 = time.time()
    if mode == 'points':
        rows = []
        for Mt, H in ((172.60, 1.2923e11), (172.87, 1.2923e11), (172.60, 1.4908e11), (172.87, 1.4908e11), (172.87, 1.0724e11)):
            d = D60(Mt, H, cHs=(0.5, 2.0)); rows.append(d)
            print(d, '%.0fs' % (time.time() - t0), flush=True)
        save('points', rows)
    elif mode == 'mtcrit':
        H = float(sys.argv[2])
        hist = []
        def f(Mt):
            d = D60(Mt, H); hist.append(d); print(d, '%.0fs' % (time.time() - t0), flush=True)
            return d['D60']
        lo = float(sys.argv[3]) if len(sys.argv) > 3 else 172.60
        hi = float(sys.argv[4]) if len(sys.argv) > 4 else 172.87
        Mc = brentq(f, lo, hi, xtol=4e-3)
        save('mtcrit_H%.4e' % H, dict(H=H, Mt_crit=Mc, evals=hist))
        print('M_t,crit(H = %.4e) = %.3f GeV' % (H, Mc))


if __name__ == '__main__':
    main()
