#!/usr/bin/env python3
"""Skeptic of inf_higgs: D with the track's S1 hard-gas induced mass, on my own mode runs, with the field grid
hbar in [h_max, 20 h_max] (the track's) and extended to 1000 h_max.

Why: at M_t >~ 173.5 GeV the barrier h_max (2e10 GeV at 173.45, 1.8e9 at 174.6) is far below the end-of-inflation
spread (sigma ~ 0.1-0.2 H_i ~ 2e10 GeV), so regions with hbar >> 20 h_max are abundant (x << 3N).  The track's
'-inf' entries (no hbar in its grid rolls) are then a truncation of the grid, not an absence of rolling regions.
  python3 -B sk_plasmaD.py            -> results/sk_plasmaD.json
"""
import sys
sys.dont_write_bytecode = True
import json

import numpy as np

import sk_common as C
import sk_bench as B
import sk_aterm as AT

MT_EXT = [172.57, 172.86, 173.15, 173.45, 173.75, 174.05, 174.35, 174.6]


def main():
    P = json.load(open(AT.PLASMA_FULL))
    out = {}
    for name in ('HT6t', 'T1', 'E1'):
        ev = B.load(B.key(name, 'MD', 0))
        ev['Rbar'] = C.period_mean(ev, 1)
        tau = ev['rows'][:, 0]
        o = {}
        for row, r in P[name]['rows'].items():
            c = r['curve']
            mp2 = np.interp(tau, c['t'], c['m2'], left=0.0, right=c['m2'][-1])
            res = {}
            for hm in (20.0, 1000.0):
                Ds, hb = [], []
                for Mt in MT_EXT:
                    b = C.danger(ev, C.sm_in(Mt), ev['Rbar'], Ns=(60.0,), mp2=mp2, hmult=hm, nh=(200 if hm == 20 else 500))[60.0]
                    Ds.append(b['D']); hb.append(b.get('hbar_over_hmax'))
                fin = [d if np.isfinite(d) else -1e9 for d in Ds]
                res[f'hmult={hm:g}'] = dict(D60=Ds, hbar_over_hmax=hb, Mt_crit=C.crit(MT_EXT, fin))
            o[row] = dict(y=r['y'], MN_over_mI=r['MN_over_mI'], **res)
            print(name, row, 'track grid', [round(x, 1) for x in res['hmult=20']['D60']],
                  '| extended', [round(x, 1) for x in res['hmult=1000']['D60']], 'crit', res['hmult=1000']['Mt_crit'], flush=True)
        out[name] = o
        C.jdump(out, C.RES / 'sk_plasmaD.json')


if __name__ == '__main__':
    main()
