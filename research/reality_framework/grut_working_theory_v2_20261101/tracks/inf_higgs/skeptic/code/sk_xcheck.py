#!/usr/bin/env python3
"""Skeptic of inf_higgs: cross-checks on the TRACK's cached mode runs (scratchpad, read-only).
  (a) my danger function on the track's HT6t MD run with the track's S1 plasma curves (does my criterion code reproduce
      the track's plasma D?  then: my mode run vs the track's mode run);
  (b) the track's R-even portal runs (lam_min, 0.1 lam_min) with the field grid extended from 20 to 1000 h_max.
  -> results/sk_xcheck.json"""
import sys
sys.dont_write_bytecode = True
import json

import numpy as np

import sk_common as C
import sk_bench as B
import sk_aterm as AT

RUNS = C.SP.parent / 'inf_higgs' / 'runs'


def trackrun(fn):
    z = np.load(RUNS / fn, allow_pickle=True)
    res = json.loads(str(z['res']))
    ev = dict(rows=z['rows'][:, :4], spec=z['spec'], ks=z['ks'], wq=z['wq'], tau_end_inf=res['tau_endinf'], xi=res['xi'],
              H_GeV=1.29e11)
    return ev, z


def main():
    out = {}
    P = json.load(open(AT.PLASMA_FULL))
    evT, z = trackrun('HT6t_MD_iH0_xi0.0877278_Nb5.0_n260_KNone_lam0_A0_mt172.57_te80.0.npz')
    ev = B.load(B.key('HT6t', 'MD', 0)); ev['Rbar'] = C.period_mean(ev, 1)
    a = {}
    for row in ('MN=0.25 y=3e-05', 'MN=0.0625 y=3e-05'):
        c = P['HT6t']['rows'][row]['curve']
        for Mt in (172.86, 173.45):
            smi = C.sm_in(Mt)
            mpT = np.interp(evT['rows'][:, 0], c['t'], c['m2'], left=0.0, right=c['m2'][-1])
            mp = np.interp(ev['rows'][:, 0], c['t'], c['m2'], left=0.0, right=c['m2'][-1])
            dT = C.danger(evT, smi, z['Rb'], Ns=(60.0,), mp2=mpT)[60.0]['D']
            dM = C.danger(ev, smi, ev['Rbar'], Ns=(60.0,), mp2=mp)[60.0]['D']
            a[f'{row} {Mt}'] = dict(my_criterion_on_track_run=dT, my_run=dM)
            print('plasma', row, Mt, 'track run %.1f  my run %.1f' % (dT, dM), flush=True)
    out['plasma_D_track_run_vs_mine'] = a
    b = {}
    for lam in ('3.36e-11', '3.36e-12'):
        evP, zP = trackrun(f'HT6t_MD_iH0_xi0.0877278_Nb5.0_n260_KNone_lam{lam}_A0_mt172.57_te80.0.npz')
        for Mt in (172.86, 173.75, 174.6):
            smi = C.sm_in(Mt)
            d20 = C.danger(evP, smi, zP['Rb'], Ns=(60.0,), m2xbar=zP['m2xbar'], hmult=20.0)[60.0]['D']
            d1k = C.danger(evP, smi, zP['Rb'], Ns=(60.0,), m2xbar=zP['m2xbar'], hmult=1000.0, nh=500)[60.0]['D']
            b[f'lam={lam} Mt={Mt}'] = dict(grid20=d20, grid1000=d1k)
            print('portal', lam, Mt, 'grid 20: %.0f  grid 1000: %.0f' % (d20, d1k), flush=True)
    out['portal_grid'] = b
    C.jdump(out, C.RES / 'sk_xcheck.json')


if __name__ == '__main__':
    main()
