#!/usr/bin/env python3
"""Skeptic of inf_higgs: post-inflation effect of an R-odd A phi_I |H|^2 on the asymmetric benchmarks.

The track bounds A only during inflation (Hawking-Moss, |A| <= 12 / 34 / 35 TeV at +1 sigma).  On HT6t and E1 the
oscillation is asymmetric, so the period mean <phi_I> != 0 and, since <phi_I> ~ Phi^2 ~ rho ~ H^2 in the matter era,
the A-term leaves a DC Higgs mass A <phi_I> = kappa_A H^2 with kappa_A constant, of the SAME sign as A phi_I during
inflation (phi_I < 0 throughout for HT6t, > 0 for E1).  On the destabilising branch it subtracts from xi <R> = 3 xi H^2.
  python3 -B sk_aterm.py kappa          # kappa_A per TeV from the background (period means)
  python3 -B sk_aterm.py run <name> <A_TeV>   # mode run (A in modes and criterion), D(M_t) -> results/sk_aterm.json
A > 0 on HT6t and A < 0 on E1 is the destabilising sign (A phi_I < 0 during inflation).
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np

import sk_common as C
import sk_bench as B

MTS = [172.57, 172.70, 172.80, 172.86, 172.90, 172.95, 173.0]


def cmd_run(name, A_TeV):
    pot = C.Pot(name)
    Acode = A_TeV * 1e3 / pot.H_GeV
    t0 = time.time()
    ev = B.run(name, 'MD', 0, Aterm=Acode, tag='', save=True)
    g = B.Dgrid(ev, Ns=(60.0,), Mts=MTS)
    post = ev['rows'][:, 0] > ev['tau_end_inf'] + 0.5
    kap = ev['m2xbar'][post] / ev['rows'][post, 2] ** 2
    o = dict(name=name, A_TeV=A_TeV, A_code=Acode, D60={str(m): g[m][60.0]['D'] for m in MTS},
             tau_c={str(m): g[m][60.0].get('tau_c') for m in MTS},
             Mt_crit_60=C.crit(MTS, [g[m][60.0]['D'] for m in MTS]),
             kappa_A_median=float(np.median(kap)), kappa_A_range=[float(np.min(kap)), float(np.max(kap))],
             m2_extra_at_S1=float(Acode * (pot.rec['phiI_S1_MP']) / pot.eta), seconds=time.time() - t0)
    path = C.RES / 'sk_aterm.json'
    try:
        allr = json.load(open(path))
    except Exception:
        allr = {}
    allr[f'{name}|A={A_TeV}TeV'] = o
    C.jdump(allr, path)
    print(name, 'A %+.1f TeV' % A_TeV, 'kappa %.4f' % o['kappa_A_median'], 'm2(S1) %.3f' % o['m2_extra_at_S1'],
          {m: round(o['D60'][str(m)], 1) for m in MTS}, 'crit', o['Mt_crit_60'], '%.0fs' % o['seconds'], flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        cmd_run(sys.argv[2], float(sys.argv[3]))


MT_EXT = [172.57, 172.86, 173.15, 173.45, 173.75, 174.05, 174.35, 174.6]
PLASMA = C.TRACKS / 'inf_higgs' / 'results' / 'plasma.json'          # track's thinned curves (read-only)
PLASMA_FULL = C.SP.parent / 'inf_higgs' / 'plasma_full.json'         # track's full curves (scratchpad, read-only)


def cmd_plasma(name, A_TeV, row):
    """D with both the A-term (mode run) and the track's S1 hard-gas induced mass m_ind^2(tau) of plasma row `row`"""
    P = json.load(open(PLASMA_FULL if PLASMA_FULL.exists() else PLASMA))
    c = P[name]['rows'][row]['curve']
    k = B.key(name, 'MD', 0, None, A_TeV * 1e3 / C.Pot(name).H_GeV, 0.0, '')
    ev = B.load(k)
    ev['Rbar'] = C.period_mean(ev, 1)
    ev['m2xbar'] = C.period_mean(ev, 5)
    tau = ev['rows'][:, 0]
    mp2 = np.interp(tau, c['t'], c['m2'], left=0.0, right=c['m2'][-1])
    Ds = []
    for Mt in MT_EXT:
        b = C.danger(ev, C.sm_in(Mt), ev['Rbar'], Ns=(60.0,), m2xbar=ev['m2xbar'], mp2=mp2, hmult=1000.0, nh=500)
        Ds.append(b[60.0]['D'])
    o = dict(name=name, A_TeV=A_TeV, plasma_row=row, Mt=MT_EXT, D60=Ds, Mt_crit=C.crit(MT_EXT, Ds),
             A_ind_track_GeV=P[name]['rows'][row]['A_ind_GeV'])
    path = C.RES / 'sk_aterm_plasma.json'
    try:
        allr = json.load(open(path))
    except Exception:
        allr = {}
    allr[f'{name}|A={A_TeV}|{row}'] = o
    C.jdump(allr, path)
    print(name, A_TeV, row, [round(x, 1) for x in Ds], 'crit', o['Mt_crit'], flush=True)


if __name__ == '__main__' and sys.argv[1] == 'plasma':
    cmd_plasma(sys.argv[2], float(sys.argv[3]), sys.argv[4])
