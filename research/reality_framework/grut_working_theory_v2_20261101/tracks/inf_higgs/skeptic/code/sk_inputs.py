#!/usr/bin/env python3
"""Skeptic of inf_higgs: D at the newer PDG inputs (M_t = 172.60 +- 0.27, M_h = 125.13, as higgs_reheat used) and the
N_*, n_s, A_s of the matter-era histories (T_RH = 1e11 GeV) from inf_potential's tables.  -> results/sk_inputs.json"""
import sys, math, json

sys.dont_write_bytecode = True
import numpy as np
import sk_common as C, sk_bench as B
o = {}
Mts = [172.60, 172.70, 172.80, 172.87, 172.90, 172.95, 173.0]
for hist in ('MD', 'G05'):
    ev = B.load(B.key('HT6t', hist, 0)); ev['Rbar'] = C.period_mean(ev, 1)
    Ds = [C.danger(ev, C.sm_in(m, 125.13), ev['Rbar'], Ns=(60.0,))[60.0]['D'] for m in Mts]
    cr = C.crit(Mts, Ds)
    o[hist] = dict(Mt=Mts, D60=Ds, Mt_crit=cr, sigma_PDG2026=(cr - 172.60) / 0.27 if cr else None)
    print(hist, [round(d, 1) for d in Ds], 'crit', cr, 'sigma(172.60+-0.27)', o[hist]['sigma_PDG2026'])
# N_* and n_s of the MD histories (T_RH = 1e11) from inf_potential's table (log interpolation between 1e10 and 1e12)
for n in ('HT6t', 'T1', 'E1'):
    r = json.load(open(C.TRACKS / 'inf_potential' / 'results' / f'bench_{n}.json'))['by_Hi'][0]
    a = {tuple(x['history'][1:]) if x['history'][0] == 'wT' else x['history'][0]: x for x in r['Nstar_alternatives']}
    lo, hi = a[(0.0, 1e10)], a[(0.0, 1e12)]
    f = 0.5
    o[n] = dict(Nstar=lo['Nstar'] + f * (hi['Nstar'] - lo['Nstar']), ns=lo['ns'] + f * (hi['ns'] - lo['ns']),
                As_fixed_shape=lo['As_at_fixed_shape'] + f * (hi['As_at_fixed_shape'] - lo['As_at_fixed_shape']))
    o[n]['ns_pull_Planck2018'] = (o[n]['ns'] - 0.9649) / 0.0042
    o[n]['ns_pull_PACT_LB'] = (o[n]['ns'] - 0.9743) / 0.0034
    o[n]['As_shift'] = o[n]['As_fixed_shape'] / 2.1e-9 - 1
    print(n, {k: round(v, 4) if abs(v) > 1e-3 else '%.3e' % v for k, v in o[n].items()})
C.jdump(o, C.RES / 'sk_inputs.json')
