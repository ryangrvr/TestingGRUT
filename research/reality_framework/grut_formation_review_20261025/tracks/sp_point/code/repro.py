"""sp_point: independent reproduction of Codex's leading point (373d37e handoff targets), own code (sprel.c).

Usage: python3 repro.py ../results/repro.json
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import spmodel as m

H = m.H_STAR
D0 = np.pi / 2
KEYS = ['epsX', 'NX', 'EX', 'JY', 'JYeig', 'YX_over_YQ', 'eps_of_means', 'pairs_per_charge', 'NY', 'Neff',
        'Neff_std', 'epsX_std', 'NX_std', 'deps', 'dNX', 'EYexc', 'EYtot']


def job(spec):
    lab, delta, meth, c, aend, wins, tang = spec
    z = m.initial_z(H, delta)
    t0 = time.time()
    r = m.run(H, z, a_end=aend, method=meth, cpar=c, windows=wins, ns=300,
              tangent=m.delta_tangent(H, delta) if tang else None)
    r['label'] = lab
    r['runtime_s'] = time.time() - t0
    r['windows'] = {k: {q: v.get(q) for q in KEYS} for k, v in r['windows'].items()}
    return r


def main(out):
    W = (1, 100, 300, 600, 1000)
    specs = [('Y8_c0.0125', D0, 8, 0.0125, 1000, W, False),
             ('Y8_c0.00625', D0, 8, 0.00625, 1000, W, False),
             ('Y8_c0.003125', D0, 8, 0.003125, 1000, W, False),
             ('DP5_1e-12', D0, 5, 1e-12, 1000, W, False),
             ('DP5_1e-13', D0, 5, 1e-13, 1000, W, False),
             ('DP5_1e-14', D0, 5, 1e-14, 1000, W, False),
             ('Y8_c0.00625_shift1e-9', D0 + 1e-9, 8, 0.00625, 1000, W, False),
             ('DP5_1e-14_shift1e-9', D0 + 1e-9, 5, 1e-14, 1000, W, False),
             ('Y8_c0.00625_tangent', D0, 8, 0.00625, 1000, W, True),
             ('Y8_c0.0125_a3000', D0, 8, 0.0125, 3000, (1, 300, 600, 1000, 2000, 3000), True),
             ('DP5_1e-13_a3000', D0, 5, 1e-13, 3000, (1, 300, 600, 1000, 2000, 3000), False)]
    with ProcessPoolExecutor(max_workers=4) as pool:
        res = list(pool.map(job, specs))
    summary = {}
    for r in res:
        summary[r['label']] = dict(controls=r['controls'], runtime_s=r['runtime_s'], windows=r['windows'])
        w = r['windows']['600-1000']
        print(r['label'], '%.1fs' % r['runtime_s'], 'eps=%.10f NX=%.7f EX=%.7f JY=%.7f Y/YQ=%.6f' %
              (w['epsX'], w['NX'], w['EX'], w['JY'], w['YX_over_YQ']), r['controls'], flush=True)
    b = summary['Y8_c0.00625']['windows']['600-1000']
    s = summary['Y8_c0.00625_shift1e-9']['windows']['600-1000']
    b14 = summary['DP5_1e-14']['windows']['600-1000']
    s14 = summary['DP5_1e-14_shift1e-9']['windows']['600-1000']
    tang = summary['Y8_c0.00625_tangent']['windows']['600-1000']
    a3 = summary['Y8_c0.0125_a3000']['windows']
    d3 = summary['DP5_1e-13_a3000']['windows']
    derived = dict(
        shift_1e9_delta_eps_Y8=s['epsX'] - b['epsX'],
        shift_1e9_delta_eps_DP5_1e14=s14['epsX'] - b14['epsX'],
        tangent_window_mean_deps_ddelta=tang['deps'],
        tangent_prediction_for_1e9=tang['deps'] * 1e-9,
        method_spread_eps=max(v['windows']['600-1000']['epsX'] for k, v in summary.items() if 'shift' not in k and 'a3000' not in k)
        - min(v['windows']['600-1000']['epsX'] for k, v in summary.items() if 'shift' not in k and 'a3000' not in k),
        DP5_1e12_vs_1e14_eps=summary['DP5_1e-12']['windows']['600-1000']['epsX'] - b14['epsX'],
        late_window_shift_Y8=a3['2000-3000']['epsX'] - a3['600-1000']['epsX'],
        late_window_shift_DP5=d3['2000-3000']['epsX'] - d3['600-1000']['epsX'],
        late_window_JY_shift_rel=a3['2000-3000']['JY'] / a3['600-1000']['JY'] - 1,
        codex_targets=dict(epsX=-0.24066124, NX=-238.64868, EX=741.64988, JY=13.26903, YX_over_YQ=-65.9823,
                           shift_1e9=2.6e-7, rtol_2e12_vs_2e14=6.7e-9, late_mean_shift=0.00046))
    print(json.dumps(derived, indent=1))
    json.dump(dict(point=dict(H_GeV=H, cH=1.0, delta=D0, amplitudes=list(m.minimum(H, 1.0)),
                              masses=m.masses(H, 1.0)),
                   runs=summary, derived=derived,
                   scope='Own C integrators (Yoshida-8 composition; Dormand-Prince 5(4)); Codex code not run.'),
              open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
