#!/usr/bin/env python3
"""ls_threshold driver, own code: ke mode with a thermodynamically consistent lump phase.

The benchmark's lump proxy puts condensed charge at chemical potential mu_c = omega(Q) = 0.861195 m_X but books
its energy at E/Q = 0.923316 m_X.  For a zero-entropy condensate the marginal energy per charge must equal its
chemical potential, otherwise every adiabatic condensation/evaporation cycle creates (E_c - mu_c) per charge
(seen as a linear-in-step energy-work residual in ke mode).  Here E_c = mu_c, bracketed by
mu_c = 0.861195 (lumps grow: omega) and mu_c = 0.923316 (new Q = 200000 lumps form: E/Q).
Also: benchmark (proxy) thresholds with each consistent choice.

Usage: OMP_NUM_THREADS=1 python3 -B run_ke_consistent.py <out.json>
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import lsrelic as LR
from run_proxy import C


def ke(a):
    cand, R, ym, uc, s_ch, N = a
    LR.UC = uc; LR.EC = uc
    r = LR.ke_run(C[cand]['eps'], R, ym, N=N, antilumps=True, s_channel=s_ch)
    r.update(kind='ke', candidate=cand, mu_c=uc, E_c=uc)
    if R not in (C[cand]['R_lower'],):
        r.pop('history')
    return r


def thr(a):
    cand, ym, uc = a
    LR.UC = uc; LR.EC = uc
    r = LR.threshold(C[cand]['eps'], ym, 20.0, 300.0, 11, N=400)
    r.update(kind='threshold', candidate=cand, mu_c=uc, E_c=uc)
    return r


if __name__ == '__main__':
    W, A = 0.8611951308217225, 0.9233161353726057
    c10, c17 = C[10], C[17]
    jobs = [(ke, (10, c10['R_lower'], c10['ym'], W, False, 600)), (ke, (17, c17['R_lower'], c17['ym'], W, False, 600)),
            (ke, (10, c10['R_lower'], c10['ym'], A, False, 600)), (ke, (17, c17['R_lower'], c17['ym'], A, False, 600)),
            (ke, (10, c10['R_lower'], c10['ym'], W, True, 600)), (ke, (17, c17['R_lower'], c17['ym'], W, True, 600)),
            (ke, (10, c10['R_lower'], c10['ym'], W, False, 1200)),
            (ke, (10, c10['R_Y'], c10['ym'], W, False, 600)), (ke, (17, c17['R_Y'], c17['ym'], W, False, 600)),
            (thr, (10, c10['ym'], W)), (thr, (10, c10['ym'], A))]
    jobs += [(ke, (10, R, c10['ym'], W, False, 500)) for R in (0.0, 3.0, 10.0, 30.0, 100.0)]
    rows = []
    with Pool(4, maxtasksperchild=1) as p:   # fresh worker per job: globals are modified
        res = [(f.__name__, a, p.apply_async(f, (a,))) for f, a in jobs]
        for name, a, r in res:
            try:
                rows.append(r.get())
            except Exception as exc:          # recorded, not hidden
                rows.append(dict(kind='failed', job=name, args=list(a), status=repr(exc)))
            Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
    for r in rows:
        if r['kind'] == 'ke':
            print('ke', r['candidate'], 'R=%.3f' % r['R'], 'ym %.3f' % r['ymult'], 'mu_c %.4f' % r['mu_c'],
                  's' if r['s_channel'] else '-', r['N'], 'f=%.4g' % r['f'],
                  'ew %.2e %.2e' % (r['energy_work_max_step_residual'], r['energy_work_cumulative_residual']),
                  'fallbacks', r['root_fallbacks'])
        elif r['kind'] == 'threshold':
            print('R_L proxy consistent', r['candidate'], r['mu_c'], r['bracket'])
        else:
            print('FAILED', r)
