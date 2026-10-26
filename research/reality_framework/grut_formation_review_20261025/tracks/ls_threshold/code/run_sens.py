#!/usr/bin/env python3
"""ls_threshold driver, own code: sensitivity of R_L and of the post-freeze-out erasure to the benchmark's
quoted inputs (conversion sigma v, 4->2 coefficients, lump binding mu_c), and an R scan for re-entrant windows.

Usage: OMP_NUM_THREADS=1 python3 -B run_sens.py <out.json>
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import lsrelic as LR
from run_proxy import C


def scan(a):
    cand, ym, R = a
    r = LR.proxy_run(C[cand]['eps'], R, ym, N=400)
    return dict(kind='scan', candidate=cand, ymult=ym, R=R, f=r['f'], sigma=r['sigma'], max_cond=r['max_cond'])


def thr(a):
    cand, ym, sv_fac, a_fac = a
    LR.SV *= sv_fac; LR.AA *= a_fac; LR.AB *= a_fac
    r = LR.threshold(C[cand]['eps'], ym, 40.0, 400.0, 11, N=400)
    r.update(kind='threshold', candidate=cand, sv_factor=sv_fac, a42_factor=a_fac)
    return r


def ke_mu(a):
    cand, uc = a
    LR.EC = LR.EC / LR.UC * uc; LR.UC = uc
    r = LR.ke_run(C[cand]['eps'], C[cand]['R_lower'], C[cand]['ym'], N=600, antilumps=True)
    r.update(kind='ke_mu_c', candidate=cand, mu_c=uc, E_c=LR.EC)
    return r


def sv_runs(a):
    mode, cand, fac, endT = a
    LR.SV *= fac
    if mode == 'ke':
        r = LR.ke_run(C[cand]['eps'], C[cand]['R_lower'], C[cand]['ym'], N=600, antilumps=True, endT=endT)
    else:
        r = LR.proxy_run(C[cand]['eps'], C[cand]['R_lower'], C[cand]['ym'], N=600, endT=endT)
    r.update(kind='sv_' + mode, candidate=cand, sv_factor=fac)
    return r


if __name__ == '__main__':
    jobs = [(thr, (10, C[10]['ym'], 72.0, 1.0)), (thr, (10, C[10]['ym'], 0.1, 1.0)),
            (thr, (10, C[10]['ym'], 1.0, 10.0)), (thr, (10, C[10]['ym'], 1.0, 0.1)),
            (thr, (17, C[17]['ym'], 72.0, 1.0)),
            (ke_mu, (10, 0.6)), (ke_mu, (10, 0.3)),
            (sv_runs, ('ke', 10, 100.0, 1e-8)), (sv_runs, ('proxy', 10, 100.0, 1e-9)), (sv_runs, ('proxy', 10, 1.0, 1e-9))]
    jobs += [(scan, (10, C[10]['ym'], R)) for R in (60, 80, 95, 100, 106, 113.069, 130, 160, 199.404, 240, 300, 400)]
    rows = []
    with Pool(4, maxtasksperchild=1) as p:   # fresh worker per job: globals are modified
        res = [(f.__name__, a, p.apply_async(f, (a,))) for f, a in jobs]
        for name, a, r in res:
            try:
                rows.append(r.get())
            except Exception as exc:          # recorded, not hidden
                rows.append(dict(kind='failed', job=name, args=list(a), status=repr(exc)))
            Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
    Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
    for r in rows:
        if r['kind'] == 'failed':
            print('FAILED', r)
        elif r['kind'] == 'scan':
            print('scan', r['R'], 'f=%.4g' % r['f'], 'maxcond %.3g' % r['max_cond'])
        elif r['kind'] == 'threshold':
            print('R_L', r['candidate'], 'sv x', r['sv_factor'], 'A42 x', r['a42_factor'], r['bracket'])
        elif r['kind'].startswith('sv_'):
            print(r['kind'], 'sv x', r['sv_factor'], 'endT', r['endT'], 'f=%.4g' % r['f'])
        else:
            print('ke mu_c', r['mu_c'], 'f=%.4g' % r['f'], 'p=%.4g' % r['p_final'])
