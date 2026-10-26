#!/usr/bin/env python3
"""ls_threshold driver, own code: premise probe with the two-potential, energy-conserving hidden state (ke mode)
and the on-shell light-scalar route, at the conservative (no cascade) ledgers and raw yields of b448f25.

Usage: OMP_NUM_THREADS=1 python3 -B run_ke.py <out.json>
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import lsrelic as LR
from run_proxy import C


def job(a):
    cand, name, ym, anti, s_ch, N = a
    r = LR.ke_run(C[cand]['eps'], C[cand][name], ym, N=N, antilumps=anti, s_channel=s_ch)
    r.update(candidate=cand, ledger=name)
    return r


if __name__ == '__main__':
    jobs = [(10, 'R_lower', C[10]['ym'], True, True, 600), (17, 'R_lower', C[17]['ym'], True, True, 600),
            (17, 'R_lower', C[17]['ym'], True, False, 600), (10, 'R_lower', C[10]['ym'], False, False, 600),
            (10, 'R_lower', C[10]['ym'], True, False, 1200), (10, 'R_Y', C[10]['ym'], True, False, 600),
            (17, 'R_Y', C[17]['ym'], True, False, 600), (10, 'R_lower', 1.0, True, False, 600)]
    with Pool(4) as p:
        rows = p.map(job, jobs)
    Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
    for r in rows:
        print(r['candidate'], r['ledger'], round(r['ymult'], 4), 'anti' if r['antilumps'] else 'noanti',
              's' if r['s_channel'] else '-', r['N'], 'p=%.4g f=%.4g sigma %.3g->%.3g' % (r['p_final'], r['f'], r['sigma0'], r['sigma_final']),
              'runtime %.0f' % r['runtime'])
