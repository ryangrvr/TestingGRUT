#!/usr/bin/env python3
"""ls_threshold driver, own code: ke mode at small injection (where the f = 0.01 crossing lies) and the
energy-work diagnostics of ke mode at the conservative ledger.  Candidate 10, raw yield.
Usage: OMP_NUM_THREADS=1 python3 -B run_ke_low.py <out.json>
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import lsrelic as LR
from run_proxy import C


def job(a):
    R, N = a
    try:
        r = LR.ke_run(C[10]['eps'], R, C[10]['ym'], N=N, antilumps=True)
    except Exception as exc:
        return dict(R=R, N=N, status=repr(exc))
    if R not in (10.0, C[10]['R_lower']):
        r.pop('history')
    return r


if __name__ == '__main__':
    jobs = [(0.0, 500), (2.0, 500), (5.0, 500), (10.0, 500), (C[10]['R_lower'], 600), (5.0, 1000)]
    with Pool(4) as p:
        rows = p.map(job, jobs)
    Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1, default=float))
    for r in rows:
        print('R', r['R'], 'N', r['N'], 'f=%s' % r.get('f'), 'ew', r.get('energy_work_max_step_residual'),
              r.get('energy_work_cumulative_residual'), r.get('status', ''))
