#!/usr/bin/env python3
"""kin_thermo (Claude): the two old case-L proxies evaluated at the heavy spot point for R = 1-90 (task item 4).

REUSE, stated: this driver imports our own earlier track code tracks/ls_threshold/code/lsrelic.py (and its
lsthermo.py) UNCHANGED, with no bytecode written into that directory.  Nothing here is new physics; it is the
reference curve f_proxy(R) the kinetic solver has to be compared with.

  proxy : single chemical potential, isentropic hidden state, lumps at mu_c = 0.861195 m (E_c = E/Q = 0.923316 m
          booked, as the benchmark specified), relaxation Gamma = 4->2 + conversion (+ s route).
  ke    : two chemical potentials, every reaction conserves energy, lumps/anti-lumps at mu_c with the consistent
          E_c = mu_c (ls_threshold run_ke_consistent choice), N = 500 steps, optional s route.
Point: |eps| = 0.14571744282991028 (INPUTS.json 225028b), yield 1.00 Y_Q, g* coarse table of lsrelic.
Usage: OMP_NUM_THREADS=1 python3 -B old_proxies.py <out.json> [R1,R2,...]   (~1.5 h on 4 cores; resumable: finished
       rows already in <out.json> are kept and their jobs skipped; the optional comma-separated R list restricts the
       run to those R values, for finishing an interrupted run in short foreground batches)
"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'ls_threshold' / 'code'))
import lsrelic as LR  # noqa: E402  (our ls_threshold code, unchanged)

EPS = 0.14571744282991028
YM = 1.0
W, A = 0.8611951308217225, 0.9233161353726057
R_GRID = [1.0, 1.14445, 2.0, 3.0, 5.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 40.5626, 41.707,
          48.499, 49.643, 60.0, 75.0, 90.0]


def proxy(a):
    R, s_ch = a
    t = time.time()
    try:
        r = LR.proxy_run(EPS, R, YM, N=800, s_channel=s_ch)
        r.update(kind='proxy', mu_c=W, E_c=A, seconds=time.time() - t)
    except Exception as exc:  # recorded, not hidden
        r = dict(kind='proxy', R=R, s_channel=s_ch, status='failure %r' % exc)
    return r


def ke(a):
    R, uc, s_ch, N = a
    LR.UC = uc
    LR.EC = uc
    t = time.time()
    try:
        r = LR.ke_run(EPS, R, YM, N=N, antilumps=True, s_channel=s_ch)
        r.pop('history', None)
        r.update(kind='ke', mu_c=uc, E_c=uc, seconds=time.time() - t)
    except Exception as exc:
        r = dict(kind='ke', R=R, mu_c=uc, s_channel=s_ch, N=N, status='failure %r' % exc)
    return r


def thr(a):
    lo, hi, it, s_ch = a
    r = LR.threshold(EPS, YM, lo, hi, it, N=400, s_channel=s_ch)
    r.update(kind='proxy_threshold', s_channel=s_ch)
    return r


def _key(r):
    """Identity of a finished row, for resuming an interrupted run."""
    k = r.get('kind')
    if k == 'proxy_threshold':
        return (k, bool(r.get('s_channel')))
    if k == 'ke':
        return (k, round(float(r['R']), 6), round(float(r['mu_c']), 6), bool(r.get('s_channel')))
    if k == 'proxy':
        return (k, round(float(r['R']), 6), bool(r.get('s_channel')))
    return None


def _job_key(f, a):
    if f is thr:
        return ('proxy_threshold', bool(a[3]))
    if f is ke:
        return ('ke', round(float(a[0]), 6), round(float(a[1]), 6), bool(a[2]))
    return ('proxy', round(float(a[0]), 6), bool(a[1]))


if __name__ == '__main__':
    out = Path(sys.argv[1])
    jobs = [(thr, (40.0, 200.0, 10, False)), (thr, (40.0, 200.0, 10, True))]
    jobs += [(ke, (R, W, False, 500)) for R in R_GRID]
    jobs += [(ke, (R, W, True, 500)) for R in R_GRID]
    jobs += [(ke, (R, A, False, 500)) for R in (3.0, 5.0, 7.0, 8.0, 10.0, 15.0)]
    jobs += [(proxy, (R, s)) for s in (False, True) for R in R_GRID]
    if len(sys.argv) > 2:   # optional R subset (threshold jobs excluded)
        want = {round(float(x), 6) for x in sys.argv[2].split(',')}
        jobs = [(f, a) for f, a in jobs if f is not thr and round(float(a[0]), 6) in want]
    rows = []
    # resume: keep every finished (non-failed) row already in the output file and skip its job
    if out.exists():
        try:
            rows = [r for r in json.loads(out.read_text())['rows'] if 'status' not in r]
        except Exception as exc:
            print('could not read existing output, starting afresh:', exc)
            rows = []
        done = {_key(r) for r in rows}
        jobs = [(f, a) for f, a in jobs if _job_key(f, a) not in done]
        print('resuming: %d rows kept, %d jobs to run' % (len(rows), len(jobs)), flush=True)
    with Pool(4, maxtasksperchild=1) as p:   # fresh worker per job: lsrelic globals are modified by ke()
        res = [(f.__name__, a, p.apply_async(f, (a,))) for f, a in jobs]
        for name, a, r in res:
            try:
                rows.append(r.get())
            except Exception as exc:
                rows.append(dict(kind='failed', job=name, args=list(a), status=repr(exc)))
            out.write_text(json.dumps(dict(method=__doc__, eps=EPS, ymult=YM, R_grid=R_GRID, rows=rows),
                                      indent=1, default=float))
            x = rows[-1]
            print(x.get('kind'), {k: x.get(k) for k in ('R', 'mu_c', 's_channel', 'f', 'p_final', 'root_fallbacks',
                                                        'energy_work_cumulative_residual', 'bracket', 'seconds')},
                  flush=True)
