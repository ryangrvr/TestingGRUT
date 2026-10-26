#!/usr/bin/env python3
"""ls_critic (Claude): cross-track check.  Does the ls_scaling heavy-mediator window survive the ls_threshold
energy-conserving (ke) premise, and is ls_scaling's interpolated benchmark R_L = 88 right?

ls_scaling judged its heavy points (eps = -0.1457, yield 1.00 Y_Q after the scaling symmetry) against the benchmark
proxy R_L ~ 88 read from a coarse sp_critic table.  ls_threshold found that an energy-conserving hidden state (ke)
moves the case-L crossing at candidate 10 (raw density) from R ~ 103 to R ~ 10.  This driver evaluates, at the heavy
point's eps and normalised yield:
  * ke f(R) at the m = M_Y/H scan values R = 0.25 (m=40), 1.1445 (m=30, base point), 8.9 (m=20), 34.2 (m=15),
    plus R = 3, 5, 15 to locate the ke crossing; mu_c = E_c = 0.861195 m_X (lumps grow) and 0.923316 (new lumps);
    the on-shell s route on/off;
  * the benchmark proxy threshold by direct bisection (replaces the log-interpolated table value).
Code: imports tracks/ls_threshold/code/lsrelic.py unchanged (our own code); no bytecode is written there.
Usage: OMP_NUM_THREADS=1 python3 -B heavy_threshold.py <out.json>     (~20-30 min on 4 cores)
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'ls_threshold' / 'code'))
import lsrelic as LR  # noqa: E402

EPS_HEAVY = 0.14571744147282          # |eps| at the ls_scaling base point (analytic.json scaling_symmetry) [Q, ours]
YM = 1.0                              # yield normalised to Y_Q by the scaling symmetry
W, A = 0.8611951308217225, 0.9233161353726057
R_SCAN = {'m40': 0.25, 'm30_base': 1.1444536115495259, 'm20': 8.9, 'm15': 34.22137642047465}


def ke(a):
    R, uc, s_ch, N = a
    LR.UC = uc; LR.EC = uc
    r = LR.ke_run(EPS_HEAVY, R, YM, N=N, antilumps=True, s_channel=s_ch)
    r.pop('history')
    r.update(kind='ke', mu_c=uc, E_c=uc)
    return r


def proxy_thr(a):
    lo, hi, it = a
    r = LR.threshold(EPS_HEAVY, YM, lo, hi, it, N=400)
    r['kind'] = 'proxy_threshold'
    return r


if __name__ == '__main__':
    jobs = [(proxy_thr, (40.0, 200.0, 9))]
    jobs += [(ke, (R, W, False, 500)) for R in (0.25, R_SCAN['m30_base'], 3.0, 5.0, 8.9, 15.0, R_SCAN['m15'])]
    jobs += [(ke, (R, W, True, 500)) for R in (8.9, 15.0, R_SCAN['m15'])]
    jobs += [(ke, (R, A, False, 500)) for R in (R_SCAN['m30_base'], 3.0, 8.9)]
    rows = []
    with Pool(4, maxtasksperchild=1) as p:     # fresh worker per job: lsrelic globals are modified
        res = [(f.__name__, a, p.apply_async(f, (a,))) for f, a in jobs]
        for name, a, r in res:
            try:
                rows.append(r.get())
            except Exception as exc:            # recorded, not hidden
                rows.append(dict(kind='failed', job=name, args=list(a), status=repr(exc)))
            Path(sys.argv[1]).write_text(json.dumps(dict(method=__doc__, eps=EPS_HEAVY, ymult=YM, R_scan=R_SCAN,
                                                         rows=rows), indent=1, default=float))
            x = rows[-1]
            if x['kind'] == 'ke':
                print('ke R=%.4g mu_c=%.4f s=%s f=%.4g p=%.4g fallbacks=%d ew=%.2e' % (
                    x['R'], x['mu_c'], x['s_channel'], x['f'], x['p_final'], x['root_fallbacks'],
                    x['energy_work_cumulative_residual']), flush=True)
            else:
                print(x, flush=True)
