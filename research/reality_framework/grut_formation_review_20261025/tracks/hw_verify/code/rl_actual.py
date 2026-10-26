#!/usr/bin/env python3
"""hw_verify (Claude, own driver; reuses our sp_yledger/gp_leftovers benchmark code, as the brief permits):
case-L f(R) at the heavy spot point's actual eps, at the normalised yield (ymult = 1) and at the raw yield of the
ladder point at H_i = 4.3487e11 GeV (ymult = raw Y/Y_Q), to bracket R_L "at the actual eps and yield".
Usage: python3 -B rl_actual.py <out.json> eps ymult1:R,R,.. ymult2:R,R,..   (~1-4 min per evaluation, 4 cores)
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
sys.dont_write_bytecode = True
YL = Path(__file__).resolve().parents[2] / 'sp_yledger' / 'code'
sys.path.insert(0, str(YL))
import thermal_relic as TR  # noqa: E402


def job(a):
    eps, ym, R = a
    try:
        r = TR.run_inj(eps, 'L', R, ym)
        return dict(eps=eps, ymult=ym, R=R, f=float(r['f_extra']), sigma=float(r['sigma']),
                    T_h0_over_m=float(r['T_h0_over_m']), cond_frac0=float(r['cond_frac0']))
    except Exception as exc:
        return dict(eps=eps, ymult=ym, R=R, status=repr(exc))


if __name__ == '__main__':
    out = Path(sys.argv[1]); eps = float(sys.argv[2])
    jobs = []
    for spec in sys.argv[3:]:
        ym, rs = spec.split(':')
        jobs += [(eps, float(ym), float(R)) for R in rs.split(',')]
    jobs.sort(key=lambda j: -j[2])
    rows = []
    with Pool(4) as p:
        for r in p.imap_unordered(job, jobs):
            rows.append(r); print(r, flush=True)
            out.write_text(json.dumps(dict(method=__doc__, rows=rows), indent=1))
